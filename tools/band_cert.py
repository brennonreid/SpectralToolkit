#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
band_cert.py — rigorous band certificate using interval arithmetic.

Purpose:
  Given a window configuration (Gaussian notch) and a band specification,
  compute rigorous lower/upper bounds for |W(f)| on each band using
  interval subdivision, and certify a global band margin.

CLI (normalized v2.1):
  --window-config   : path to window JSON (from window_gen.py)
  --bands           : path to bands JSON
  --out             : output JSON path
  --dps             : decimal precision for mpmath
  --tol             : subdivision tolerance for band_min_bounds
  --max-parts       : maximum number of sub-intervals
  --tqdm            : show per-band progress (optional)

JSON (normalized v2.1 pattern):
  {
    "kind": "band_cert",
    "created_utc": "<ISO8601 UTC>",
    "inputs": {
      "window_config_path": "...",
      "bands_path": "...",
      "mode": "gauss",
      "sigma": "<string>",
      "k0": "<string>",
      "dps": <int>
    },
    "numbers": {
      "band_margin_lo": "<string>",
      "band_margin_hi": "<string>",
      "bands_count": <int>
    },
    "band_cert": {
      "band_margin": {
        "lo": "<string>",
        "hi": "<string>"
      },
      "per_band": [
        {
          "label": "...",
          "left": "<string>",
          "right": "<string>",
          "min_abs_lo": "<string>",
          "min_abs_hi": "<string>"
        },
        ...
      ],
      "status": "PASS" | "FAIL",
      "critical_band": {
        "left": "<string>",
        "right": "<string>"
      } | null
    },
    "PASS": true | false,
    "meta": {
      "tool": "band_cert",
      "dps": <int>,
      "sha256": "<digest>"
    }
  }

All real-valued numeric fields are serialized as strings.
"""

import sys
import json
import argparse
import io
import pathlib
import hashlib
from typing import List, Tuple, Dict
from datetime import datetime, timezone

try:
    from mpmath import mp, iv
except Exception as e:
    sys.exit(
        "[error] mpmath with interval arithmetic is required (mp, iv). "
        f"Please install a recent mpmath. Original error: {e!r}"
    )

try:
    from tqdm import tqdm as _tqdm
except Exception:
    _tqdm = None


# --- Generic utilities -----------------------------------------------------


def fail(msg: str):
    sys.exit(f"[error] {msg}")


def read_json(path: str):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        fail(f"failed to read JSON {path}: {e}")


def to_mpf(x):
    return mp.mpf(str(x))


def mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def now_utc_iso() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def write_json(obj: Dict, path: str):
    """
    Write JSON with deterministic sha256 stored under meta.sha256.

    The sha256 is computed over the canonical serialized form
    (sorted keys, UTF-8, indent=2).
    """
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    if "meta" not in obj or not isinstance(obj["meta"], dict):
        obj["meta"] = {}

    js = json.dumps(obj, sort_keys=True, ensure_ascii=False, indent=2)
    sha = hashlib.sha256(js.encode("utf-8")).hexdigest()
    obj["meta"]["sha256"] = sha

    tmp = str(p) + ".tmp"
    with io.open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, sort_keys=True)
        f.write("\n")
    pathlib.Path(tmp).replace(p)
    print(f"[ok] band_cert -> {path} sha256={sha}")


# --- Window model (canonical sigma/k0 only) --------------------------------


def make_window(window_js: Dict):
    """
    Expect a canonical window configuration of the form produced by window_gen:

      {
        "kind": "window",
        "mode": "gauss",
        "sigma": "...",
        "k0": "...",
        "window": {
          "mode": "gauss",
          "sigma": "...",
          "k0": "..."
        },
        ...
      }

    Either the top-level or the nested 'window' block may be present; both
    use canonical field names 'sigma' and 'k0'.
    """
    if not isinstance(window_js, dict):
        fail("window-config JSON must be an object")

    src = window_js.get("window") if isinstance(window_js.get("window"), dict) else window_js

    mode = src.get("mode", "gauss")
    if mode != "gauss":
        fail(f"unsupported window mode {mode!r}; only 'gauss' is implemented")

    if "sigma" not in src or "k0" not in src:
        fail("window-config must contain 'sigma' and 'k0'")

    try:
        sigma = mp.mpf(str(src["sigma"]))
        k0 = mp.mpf(str(src["k0"]))
    except Exception:
        fail("window-config sigma/k0 not numeric")

    if sigma <= 0:
        fail("window sigma must be > 0")
    if k0 <= 0:
        fail("window k0 must be > 0")

    def W_abs_iv_on(a: mp.mpf, b: mp.mpf):
        I = iv.mpf([a, b])
        g = iv.exp(-(I / sigma) ** 2)  # Gaussian envelope
        n = 1 - iv.exp(-(I / k0) ** 2)  # multiplicative notch
        w = g * n
        lo = mp.mpf(w.a if w.a is not None else 0)
        hi = mp.mpf(w.b if w.b is not None else 0)
        if lo < 0:
            lo = mp.mpf("0")
        if hi < 0:
            hi = mp.mpf("0")
        return lo, hi

    return W_abs_iv_on, sigma, k0, mode


# --- Band parsing ----------------------------------------------------------


def parse_bands_generic(bands_js) -> List[Tuple[mp.mpf, mp.mpf, str]]:
    """
    Accepts multiple shapes:
      1) list: [ {left,right,label?}, ... ]
      2) {"bands": [ ... ]}
      3) {"bands": { name: {left,right,...}, ... }}
      4) {"named_grids": { name: {left,right,...}, ... }}
    """

    def coerce_list(obj):
        out = []
        for i, b in enumerate(obj):
            if "left" not in b or "right" not in b:
                fail(f"band {i} missing 'left' or 'right'")
            L = to_mpf(b["left"])
            R = to_mpf(b["right"])
            if not (L < R):
                fail(f"band {i} has non-increasing interval [{L},{R}]")
            label = str(b.get("label", f"band_{i}"))
            out.append((L, R, label))
        return out

    if isinstance(bands_js, list):
        return coerce_list(bands_js)

    if isinstance(bands_js, dict):
        if "bands" in bands_js and isinstance(bands_js["bands"], list):
            return coerce_list(bands_js["bands"])
        for key in ("bands", "named_grids"):
            obj = bands_js.get(key)
            if isinstance(obj, dict) and obj:
                out = []
                for j, (name, spec) in enumerate(obj.items()):
                    if "left" not in spec or "right" not in spec:
                        fail(f"{key}['{name}'] missing left/right")
                    L = to_mpf(spec["left"])
                    R = to_mpf(spec["right"])
                    if not (L < R):
                        fail(f"{key}['{name}'] non-increasing [{L},{R}]")
                    out.append((L, R, str(name)))
                return out

    fail("bands JSON has no recognizable band list")


# --- Interval refinement ---------------------------------------------------


def band_min_bounds(
    W_abs_iv_on, L: mp.mpf, R: mp.mpf, max_parts=16384, tol=mp.mpf("1e-30")
):
    """
    Returns (lo, hi) bounds for min_{f in [L,R]} |W(f)| using interval subdivision.
    """
    from heapq import heappush, heappop

    pq = []

    def push(a, b):
        lo, hi = W_abs_iv_on(a, b)
        heappush(pq, (float(lo), a, b, lo, hi))

    push(L, R)
    best_lo = mp.inf
    best_hi = mp.inf
    parts = 1

    while pq and parts <= max_parts:
        _, a, b, lo, hi = heappop(pq)
        if lo < best_lo:
            best_lo = lo
        if hi < best_hi:
            best_hi = hi
        if hi - lo <= tol:
            pass
        else:
            mid = (a + b) / 2
            push(a, mid)
            push(mid, b)
            parts += 1
        if best_hi - best_lo <= tol:
            break

    return best_lo, best_hi


# --- Main ------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Rigorous band certificate using interval arithmetic."
    )
    ap.add_argument(
        "--window-config",
        required=True,
        help="Window configuration JSON (from window_gen.py).",
    )
    ap.add_argument(
        "--bands",
        required=True,
        help="Bands specification JSON.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output band_cert JSON path.",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=220,
        help="Decimal precision for mpmath.",
    )
    ap.add_argument(
        "--tol",
        type=str,
        default="1e-30",
        help="Subdivision tolerance for band_min_bounds.",
    )
    ap.add_argument(
        "--max-parts",
        type=int,
        default=16384,
        help="Maximum number of sub-intervals per band.",
    )
    ap.add_argument(
        "--tqdm",
        action="store_true",
        help="Show per-band progress bar.",
    )
    args = ap.parse_args()

    mp.dps = int(args.dps)
    try:
        _ = iv
    except Exception:
        fail("mpmath interval arithmetic unavailable (iv)")

    window_js = read_json(args.window_config)
    bands_js = read_json(args.bands)

    W_abs_iv_on, sigma, k0, mode = make_window(window_js)
    bands = parse_bands_generic(bands_js)
    tol = mp.mpf(str(args.tol))

    per_band = []
    glo_lo = mp.inf
    glo_hi = mp.inf

    iterator = enumerate(bands)
    if args.tqdm and _tqdm is not None:
        iterator = _tqdm(
            iterator,
            total=len(bands),
            desc="[band_cert] bands",
            leave=False,
        )

    for idx, (L, R, label) in iterator:
        lo, hi = band_min_bounds(
            W_abs_iv_on, L, R, max_parts=int(args.max_parts), tol=tol
        )
        per_band.append(
            {
                "label": label,
                "left": mp_str(L),
                "right": mp_str(R),
                "min_abs_lo": mp_str(lo),
                "min_abs_hi": mp_str(hi),
            }
        )
        if lo < glo_lo:
            glo_lo = lo
        if hi < glo_hi:
            glo_hi = hi

    PASS_bool = bool(glo_lo > 0)
    PASS_str = "PASS" if PASS_bool else "FAIL"

    # Optional critical band mirror (by label)
    critical = next(
        (b for b in per_band if b.get("label") == "critical"),
        None,
    )
    if critical is not None:
        critical_band = {
            "left": critical["left"],
            "right": critical["right"],
        }
    else:
        critical_band = None

    payload = {
        "kind": "band_cert",
        "created_utc": now_utc_iso(),
        "inputs": {
            "window_config_path": args.window_config,
            "bands_path": args.bands,
            "mode": mode,
            "sigma": mp_str(sigma),
            "k0": mp_str(k0),
            "dps": int(args.dps),
        },
        "numbers": {
            "band_margin_lo": mp_str(glo_lo),
            "band_margin_hi": mp_str(glo_hi),
            "bands_count": len(per_band),
        },
        "band_cert": {
            "band_margin": {
                "lo": mp_str(glo_lo),
                "hi": mp_str(glo_hi),
            },
            "per_band": per_band,
            "status": PASS_str,
            "critical_band": critical_band,
        },
        "PASS": PASS_bool,
        "meta": {
            "tool": "band_cert",
            "dps": int(args.dps),
        },
    }

    write_json(payload, args.out)
    print(
        f"[info] band_margin_lo(cert)={mp_str(glo_lo)}  "
        f"bands={len(per_band)}  status={PASS_str}"
    )


if __name__ == "__main__":
    main()
