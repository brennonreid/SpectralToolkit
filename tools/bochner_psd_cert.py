#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bochner_psd_cert.py — Bochner PSD certificate for Gaussian-notch window (v2.1 normalized)

Purpose:
  Certify that the Gaussian-notch window has a nonnegative spectral density
  h_hat(t), so the associated kernel is positive-definite by Bochner's theorem.

  Frequency-domain form:
    h_hat(t) = exp(-(t/sigma)^2) * (1 - exp(-(t/k0)^2))

  For all real t, exp(-(t/sigma)^2) >= 0 and 1 - exp(-(t/k0)^2) >= 0 when
  sigma > 0 and k0 > 0, so h_hat(t) >= 0 on R. This is the analytic PSD
  certificate; a numeric sweep is included only as telemetry.

CLI (v2.1 normalized):
  --window         : path to window.json (Gaussian-notch config)
  --out            : output JSON path (e.g. PROOF_PACKET/weil_psd_bochner.json)
  --dps            : decimal precision for mpmath (default 200)
  --sweep-T        : half-width of numeric sweep in t (default 200.0)
  --sweep-steps    : number of sweep sample points (default 40001)

JSON (v2.1 normalized):
  kind = "weil_psd_bochner"

  inputs {
    window_path,
    dps,
    sweep_T,
    sweep_steps
  }

  weil_psd_bochner {
    PSD_verified,       # analytic PSD flag (primary certificate)
    PASS,               # alias for PSD_verified for backward compatibility
    reason,             # short justification string
    certificate {
      type,             # "Bochner"
      mode,             # window mode (e.g. "gauss")
      sigma,            # Gaussian width (stringified mp.mpf)
      notch_k0          # notch location k0 (stringified mp.mpf)
    }
    eval {
      min_hat_h_sample, # numeric min of h_hat(t) on sweep
      t_at_min,         # t where min was observed
      sweep_T,
      sweep_steps
    }
    window_echo         # raw window.json echo for traceability
  }

  meta {
    tool,
    dps,
    created_utc,
    sha256
  }

Notes:
  - Analytic PSD certificate is the only thing needed for the proof budget:
    PSD_verified is True iff sigma > 0 and k0 > 0.
  - Numeric sweep on [-T, T] is for telemetry / sanity checking only.
"""

import os
import json
import argparse
from mpmath import mp
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def now_utc_iso() -> str:
    """Return current UTC time as an ISO-8601 string with trailing 'Z'."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def set_precision(dps: int) -> None:
    """Set global mpmath precision."""
    mp.dps = int(dps)


def load_json(path: str):
    """Load JSON from disk."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def coalesce(d, *keys, default=None):
    """
    Search (possibly nested) dict d for the first present key path.
    Each key can be a string 'a' or a tuple path like ('window','sigma').
    """
    for k in keys:
        if isinstance(k, tuple):
            cur = d
            ok = True
            for part in k:
                if isinstance(cur, dict) and part in cur:
                    cur = cur[part]
                else:
                    ok = False
                    break
            if ok:
                return cur
        else:
            if isinstance(d, dict) and k in d:
                return d[k]
    return default


def _sha256_bytes(b: bytes) -> str:
    """Compute SHA-256 hex digest of a bytes object."""
    import hashlib
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def _write_json(path: str, obj: dict) -> str:
    """
    Write JSON to path (UTF-8, pretty-printed) and return SHA-256 digest
    of the on-disk bytes.
    """
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    s = json.dumps(obj, indent=2, ensure_ascii=False).encode("utf-8")
    with open(path, "wb") as f:
        f.write(s)
    return _sha256_bytes(s)


# ---------------------------------------------------------------------------
# Window parsing and PSD logic
# ---------------------------------------------------------------------------

def parse_window(path: str):
    """
    Parse Gaussian-notch window JSON from multiple possible schemas, e.g.:

      { "mode":"gauss","sigma":5.6,"notch_k0":0.25 }
      { "window": {"mode":"gauss","sigma":5.6,"notch_k0":0.25} }
      { "mode":"gauss","gauss_sigma":5.6,"notch":{"k0":0.25} }
      { "window":{"mode":"gauss","params":{"sigma":5.6,"notch_k0":0.25}} }

    Returns (mode, sigma, k0, raw_json).
    """
    js = load_json(path)

    mode = coalesce(
        js,
        "mode",
        ("window", "mode"),
        ("window", "type"),
        default="gauss",
    )

    sigma = coalesce(
        js,
        "sigma",
        "gauss_sigma",
        ("window", "sigma"),
        ("window", "gauss_sigma"),
        ("window", "params", "sigma"),
    )

    k0 = coalesce(
        js,
        "notch_k0",
        "k0",
        ("window", "notch_k0"),
        ("window", "params", "notch_k0"),
        ("notch", "k0"),
        ("window", "notch", "k0"),
    )

    if sigma is None:
        raise KeyError(
            "window.json missing sigma "
            "(tried: 'sigma', 'gauss_sigma', window.sigma, window.params.sigma)"
        )
    if k0 is None:
        raise KeyError(
            "window.json missing notch k0 "
            "(tried: 'notch_k0', 'k0', window.notch_k0, window.params.notch_k0)"
        )

    sigma = mp.mpf(str(sigma))
    k0 = mp.mpf(str(k0))
    return mode, sigma, k0, js


def hhat(t, sigma, k0):
    """
    Frequency-domain kernel:

      h_hat(t) = exp(-(t/sigma)^2) * (1 - exp(-(t/k0)^2))
    """
    return mp.e ** (-(t / sigma) ** 2) * (1 - mp.e ** (-(t / k0) ** 2))


def analytic_psd_holds(sigma, k0) -> bool:
    """
    Analytic Bochner certificate:

      - exp(-(t/sigma)^2) >= 0 for all real t
      - 1 - exp(-(t/k0)^2) >= 0 for all real t when k0 > 0

    With sigma > 0 and k0 > 0, h_hat(t) >= 0 on R, so the kernel is PSD.
    """
    return (sigma > 0) and (k0 > 0)


def numeric_sweep(sigma, k0, T=100.0, steps=20001):
    """
    Numerically sample h_hat(t) on [-T, T] to find a minimum value.
    This is telemetry only and not used in the formal inequality budget.
    """
    T = mp.mpf(T)
    steps = int(steps)
    if steps < 2:
        steps = 2

    mn = mp.inf
    argmin = mp.ninf

    for i in range(steps):
        t = -T + (2 * T) * mp.mpf(i) / (steps - 1)
        val = hhat(t, sigma, k0)
        if val < mn:
            mn = val
            argmin = t

    return mn, argmin


# ---------------------------------------------------------------------------
# Main CLI entry point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Bochner PSD certificate for Gaussian-notch window"
    )
    ap.add_argument("--window", required=True, help="path to window.json")
    ap.add_argument("--out", required=True, help="output JSON certificate path")
    ap.add_argument(
        "--dps", type=int, default=200, help="decimal precision for mpmath"
    )
    ap.add_argument(
        "--sweep-T",
        dest="sweep_T",
        type=float,
        default=200.0,
        help="numeric sweep half-width",
    )
    ap.add_argument(
        "--sweep-steps",
        type=int,
        default=40001,
        help="number of sample points in numeric sweep",
    )
    args = ap.parse_args()

    set_precision(args.dps)

    mode, sigma, k0, wraw = parse_window(args.window)

    # Analytic Bochner PSD check
    passed = analytic_psd_holds(sigma, k0)

    # Numeric telemetry sweep
    min_val, t_at_min = numeric_sweep(
        sigma, k0, T=args.sweep_T, steps=args.sweep_steps
    )

    payload = {
        "kind": "weil_psd_bochner",
        "inputs": {
            "window_path": args.window,
            "dps": str(mp.dps),
            "sweep_T": float(args.sweep_T),
            "sweep_steps": int(args.sweep_steps),
        },
        "weil_psd_bochner": {
            # STP expects PSD_verified; PASS is kept as an alias.
            "PSD_verified": bool(passed),
            "PASS": bool(passed),
            "reason": (
                "nonnegative_spectral_density: "
                "exp(-(t/sigma)^2) * (1 - exp(-(t/k0)^2)) >= 0 for all real t; "
                "Bochner implies PSD"
            ),
            "certificate": {
                "type": "Bochner",
                "mode": str(mode),
                "sigma": mp.nstr(sigma, n=mp.dps),
                "notch_k0": mp.nstr(k0, n=mp.dps),
            },
            "eval": {
                "min_hat_h_sample": mp.nstr(min_val, n=mp.dps),
                "t_at_min": mp.nstr(t_at_min, n=mp.dps),
                "sweep_T": float(args.sweep_T),
                "sweep_steps": int(args.sweep_steps),
            },
            "window_echo": wraw,
        },
        "meta": {
            "tool": "bochner_psd_cert",
            "dps": str(mp.dps),
            "created_utc": now_utc_iso(),
            "sha256": "",
        },
    }

    # First write to compute SHA-256, then embed it and rewrite.
    sha = _write_json(args.out, payload)
    payload["meta"]["sha256"] = sha
    sha = _write_json(args.out, payload)

    status = "PASS=True" if passed else "PASS=False"
    print(f"[bochner_psd_cert] {status} wrote {args.out} (sha256={sha})")


if __name__ == "__main__":
    main()
