#!/usr/bin/env python3
"""
Generate a Gaussian (notched) window configuration JSON for the kernel.

CLI (normalized v2.1):
  --mode    : window type (currently "gauss")
  --sigma   : Gaussian width parameter (>0), real, stored as string in JSON
  --k0      : notch parameter k0 (>0), real, stored as string in JSON
  --dps     : decimal precision for mpmath
  --out     : output JSON path

JSON (normalized v2.1):
  kind  = "window"
  mode
  sigma   (string)
  k0      (string)
  window {
    mode
    sigma
    k0
  }
  meta {
    tool        = "window_gen"
    dps         = <int>
    created_utc = "<ISO8601 UTC>"
    sha256      = "<digest>"   # added after serialization
  }

All real-valued numeric fields are emitted as strings to preserve precision.
"""

import argparse
import json
import hashlib
from mpmath import mp
from datetime import datetime, timezone
from pathlib import Path


def now_utc_iso() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def write_json(path: Path, payload: dict) -> None:
    """
    Write JSON to path (Windows-safe newlines), computing sha256 over the
    canonical serialized form and storing it in meta.sha256.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    # Serialize once to compute sha256.
    js = json.dumps(payload, sort_keys=True, ensure_ascii=False, indent=2)
    sha = hashlib.sha256(js.encode("utf-8")).hexdigest()
    payload["meta"]["sha256"] = sha

    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, sort_keys=True, ensure_ascii=False, indent=2)
        f.write("\n")
    Path(tmp).replace(path)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Generate Gaussian (notched) window configuration JSON."
    )
    ap.add_argument(
        "--mode",
        choices=["gauss"],
        required=True,
        help="Window type (currently only 'gauss').",
    )
    ap.add_argument(
        "--sigma",
        type=str,
        required=True,
        help="Gaussian width parameter sigma (>0).",
    )
    ap.add_argument(
        "--k0",
        type=str,
        required=True,
        help="Notch parameter k0 (>0).",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=220,
        help="Decimal precision for mpmath / printing.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path.",
    )
    args = ap.parse_args()

    mp.dps = int(args.dps)

    sigma = mp.mpf(args.sigma)
    k0 = mp.mpf(args.k0)
    if sigma <= 0 or k0 <= 0:
        raise ValueError("sigma and k0 must be positive")

    sigma_s = mp_str(sigma)
    k0_s = mp_str(k0)

    payload = {
        "kind": "window",
        "mode": args.mode,
        "sigma": sigma_s,
        "k0": k0_s,
        "window": {
            "mode": args.mode,
            "sigma": sigma_s,
            "k0": k0_s,
        },
        "meta": {
            "tool": "window_gen",
            "dps": int(args.dps),
            "created_utc": now_utc_iso(),
        },
    }

    out_path = Path(args.out)
    write_json(out_path, payload)
    print(f"[window_gen] -> {out_path}  mode={args.mode} sigma={sigma_s} k0={k0_s}")
    print(f"[info] sha256={payload['meta']['sha256']}")


if __name__ == "__main__":
    main()
