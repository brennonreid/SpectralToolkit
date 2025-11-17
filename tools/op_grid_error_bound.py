#!/usr/bin/env python3
# tools/op_grid_error_bound.py
"""
op_grid_error_bound.py

Conservative quadrature/grid error bound for a continuum operator on [a, b],
using a trapezoidal-rule style estimate with an a priori bound on |f''(x)|.
"""

import argparse
import json
import hashlib
import datetime as _dt
import os
from mpmath import mp

try:
    # Preferred: shared numeric helpers
    from lib.io_num import set_precision, mp_str, outward_hi
except Exception:
    # Fallbacks if lib.io_num is not available
    def set_precision(dps: int) -> None:
        mp.dps = int(dps)

    def mp_str(x):
        return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)

    def outward_hi(x):
        try:
            return mp.nextabove(mp.mpf(x))
        except Exception:
            return mp.mpf(x) * (1 + mp.mpf("1e-30"))


def _sha256_bytes(b: bytes) -> str:
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def _write_json(path: str, obj: dict) -> str:
    """Write JSON with Windows-safe newlines and return sha256 of the payload."""
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    blob = json.dumps(obj, indent=2, sort_keys=True).encode("utf-8")
    with open(path, "wb") as f:
        f.write(blob)
    return _sha256_bytes(blob)


def main():
    """
    CLI (normalized v2.1):

      --a             : left endpoint of interval [a, b]
      --b             : right endpoint of interval [a, b]
      --grid-points   : number of grid points (>= 2)
      --M2-hi         : upper bound on |f''(x)| over [a, b]
      --cap-hi        : optional hard cap override for the bound
      --dps           : decimal precision for mpmath
      --out           : output JSON path

    JSON (normalized v2.1):

      {
        "kind": "grid_error_bound",
        "inputs": {
          "a": "<string>",
          "b": "<string>",
          "grid_points": <int>,
          "M2_hi": "<string>",
          "cap_hi": "<string or null>"
        },
        "grid_error_bound": {
          "bound_hi": "<string>"
        },
        "meta": {
          "tool": "op_grid_error_bound",
          "dps": <int>,
          "created_utc": "<ISO8601 UTC>",
          "sha256": "<digest>"
        }
      }
    """
    ap = argparse.ArgumentParser(
        description="Conservative grid/quadrature error bound (trap rule style)"
    )
    ap.add_argument("--a", type=str, default="0")
    ap.add_argument("--b", type=str, default="1")
    ap.add_argument(
        "--grid-points",
        type=int,
        default=6000,
        help="Number of grid points (>= 2).",
    )
    ap.add_argument(
        "--M2-hi",
        type=str,
        default="1e-3",
        help="Upper bound on |f''(x)| over [a,b].",
    )
    ap.add_argument(
        "--cap-hi",
        type=str,
        default=None,
        help="Optional hard cap override for the bound.",
    )
    ap.add_argument("--dps", type=int, default=220)
    ap.add_argument("--out", required=True)

    args = ap.parse_args()
    set_precision(args.dps)

    a = mp.mpf(args.a)
    b = mp.mpf(args.b)
    grid_points = int(args.grid_points)

    if grid_points < 2:
        raise SystemExit("grid_points must be >= 2")

    N = grid_points - 1
    M2_hi = mp.mpf(args.M2_hi)
    width = b - a

    # Trap-rule style error bound:
    #   |error| <= (width * M2_hi) / (12 * N^2)
    bound = (width * M2_hi) / (12 * (N ** 2))
    bound = outward_hi(bound)

    if args.cap_hi is not None:
        cap_hi = mp.mpf(args.cap_hi)
        if cap_hi < bound:
            bound = cap_hi

    inputs = {
        "a": mp_str(a),
        "b": mp_str(b),
        "grid_points": grid_points,
        "M2_hi": mp_str(M2_hi),
        "cap_hi": mp_str(cap_hi) if args.cap_hi is not None else None,
    }

    grid_error_block = {
        "bound_hi": mp_str(bound),
    }

    payload = {
        "kind": "grid_error_bound",
        "inputs": inputs,
        "grid_error_bound": grid_error_block,
        "meta": {
            "tool": "op_grid_error_bound",
            "dps": int(mp.dps),
            "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        },
    }

    sha = _write_json(args.out, payload)
    payload["meta"]["sha256"] = sha
    _write_json(args.out, payload)

    print(
        f"[ok] grid_error_bound -> {args.out}  "
        f"bound_hi={mp.nstr(bound, 10)}"
    )


if __name__ == "__main__":
    main()

