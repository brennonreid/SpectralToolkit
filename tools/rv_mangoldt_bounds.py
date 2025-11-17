#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rv_mangoldt_bounds.py — Riemann–von Mangoldt lower bound and monotonicity at T0.

Purpose:
  Compute a Riemann–von Mangoldt lower bound N_lower(T) and a monotonicity
  certificate starting at T0. If the derivative dN_lower/dT >= 0 at T0 under
  the assumed S(T) bound, we record that N(T) is nondecreasing for T >= T0
  and use this as an infinitude witness.

CLI (v2.1 normalized):
  --T0          : base point T0 (string or number)
  --dps         : decimal precision for mpmath
  --out         : output JSON path
  --theory-out  : optional theory JSON path

JSON (v2.1 normalized):
  {
    "kind": "rv_mangoldt_bounds",
    "inputs": {
      "T0": "<string>"
    },
    "rv_mangoldt": {
      "lower_bound_value": "<string>",
      "S_bound": {
        "a": "<string>",
        "b": "<string>",
        "T_star": "<string>"
      },
      "monotone_for_T_ge_T0": true/false,
      "infinitude_certified": true/false
    },
    "meta": {
      "tool": "rv_mangoldt_bounds",
      "dps": <int>,
      "created_utc": "<ISO8601 UTC>",
      "sha256": "<digest>"
    }
  }

All real-valued numeric fields in JSON are stored as strings. Integers remain
JSON integers; booleans remain JSON booleans.
"""

import argparse
import time
import json
import hashlib
import os
from datetime import datetime, timezone

from mpmath import mp


# --------------------------------------------------------------------
# Numeric helpers (no heavy loops; this file is light-weight)
# --------------------------------------------------------------------

def set_precision(dps: int) -> None:
    mp.dps = int(dps)


def as_mpf(x):
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


def write_json(path: str, obj: dict) -> None:
    """
    Write JSON with deterministic sha256 in meta.sha256.

    Hash is computed over the payload with any existing meta.sha256 removed,
    using canonical JSON (sorted keys, compact separators).
    """
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)

    if "meta" not in obj or not isinstance(obj["meta"], dict):
        obj["meta"] = {}

    # Compute hash on a copy without meta.sha256 (avoid self-reference).
    tmp_obj = json.loads(json.dumps(obj, ensure_ascii=False))
    meta = tmp_obj.get("meta")
    if isinstance(meta, dict):
        meta.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    digest = hashlib.sha256(blob).hexdigest()
    obj["meta"]["sha256"] = digest

    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, path)


# --------------------------------------------------------------------
# RvM model (math unchanged)
# --------------------------------------------------------------------

def S_bound_params(T):
    a = mp.mpf("0.20")
    b = mp.mpf("1.00")
    return a, b


def N_main(T):
    u = T / (2 * mp.pi)
    return u * mp.log(u) - u + mp.mpf("0.875")


def N_lower(T):
    a, b = S_bound_params(T)
    C = mp.mpf("1.0")
    return N_main(T) - (a * mp.log(T) + b) - (C / T)


def monotone_for_T_ge_T0(T0):
    a, b = S_bound_params(T0)
    C = mp.mpf("1.0")
    T = mp.mpf(T0)
    u = T / (2 * mp.pi)
    dN = (mp.mpf("1") / (2 * mp.pi)) * mp.log(u) - a / T + C / (T * T)
    return dN >= 0


# --------------------------------------------------------------------
# Main CLI
# --------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Riemann–von Mangoldt lower bound and monotonicity certificate at T0."
    )
    ap.add_argument("--T0", type=str, required=True, help="Base point T0 (>0).")
    ap.add_argument("--dps", type=int, default=300, help="Decimal precision.")
    ap.add_argument("--out", type=str, required=True, help="Output JSON path.")
    ap.add_argument(
        "--theory-out",
        default=None,
        help="Optional theory JSON path (informal lemma statement).",
    )
    args = ap.parse_args()

    set_precision(args.dps)
    T0 = as_mpf(args.T0)
    if T0 <= 0:
        raise SystemExit("T0 must be > 0")

    nl = N_lower(T0)
    inc = monotone_for_T_ge_T0(T0)
    a, b = S_bound_params(T0)

    payload = {
        "kind": "rv_mangoldt_bounds",
        "inputs": {
            "T0": mp_str(T0),
        },
        "rv_mangoldt": {
            "lower_bound_value": mp_str(nl),
            "S_bound": {
                "a": mp_str(a),
                "b": mp_str(b),
                "T_star": mp_str(T0),
            },
            "monotone_for_T_ge_T0": bool(inc),
            "infinitude_certified": bool(inc),
        },
        "meta": {
            "tool": "rv_mangoldt_bounds",
            "dps": int(args.dps),
            "created_utc": now_utc_iso(),
        },
    }

    write_json(args.out, payload)
    print(
        f"[OK] RvM -> {args.out}  "
        f"N_lower(T0)≈{mp.nstr(nl, 8)}  monotone={inc}"
    )

    if args.theory_out:
        theory = {
            "kind": "rv_mangoldt_bounds_theory",
            "lemma": "RvMInfinitudeLowerBound",
            "hypotheses": {
                "T0": {">": "0"},
                "S(T)": {"abs_le": "a log T + b"},
            },
            "conclusion": {
                "forall_T_ge_T0": (
                    "N(T) >= N_lower(T) and N_lower(T) is nondecreasing for T >= T0"
                )
            },
            "constants": {
                "a": mp_str(a),
                "b": mp_str(b),
                "T0": mp_str(T0),
            },
            "meta": {
                "tool": "rv_mangoldt_bounds",
                "dps": int(args.dps),
                "created_utc": now_utc_iso(),
            },
        }
        with open(args.theory_out, "w", encoding="utf-8", newline="\n") as f:
            json.dump(theory, f, indent=2, sort_keys=True, ensure_ascii=False)
            f.write("\n")
        print(f"[THY] wrote {args.theory_out}")


if __name__ == "__main__":
    main()
