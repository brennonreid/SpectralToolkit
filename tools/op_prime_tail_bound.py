#!/usr/bin/env python3
"""
op_prime_tail_bound.py — analytic prime tail bound beyond x0.

Model (conservative, math unchanged):
    C_tail_hi = 2 * (K + 1) * A_prime
    and optionally C_tail_hi /= log(x0) if --scale-by-log is set.

JSON (v2.1 normalized):
  kind = "prime_tail_bound"
  inputs {
    x0,
    A_prime,
    K,
    scale_by_log
  }
  prime_tail_bound {
    C_tail_hi,
    model
  }
  meta {
    tool        = "op_prime_tail_bound",
    dps,
    created_utc,
    sha256
  }
"""

import argparse
import json
import hashlib
import datetime as _dt
import os

from mpmath import mp


# ---------------------------------------------------------------------
# Precision / JSON helpers (interface-level only)
# ---------------------------------------------------------------------

def set_precision(dps: int) -> None:
    mp.dps = int(dps)


def mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def outward_hi(x):
    """Outward rounding upward for a conservative upper bound."""
    try:
        return mp.nextabove(mp.mpf(x))
    except Exception:
        return mp.mpf(x) * (1 + mp.mpf("1e-30"))


def now_utc_iso() -> str:
    return (
        _dt.datetime.now(_dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def write_json(path: str, payload: dict, dps: int) -> None:
    """
    Write JSON with canonical meta block:

      meta.tool        = "op_prime_tail_bound"
      meta.dps         = dps
      meta.created_utc = ISO8601 UTC
      meta.sha256      = hash over payload without meta.sha256

    Hash is computed using sorted keys and compact separators.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    if "meta" not in payload or not isinstance(payload["meta"], dict):
        payload["meta"] = {}
    meta = payload["meta"]
    meta["tool"] = "op_prime_tail_bound"
    meta["dps"] = int(dps)
    meta["created_utc"] = now_utc_iso()

    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    meta2 = tmp_obj.get("meta")
    if isinstance(meta2, dict):
        meta2.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    sha = hashlib.sha256(blob).hexdigest()
    meta["sha256"] = sha

    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------------
# Main (math unchanged, interface normalized)
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Analytic prime tail bound constant.")
    ap.add_argument(
        "--x0",
        type=str,
        default="1e6",
        help="Cutover x0 for explicit vs tail (string, converted via mpmath).",
    )
    ap.add_argument(
        "--A-prime",
        dest="A_prime",
        type=str,
        default="1.2762",
        help="A_prime (prime-sum coefficient).",
    )
    ap.add_argument(
        "--K",
        type=int,
        default=3,
        help="Auxiliary index K used in the tail model.",
    )
    ap.add_argument(
        "--scale-by-log",
        action="store_true",
        help="If set, scales C_tail_hi by 1/log(x0) (conservative).",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path.",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=220,
        help="Decimal precision for mpmath.",
    )
    args = ap.parse_args()

    set_precision(args.dps)

    x0 = mp.mpf(args.x0)
    A = mp.mpf(args.A_prime)
    K = int(args.K)

    # Math model (unchanged):
    #   C_tail_hi = 2 * (K + 1) * A_prime [* 1/log(x0) if scale_by_log]
    C = 2 * (K + 1) * A
    if args.scale_by_log:
        C = C / mp.log(x0)
    C = outward_hi(C)

    payload = {
        "kind": "prime_tail_bound",
        "inputs": {
            "x0": mp_str(x0),
            "A_prime": mp_str(A),
            "K": K,
            "scale_by_log": bool(args.scale_by_log),
        },
        "prime_tail_bound": {
            "C_tail_hi": mp_str(C),
            "model": "C_tail_hi = 2*(K+1)*A_prime [* 1/log(x0) if scale_by_log].",
        },
        "meta": {},
    }

    write_json(args.out, payload, args.dps)
    print(f"[ok] prime_tail_bound -> {args.out}  C_tail_hi={mp_str(C)}")


if __name__ == "__main__":
    main()
