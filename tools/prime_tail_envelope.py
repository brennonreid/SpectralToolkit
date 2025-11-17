#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prime_tail_envelope.py — prime tail envelope at T0 using A_prime and K.

Purpose:
  Provide a conservative analytic envelope for the prime tail at a cutoff T0.
  The envelope is modeled as a decaying function of (sigma, k0, T0, A_prime, K)
  and used as a scalar prime_tail_norm in downstream rollups.

CLI (v2.1 normalized):
  --T0          : cutoff T0 (>0)
  --sigma       : Gaussian width sigma (>0)
  --k0          : notch parameter k0
  --A-prime     : A_prime parameter for the prime tail model
  --K           : integer parameter shaping the envelope
  --dps         : decimal precision for mpmath
  --out         : output JSON path
  --theory-out  : optional theory JSON path
  --env-T0-hi   : optional explicit override for env_T0_hi

JSON (v2.1 normalized):
  {
    "kind": "prime_tail_envelope",
    "inputs": {
      "T0": "<string>",
      "sigma": "<string>",
      "k0": "<string>",
      "A_prime": "<string>",
      "K": <int>
    },
    "prime_tail": {
      "env_T0_hi": "<string>",
      "norm": "<string>"
    },
    "meta": {
      "tool": "prime_tail_envelope",
      "dps": <int>,
      "created_utc": "<ISO8601 UTC>",
      "sha256": "<digest>"
    }
  }

All real-valued numeric fields are stored as strings.
Integers remain JSON integers; booleans remain JSON booleans.
"""

import argparse
import json
import os
import hashlib
from datetime import datetime, timezone

from mpmath import mp


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def set_prec(dps: int) -> None:
    mp.dps = int(dps)


def mpstr(x) -> str:
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


# ---------------------------------------------------------------------
# Prime tail envelope model (math unchanged)
# ---------------------------------------------------------------------

def conservative_prime_tail_env_T0(sigma, k0, T0, A_prime, K):
    """
    Conservative envelope for the prime tail at T0.

    Math is preserved from the original implementation:
      x    = sigma * k0 * T0
      base = exp(-(x^2)/4)
      env  = A_prime * base / (1 + K)
    """
    x = mp.mpf(sigma) * mp.mpf(k0) * mp.mpf(T0)
    base = mp.e ** (-(x ** 2) / 4)
    env = mp.mpf(A_prime) * base / (1 + mp.mpf(K))
    return env


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Prime tail envelope at T0 using A_prime and K."
    )
    ap.add_argument(
        "--T0",
        type=str,
        required=True,
        help="Cutoff T0 (>0).",
    )
    ap.add_argument(
        "--sigma",
        type=str,
        required=True,
        help="Gaussian width sigma (>0).",
    )
    ap.add_argument(
        "--k0",
        type=str,
        required=True,
        help="Notch parameter k0.",
    )
    ap.add_argument(
        "--A-prime",
        dest="A_prime",
        type=str,
        default=(
            "0.0076384524109054769957964191958869400469835723173758547544266337342484501729838427610"
            "5585333349922026737375344819040959339805474951892698656605316742311018265496226235760840"
            "5931698972986515621555230774697452540257558441032535745760914086692541001577311615455938"
            "56661222382687673426384518037616762448945"
        ),
        help="A_prime parameter for the prime tail model.",
    )
    ap.add_argument(
        "--K",
        type=int,
        default=3,
        help="Auxiliary integer parameter shaping the envelope.",
    )
    ap.add_argument(
        "--env-T0-hi",
        dest="env_T0_hi",
        type=str,
        default=None,
        help="Explicit value for env_T0_hi (overrides derived envelope).",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=200,
        help="Decimal precision for mpmath.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path (prime_tail_envelope).",
    )
    ap.add_argument(
        "--theory-out",
        default=None,
        help="Optional sidecar theory JSON path.",
    )
    args = ap.parse_args()

    set_prec(args.dps)

    T0 = mp.mpf(args.T0)
    sigma = mp.mpf(args.sigma)
    k0 = mp.mpf(args.k0)
    A_prime = mp.mpf(args.A_prime)
    K = int(args.K)

    if T0 <= 0:
        raise SystemExit("T0 must be > 0")
    if sigma <= 0:
        raise SystemExit("sigma must be > 0")

    if args.env_T0_hi is not None:
        env_T0 = mp.mpf(args.env_T0_hi)
    else:
        env_T0 = conservative_prime_tail_env_T0(sigma, k0, T0, A_prime, K)

    payload = {
        "kind": "prime_tail_envelope",
        "inputs": {
            "T0": mpstr(T0),
            "sigma": mpstr(sigma),
            "k0": mpstr(k0),
            "A_prime": mpstr(A_prime),
            "K": K,
        },
        "prime_tail": {
            "env_T0_hi": mpstr(env_T0),
            "norm": mpstr(env_T0),
        },
        "meta": {
            "tool": "prime_tail_envelope",
            "dps": int(mp.dps),
            "created_utc": now_utc_iso(),
        },
    }

    write_json(args.out, payload)
    print(
        f"[ok] prime_tail_envelope -> {args.out}  "
        f"T0={mpstr(T0)}  env_T0_hi={mpstr(env_T0)}"
    )

    if args.theory_out:
        theory = {
            "kind": "prime_tail_envelope_theory",
            "inputs": {
                "T0": mpstr(T0),
                "sigma": mpstr(sigma),
                "k0": mpstr(k0),
                "A_prime": mpstr(A_prime),
                "K": K,
            },
            "theory": {
                "lemma": "PrimeTailEnvelope",
                "statement": (
                    "For T >= T0, the prime tail contribution is bounded above by "
                    "an envelope depending on (sigma, k0, A_prime, K). This module "
                    "records the value at T0 as env_T0_hi and uses it as a scalar "
                    "prime_tail_norm in downstream rollups."
                ),
            },
            "meta": {
                "tool": "prime_tail_envelope",
                "dps": int(mp.dps),
                "created_utc": now_utc_iso(),
            },
        }
        write_json(args.theory_out, theory)
        print(f"[ok] prime_tail_envelope theory -> {args.theory_out}")


if __name__ == "__main__":
    main()
