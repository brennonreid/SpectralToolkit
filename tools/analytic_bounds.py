#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analytic_bounds.py — Analytic helper bounds (v2.1 normalized)

Purpose
-------
Emit:
  - gamma-tail constants c1, c2 as sigma-dependent envelopes
  - prime-tail constant C(K, A') = 2*(K+1)*A'
  - effective epsilon lower bound eps_eff_lo (alias of C)
  - grid error upper bound grid_error_hi (from grid_error_bound.json)

CLI (v2.1)
----------
  --sigma       : Gaussian width parameter (string, high-precision)
  --A-prime     : prime-tail constant A' (string)
  --K           : integer K in C(K, A')
  --dps         : decimal precision for mpmath
  --out         : output JSON path

JSON (v2.1)
-----------
  kind = "analytic_bounds"
  inputs {
    sigma,
    A_prime,
    K,
    dps
  }
  analytic_bounds {
    gamma {
      c1,
      c2,
      C,
      a
    }
    gamma_tail {         # alias of gamma
      c1,
      c2,
      C,
      a
    }
    prime {
      C,
      K,
      A_prime,
      a
    }
    prime_tail {         # alias of prime
      C,
      K,
      A_prime,
      a
    }
    eps_eff_lo,
    grid_error_hi
  }
  gamma_tail     (top-level alias)
  prime_tail     (top-level alias)
  eps_eff_lo     (top-level alias)
  grid_error_hi  (top-level alias)
  meta {
    tool,
    dps,
    created_utc,
    sha256
  }
"""

import argparse
import json
import hashlib
import os
from datetime import datetime, timezone

from mpmath import mp


# ---------- helpers ----------

def set_precision(dps: int) -> None:
    mp.dps = int(dps)


def mp_str(x):
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def sha256_bytes(b: bytes) -> str:
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def write_json_with_sha(path: str, payload: dict) -> str:
    raw = json.dumps(payload, indent=2, ensure_ascii=False).encode("utf-8")
    with open(path, "wb") as f:
        f.write(raw)
    return sha256_bytes(raw)


# ---------- analytic bounds (math unchanged) ----------

def gamma_tail_constants(sigma):
    """
    Simple monotone decreasing envelopes ~ 1 + 1/sigma, 1 + 1/sigma^2.
    """
    sigma = mp.mpf(sigma)
    return (mp.mpf("1") + 1 / sigma, mp.mpf("1") + 1 / (sigma * sigma))


def prime_tail_constant(A_prime, K=3):
    """
    C(K, A') = 2 * (K + 1) * A'
    """
    return mp.mpf("2") * (K + 1) * mp.mpf(A_prime)


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Emit analytic helper bounds.")
    ap.add_argument("--sigma", type=str, default="6.0")
    ap.add_argument("--A-prime", dest="A_prime", type=str, default="1.2762")
    ap.add_argument("--K", type=int, default=3)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dps", type=int, default=220)
    args = ap.parse_args()

    set_precision(args.dps)

    # Core analytic constants
    c1, c2 = gamma_tail_constants(args.sigma)
    C = prime_tail_constant(args.A_prime, args.K)

    # Pull grid_error_hi from grid_error_bound.json in the same directory as --out
    out_dir = os.path.dirname(os.path.abspath(args.out))
    geb_path = os.path.join(out_dir, "grid_error_bound.json")

    try:
        with open(geb_path, "r", encoding="utf-8") as f:
            geb = json.load(f)
        # Expected schema from op_grid_error_bound.py:
        # { "grid_error_bound": { "bound_hi": "..." , ... }, ... }
        grid_error_hi = mp.mpf(geb["grid_error_bound"]["bound_hi"])
    except FileNotFoundError:
        raise SystemExit(
            f"[analytic_bounds] ERROR: grid_error_bound.json not found at {geb_path}"
        )
    except KeyError as e:
        raise SystemExit(
            f"[analytic_bounds] ERROR: malformed grid_error_bound.json, missing key {e!r}"
        )

    # a=1.0 is the default “amplitude” used downstream (matches analytic_tail_fit output)
    a_val = mp.mpf("1")

    gamma_block = {
        "c1": mp_str(c1),
        "c2": mp_str(c2),
        # Treat c1 as the effective C used for gamma tails (matches the 1 + 1/sigma shape)
        "C": mp_str(c1),
        "a": mp_str(a_val),
    }

    prime_block = {
        "C": mp_str(C),
        "K": int(args.K),
        "A_prime": mp_str(args.A_prime),
        "a": mp_str(a_val),
    }

    payload = {
        "kind": "analytic_bounds",
        "inputs": {
            "sigma": args.sigma,
            "A_prime": args.A_prime,
            "K": int(args.K),
            "dps": str(args.dps),
        },
        "analytic_bounds": {
            # Primary gamma block
            "gamma": gamma_block,
            # Alias block for tools expecting "gamma_tail" under analytic_bounds
            "gamma_tail": gamma_block,
            # Primary prime-tail block
            "prime": prime_block,
            # Alias block for tools expecting "prime_tail" under analytic_bounds
            "prime_tail": prime_block,
            # Effective epsilon lower bound (alias of C)
            "eps_eff_lo": mp_str(C),
            # Grid error upper bound re-exported from op_grid_error_bound.py
            "grid_error_hi": mp_str(grid_error_hi),
        },
        # Top-level aliases for tools that read these directly
        "gamma_tail": gamma_block,
        "prime_tail": prime_block,
        "eps_eff_lo": mp_str(C),
        "grid_error_hi": mp_str(grid_error_hi),
        "meta": {
            "tool": "analytic_bounds",
            "dps": str(args.dps),
            "created_utc": datetime.now(timezone.utc)
            .isoformat()
            .replace("+00:00", "Z"),
            "sha256": "",
        },
    }

    sha = write_json_with_sha(args.out, payload)
    payload["meta"]["sha256"] = sha
    write_json_with_sha(args.out, payload)

    print(f"[ok] analytic_bounds -> {args.out}")
    print(
        f"[analytic_bounds] sigma={args.sigma}  A_prime={args.A_prime}  "
        f"K={args.K}  C={mp_str(C)}  grid_error_hi={mp_str(grid_error_hi)}"
    )


if __name__ == "__main__":
    main()
