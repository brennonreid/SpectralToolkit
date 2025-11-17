#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
continuum_operator_rollup.py — continuum operator inequality (v2.1 normalized)

Purpose
-------
Combine:
  - band margin (band_cert)
  - prime block cap (prime_block_norm)
  - prime tail norm (prime_tail_envelope)
  - gamma tail envelope at T0 (gamma_tails)
  - grid error bound (op_grid_error_bound)

into a single inequality:

  lhs_total = prime_block_cap + prime_tail_norm + grid_error_norm
  epsilon_eff = band_margin - gamma_tails

and certify PASS iff:

  lhs_total <= epsilon_eff.

CLI (v2.1)
----------
  --band-cert     : path to band_cert.json
  --prime-block   : path to prime_block_norm.json
  --prime-tail    : path to prime_tail_envelope.json
  --gamma-tails   : path to gamma_tails.json (core_integral_prover)
  --grid-error    : path to grid_error_bound.json (op_grid_error_bound)
  --dps           : decimal precision
  --out           : output JSON path

JSON (v2.1)
-----------
  kind = "continuum_operator_cert"
  inputs {
    band_cert_path,
    prime_block_path,
    prime_tail_path,
    gamma_tails_path,
    grid_error_path
  }
  numbers {
    band_margin,
    prime_block_cap,
    prime_tail_norm,
    gamma_tails,
    grid_error_norm,
    lhs_total,
    epsilon_eff
  }
  PASS  (boolean, lhs_total <= epsilon_eff)
  meta {
    tool        = "continuum_operator_rollup",
    dps,
    created_utc
  }
"""

import argparse
import json
import os
import time
from mpmath import mp


# ---------- basic helpers ----------

def set_prec(dps):
    mp.dps = int(dps)


def mpf(x):
    return mp.mpf(str(x))


def mpstr(x):
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def jload(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def dig(obj, path):
    """Safe nested lookup; returns None if any key is missing."""
    o = obj
    for k in path:
        if not isinstance(o, dict) or k not in o:
            return None
        o = o[k]
    return o


def mpf_from_paths(js, paths, default=None):
    """
    Try multiple key paths until one succeeds; return mp.mpf.
    If all fail and default is None, raise KeyError.
    """
    for p in paths:
        v = dig(js, p)
        if v not in (None, "", "null"):
            return mpf(v)
    if default is not None:
        return mpf(default)
    raise KeyError(f"none of the paths exist: {paths}")


# ---------- v2.1-aware extractors ----------

def get_band_margin(js):
    """
    band_cert.json (v2.1 expected):

      kind = "band_cert"
      band_cert {
        band_margin_lo,
        ...
      }

    Fallbacks kept for safety.
    """
    return mpf_from_paths(
        js,
        [
            ("band_cert", "band_margin_lo"),
            ("band_cert", "band_margin", "lo"),
            ("numbers", "band_margin"),
            ("band_margin",),
        ],
    )


def get_prime_block_cap(js):
    """
    prime_block_norm.json (v2.1 expected):

      prime_block_norm {
        used_operator_norm,
        operator_norm_cap_hi,
        method,
        sum_zeros_contrib,
        tail_bound_hi,
        cap_total_hi
      }

    Canonical cap is cap_total_hi; fall back to used_operator_norm, etc.
    """
    return mpf_from_paths(
        js,
        [
            ("prime_block_norm", "cap_total_hi"),
            ("prime_block_norm", "used_operator_norm"),
            ("numbers", "cap_total_hi"),
            ("used_operator_norm",),
            ("operator_norm_cap_hi",),
            ("operator_norm_cap",),
        ],
    )


def get_prime_tail_norm(js):
    """
    prime_tail_envelope.json (v2.1 expected):

      kind = "prime_tail_envelope"
      prime_tail {
        env_T0_hi,
        norm
      }

    Canonical scalar is prime_tail.norm.
    """
    return mpf_from_paths(
        js,
        [
            ("prime_tail", "norm"),
            ("numbers", "prime_tail_norm"),
            ("prime_tail_norm",),
        ],
    )


def get_gamma_env_T0(js):
    """
    gamma_tails.json (from core_integral_prover, v2.1):

      kind = "gamma_tails"
      gamma_tails {
        gamma_env_at_T0,
        c1,
        c2,
        tails_total
      }

    Canonical quantity for the continuum budget is gamma_env_at_T0.

    If absent, fall back to tails_total or 0.
    """
    return mpf_from_paths(
        js,
        [
            ("gamma_tails", "gamma_env_at_T0"),
            ("gamma_env_at_T0",),
            ("gamma_tails", "tails_total"),
            ("tails_total",),
        ],
        default="0",
    )


def get_grid_error(js):
    """
    grid_error_bound.json (v2.1):

      kind = "grid_error_bound"
      grid_error_bound {
        bound_hi
      }

    If missing, default to 0 (grid error disabled).
    """
    return mpf_from_paths(
        js,
        [
            ("grid_error_bound", "bound_hi"),
            ("numbers", "grid_error_norm"),
            ("grid_error_norm",),
        ],
        default="0",
    )


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(
        description="Continuum operator inequality rollup (v2.1 normalized)."
    )
    ap.add_argument("--band-cert", required=True, help="band_cert.json")
    ap.add_argument("--prime-block", required=True, help="prime_block_norm.json")
    ap.add_argument("--prime-tail", required=True, help="prime_tail_envelope.json")
    ap.add_argument(
        "--gamma-tails",
        required=False,
        default=None,
        help="gamma_tails.json (core_integral_prover)",
    )
    ap.add_argument(
        "--grid-error",
        required=False,
        default=None,
        help="grid_error_bound.json (op_grid_error_bound)",
    )
    ap.add_argument("--dps", type=int, default=200, help="mpmath precision")
    ap.add_argument("--out", required=True, help="Output JSON path.")
    args = ap.parse_args()

    set_prec(args.dps)

    band_js = jload(args.band_cert)
    pblk_js = jload(args.prime_block)
    ptail_js = jload(args.prime_tail)

    gamma_js = {}
    if args.gamma_tails and os.path.exists(args.gamma_tails):
        gamma_js = jload(args.gamma_tails)

    grid_js = {}
    if args.grid_error and os.path.exists(args.grid_error):
        grid_js = jload(args.grid_error)

    band_margin = get_band_margin(band_js)
    prime_block_cap = get_prime_block_cap(pblk_js)
    prime_tail_norm = get_prime_tail_norm(ptail_js)
    gamma_env_T0 = get_gamma_env_T0(gamma_js) if gamma_js else mp.mpf("0")
    grid_error_norm = get_grid_error(grid_js) if grid_js else mp.mpf("0")

    lhs_total = prime_block_cap + prime_tail_norm + grid_error_norm
    epsilon_eff = band_margin - gamma_env_T0
    PASS = bool(lhs_total <= epsilon_eff)

    inputs_block = {
        "band_cert_path": args.band_cert or "",
        "prime_block_path": args.prime_block or "",
        "prime_tail_path": args.prime_tail or "",
        "gamma_tails_path": args.gamma_tails or "",
        "grid_error_path": args.grid_error or "",
    }

    out = {
        "kind": "continuum_operator_cert",
        "inputs": inputs_block,
        "numbers": {
            "band_margin": mpstr(band_margin),
            "prime_block_cap": mpstr(prime_block_cap),
            "prime_tail_norm": mpstr(prime_tail_norm),
            "gamma_tails": mpstr(gamma_env_T0),
            "grid_error_norm": mpstr(grid_error_norm),
            "lhs_total": mpstr(lhs_total),
            "epsilon_eff": mpstr(epsilon_eff),
        },
        "PASS": PASS,
        "meta": {
            "tool": "continuum_operator_rollup",
            "dps": str(mp.dps),
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    }

    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    status = "PASS" if PASS else "FAIL"
    print(
        f"[{status}] continuum op -> {args.out}  "
        f"lhs_total={mpstr(lhs_total)}  epsilon_eff={mpstr(epsilon_eff)}"
    )


if __name__ == "__main__":
    main()
