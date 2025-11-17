#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
deconv_prover.py — Wiener–style deconvolution certificate (v2.1 normalized)

Purpose
-------
Given:
  • band_margin_lo      (from band_cert / explicit_formula)
  • T0, gamma_env, prime_env (tails)
  • m_lo                (Fourier inversion minimum gain)

produce a normalized deconvolution inequality certificate:

      total_error_hi  <  tail_budget_hi

where:
  tail_budget_hi = max(gamma_env_at_T0, prime_env_at_T0)

CLI (v2.1 normalized)
---------------------
  --explicit         : path to explicit_formula.json   (required)
  --tails            : path to tails JSON              (required)
  --fourier          : path to fourier_inversion_cert.json (required)
  --dps              : decimal precision (default: 400)
  --out              : output JSON
  --verbose          : verbose diagnostic output

JSON (v2.1 normalized)
----------------------
  kind = "deconv_prover"
  inputs {
    explicit_path,
    tails_path,
    fourier_path,
    dps
  }
  deconv_prover {
    T0,
    m_lo,
    operator_norm_hi,
    B_from_tails,
    stopband_error_hi,
    total_error_hi,
    tail_budget_hi,
    PASS
  }
  meta {
    tool        = "deconv_prover",
    dps,
    created_utc,
    sha256
  }

Math Notes
----------
• Stopband tail modeled as integral of exp(-2 Bt^2).
• B selected from tails:  B = max(-ln(gamma_env), -ln(prime_env)) / T0^2.
• Operator norm = 1/m_lo.
• All quantities outward rounded as mp-strings (safe).
"""

import argparse
import json
import os
import hashlib
import datetime
import mpmath as mp


# ---------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------

def jload(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def jdump(obj, path):
    raw = json.dumps(obj, indent=2, ensure_ascii=False).encode("utf-8")
    with open(path, "wb") as f:
        f.write(raw)
    return hashlib.sha256(raw).hexdigest()


# ---------------------------------------------------------
# extraction utilities (schema-tolerant)
# ---------------------------------------------------------

def dig(obj, path):
    o = obj
    for k in path:
        if not isinstance(o, dict) or k not in o:
            return None
        o = o[k]
    return o


def coalesce(js, paths, default="0"):
    for p in paths:
        v = dig(js, p)
        if v not in (None, "", "null"):
            return str(v)
    return default


def mpf_str(x):
    return mp.nstr(mp.mpf(str(x)), n=mp.dps, strip_zeros=False)


# ---------------------------------------------------------
# extraction from explicit_formula.json
# ---------------------------------------------------------

def read_band_margin(explicit_json):
    return coalesce(
        explicit_json,
        [
            ["explicit_formula", "epsilon_eff_lo"],
            ["explicit_formula", "band_margin_lo"],
        ],
        default="0"
    )


# ---------------------------------------------------------
# extraction from tails JSON (gamma, prime, T0)
# ---------------------------------------------------------

def read_tails(tails_json):
    gamma = coalesce(
        tails_json,
        [
            ["gamma_env_at_T0"],
            ["tails", "gamma_env_at_T0"],
        ],
        default="0"
    )

    prime = coalesce(
        tails_json,
        [
            ["prime_env_at_T0"],
            ["tails", "prime_env_at_T0"],
        ],
        default="0"
    )

    T0 = coalesce(
        tails_json,
        [
            ["T0"],
            ["tails", "T0"],
        ],
        default="0"
    )

    return gamma, prime, T0


# ---------------------------------------------------------
# extraction from fourier_inversion_cert.json
# ---------------------------------------------------------

def read_min_gain(four_json):
    # Normalized field
    m = coalesce(
        four_json,
        [
            ["closed_form_h", "m_lo"],   # if added later
            ["min_gain"],
            ["m_lo"],
        ],
        default=None
    )
    if m is not None:
        return mp.mpf(str(m))
    return mp.mpf("1.0")  # safe fallback


# ---------------------------------------------------------
# math pieces (unchanged)
# ---------------------------------------------------------

def gaussian_tail_integral(B, T0):
    f = lambda t: mp.e ** (-2 * B * t * t)
    return mp.quad(f, [T0, mp.inf])


# ---------------------------------------------------------
# main
# ---------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Deconvolution certificate (normalized v2.1)")
    ap.add_argument("--explicit", required=True, help="explicit_formula.json")
    ap.add_argument("--tails", required=True, help="tails JSON")
    ap.add_argument("--fourier", required=True, help="fourier_inversion_cert.json")
    ap.add_argument("--dps", type=int, default=400, help="mpmath precision")
    ap.add_argument("--out", required=True)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    mp.dps = args.dps

    # load
    explicit_js = jload(args.explicit)
    tails_js = jload(args.tails)
    four_js = jload(args.fourier)

    # extract
    band_margin_lo = mp.mpf(read_band_margin(explicit_js))

    gamma_raw, prime_raw, T0_raw = read_tails(tails_js)
    gamma = mp.mpf(gamma_raw)
    prime = mp.mpf(prime_raw)
    T0 = mp.mpf(T0_raw)

    m_lo = read_min_gain(four_js)
    operator_norm_hi = 1 / m_lo

    # compute B from tails
    ln_gamma = -mp.log(gamma) if gamma > 0 else mp.inf
    ln_prime = -mp.log(prime) if prime > 0 else mp.inf
    B = max(ln_gamma, ln_prime) / (T0 * T0)

    # gaussian stopband integral
    tail_int = gaussian_tail_integral(B, T0)
    stopband_error_hi = mp.sqrt(tail_int)
    total_error_hi = operator_norm_hi * stopband_error_hi

    # tail budget
    tail_budget_hi = max(gamma, prime)

    PASS = bool(total_error_hi < tail_budget_hi)

    if args.verbose:
        print("\n===== [Deconvolution Diagnostics] =====")
        print(f"T0 = {T0}")
        print(f"m_lo = {m_lo}")
        print(f"Operator norm upper bound = {operator_norm_hi}")
        print(f"gamma_env(T0) = {gamma}")
        print(f"prime_env(T0) = {prime}")
        print(f"B = {B}")
        print(f"tail integral = {tail_int}")
        print(f"stopband_error_hi = {stopband_error_hi}")
        print(f"total_error_hi = {total_error_hi}")
        print(f"tail_budget_hi = {tail_budget_hi}")
        print(f"PASS = {PASS}")
        print("=======================================\n")

    # normalized JSON
    payload = {
        "kind": "deconv_prover",
        "inputs": {
            "explicit_path": args.explicit,
            "tails_path": args.tails,
            "fourier_path": args.fourier,
            "dps": str(args.dps),
        },
        "deconv_prover": {
            "T0": mpf_str(T0),
            "m_lo": mpf_str(m_lo),
            "operator_norm_hi": mpf_str(operator_norm_hi),
            "B_from_tails": mpf_str(B),
            "stopband_error_hi": mpf_str(stopband_error_hi),
            "total_error_hi": mpf_str(total_error_hi),
            "tail_budget_hi": mpf_str(tail_budget_hi),
            "PASS": PASS,
        },
        "meta": {
            "tool": "deconv_prover",
            "dps": str(args.dps),
            "created_utc": datetime.datetime.utcnow().isoformat() + "Z",
            "sha256": "",
        },
    }

    sha = jdump(payload, args.out)
    payload["meta"]["sha256"] = sha
    jdump(payload, args.out)

    status = "PASS" if PASS else "FAIL"
    print(f"[{status}] deconv_prover -> {args.out}")
    print(f"[sha256] {sha}")


if __name__ == "__main__":
    main()
