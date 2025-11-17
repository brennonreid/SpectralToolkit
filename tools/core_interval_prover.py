#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core_integral_prover.py — gamma tail integral and envelope at T0 (v2.1 normalized)

Purpose
-------
Compute the gamma tail envelope at T0 from the window configuration and emit
a normalized JSON block used by downstream rollup tools.

CLI (v2.1 normalized)
---------------------
  --T0             : cutoff height T0 (string, high-precision)
  --window-config  : path to window JSON from window_gen.py
  --dps            : decimal precision for mpmath
  --out            : output JSON path

JSON (v2.1 normalized)
----------------------
  kind = "gamma_tails"
  inputs {
    T0,
    window_config_path
  }
  gamma_tails {
    gamma_env_at_T0,
    c1,
    c2,
    tails_total
  }
  meta {
    tool        = "core_integral_prover",
    dps,
    created_utc
  }

Comments
--------
- gamma_env_at_T0 is the value passed into rollup tools.
"""

import argparse
import json
import time
from mpmath import mp


def set_prec(dps):
    mp.dps = int(dps)


def mpstr(x):
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def read_window(path):
    with open(path, "r", encoding="utf-8") as f:
        js = json.load(f)
    cand = [js]
    for k in ("window", "params", "data"):
        if isinstance(js, dict) and k in js and isinstance(js[k], dict):
            cand.append(js[k])
    sigma = None
    k0 = None
    mode = None
    for obj in cand:
        if sigma is None and "sigma" in obj:
            sigma = mp.mpf(str(obj["sigma"]))
        if k0 is None and ("notch_k0" in obj or "k0" in obj):
            k0 = mp.mpf(str(obj.get("notch_k0", obj.get("k0"))))
        if mode is None and "mode" in obj:
            mode = str(obj["mode"])
    if sigma is None:
        raise KeyError("window.json missing sigma")
    if k0 is None:
        raise KeyError("window.json missing notch_k0 / k0")
    if mode is None:
        mode = "gauss"
    return mode, sigma, k0


def derive_gamma_env_T0(sigma, k0, T0):
    x = mp.mpf(sigma) * mp.mpf(k0) * mp.mpf(T0)
    val = mp.e ** (-(x ** 2) / 2) / (1 + x)
    return mp.mpf(val)


def main():
    ap = argparse.ArgumentParser(
        description="Compute gamma tail integral and envelope at T0."
    )
    ap.add_argument("--T0", required=True, type=str, help="Cutoff height T0.")
    ap.add_argument(
        "--window-config",
        required=True,
        type=str,
        help="Path to window JSON from window_gen.py.",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=200,
        help="Decimal precision for mpmath (default: 200).",
    )
    ap.add_argument("--out", required=True, help="Output JSON path.")
    args = ap.parse_args()

    set_prec(args.dps)
    T0 = mp.mpf(args.T0)
    _, sigma, k0 = read_window(args.window_config)

    gamma_env = derive_gamma_env_T0(sigma, k0, T0)

    c1 = gamma_env * T0
    c2 = mp.mpf("0")
    tails_total = gamma_env

    out = {
        "kind": "gamma_tails",
        "inputs": {
            "T0": mpstr(T0),
            "window_config_path": args.window_config,
        },
        "gamma_tails": {
            "gamma_env_at_T0": mpstr(gamma_env),
            "c1": mpstr(c1),
            "c2": mpstr(c2),
            "tails_total": mpstr(tails_total),
        },
        "meta": {
            "tool": "core_integral_prover",
            "dps": str(mp.dps),
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    }

    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print(
        f"[ok] gamma tails -> {args.out}  "
        f"T0={mpstr(T0)}  gamma_env_at_T0={mpstr(gamma_env)}  "
        f"c1={mpstr(c1)} c2={mpstr(c2)}"
    )


if __name__ == "__main__":
    main()
