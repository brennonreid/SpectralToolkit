#!/usr/bin/env python3
# tools/analytic_tail_fit.py

import argparse, json, os, sys, time
from typing import Any, Dict
from mpmath import mp, mpf


def nows() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def dig(obj: Any, path):
    cur = obj
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def coalesce(obj, paths, default=None):
    for p in paths:
        v = dig(obj, p)
        if v not in (None, "", "null"):
            return v
    return default


def mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def main():
    ap = argparse.ArgumentParser(
        description="Fit analytic 1/T^a tail models to numeric gamma/prime tails."
    )
    ap.add_argument(
        "--packet-dir",
        dest="packet_dir",
        required=True,
        help="Path to PROOF_PACKET directory.",
    )
    ap.add_argument(
        "--Ap",
        type=str,
        default="1.0",
        help="Exponent a_p for prime tail model.",
    )
    ap.add_argument(
        "--Ag",
        type=str,
        default="1.0",
        help="Exponent a_g for gamma tail model.",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="Output analytic_tail_fit.json (default: <packet-dir>/analytic_tail_fit.json)",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=220,
        help="Decimal precision for mpmath.",
    )
    args = ap.parse_args()

    mp.dps = args.dps

    pkt = args.packet_dir
    out_path = args.out or os.path.join(pkt, "analytic_tail_fit.json")

    # Required files (relative to packet-dir)
    p_cont = os.path.join(pkt, "continuum_operator_cert.json")
    p_grid = os.path.join(pkt, "grid_error_bound.json")
    p_pt = os.path.join(pkt, "prime_tail_envelope.json")

    # Gamma tails: prefer v2.1 gamma_tails.json, fall back to legacy gamma_tail.json
    p_gt_v21 = os.path.join(pkt, "gamma_tails.json")
    p_gt_legacy = os.path.join(pkt, "gamma_tail.json")
    if os.path.exists(p_gt_v21):
        p_gt = p_gt_v21
    else:
        p_gt = p_gt_legacy

    try:
        cont = load_json(p_cont)
        grid = load_json(p_grid)
        pt = load_json(p_pt)
        gt = load_json(p_gt)
    except FileNotFoundError as e:
        print(f"[analytic_tail_fit] ERROR: missing file: {e}", file=sys.stderr)
        sys.exit(1)

    # eps_eff (lower bound), prefer continuum_operator_cert.numbers.epsilon_eff (or eps_eff)
    eps_eff_s = coalesce(
        cont,
        [["numbers", "epsilon_eff"], ["numbers", "eps_eff"]],
        None,
    )
    if eps_eff_s is None:
        print(
            "[analytic_tail_fit] ERROR: could not locate epsilon_eff in "
            "continuum_operator_cert.json (numbers.epsilon_eff / numbers.eps_eff)",
            file=sys.stderr,
        )
        sys.exit(1)
    eps_eff_lo = mp.mpf(str(eps_eff_s))

    # grid_error upper bound (constant)
    grid_hi_s = coalesce(
        grid,
        [["grid_error_bound", "bound_hi"], ["numbers", "grid_error_norm"]],
        None,
    )
    if grid_hi_s is None:
        print(
            "[analytic_tail_fit] ERROR: could not locate grid bound "
            "(grid_error_bound.bound_hi or numbers.grid_error_norm).",
            file=sys.stderr,
        )
        sys.exit(1)
    grid_hi = mp.mpf(str(grid_hi_s))

    # Prime tail T0 and env_T0_hi
    T0_pt_s = coalesce(
        pt,
        [["inputs", "T0"], ["T0"], ["prime_tail", "T0"]],
        None,
    )
    env_pt_s = coalesce(
        pt,
        [
            ["prime_tail", "env_T0_hi"],
            ["prime_tail_envelope", "env_T0_hi"],
            ["numbers", "prime_tail_norm"],
        ],
        None,
    )
    if T0_pt_s is None or env_pt_s is None:
        print(
            "[analytic_tail_fit] ERROR: could not locate prime T0 and/or env_T0_hi "
            "in prime_tail_envelope.json.",
            file=sys.stderr,
        )
        sys.exit(1)
    T0_pt = mp.mpf(str(T0_pt_s))
    env_pt = mp.mpf(str(env_pt_s))

    # Gamma tail T0 and gamma_env_at_T0 (support v2.1 gamma_tails + legacy gamma_tail)
    T0_gt_s = coalesce(
        gt,
        [["inputs", "T0"], ["T0"], ["gamma_tail", "T0"]],
        None,
    )
    env_gt_s = coalesce(
        gt,
        [
            ["gamma_tails", "gamma_env_at_T0"],
            ["gamma_tail", "gamma_env_at_T0"],
            ["gamma_env_at_T0"],
        ],
        None,
    )
    if T0_gt_s is None or env_gt_s is None:
        print(
            "[analytic_tail_fit] ERROR: could not locate gamma T0 and/or gamma_env_at_T0 "
            "in gamma_tails/gamma_tail JSON.",
            file=sys.stderr,
        )
        sys.exit(1)
    T0_gt = mp.mpf(str(T0_gt_s))
    env_gt = mp.mpf(str(env_gt_s))

    # Exponents
    a_p = mp.mpf(str(args.Ap))
    a_g = mp.mpf(str(args.Ag))
    if a_p <= 0 or a_g <= 0:
        print(
            "[analytic_tail_fit] ERROR: exponents Ap, Ag must be positive.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Derive conservative C so that C/T^a >= env_T0_hi at T0
    # => C = env_T0_hi * T0^a
    C_prime = env_pt * (T0_pt ** a_p)
    C_gamma = env_gt * (T0_gt ** a_g)

    out = {
        "kind": "analytic_tail_fit",
        "inputs": {
            "packet_dir": pkt,
            "Ap": str(args.Ap),
            "Ag": str(args.Ag),
            "dps": str(args.dps),
        },
        "bounds": {
            "eps_eff_lo": mp_str(eps_eff_lo),
            "grid_error_hi": mp_str(grid_hi),
            "prime_tail": {
                "C": mp_str(C_prime),
                "a": mp_str(a_p),
                "T0": mp_str(T0_pt),
                "env_T0_hi": mp_str(env_pt),
            },
            "gamma_tail": {
                "C": mp_str(C_gamma),
                "a": mp_str(a_g),
                "T0": mp_str(T0_gt),
                "env_T0_hi": mp_str(env_gt),
            },
        },
        "meta": {
            "tool": "analytic_tail_fit",
            "dps": str(args.dps),
            "created_utc": nows(),
        },
    }

    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print(f"[analytic_tail_fit] wrote {out_path}")
    print(f"  eps_eff_lo   = {eps_eff_lo}")
    print(f"  grid_error_hi= {grid_hi}")
    print(
        f"  prime: C={C_prime}  a={a_p}  "
        f"(from env_T0_hi={env_pt} @ T0={T0_pt})"
    )
    print(
        f"  gamma: C={C_gamma}  a={a_g}  "
        f"(from env_T0_hi={env_gt} @ T0={T0_gt})"
    )


if __name__ == "__main__":
    main()
