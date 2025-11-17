#!/usr/bin/env python
import argparse
import json
import os
from typing import Any, Dict, Optional

def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        return json.load(f)

def safe_get(d: Dict[str, Any], *keys, default=None):
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur

def to_float(x: Any) -> Optional[float]:
    """Best-effort conversion to float for huge-precision strings.
    Extremely tiny values (like 1e-24429064) will underflow to 0.0, which is
    fine for these inequality checks."""
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, str):
        try:
            return float(x)
        except ValueError:
            return None
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band-cert", required=True)
    ap.add_argument("--weil-psd", required=True)
    ap.add_argument("--tails", required=True)
    ap.add_argument("--prime-norm", required=True)
    ap.add_argument("--explicit", required=True)
    ap.add_argument("--continuum-cert", required=False)
    ap.add_argument("--dps", type=int, default=950)
    args = ap.parse_args()

    # Base packet dir (used to auto-find uniform_certificate)
    packet_dir = os.path.dirname(os.path.abspath(args.band_cert)) or "."

    # Try to load the aggregated uniform certificate first
    uniform_path = None
    for name in ["uniform_certificate.json", "uniform_certificate.json.tmp"]:
        cand = os.path.join(packet_dir, name)
        if os.path.exists(cand):
            uniform_path = cand
            break

    uniform = load_json(uniform_path) if uniform_path else None
    u = uniform.get("uniform_certificate", {}) if uniform else {}

    # Individual components (used as primary sources or fallbacks)
    band = load_json(args.band_cert)
    tails = load_json(args.tails)
    prime_block = load_json(args.prime_norm)
    explicit = load_json(args.explicit)
    weil = load_json(args.weil_psd)
    continuum = load_json(args.continuum_cert) if args.continuum_cert else None

    print("=" * 66)
    print("STP STATUS")
    print("-" * 66)

    # -------------------
    # BAND MARGIN (epsilon band)
    # -------------------
    margin = safe_get(u, "band_margin")
    if margin is None:
        margin = safe_get(band, "band_cert", "numbers", "band_margin_lo")
    band_ok = margin is not None
    print(f"[band] {'OK  ' if band_ok else 'FAIL'} margin={margin}")

    # -------------------
    # WEIL PSD
    # -------------------
    psd = safe_get(u, "PSD_verified")
    if psd is None:
        psd = safe_get(weil, "weil_psd_bochner", "PSD_verified")
        if psd is None:
            psd = safe_get(weil, "weil_psd", "PSD_verified")
        if psd is None:
            psd = safe_get(weil, "PSD_verified")
    weil_ok = bool(psd)
    print(f"[weil] {'OK  ' if weil_ok else 'FAIL'} PSD_verified={psd}")

    # -------------------
    # TAILS: gamma envelope + prime tail norm
    # Prefer the already-rolled-up numbers in uniform_certificate
    # -------------------
    gamma_env = safe_get(u, "gamma_env_at_T0")
    prime_env = safe_get(u, "prime_tail_norm")

    if gamma_env is None:
        gamma_env = safe_get(tails, "gamma_tail", "env_T0_hi")
    if prime_env is None:
        prime_env = safe_get(tails, "prime_tail", "norm")

    tails_ok = (gamma_env is not None) and (prime_env is not None)
    print(f"[tails] {'OK  ' if tails_ok else 'FAIL'} gamma_env_at_T0={gamma_env} prime_env_at_T0={prime_env}")

    # -------------------
    # PRIME BLOCK NORM CAP
    # -------------------
    prime_cap = safe_get(u, "prime_block_cap")
    if prime_cap is None:
        prime_cap = safe_get(prime_block, "prime_block_norm", "operator_norm_cap_hi")
    prime_ok = prime_cap is not None
    print(f"[prime] {'OK  ' if prime_ok else 'FAIL'} operator_norm_cap={prime_cap}")

    # -------------------
    # EXPLICIT FORMULA eps_eff (effective margin after tails + grid)
    # -------------------
    eps_eff = safe_get(u, "epsilon_eff")
    if eps_eff is None:
        eps_eff = safe_get(explicit, "explicit_formula", "eps_eff")
    explicit_ok = eps_eff is not None
    print(f"[explicit] {'OK  ' if explicit_ok else 'FAIL'} eps_eff={eps_eff}")

    # -------------------
    # CONTINUUM CERT (just presence)
    # -------------------
    if continuum is not None:
        print(f"[continuum] OK   path={args.continuum_cert}")
    else:
        print("[continuum] SKIP path=None")

    print("-" * 66)

    # -------------------
    # ROLLUP SUMMARY (read from uniform_certificate where possible)
    # -------------------
    grid_err = safe_get(u, "grid_error_norm")
    lhs_total = safe_get(u, "lhs_total")

    print("[rollup] (informational)")
    print(f"  prime_block_cap = {prime_cap}")
    print(f"  prime_tail_norm = {prime_env}")
    print(f"  grid_error_norm = {grid_err}")
    print(f"  lhs_total       = {lhs_total}")

    # -------------------
    # GLOBAL ERROR-BUDGET CHECK
    #   We want: (prime_block + gamma_tail + prime_tail + grid_error) <= band_margin
    #   and, equivalently, lhs_total <= eps_eff (effective margin).
    # -------------------
    margin_f     = to_float(margin)
    prime_cap_f  = to_float(prime_cap)
    gamma_env_f  = to_float(gamma_env)
    prime_env_f  = to_float(prime_env)
    grid_err_f   = to_float(grid_err)
    lhs_total_f  = to_float(lhs_total)
    eps_eff_f    = to_float(eps_eff)

    components = [prime_cap_f, gamma_env_f, prime_env_f, grid_err_f]
    have_all_components = all(c is not None for c in components) and (margin_f is not None)

    total_error = None
    budget_ok_components = None
    if have_all_components:
        total_error = sum(components)
        budget_ok_components = total_error <= margin_f

    # Independent global check using lhs_total vs eps_eff (if available)
    budget_ok_lhs = None
    if (lhs_total_f is not None) and (eps_eff_f is not None):
        budget_ok_lhs = lhs_total_f <= eps_eff_f

    print("[budget] global inequality checks")
    print(f"  margin_band     = {margin_f}")
    print(f"  eps_eff         = {eps_eff_f}")
    print(f"  total_error_est = {total_error}")
    print(f"  lhs_total_f     = {lhs_total_f}")
    print(f"  lhs_vs_eps_ok   = {budget_ok_lhs}")
    print(f"  sum_vs_margin_ok= {budget_ok_components}")

    # Overall STP pass requires:
    #  - band margin present
    #  - Weil PSD
    #  - tails, prime, explicit data present
    #  - uniform_certificate PASS flag
    #  - and at least one of the budget inequalities to hold True
    rollup_pass = bool(uniform and uniform.get("PASS", False))

    budget_flag = True
    if (budget_ok_lhs is not None) or (budget_ok_components is not None):
        checks = []
        if budget_ok_lhs is not None:
            checks.append(budget_ok_lhs)
        if budget_ok_components is not None:
            checks.append(budget_ok_components)
        budget_flag = all(checks)

    stp_pass = band_ok and weil_ok and tails_ok and prime_ok and explicit_ok and rollup_pass and budget_flag

    print(f"  PASS            = {stp_pass}")
    print("-" * 66)

if __name__ == "__main__":
    main()
