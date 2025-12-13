#!/usr/bin/env python3
"""
quartet_perturbation_scan.py

Prototype module for scanning the spectral quartet penalty for a given window
and computing the contradiction threshold delta_c via the Spectral Quartet Lemma.

IMPORTANT:
- This is a structural / architectural skeleton.
- You MUST implement the actual kernel H_theta(s; T) that matches your explicit
  formula. The placeholder below raises NotImplementedError.
- The estimates for C0 and R_max are numerical, not interval-rigorous.
  For production RH work, replace with outward-rounded bounds.
"""

import argparse
import json
import math
import hashlib
from datetime import datetime

from mpmath import mp, mpc, re as mp_re

# ----------------------------------------------------------------------
# 1. Window / kernel hook
# ----------------------------------------------------------------------

class Kernel:
    """
    Abstract interface for the explicit-formula kernel H_theta(s; T) corresponding
    to your test function h_theta(t).

    You must implement eval_H and (preferably) eval_H_second_derivative to match
    the exact kernel used in your Q functional.

    NOTE: The default implementation raises NotImplementedError so you cannot
    accidentally treat this as the real RH kernel.
    """
    def __init__(self, window_config: dict):
        self.window_config = window_config
        # TODO: parse sigma, k0, notch params, etc. from window_config as needed.

    def eval_H(self, s: mpc, T: mp.mpf) -> mpc:
        """
        H_theta(s; T) as used in your explicit formula.

        Parameters:
            s : complex (mp.mpc) argument, typically 1/2 + i*gamma or 1/2 +/- delta + i*gamma
            T : real height parameter T (mp.mpf)

        Returns:
            complex (mp.mpc) value H_theta(s; T)

        You MUST implement this using your actual test function.
        """
        raise NotImplementedError("Kernel.eval_H must be implemented for your window/test function.")

    def eval_H_second_derivative(self, s: mpc, T: mp.mpf) -> mpc:
        """
        Second derivative of H_theta(s; T) with respect to the real part of s
        (i.e., along the sigma direction) evaluated at s.

        This is used to approximate the curvature coefficient C(γ, T).

        Default: numerical second derivative via mp.diff in the real direction.
        For RH-grade rigor, replace with an analytic formula and outward rounding.
        """
        # Numerical second derivative in the sigma direction:
        # H''(s) ~ d^2/dsigma^2 H(sigma + i*t, T)
        h = mp.mpf('1e-6')  # step size for finite difference (tunable)
        sigma = mp.re(s)
        t = mp.im(s)

        def H_sigma(sig):
            return self.eval_H(mpc(sig, t), T)

        # Central finite difference for second derivative:
        # H''(sigma) ≈ (H(sigma + h) - 2 H(sigma) + H(sigma - h)) / h^2
        H_plus = H_sigma(sigma + h)
        H_0 = H_sigma(sigma)
        H_minus = H_sigma(sigma - h)
        return (H_plus - 2 * H_0 + H_minus) / (h * h)


# ----------------------------------------------------------------------
# 2. Quartet penalty and curvature
# ----------------------------------------------------------------------

def delta_Q_quartet(kernel: Kernel, delta: mp.mpf, gamma: mp.mpf, T: mp.mpf) -> mp.mpf:
    """
    Quartet contribution ΔQ_theta(delta, gamma; T) for a single hypothetical quartet
    at off-line distance delta and height gamma.

    ΔQ_theta(delta, gamma; T) = -4 * Re[ H(1/2 + delta + i*gamma; T) + H(1/2 - delta + i*gamma; T) ]

    Sign convention: negative ΔQ pushes Q downward.
    """
    s1 = mpc(mp.mpf('0.5') + delta, gamma)
    s2 = mpc(mp.mpf('0.5') - delta, gamma)
    H1 = kernel.eval_H(s1, T)
    H2 = kernel.eval_H(s2, T)
    return -4 * mp_re(H1 + H2)


def curvature_C_gamma_T(kernel: Kernel, gamma: mp.mpf, T: mp.mpf) -> mp.mpf:
    """
    Curvature coefficient C(γ, T) from the second derivative of H at s = 1/2 + i*γ:

        C(γ,T) = 4 * Re( H''(1/2 + iγ; T) ).

    NOTE: This assumes eval_H_second_derivative computes second derivative in the
    sigma direction (real part).
    """
    s0 = mpc(mp.mpf('0.5'), gamma)
    Hpp = kernel.eval_H_second_derivative(s0, T)
    return 4 * mp_re(Hpp)


# ----------------------------------------------------------------------
# 3. Numerical estimates for C0 and R_max over grids
# ----------------------------------------------------------------------

def estimate_C0(kernel: Kernel, gamma_min, gamma_max, gamma_steps,
                T_min, T_max, T_steps):
    """
    Estimate a numerical lower bound C0 ≈ min_{γ,T} C(γ,T) over the given grids.

    This is NOT rigorous; it just takes the minimum over sample points.

    Returns:
        C0_est (mp.mpf), C_grid (list of dicts) for debugging / inspection.
    """
    C0 = None
    C_grid = []
    for j in range(gamma_steps):
        gamma = gamma_min + (gamma_max - gamma_min) * j / (gamma_steps - 1 if gamma_steps > 1 else 1)
        gamma = mp.mpf(gamma)
        for k in range(T_steps):
            T = T_min + (T_max - T_min) * k / (T_steps - 1 if T_steps > 1 else 1)
            T = mp.mpf(T)
            C_val = curvature_C_gamma_T(kernel, gamma, T)
            C_grid.append({
                "gamma": str(gamma),
                "T": str(T),
                "C_gamma_T": str(C_val),
            })
            if C0 is None or C_val < C0:
                C0 = C_val
    return C0, C_grid


def estimate_Rmax(kernel: Kernel, C0: mp.mpf, delta_max_radius: mp.mpf,
                  gamma_min, gamma_max, gamma_steps,
                  T_min, T_max, T_steps,
                  delta_samples: int):
    """
    Numerical estimate of R_max for the expansion:

        ΔQ(delta) <= -C0 * delta^2 + R_max * delta^4

    We approximate R_max as:

        R_max >= sup_{|delta|<=delta_max_radius, γ,T}
                    (ΔQ(delta) + C0 delta^2) / delta^4

    using a discrete grid in delta, gamma, T.

    Returns:
        R_max_est (mp.mpf)
    """
    R_max = mp.mpf('0')
    # Sample symmetric deltas in (0, delta_max_radius]
    for i in range(1, delta_samples + 1):
        delta = delta_max_radius * i / delta_samples
        delta = mp.mpf(delta)
        for j in range(gamma_steps):
            gamma = gamma_min + (gamma_max - gamma_min) * j / (gamma_steps - 1 if gamma_steps > 1 else 1)
            gamma = mp.mpf(gamma)
            for k in range(T_steps):
                T = T_min + (T_max - T_min) * k / (T_steps - 1 if T_steps > 1 else 1)
                T = mp.mpf(T)
                dQ = delta_Q_quartet(kernel, delta, gamma, T)
                # rearrange: ΔQ(delta) <= -C0 delta^2 + R_max delta^4
                # => R_max >= (ΔQ(delta) + C0 delta^2) / delta^4
                numer = dQ + C0 * delta * delta
                denom = delta**4
                R_candidate = numer / denom
                if R_candidate > R_max:
                    R_max = R_candidate
    # Ensure nonnegative
    if R_max < 0:
        R_max = mp.mpf('0')
    return R_max


# ----------------------------------------------------------------------
# 4. Contradiction threshold computation (from the lemma)
# ----------------------------------------------------------------------

def compute_delta_c(C0: mp.mpf, R_max: mp.mpf, m0: mp.mpf, delta_valid: mp.mpf):
    """
    Given C0, R_max, margin m0 and validity radius delta_valid (delta_1),
    compute the contradiction threshold delta_c and the upper bound sqrt(x_max)
    according to the Spectral Quartet Exclusion Lemma.

    Returns:
        {
          "discriminant": str(...),
          "has_leverage": bool,
          "delta_c": str or None,
          "sqrt_x_max": str or None,
          "valid_region_min": str or None,
          "valid_region_max": str or None
        }
    """
    two = mp.mpf('2')
    four = mp.mpf('4')

    disc = C0*C0 - four*R_max*m0

    result = {
        "discriminant": str(disc),
        "has_leverage": False,
        "delta_c": None,
        "sqrt_x_max": None,
        "valid_region_min": None,
        "valid_region_max": None
    }

    if disc < 0:
        # No real roots: this band is unresolved with current bounds
        return result

    sqrt_disc = mp.sqrt(disc)
    x_min = (C0 - sqrt_disc) / (two * R_max)
    x_max = (C0 + sqrt_disc) / (two * R_max)

    # Sanity: x_min, x_max >= 0
    if x_min < 0:
        x_min = mp.mpf('0')
    if x_max < 0:
        x_max = mp.mpf('0')

    delta_c = mp.sqrt(x_min)
    sqrt_x_max = mp.sqrt(x_max)

    # Valid region is |delta| in [delta_c, min(sqrt_x_max, delta_valid)]
    valid_max = sqrt_x_max if sqrt_x_max < delta_valid else delta_valid

    result["has_leverage"] = delta_c <= delta_valid
    result["delta_c"] = str(delta_c)
    result["sqrt_x_max"] = str(sqrt_x_max)
    result["valid_region_min"] = str(delta_c)
    result["valid_region_max"] = str(valid_max)
    return result


# ----------------------------------------------------------------------
# 5. CLI and JSON output
# ----------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Prototype spectral quartet perturbation scan."
    )
    p.add_argument("--window-config", type=str, required=True,
                   help="Path to window.json (test function parameters).")
    p.add_argument("--out", type=str, required=True,
                   help="Output JSON file for quartet scan results.")

    p.add_argument("--dps", type=int, default=200,
                   help="mpmath decimal precision.")
    p.add_argument("--gamma-min", type=str, required=True,
                   help="Minimum gamma (imaginary part) in the band.")
    p.add_argument("--gamma-max", type=str, required=True,
                   help="Maximum gamma in the band.")
    p.add_argument("--gamma-steps", type=int, default=5,
                   help="Number of gamma samples in the band.")

    p.add_argument("--T-min", type=str, required=True,
                   help="Minimum T in the tested range.")
    p.add_argument("--T-max", type=str, required=True,
                   help="Maximum T in the tested range.")
    p.add_argument("--T-steps", type=int, default=5,
                   help="Number of T samples in the range.")

    p.add_argument("--delta-valid", type=str, required=True,
                   help="Validity radius delta_1 for the Taylor bound.")
    p.add_argument("--delta-samples", type=int, default=5,
                   help="Number of delta samples in (0, delta_valid] for R_max estimate.")

    p.add_argument("--margin-m0", type=str, required=True,
                   help="Certified positive margin m0 for Q_nooff in this band.")

    return p.parse_args()


def main():
    args = parse_args()

    mp.mp.dps = args.dps

    # Load window config (you already know the schema; we do not assume specifics)
    with open(args.window_config, "r", encoding="utf-8") as f:
        window_config = json.load(f)

    # Instantiate kernel (you must implement eval_H)
    kernel = Kernel(window_config)

    gamma_min = mp.mpf(args.gamma_min)
    gamma_max = mp.mpf(args.gamma_max)
    T_min = mp.mpf(args.T_min)
    T_max = mp.mpf(args.T_max)
    delta_valid = mp.mpf(args.delta_valid)
    m0 = mp.mpf(args.margin_m0)

    # Estimate C0 over (gamma, T) grid
    C0_est, C_grid = estimate_C0(
        kernel,
        gamma_min, gamma_max, args.gamma_steps,
        T_min, T_max, args.T_steps
    )

    # Estimate R_max using delta in (0, delta_valid]
    R_max_est = estimate_Rmax(
        kernel,
        C0_est,
        delta_valid,
        gamma_min, gamma_max, args.gamma_steps,
        T_min, T_max, args.T_steps,
        args.delta_samples
    )

    # Compute contradiction threshold delta_c, etc.
    threshold_info = compute_delta_c(C0_est, R_max_est, m0, delta_valid)

    # Build result JSON
    result = {
        "tool": "quartet_perturbation_scan",
        "version": "0.1-prototype",
        "meta": {
            "created_utc": datetime.utcnow().isoformat() + "Z",
            "dps": args.dps,
        },
        "window_config": window_config,
        "inputs": {
            "gamma_min": args.gamma_min,
            "gamma_max": args.gamma_max,
            "gamma_steps": args.gamma_steps,
            "T_min": args.T_min,
            "T_max": args.T_max,
            "T_steps": args.T_steps,
            "delta_valid": args.delta_valid,
            "delta_samples": args.delta_samples,
            "margin_m0": args.margin_m0,
        },
        "estimates": {
            "C0": str(C0_est),
            "R_max": str(R_max_est),
        },
        "threshold": threshold_info,
        "debug": {
            # beware: can be large; trim or drop in production
            "C_grid": C_grid
        }
    }

    # SHA256 of the JSON (canonicalized)
    canonical = json.dumps(result, sort_keys=True).encode("utf-8")
    sha = hashlib.sha256(canonical).hexdigest()
    result["meta"]["sha256"] = sha

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, sort_keys=True)

    print(f"[quartet_perturbation_scan] -> {args.out}")
    print(f"[quartet_perturbation_scan] sha256={sha}")


if __name__ == "__main__":
    main()
