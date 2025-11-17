#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fourier_inversion_cert.py — closed-form Fourier inversion certificate.

Purpose:
  Certify the closed-form h(x) associated with the Gaussian notch kernel and
  compute bounds used by downstream tools (e.g. deconv_prover).

CLI (v2.1 normalized):
  --window-config   : path to window JSON from window_gen.py
  --sigma           : Gaussian width (alternative to --window-config)
  --k0              : notch parameter (alternative to --window-config)
  --xmax            : half-width for optional sup probe
  --simpson-n       : retained for compatibility; unused in closed-form mode
  --dps             : decimal precision
  --out             : path to output JSON
  --probe           : enable numeric probe for sup|h|
"""

import argparse
import datetime
import json
import sys

from mpmath import mp


# ---------------- I/O helpers ----------------
def jload(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def jdump(obj, path):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


# ------------- Window parameters -------------
def read_window(args):
    """
    Resolve (sigma, k0) either from a window-config JSON or from raw CLI
    parameters. Exactly one mode must be used per run.

    Returns:
      mode_tag        : "window" or "raw_params"
      window_mode     : window mode string (e.g. "gauss")
      sigma_mp        : mpmath mpf
      k0_mp           : mpmath mpf
      window_path     : path to window JSON or None
    """
    if args.window_config and (args.sigma is not None or args.k0 is not None):
        raise ValueError(
            "Use either --window-config or (--sigma and --k0), but not both."
        )

    if args.window_config:
        js = jload(args.window_config)

        # Prefer a dedicated window block if present, otherwise fall back
        window = js.get("window", js)

        sigma = window.get("sigma")
        k0 = window.get("k0")
        window_mode = window.get("mode", js.get("mode", "gauss"))

        if sigma is None or k0 is None:
            raise KeyError("window JSON missing canonical 'sigma' or 'k0' fields")

        return (
            "window",
            str(window_mode),
            mp.mpf(str(sigma)),
            mp.mpf(str(k0)),
            args.window_config,
        )

    # Raw parameter mode
    if args.sigma is None or args.k0 is None:
        raise ValueError(
            "You must provide either --window-config or both --sigma and --k0."
        )

    return (
        "raw_params",
        "gauss",
        mp.mpf(str(args.sigma)),
        mp.mpf(str(args.k0)),
        None,
    )


# Closed-form mapping for Gaussian notch window
def gauss_minus_gauss_closed_form(sigma, k0):
    """
    Map (sigma, k0) to parameters of

        h(x) = A1 * exp(-a1 * x^2) - A2 * exp(-a2 * x^2),

    matching the earlier closed-form choice used in the toolkit.
    """
    alpha1 = 1 / (sigma * sigma)
    alpha2 = 1 / (sigma * sigma) + 1 / (k0 * k0)

    A1 = sigma
    a1 = (mp.pi * sigma) ** 2

    A2 = 1 / mp.sqrt(alpha2)
    a2 = (mp.pi ** 2) / alpha2

    return A1, a1, A2, a2


# -------- Closed-form bounds (no Simpson quadrature) --------
def L1_gpp_of_gaussian(a):
    """
    L1 norm of the second derivative for g(x) = exp(-a x^2), a > 0.

    Uses the exact piecewise integral of |(4 a^2 x^2 - 2 a) e^{-a x^2}|.
    """
    a = mp.mpf(a)
    x0 = mp.sqrt(1 / (2 * a))

    def I0(u):
        return mp.mpf("0.5") * mp.sqrt(mp.pi / a) * mp.erf(mp.sqrt(a) * u)

    def I2(u):
        return (mp.sqrt(mp.pi) * mp.erf(mp.sqrt(a) * u)) / (
            4 * a ** (mp.mpf("1.5"))
        ) - (u * mp.e ** (-a * u * u)) / (2 * a)

    # [0, x0]
    J1 = (2 * a) * I0(x0) - (4 * a * a) * I2(x0)

    # [x0, +inf)
    I0_inf = mp.mpf("0.5") * mp.sqrt(mp.pi / a)
    I2_inf = mp.sqrt(mp.pi) / (4 * a ** (mp.mpf("1.5")))
    J2 = (4 * a * a) * (I2_inf - I2(x0)) - (2 * a) * (I0_inf - I0(x0))

    # Even integrand
    return 2 * (J1 + J2)


def sup_h_hi_triangle(A1, a1, A2, a2):
    """
    Crude but safe bound: sup_x |h(x)| <= |A1| + |A2|.
    """
    return abs(A1) + abs(A2)


def sup_hpp_hi_triangle(A1, a1, A2, a2):
    """
    For g(x) = exp(-a x^2),
      g''(x) = (4 a^2 x^2 - 2 a) e^{-a x^2},
    and |g''(x)| has maximum 2 a at x = 0.
    """
    return 2 * abs(A1) * a1 + 2 * abs(A2) * a2


def sup_probe_h(A1, a1, A2, a2, xmax=mp.mpf("6.0"), points=2049):
    """
    Optional numeric probe of sup_x |h(x)| on [-xmax, xmax].
    """
    supv = mp.mpf("0")
    supx = mp.mpf("0")
    for i in range(points):
        x = -xmax + (2 * xmax) * i / (points - 1)
        v = abs(A1 * mp.e ** (-a1 * x * x) - A2 * mp.e ** (-a2 * x * x))
        if v > supv:
            supv, supx = v, x
    return supv, supx


# ------------- Main entrypoint -------------
def main():
    ap = argparse.ArgumentParser(
        description="Fourier inversion certificate for Gaussian notch kernel "
        "(closed-form, no quadrature)."
    )
    ap.add_argument(
        "--window-config",
        dest="window_config",
        help="Path to window JSON (from window_gen.py).",
    )
    ap.add_argument(
        "--sigma",
        type=str,
        help="Gaussian width (if not using --window-config).",
    )
    ap.add_argument(
        "--k0",
        type=str,
        help="Notch parameter (if not using --window-config).",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Path to output JSON.",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=220,
        help="Decimal precision.",
    )
    ap.add_argument(
        "--xmax",
        type=str,
        default="6.0",
        help="Half-width for sup probe (used only when --probe is set).",
    )
    ap.add_argument(
        "--simpson-n",
        type=int,
        default=0,
        help="Retained for compatibility; unused in closed-form mode.",
    )
    ap.add_argument(
        "--probe",
        action="store_true",
        help="Enable numeric probe for sup|h| on [-xmax, xmax].",
    )

    args = ap.parse_args()
    mp.dps = args.dps

    mode_tag, window_mode, sigma, k0, window_path = read_window(args)

    if window_mode != "gauss":
        print(
            "[fourier_inversion_cert] warning: window mode is not 'gauss'; "
            "using Gaussian closed-form formulas anyway.",
            file=sys.stderr,
        )

    # Closed-form parameters for h(x)
    A1, a1, A2, a2 = gauss_minus_gauss_closed_form(sigma, k0)

    # h(0) = A1 - A2 and epsilon_eff is stored for downstream rollups
    h0 = A1 - A2
    epsilon_eff = mp.mpf("0.5") * h0

    # L1 norm of h'' as a combination of Gaussian pieces
    L1_h2_hi = abs(A1) * L1_gpp_of_gaussian(a1) + abs(A2) * L1_gpp_of_gaussian(a2)

    # Triangle-type bounds for sup|h| and sup|h''|
    sup_h_tri = sup_h_hi_triangle(A1, a1, A2, a2)
    sup_h2_tri = sup_hpp_hi_triangle(A1, a1, A2, a2)

    if args.probe:
        sup_h_probe, where_h = sup_probe_h(
            A1, a1, A2, a2, xmax=mp.mpf(args.xmax)
        )
        sup_h = max(sup_h_tri, sup_h_probe)
        where_h2 = mp.mpf("0.0")
    else:
        sup_h = sup_h_tri
        where_h = mp.mpf("0.0")
        where_h2 = mp.mpf("0.0")

    inputs_block = {
        "mode": mode_tag,
        "sigma": str(sigma),
        "k0": str(k0),
    }
    if window_path is not None:
        inputs_block["window_config_path"] = window_path

    out = {
        "kind": "fourier_inversion_cert",
        "inputs": inputs_block,
        "closed_form_h": {
            "A1": str(A1),
            "a1": str(a1),
            "A2": str(A2),
            "a2": str(a2),
            "h0": str(h0),
            "epsilon_eff": str(epsilon_eff),
        },
        "bounds": {
            "L1_h2_hi": str(L1_h2_hi),
            "sup_h_hi": str(sup_h),
            "sup_h2_hi": str(sup_h2_tri),
            "where_sup_h": str(where_h),
            "where_sup_h2": str(where_h2),
            "xmax_used": str(mp.mpf(args.xmax)),
            "simpson_n": int(args.simpson_n),
        },
        "meta": {
            "tool": "fourier_inversion_cert",
            "dps": str(mp.dps),
            "created_utc": datetime.datetime.utcnow().isoformat() + "Z",
            "fourier_convention": (
                "F1: f^(xi)=∫ f(x) e^{-2π i x xi} dx; inverse uses +2π"
            ),
        },
    }

    jdump(out, args.out)

    print(
        "[fourier_inversion_cert] wrote {}  eps_eff={}  A1={}  A2={}".format(
            args.out, epsilon_eff, A1, A2
        )
    )
    print(
        "[fourier_inversion_cert] L1_h2_hi={}  sup|h|<={}  sup|h''|<={}".format(
            L1_h2_hi, sup_h, sup_h2_tri
        )
    )


if __name__ == "__main__":
    main()
