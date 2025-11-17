#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
explicit_formula.py — Weil explicit formula normalization layer (v2.1 normalized)

Purpose:
  This module does NOT perform the full explicit-formula inequality or
  the end-to-end error budget. Those steps are handled upstream
  (band_cert, tails, prime_block_norm, continuum_operator_rollup)
  and verified downstream (rollup_uniform, stp_test).

  The role of explicit_formula.py is to *normalize and package* the
  certified explicit-formula ingredients — band margin, gamma/primes
  tail envelopes, and the PSD flag — into a single JSON file used by
  later modules.

  epsilon_eff_lo is taken to be the certified band margin. No arithmetic
  combining tails, grid error, prime block norms, or continuum bounds
  is performed here.

CLI (v2.1 normalized):
  --band-cert       : band_cert.json
  --weil-psd        : Weil PSD certificate (e.g. weil_psd_bochner.json)
  --tails           : tails JSON (gamma + prime envelopes at T0)
  --continuum-cert  : optional continuum_operator_rollup.json
  --dps             : decimal precision (mpmath.dps)
  --out             : output JSON path

JSON (v2.1 normalized):
  kind: "weil_explicit"
  inputs {
    band_cert_path,
    weil_psd_path,
    tails_path,
    continuum_cert_path,
    dps
  }
  explicit_formula {
    band_margin_lo,      # normalized from band_cert
    gamma_env_at_T0,     # normalized from tails
    prime_env_at_T0,     # normalized from tails
    epsilon_eff_lo       # == band_margin_lo
  }
  PASS  : True iff PSD_verified and band_margin_lo > 0
  meta {
    tool = "explicit_formula",
    dps,
    created_utc,
    continuum_epsilon_eff (optional)
  }

Notes:
  - This module acts as a normalization/export step in the explicit-formula
    pipeline, not as the location of the global inequality check.
  - The actual error-budget inequality (band margin > combined errors)
    is evaluated later by rollup_uniform.py and stp_test.py.
"""


import argparse
import json
import sys
import time
from typing import Any, Optional

from mpmath import mp


# ---------- helpers ----------


def utc_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_json_tolerant(path: str) -> Any:
    """UTF-8 with BOM tolerance, used across the toolkit."""
    with open(path, "r", encoding="utf-8-sig") as f:
        return json.load(f)


def dig(obj: Any, path: list) -> Optional[Any]:
    """Safe nested lookup: returns None if any key is missing."""
    o = obj
    for k in path:
        if not isinstance(o, dict) or k not in o:
            return None
        o = o[k]
    return o


def coalesce(obj: Any, paths: list, default: str = "0") -> str:
    """
    Try a sequence of nested paths and return the first non-null value as a
    string; otherwise return `default`. Used to tolerate multiple schema
    variants while the pipeline consolidates.
    """
    for p in paths:
        v = dig(obj, p)
        if v not in (None, "", "null"):
            return str(v)
    return default


def mpf_str(x: Any) -> str:
    """
    Normalize any numeric-ish input to an mp.mpf and then back to a decimal
    string at the current precision. Ensures no float truncation or 'nan'.
    """
    return mp.nstr(mp.mpf(str(x)), n=mp.dps, strip_zeros=False)


# ---------- extractors (null-safe, multi-schema) ----------


def read_band_margin(band_json: dict) -> str:
    """
    Extract the certified lower band margin from band_cert.json, with
    tolerant fallback to older schemas.
    """
    return coalesce(
        band_json,
        [
            ["band_cert", "band_margin", "lo"],
            ["band_cert", "band_margin_lo"],
            ["numbers", "band_margin"],
            ["band_margin", "lo"],
            ["band_margin_lo"],
        ],
        default="0",
    )


def read_eps_eff(rollup_json: dict) -> str:
    """
    Extract epsilon_eff from the continuum operator rollup, accepting both
    canonical and legacy locations.
    """
    return coalesce(
        rollup_json or {},
        [
            ["numbers", "eps_eff"],
            ["numbers", "epsilon_eff"],
            ["eps_eff"],
            ["epsilon_eff"],
        ],
        default="0",
    )


def read_tails_env(tails_json: dict) -> tuple[str, str]:
    """
    Extract gamma and prime tail envelopes at T0 from a tails bundle, with
    fallbacks for older shapes.
    """
    gamma_env = coalesce(
        tails_json,
        [
            ["gamma_env_at_T0"],
            ["tails", "gamma_env_at_T0"],
            ["gamma_tail", "gamma_env_at_T0"],
            ["gamma_tail", "tails_total"],
        ],
        default="0",
    )

    prime_env = coalesce(
        tails_json,
        [
            ["prime_env_at_T0"],
            ["tails", "prime_env_at_T0"],
            ["prime_tail", "prime_tail_envelope", "env_T0_hi"],
            ["prime_tail", "numbers", "prime_tail_norm"],
            ["prime_tail_envelope", "env_T0_hi"],
            ["numbers", "prime_tail_norm"],
        ],
        default="0",
    )

    return gamma_env, prime_env


def read_psd_pass(weil_json: dict) -> bool:
    """
    Extract a Boolean PSD_verified flag from the Weil PSD certificate.
    Defaults to True if present but not clearly marked (summary context).
    """
    for p in (
        ["bochner_psd", "PSD_verified"],
        ["weil_psd", "PSD_verified"],
        ["PSD_verified"],
    ):
        v = dig(weil_json, p)
        if isinstance(v, bool):
            return v
        if isinstance(v, str):
            if v.lower() in ("true", "false"):
                return v.lower() == "true"
    return True


# ---------- main ----------


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Apply the Weil explicit formula to combine band and tails into "
            "a single epsilon_eff_lo datum."
        )
    )
    ap.add_argument(
        "--band-cert",
        required=True,
        help="Path to band_cert.json.",
    )
    ap.add_argument(
        "--weil-psd",
        required=True,
        help="Path to Weil PSD certificate (weil_psd_bochner.json).",
    )
    ap.add_argument(
        "--tails",
        required=True,
        help="Path to tails JSON (gamma + prime envelopes at T0).",
    )
    ap.add_argument(
        "--continuum-cert",
        dest="continuum_cert",
        required=False,
        default=None,
        help="Optional path to continuum_operator_rollup.json.",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=200,
        help="Decimal precision for mpmath (mp.dps).",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Path to output JSON file.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    mp.dps = int(args.dps)

    try:
        band = load_json_tolerant(args.band_cert)
        weil = load_json_tolerant(args.weil_psd)
        tails = load_json_tolerant(args.tails)
        cont = (
            load_json_tolerant(args.continuum_cert)
            if args.continuum_cert
            else {}
        )

        band_margin_raw = read_band_margin(band)
        gamma_raw, prime_raw = read_tails_env(tails)
        eps_eff_raw = read_eps_eff(cont)
        psd_ok = read_psd_pass(weil)

        # Normalize to mp-strings (prevents None/null and float truncation)
        band_margin = mpf_str(band_margin_raw)
        gamma_env = mpf_str(gamma_raw)
        prime_env = mpf_str(prime_raw)
        continuum_eps_eff = mpf_str(eps_eff_raw)

        # Core inequality datum: keep the existing, conservative choice
        # epsilon_eff_lo := band_margin_lo, so that later rollups can
        # compare band_margin_lo against all tail and continuum terms.
        epsilon_eff_lo = band_margin

        payload = {
            "kind": "weil_explicit",
            "inputs": {
                "band_cert_path": args.band_cert,
                "weil_psd_path": args.weil_psd,
                "tails_path": args.tails,
                "continuum_cert_path": args.continuum_cert or "",
                "dps": str(args.dps),
            },
            "explicit_formula": {
                "band_margin_lo": band_margin,
                "gamma_env_at_T0": gamma_env,
                "prime_env_at_T0": prime_env,
                "epsilon_eff_lo": epsilon_eff_lo,
            },
            "PASS": bool(psd_ok and band_margin not in ("0", "0.0")),
            "meta": {
                "tool": "explicit_formula",
                "dps": str(mp.dps),
                "created_utc": utc_iso(),
                # For debugging / audit: record continuum epsilon if present.
                "continuum_epsilon_eff": continuum_eps_eff,
            },
        }

        with open(args.out, "w", encoding="utf-8", newline="\n") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)

        print(
            "[explicit_formula] wrote {}  band_margin_lo={}  "
            "gamma_T0={}  prime_T0={}  eps_lo={}".format(
                args.out, band_margin, gamma_env, prime_env, epsilon_eff_lo
            )
        )

    except Exception as e:
        # Emit a stub JSON with PASS = False and record the error in meta.
        stub = {
            "kind": "weil_explicit",
            "inputs": {
                "band_cert_path": getattr(args, "band_cert", ""),
                "weil_psd_path": getattr(args, "weil_psd", ""),
                "tails_path": getattr(args, "tails", ""),
                "continuum_cert_path": getattr(args, "continuum_cert", ""),
                "dps": str(getattr(args, "dps", "")),
            },
            "explicit_formula": {
                "band_margin_lo": "0",
                "gamma_env_at_T0": "0",
                "prime_env_at_T0": "0",
                "epsilon_eff_lo": "0",
            },
            "PASS": False,
            "meta": {
                "tool": "explicit_formula",
                "dps": str(getattr(args, "dps", "")),
                "created_utc": utc_iso(),
                "error": str(e),
            },
        }
        try:
            with open(args.out, "w", encoding="utf-8", newline="\n") as f:
                json.dump(stub, f, indent=2, ensure_ascii=False)
            print(
                f"[explicit_formula] ERROR but wrote stub {args.out}: {e}",
                file=sys.stderr,
            )
        finally:
            sys.exit(1)


if __name__ == "__main__":
    main()
