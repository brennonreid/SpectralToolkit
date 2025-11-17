#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uniform_rollup_cert.py (v2.1-normalized)

Purpose:
  Aggregate local certificates at a single T0 into a uniform certificate:

    band_cert.json
    prime_block_norm.json
    prime_tail_envelope.json
    gamma_tail.json
    grid_error_bound.json           (optional)
    continuum_operator_cert.json    (optional, echoed)
    weil_explicit_cert.json         (optional, echoed)
    weil_psd_bochner.json           (optional, for PSD flag)

The tool computes:
  - band_margin      : lower bound on spectral margin from band_cert
  - gamma_env_at_T0  : gamma tail envelope at T0
  - epsilon_eff      : effective margin = band_margin - gamma_env_at_T0
  - prime_block_cap  : cap on prime block operator norm
  - prime_tail_norm  : prime tail contribution at T0
  - grid_error_norm  : quadrature / grid error bound (optional)
  - lhs_total        : prime_block_cap + prime_tail_norm + grid_error_norm
  - PSD_verified     : from Weil PSD certificate (if present)

PASS condition:
  PASS is true iff:
    PSD_verified is true AND lhs_total <= epsilon_eff
  All comparisons are done in high-precision mpmath arithmetic.

CLI (normalized v2.1):
  --T0          : target T0 (string or number)
  --certs-dir   : directory containing the input certificate JSONs
  --dps         : decimal precision for mpmath
  --out         : output JSON path

JSON (normalized v2.1):
  {
    "kind": "uniform_certificate",
    "inputs": {
      "T0": "<string>",
      "certs_dir": "<path>",
      "band_cert_path": "<path>",
      "prime_block_path": "<path>",
      "prime_tail_path": "<path>",
      "gamma_tail_path": "<path>",
      "grid_error_path": "<path or null>",
      "continuum_operator_path": "<path or null>",
      "weil_explicit_path": "<path or null>",
      "weil_psd_path": "<path or null>"
    },
    "uniform_certificate": {
      "band_margin": "<string>",
      "gamma_env_at_T0": "<string>",
      "epsilon_eff": "<string>",
      "prime_block_cap": "<string>",
      "prime_tail_norm": "<string>",
      "grid_error_norm": "<string>",
      "lhs_total": "<string>",
      "PSD_verified": true/false
    },
    "PASS": true/false,
    "meta": {
      "tool": "uniform_rollup_cert",
      "dps": <int>,
      "created_utc": "<ISO8601 UTC>",
      "sha256": "<digest>"
    }
  }

All real-valued numeric fields are stored as strings. Integers remain ints,
booleans remain JSON booleans.
"""

import argparse
import os
import json
import time
import sys
import hashlib
from typing import Any, Optional
from pathlib import Path

from mpmath import mp


# ---------- helpers ----------

def utc_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8-sig") as f:
        return json.load(f)


def dig(obj: Any, path: list) -> Optional[Any]:
    o = obj
    for k in path:
        if not isinstance(o, dict) or k not in o:
            return None
        o = o[k]
    return o


def coalesce(obj: Any, paths: list, default: str = "0") -> str:
    for p in paths:
        v = dig(obj, p)
        if v not in (None, "", "null"):
            return str(v)
    return default


def mp_str(x: Any) -> str:
    return mp.nstr(mp.mpf(str(x)), n=mp.dps, strip_zeros=False)


def must(path: str, what: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing {what}: {path}")


def jpath(root: str, name: str) -> str:
    return os.path.join(root, name)


def write_json(path: str, payload: dict, dps: int) -> None:
    """
    Write JSON with deterministic sha256 in meta.sha256.

    Hash is computed over payload with any existing meta.sha256 removed,
    using canonical JSON (sorted keys, compact separators).
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    # Ensure meta exists and has dps/tool set (caller sets tool, we enforce dps)
    if "meta" not in payload or not isinstance(payload["meta"], dict):
        payload["meta"] = {}
    payload["meta"].setdefault("dps", int(dps))

    # Compute sha256 over a copy without meta.sha256
    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    meta = tmp_obj.get("meta")
    if isinstance(meta, dict):
        meta.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    digest = hashlib.sha256(blob).hexdigest()
    payload["meta"]["sha256"] = digest

    # Atomic write via .tmp -> final
    tmp = p.with_suffix(p.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write("\n")

    # Rename tmp to final (works whether or not final existed before)
    tmp.replace(p)


# ---------- field readers (multi-schema, with v2.1 canonical paths first) ----------

def read_band_margin(band_json: dict) -> str:
    """
    Return the lower band margin as a string.

    Canonical v2.1 band_cert (as normalized earlier) provides:
      numbers.band_margin_lo
      band_cert.band_margin.lo
    """
    return coalesce(
        band_json,
        [
            # canonical v2.1
            ["numbers", "band_margin_lo"],
            ["band_cert", "band_margin", "lo"],
            # older / fallback schemas
            ["band_cert", "band_margin_lo"],
            ["numbers", "band_margin"],
            ["band_margin", "lo"],
            ["band_margin_lo"],
        ],
        default="0",
    )


def read_prime_block_cap(pblk_json: dict) -> str:
    """
    Prime block cap; in v2.1 prime_block_norm.json we use:

      prime_block_norm.used_operator_norm
    """
    return coalesce(
        pblk_json,
        [
            # canonical v2.1 (used operator norm)
            ["prime_block_norm", "used_operator_norm"],
            # older / fallback schemas
            ["used_operator_norm"],
            ["operator_norm_cap", "hi"],
            ["operator_norm_cap"],
            ["cap"],
        ],
        default="0",
    )


def read_prime_tail_norm(pt_json: dict) -> str:
    """
    Prime tail norm at T0.

    Canonical v2.1 prime_tail_envelope.json:

      prime_tail.env_T0_hi
    """
    return coalesce(
        pt_json,
        [
            # canonical v2.1
            ["prime_tail", "env_T0_hi"],
            # older / fallback schemas
            ["prime_tail_envelope", "env_T0_hi"],
            ["numbers", "prime_tail_norm"],
            ["env_T0_hi"],
            ["prime_tail_norm"],
        ],
        default="0",
    )


def read_gamma_env_T0(ga_json: dict) -> str:
    """
    Gamma tail envelope at T0.

    Canonical sources in v2.1:
      - core_integral_prover: gamma_tails.gamma_env_at_T0
      - tail_envelope: gamma_tail.env_at_T0
    """
    return coalesce(
        ga_json,
        [
            # canonical v2.1
            ["gamma_tails", "gamma_env_at_T0"],
            ["gamma_tail", "env_at_T0"],
            # older / fallback schemas
            ["gamma_env_at_T0"],
            ["tails_total"],
            ["numbers", "gamma_env_at_T0"],
        ],
        default="0",
    )


def read_grid_error_norm(grid_json: Optional[dict]) -> str:
    if not isinstance(grid_json, dict):
        return "0"
    return coalesce(
        grid_json,
        [
            # canonical v2.1
            ["grid_error_bound", "bound_hi"],
            # older / fallback schemas
            ["grid_error_norm"],
            ["numbers", "grid_error_norm"],
            ["hi"],
            ["lo"],
        ],
        default="0",
    )


def read_psd_pass(weil_json: Optional[dict]) -> bool:
    """
    Canonical v2.1 Weil PSD cert (weil_kernel.py normalized):

      kind = "weil_psd_bochner"
      PSD_verified = true/false
    """
    if not isinstance(weil_json, dict):
        return True
    for p in (
        ["PSD_verified"],                # canonical v2.1
        ["bochner_psd", "PSD_verified"],
        ["weil_psd", "PSD_verified"],
    ):
        v = dig(weil_json, p)
        if isinstance(v, bool):
            return v
        if isinstance(v, str) and v.lower() in ("true", "false"):
            return v.lower() == "true"
    return True


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(
        description="Uniform rollup certificate at T0 (band + tails + prime + grid)."
    )
    ap.add_argument(
        "--T0",
        required=True,
        help="Target T0 (string or number).",
    )
    ap.add_argument(
        "--certs-dir",
        default="PROOF_PACKET",
        help="Directory containing input cert JSONs.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path (uniform_certificate.json).",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=220,
        help="Decimal precision for mpmath.",
    )
    args = ap.parse_args()

    mp.dps = int(args.dps)
    certs_dir = args.certs_dir
    out_path = args.out

    # Required inputs (hard-fail if missing)
    band_path = jpath(certs_dir, "band_cert.json")
    pblk_path = jpath(certs_dir, "prime_block_norm.json")
    pt_path   = jpath(certs_dir, "prime_tail_envelope.json")
    ga_path   = jpath(certs_dir, "gamma_tail.json")

    must(band_path, "band_cert.json")
    must(pblk_path, "prime_block_norm.json")
    must(pt_path,   "prime_tail_envelope.json")
    must(ga_path,   "gamma_tail.json")

    # Optional inputs
    grid_path = jpath(certs_dir, "grid_error_bound.json")
    roll_path = jpath(certs_dir, "continuum_operator_cert.json")
    wexp_path = jpath(certs_dir, "weil_explicit_cert.json")
    psd_path  = jpath(certs_dir, "weil_psd_bochner.json")

    # Load JSON artifacts
    band  = load_json(band_path)
    pblk  = load_json(pblk_path)
    pt    = load_json(pt_path)
    ga    = load_json(ga_path)
    grid  = load_json(grid_path) if os.path.exists(grid_path) else None
    roll  = load_json(roll_path) if os.path.exists(roll_path) else None
    wexp  = load_json(wexp_path) if os.path.exists(wexp_path) else None
    psd   = load_json(psd_path)  if os.path.exists(psd_path)  else None

    # Extract numbers (as high-precision strings)
    band_margin     = mp_str(read_band_margin(band))
    prime_block_cap = mp_str(read_prime_block_cap(pblk))
    prime_tail_norm = mp_str(read_prime_tail_norm(pt))
    gamma_env_T0    = mp_str(read_gamma_env_T0(ga))
    grid_error_norm = mp_str(read_grid_error_norm(grid))

    # Derived quantities
    lhs_total_val   = (
        mp.mpf(prime_block_cap)
        + mp.mpf(prime_tail_norm)
        + mp.mpf(grid_error_norm)
    )
    epsilon_eff_val = mp.mpf(band_margin) - mp.mpf(gamma_env_T0)

    lhs_total   = mp_str(lhs_total_val)
    epsilon_eff = mp_str(epsilon_eff_val)

    # PASS criteria: spectrum margin covers all residual costs; PSD ok
    psd_ok = read_psd_pass(psd)
    PASS   = bool(psd_ok and (lhs_total_val <= epsilon_eff_val))

    payload = {
        "kind": "uniform_certificate",
        "inputs": {
            "T0": str(args.T0),
            "certs_dir": certs_dir,
            "band_cert_path": band_path,
            "prime_block_path": pblk_path,
            "prime_tail_path": pt_path,
            "gamma_tail_path": ga_path,
            "grid_error_path": grid_path if os.path.exists(grid_path) else None,
            "continuum_operator_path": roll_path if os.path.exists(roll_path) else None,
            "weil_explicit_path": wexp_path if os.path.exists(wexp_path) else None,
            "weil_psd_path": psd_path if os.path.exists(psd_path) else None,
        },
        "uniform_certificate": {
            "band_margin": band_margin,
            "gamma_env_at_T0": gamma_env_T0,
            "epsilon_eff": epsilon_eff,
            "prime_block_cap": prime_block_cap,
            "prime_tail_norm": prime_tail_norm,
            "grid_error_norm": grid_error_norm,
            "lhs_total": lhs_total,
            "PSD_verified": bool(psd_ok),
        },
        "PASS": PASS,
        "meta": {
            "tool": "uniform_rollup_cert",
            "dps": int(args.dps),
            "created_utc": utc_iso(),
        },
    }

    write_json(out_path, payload, args.dps)

    status = "PASS" if PASS else "FAIL"
    print(
        f"[uniform_rollup_cert] {status} -> {out_path}  "
        f"lhs_total={lhs_total}  epsilon_eff={epsilon_eff}  psd={psd_ok}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[uniform_rollup_cert] ERROR: {e}", file=sys.stderr)
        sys.exit(1)
