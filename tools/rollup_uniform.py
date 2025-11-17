#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rollup_uniform.py — analytic uniform roll-up at T0.

Purpose:
  Compare the analytic envelopes

      B(T0) = gamma_envelope(T0) + prime_envelope(T0)
      E(T0) = epsilon_eff_lo(T0)

  and certify a uniform inequality on [T0, ∞) under the assumption that the
  input envelopes are monotone in T. This module only compares values at T0;
  monotonicity is provided by the input envelope constructions.

CLI (v2.1 normalized):
  --T0               : base point T0 (string or number)
  --gamma-envelope   : path to gamma envelope JSON (e.g. tail_envelope.py)
  --prime-envelope   : path to prime envelope JSON (e.g. prime_tail_envelope.py)
  --explicit-formula : path to explicit_formula.json
  --dps              : decimal precision for mpmath
  --out              : output JSON path
  --theory-out       : optional theory JSON path

JSON (v2.1 normalized):
  {
    "kind": "uniform_analytic_rollup",
    "inputs": {
      "T0": "<string>",
      "gamma_envelope_path": "<path>",
      "prime_envelope_path": "<path>",
      "explicit_formula_path": "<path>"
    },
    "uniform_certificate": {
      "E_T0_lo": "<string>",
      "B_T0_hi": "<string>",
      "PASS": true/false
    },
    "PASS": true/false,
    "meta": {
      "tool": "rollup_uniform",
      "dps": <int>,
      "created_utc": "<ISO8601 UTC>",
      "sha256": "<digest>"
    }
  }

All real-valued numeric fields are stored as strings.
Integers remain JSON integers; booleans remain JSON booleans.
"""

import argparse
import json
import os
import hashlib
from pathlib import Path
from datetime import datetime, timezone

from mpmath import mp


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def set_precision(dps: int) -> None:
    mp.dps = int(dps)


def mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def now_utc_iso() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def write_json(path: str, payload: dict) -> None:
    """
    Write JSON with deterministic sha256 in meta.sha256.

    Hash is computed over the payload with any existing meta.sha256 removed,
    using canonical JSON (sorted keys, compact separators).
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    if "meta" not in payload or not isinstance(payload["meta"], dict):
        payload["meta"] = {}

    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    meta = tmp_obj.get("meta")
    if isinstance(meta, dict):
        meta.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    digest = hashlib.sha256(blob).hexdigest()
    payload["meta"]["sha256"] = digest

    tmp = str(p) + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, p)


def load_json(path: str) -> dict:
    with open(path, "r", encoding="utf-8-sig") as f:
        return json.load(f)


# ---------------------------------------------------------------------
# Value readers (canonical v2.1 shapes)
# ---------------------------------------------------------------------

def read_gamma_env_T0(js: dict) -> mp.mpf:
    """
    Gamma envelope at T0.

    Canonical v2.1 tail_envelope.py:
      kind = "gamma_tail_envelope"
      gamma_tail.env_at_T0  (string)
    """
    try:
        val = js["gamma_tail"]["env_at_T0"]
    except Exception as e:
        raise KeyError("gamma_envelope JSON must contain gamma_tail.env_at_T0") from e
    return mp.mpf(str(val))


def read_prime_env_T0(js: dict) -> mp.mpf:
    """
    Prime envelope at T0.

    Canonical v2.1 prime_tail_envelope.py:
      kind = "prime_tail_envelope"
      prime_tail.env_T0_hi  (string)
    """
    try:
        val = js["prime_tail"]["env_T0_hi"]
    except Exception as e:
        raise KeyError("prime_envelope JSON must contain prime_tail.env_T0_hi") from e
    return mp.mpf(str(val))


def read_epsilon_eff_lo(js: dict) -> mp.mpf:
    """
    Epsilon_eff lower bound at T0 from explicit_formula.

    Canonical v2.1 explicit_formula.py:
      kind = "weil_explicit"
      explicit_formula.epsilon_eff_lo  (string)
    """
    try:
        val = js["explicit_formula"]["epsilon_eff_lo"]
    except Exception as e:
        raise KeyError(
            "explicit_formula JSON must contain explicit_formula.epsilon_eff_lo"
        ) from e
    return mp.mpf(str(val))


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Analytic uniform roll-up at T0: compare B(T0) and E(T0)."
    )
    ap.add_argument(
        "--T0",
        required=True,
        help="Base point T0 (string or number).",
    )
    ap.add_argument(
        "--gamma-envelope",
        required=True,
        help="Gamma envelope JSON path (e.g. gamma_tail_envelope).",
    )
    ap.add_argument(
        "--prime-envelope",
        required=True,
        help="Prime envelope JSON path (e.g. prime_tail_envelope).",
    )
    ap.add_argument(
        "--explicit-formula",
        required=True,
        help="Explicit formula JSON path (explicit_formula.py output).",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=220,
        help="Decimal precision for mpmath.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path (uniform_analytic_rollup).",
    )
    ap.add_argument(
        "--theory-out",
        default=None,
        help="Optional theory JSON path.",
    )
    args = ap.parse_args()

    set_precision(args.dps)

    gamma_js = load_json(args.gamma_envelope)
    prime_js = load_json(args.prime_envelope)
    expl_js = load_json(args.explicit_formula)

    G = read_gamma_env_T0(gamma_js)   # gamma envelope at T0
    P = read_prime_env_T0(prime_js)   # prime envelope at T0
    E = read_epsilon_eff_lo(expl_js)  # epsilon_eff lower bound at T0

    B_hi = G + P
    E_lo = E

    PASS = bool(E_lo > B_hi)

    payload = {
        "kind": "uniform_analytic_rollup",
        "inputs": {
            "T0": str(args.T0),
            "gamma_envelope_path": args.gamma_envelope,
            "prime_envelope_path": args.prime_envelope,
            "explicit_formula_path": args.explicit_formula,
        },
        "uniform_certificate": {
            "E_T0_lo": mp_str(E_lo),
            "B_T0_hi": mp_str(B_hi),
            "PASS": PASS,
        },
        "PASS": PASS,
        "meta": {
            "tool": "rollup_uniform",
            "dps": int(args.dps),
            "created_utc": now_utc_iso(),
        },
    }

    write_json(args.out, payload)
    print(
        f"[rollup_uniform] {'PASS' if PASS else 'FAIL'} -> {args.out}  "
        f"E_T0_lo={mp_str(E_lo)}  B_T0_hi={mp_str(B_hi)}"
    )

    if args.theory_out:
        theory = {
            "kind": "uniform_analytic_rollup_theory",
            "inputs": {
                "T0": str(args.T0),
            },
            "theory": {
                "lemma": "UniformAnalyticRollup",
                "statement": (
                    "If B(T0) < E(T0) and the gamma/prime envelopes B(T) and "
                    "epsilon envelope E(T) are monotone on [T0, ∞), then "
                    "for all T ≥ T0 one has B(T) < E(T)."
                ),
                "notes": (
                    "This module checks the inequality only at T0; monotonicity "
                    "assumptions are provided by the envelope constructions."
                ),
            },
            "meta": {
                "tool": "rollup_uniform",
                "dps": int(args.dps),
                "created_utc": now_utc_iso(),
            },
        }
        write_json(args.theory_out, theory)
        print(f"[rollup_uniform] theory -> {args.theory_out}")


if __name__ == "__main__":
    main()
