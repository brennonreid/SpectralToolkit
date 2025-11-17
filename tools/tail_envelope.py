#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tail_envelope.py — analytic gamma tail envelope at T0.

Purpose:
  Analytic envelope for the gamma tail beyond T0, of the form
    G(T) <= c1(σ)/T + c2(σ)/T^2   for T >= T0,
  with parameters
    c1 = 1 + 1/σ,
    c2 = 1/2 + 1/σ^2.

CLI (normalized v2.1):
  --sigma       : Gaussian width sigma (>0)
  --T0          : base point T0 (>0)
  --dps         : decimal precision for mpmath
  --out         : output JSON path
  --theory-out  : optional theory JSON path

JSON (normalized v2.1):
  {
    "kind": "gamma_tail_envelope",
    "inputs": {
      "sigma": "<string>",
      "T0": "<string>"
    },
    "gamma_tail": {
      "c1": "<string>",
      "c2": "<string>",
      "env_at_T0": "<string>",
      "monotone": true
    },
    "meta": {
      "tool": "tail_envelope",
      "dps": <int>,
      "created_utc": "<ISO8601 UTC>",
      "sha256": "<digest>"
    }
  }

All real-valued numeric fields are stored as strings. Integers remain ints.
"""

import argparse
import json
import os
import hashlib
from datetime import datetime, timezone

from mpmath import mp


def _now_utc_iso() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def _sha256_bytes(b: bytes) -> str:
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def _write_json(path: str, obj: dict) -> None:
    """
    Write JSON with deterministic sha256 in meta.sha256.
    Hash is computed over the payload with any existing meta.sha256 removed.
    """
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)

    if "meta" not in obj or not isinstance(obj["meta"], dict):
        obj["meta"] = {}

    # Compute hash on a copy without meta.sha256 (avoid self-reference).
    tmp_obj = json.loads(json.dumps(obj, ensure_ascii=False))
    meta = tmp_obj.get("meta")
    if isinstance(meta, dict):
        meta.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    digest = _sha256_bytes(blob)
    obj["meta"]["sha256"] = digest

    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(
        description="Analytic gamma tail envelope G(T) <= c1(σ)/T + c2(σ)/T^2 for T >= T0."
    )
    ap.add_argument(
        "--sigma",
        default="6.0",
        help="Gaussian width sigma (>0).",
    )
    ap.add_argument(
        "--T0",
        default="1000000",
        help="Base point T0 (>0).",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=300,
        help="Decimal precision for mpmath.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path (gamma_tail_envelope).",
    )
    ap.add_argument(
        "--theory-out",
        default=None,
        help="Optional theory JSON path.",
    )
    args = ap.parse_args()

    mp.dps = int(args.dps)
    sigma = mp.mpf(args.sigma)
    T0 = mp.mpf(args.T0)

    if sigma <= 0:
        raise SystemExit("sigma must be > 0")
    if T0 <= 0:
        raise SystemExit("T0 must be > 0")

    # Math: unchanged from original version
    c1 = mp.mpf("1.0") + 1 / sigma
    c2 = mp.mpf("0.5") + 1 / (sigma * sigma)
    env_T0 = c1 / T0 + c2 / (T0 * T0)

    # Main payload: legacy-compatible and v2.1-friendly.
    payload = {
        "kind": "gamma_tail_envelope",
        "inputs": {
            "sigma": _mp_str(sigma),
            "T0": _mp_str(T0),
        },
        "gamma_tail": {
            "c1": _mp_str(c1),
            "c2": _mp_str(c2),
            # legacy name used by older readers
            "env_at_T0": _mp_str(env_T0),
            # new-style name expected by some aggregators
            "gamma_env_at_T0": _mp_str(env_T0),
            # make T0 visible next to the env value as well
            "T0": _mp_str(T0),
            "monotone": True,
        },
        # optional convenience mirror for tools that look for gamma_tails.*
        "gamma_tails": {
            "T0": _mp_str(T0),
            "gamma_env_at_T0": _mp_str(env_T0),
        },
        "meta": {
            "tool": "tail_envelope",
            "dps": int(mp.dps),
            "created_utc": _now_utc_iso(),
        },
    }


    _write_json(args.out, payload)
    print(f"[ok] tail_envelope -> {args.out}  env_at_T0={_mp_str(env_T0)}")

    if args.theory_out:
        theory = {
            "kind": "gamma_tail_envelope_theory",
            "theory": {
                "lemma": "GammaTailEnvelope",
                "statement": (
                    "For T >= T0, G(T) <= c1(σ)/T + c2(σ)/T^2 with "
                    "c1 = 1 + 1/σ, c2 = 1/2 + 1/σ^2; the envelope is decreasing in T."
                ),
            },
            "inputs": {
                "sigma": _mp_str(sigma),
                "T0": _mp_str(T0),
            },
            "meta": {
                "tool": "tail_envelope_theory",
                "dps": int(mp.dps),
                "created_utc": _now_utc_iso(),
            },
        }
        _write_json(args.theory_out, theory)
        print(f"[ok] tail_envelope theory -> {args.theory_out}")


if __name__ == "__main__":
    main()
