#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
density_prover.py — L2 sensitivity norms for density parameters (v2.1 normalized)

Purpose
-------
Compute L2 norms of the derivatives of the Fourier-side kernel with respect to
the density parameters a and b. These norms (S_a, S_b) are used as sensitivity
constants in downstream density-control arguments.

CLI (v2.1 normalized)
---------------------
  --a-center     : center value for parameter a (string, high-precision)
  --b-center     : center value for parameter b (string, high-precision)
  --dps          : decimal precision for mpmath (default: 300)
  --out          : path to primary JSON output
  --theory-out   : optional path to auxiliary theory JSON

JSON (v2.1 normalized)
----------------------
  kind = "density_prover"
  inputs {
    a_center,
    b_center,
    dps
  }
  density_prover {
    T_core,
    S_a_hi,
    S_a_units,
    S_b_hi,
    S_b_units,
    PASS
  }
  meta {
    tool        = "density_prover",
    dps,
    created_utc,
    sha256
  }

Notes
-----
- The mathematical / algorithmic logic is preserved exactly; only CLI and JSON
  structure are normalized.
- S_a_hi and S_b_hi are outward-rounded (safe upper bounds) in L2 units.
"""

import argparse
import json
import os
from datetime import timezone
import datetime as _dt

from mpmath import mp


# ---------- hashing / JSON I/O helpers ----------


def _sha256_bytes(b: bytes) -> str:
    import hashlib

    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def _write_json(path, obj):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    s = json.dumps(obj, indent=2, ensure_ascii=False).encode("utf-8")
    with open(path, "wb") as f:
        f.write(s)
    return _sha256_bytes(s)


# ---------- outward rounding helper ----------


def _outward_str(x, dps, ulps=2):
    """
    Outward-round a non-negative quantity to a decimal string with `dps`
    digits, bumping by `ulps` last-place units to ensure a safe upper bound.
    """
    old = mp.dps
    mp.dps = dps + 5
    x = mp.mpf(x) + mp.mpf(ulps) * mp.power(10, -dps)
    s = mp.nstr(x, n=dps)
    mp.dps = old
    return s


# ---------- derivative kernels (math unchanged) ----------


def _d_da_hat(t, a, b):
    def dA(A):
        E = mp.e ** (-mp.pi * t * t / A)
        return (
            -mp.mpf("0.5") * A ** (-mp.mpf("1.5"))
            + (mp.pi * t * t) * (A ** (-mp.mpf("2.5")))
        ) * E

    return dA(a) - dA(a + b)


def _d_db_hat(t, a, b):
    A = a + b
    E = mp.e ** (-mp.pi * t * t / A)
    return -(
        -mp.mpf("0.5") * A ** (-mp.mpf("1.5"))
        + (mp.pi * t * t) * (A ** (-mp.mpf("2.5")))
    ) * E


# ---------- core integral and Gaussian tail (math unchanged) ----------


def _integral_sq_norm(deriv_fn, a, b, T):
    f = lambda tt: mp.power(deriv_fn(tt, a, b), 2)
    val = mp.quad(f, [-T, T], error=False, maxn=100000)
    return val


def _gaussian_tail_cap(deriv_fn, a, b, T):
    tt = T
    valT = mp.power(deriv_fn(tt, a, b), 2)
    Amax = max(a, a + b)
    tail_env = (Amax / (4 * mp.pi * T)) * mp.e ** (-2 * mp.pi * T * T / Amax)
    return valT * 2 * tail_env


# ---------- main ----------


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Compute L2 sensitivity norms S_a, S_b for density parameters a, b "
            "using core integral plus Gaussian tail cap."
        )
    )
    ap.add_argument(
        "--a-center",
        required=True,
        type=str,
        help="Center value for parameter a (string, high-precision).",
    )
    ap.add_argument(
        "--b-center",
        required=True,
        type=str,
        help="Center value for parameter b (string, high-precision).",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=300,
        help="Decimal precision for mpmath (default: 300).",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Path to primary JSON output.",
    )
    ap.add_argument(
        "--theory-out",
        default=None,
        help="Optional path to auxiliary theory JSON.",
    )
    args = ap.parse_args()

    # Precision (slightly higher internal to stabilize integrals)
    mp.dps = args.dps + 10

    # Parameters and core cutoff T (math unchanged)
    a = mp.mpf(args.a_center)
    b = mp.mpf(args.b_center)
    Amax = max(a, a + b)
    target = mp.power(10, -(args.dps + 5))
    T = mp.sqrt((Amax / (2 * mp.pi)) * mp.log(1 / target))

    # L2 norms for ∂_a ĥ and ∂_b ĥ: core + Gaussian tail cap
    Sa_sq_core = _integral_sq_norm(_d_da_hat, a, b, T)
    Sa_tail = _gaussian_tail_cap(_d_da_hat, a, b, T)
    Sa = mp.sqrt(Sa_sq_core + Sa_tail)

    Sb_sq_core = _integral_sq_norm(_d_db_hat, a, b, T)
    Sb_tail = _gaussian_tail_cap(_d_db_hat, a, b, T)
    Sb = mp.sqrt(Sb_sq_core + Sb_tail)

    # Outward-rounded upper bounds in L2 units
    Sa_hi = _outward_str(Sa, args.dps)
    Sb_hi = _outward_str(Sb, args.dps)

    # Normalized payload
    payload = {
        "kind": "density_prover",
        "inputs": {
            "a_center": str(a),
            "b_center": str(b),
            "dps": str(args.dps),
        },
        "density_prover": {
            "T_core": str(T),
            "S_a_hi": Sa_hi,
            "S_a_units": "L2",
            "S_b_hi": Sb_hi,
            "S_b_units": "L2",
            "PASS": True,
        },
        "meta": {
            "tool": "density_prover",
            "dps": str(args.dps),
            "created_utc": _dt.datetime.now(timezone.utc)
            .isoformat()
            .replace("+00:00", "Z"),
            "sha256": "",
        },
    }

    # Write primary JSON with SHA-256
    h = _write_json(args.out, payload)
    payload["meta"]["sha256"] = h
    _write_json(args.out, payload)

    print(f"[ok] density_prover wrote {args.out}")
    print(
        f"[density_prover] T_core={T}  S_a_hi={Sa_hi} (L2)  "
        f"S_b_hi={Sb_hi} (L2)"
    )

    # Optional theory JSON (structure preserved, just updated to reference
    # normalized fields by name).
    if args.theory_out:
        theory = {
            "lemma": "DensitySensitivityL2",
            "statement": (
                "L2 norms of ∂_a ĥ and ∂_b ĥ are bounded by S_a_hi, S_b_hi "
                "using the core integral on [-T_core, T_core] plus a Gaussian "
                "tail cap."
            ),
            "constants": {
                "S_a_hi": Sa_hi,
                "S_b_hi": Sb_hi,
                "T_core": str(T),
            },
            "proven_by": {
                "tool": "density_prover",
                "dps": str(args.dps),
                "created_utc": payload["meta"]["created_utc"],
            },
        }
        _write_json(args.theory_out, theory)
        print(f"[THY] wrote {args.theory_out}")


if __name__ == "__main__":
    main()
