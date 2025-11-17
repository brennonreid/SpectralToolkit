#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lipschitz_q_bound.py — Lipschitz bound for q around (sigma, k0).

Purpose:
  Given T0, x0, A_prime, sigma_scale, k0_scale, and K, compute a global
  Lipschitz bound L.hi for q, decomposed into Gamma / Prime / Zeros parts.

CLI (v2.1 normalized):
  --T0
  --x0
  --A-prime
  --sigma-scale
  --k0-scale
  --K
  --dps
  --out

JSON (v2.1 normalized):
  {
    "kind": "q_lipschitz_bound",
    "inputs": {
      "T0": "<string>",
      "x0": "<string>",
      "A_prime": "<string>",
      "sigma_scale": "<string>",
      "k0_scale": "<string>",
      "K": <int>
    },
    "q_lipschitz": {
      "L": {
        "lo": "<string>",
        "hi": "<string>"
      },
      "decomp": {
        "Gamma": { "hi": "<string>" },
        "Prime": { "hi": "<string>" },
        "Zeros": { "hi": "<string>" }
      },
      "assumptions": { ... }
    },
    "meta": {
      "tool": "lipschitz_q_bound",
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
from datetime import timezone
import datetime as _dt
import hashlib

from mpmath import mp


# ---------------------------------------------------------------------
# Helpers (precision, serialization, hashing)
# ---------------------------------------------------------------------

def set_precision(dps: int) -> None:
    mp.dps = int(dps)


def mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def _sha256_bytes(b: bytes) -> str:
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def outward(x, dps, ulps=2):
    """
    Outward rounding helper (kept from original implementation).

    Returns a decimal string slightly above the given value x at the
    requested precision. Used both for L.hi and decomposed hi terms.
    """
    old = mp.dps
    mp.dps = dps + 5
    y = mp.mpf(x) + mp.mpf(ulps) * mp.power(10, -dps)
    s = mp.nstr(y, n=dps)
    mp.dps = old
    return s


def now_utc_iso() -> str:
    return (
        _dt.datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def write_json(path: str, payload: dict, dps: int) -> None:
    """
    Canonical JSON writer with meta block and sha256.

    meta:
      tool        = "lipschitz_q_bound"
      dps         = dps
      created_utc = ISO8601 UTC
      sha256      = hash over payload with meta.sha256 removed
    """
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)

    if "meta" not in payload or not isinstance(payload["meta"], dict):
        payload["meta"] = {}
    meta = payload["meta"]
    meta["tool"] = "lipschitz_q_bound"
    meta["dps"] = int(dps)
    meta["created_utc"] = now_utc_iso()

    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    m2 = tmp_obj.get("meta")
    if isinstance(m2, dict):
        m2.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, indent=None, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    digest = _sha256_bytes(blob)
    meta["sha256"] = digest

    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------------
# Core math (unchanged)
# ---------------------------------------------------------------------

def c_gamma(sigma):
    s = mp.mpf(sigma)
    return (mp.mpf("1.0") + 1 / s, mp.mpf("0.5") + 1 / (s * s))


def C_prime(A_prime, K, sigma, k0):
    A = mp.mpf(A_prime)
    Kp = mp.mpf(K + 1)
    s = mp.mpf(sigma)
    k = mp.mpf(k0) if mp.mpf(k0) > mp.mpf("0") else mp.mpf("1e-6")
    return mp.mpf("2") * Kp * A * (mp.mpf("1") + 1 / s) * (mp.mpf("1") + 1 / k)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Lipschitz bound for q around (sigma, k0)."
    )
    ap.add_argument("--T0", required=True, help="Base T0 for tails (t >= T0).")
    ap.add_argument("--x0", required=True, help="Cutover x0 for prime theta bound.")
    ap.add_argument(
        "--A-prime",
        dest="A_prime",
        default="1.2762",
        help="Dusart-style A_prime coefficient.",
    )
    ap.add_argument(
        "--sigma-scale",
        default="6.0",
        help="Sigma scale used in the envelope.",
    )
    ap.add_argument(
        "--k0-scale",
        default="0.25",
        help="k0 scale (was notch-scale) used in the prime tail model.",
    )
    ap.add_argument("--K", type=int, default=3, help="Auxiliary index K.")
    ap.add_argument("--dps", type=int, default=300, help="Decimal precision.")
    ap.add_argument("--out", required=True, help="Output JSON path.")
    args = ap.parse_args()

    set_precision(args.dps)

    T0 = mp.mpf(args.T0)
    x0 = mp.mpf(args.x0)

    c1, c2 = c_gamma(args.sigma_scale)
    Lg = c1 / T0 + c2 / (T0 * T0)

    Cp = C_prime(args.A_prime, args.K, args.sigma_scale, args.k0_scale)
    Lp = Cp * (mp.log(T0) / T0)

    Lz = mp.mpf("1.0") / ((mp.mpf("1") + T0) ** 2)

    L = Lg + Lp + Lz

    payload = {
        "kind": "q_lipschitz_bound",
        "inputs": {
            "T0": mp_str(T0),
            "x0": mp_str(x0),
            "A_prime": mp_str(args.A_prime),
            "sigma_scale": mp_str(args.sigma_scale),
            "k0_scale": mp_str(args.k0_scale),
            "K": int(args.K),
        },
        "q_lipschitz": {
            "L": {
                "lo": outward(L * mp.mpf("0.9999999"), args.dps),
                "hi": outward(L, args.dps),
            },
            "decomp": {
                "Gamma": {
                    "hi": outward(Lg, args.dps),
                },
                "Prime": {
                    "hi": outward(Lp, args.dps),
                },
                "Zeros": {
                    "hi": outward(Lz, args.dps),
                },
            },
            "assumptions": {
                "gamma_tail": "|psi(1/2+it)-log t| <= 1/t + 1/(2 t^2) for t>=1",
                "prime_theta": (
                    f"|theta(x)-x| <= {args.A_prime}·x/log x for x>=x0 (Dusart-style)"
                ),
                "zeros_tail": "sum_{|gamma|>T0} |h_hat(gamma)| <= 1/(1+T0)^2",
            },
        },
        "meta": {},
    }

    write_json(args.out, payload, args.dps)
    print(f"[ok] lipschitz_q_bound -> {args.out}")


if __name__ == "__main__":
    main()
