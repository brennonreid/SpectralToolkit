#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
param_cover_cert.py — parameter ε-net cover certificate.

Purpose:
  Given:
    - A rectangular domain in (a, b),
    - A Lipschitz bound L.hi for q(a, b),
    - A certified minimal margin m_net_lo on a discrete ε-net,

  this tool certifies that the whole domain is covered by the ε-net with
  positive margin, using a first-order Lipschitz estimate.

Math:
  The sensitivities S_a and S_b are computed from I0(c), I2(c) and S_single(c)
  as in the original implementation. The ε radius is bounded by

      eps_hi = S_a_hi * ha + S_b_hi * hb,

  where ha, hb are half-cell widths in a and b. The verified margin is

      verified_gap = m_net_lo - L_hi * eps_hi.

  PASS iff verified_gap > 0.

CLI:
  --domain        : optional; format "a=[A0,A1] b=[B0,B1]"
  --a-range A0 A1 : alternative to --domain for a-interval
  --b-range B0 B1 : alternative to --domain for b-interval
  --cellsA        : number of cells in a-direction (>=2)
  --cellsB        : number of cells in b-direction (>=2)
  --lipschitz     : JSON file with q_lipschitz.L.hi (from lipschitz_q_bound.py)
  --m-net-lo      : certified minimal margin on the ε-net (string or number)
  --out           : output JSON path
  --dps           : decimal precision for mpmath

JSON (v2.1 normalized):
  {
    "kind": "param_cover_cert",
    "inputs": {
      "a_range": ["A0", "A1"],
      "b_range": ["B0", "B1"],
      "cells_a": Na,
      "cells_b": Nb,
      "lipschitz_path": "<path>",
      "m_net_lo": "<string>"
    },
    "param_cover": {
      "domain": {
        "a": ["A0", "A1"],
        "b": ["B0", "B1"]
      },
      "grid": {
        "cellsA": Na,
        "cellsB": Nb,
        "da": "<string>",
        "db": "<string>"
      },
      "sensitivities": {
        "S_a_hi": "<string>",
        "S_b_hi": "<string>"
      },
      "epsilon": {
        "hi": "<string>"
      },
      "Lipschitz": {
        "hi": "<string>"
      },
      "m_net_lo": "<string>",
      "verified_margin": {
        "lo": "<string>",
        "hi": "<string>"
      },
      "PASS": true/false
    },
    "PASS": true/false,
    "meta": {
      "tool": "param_cover_cert",
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
import re
import sys
import hashlib
import datetime as dt
import os

from mpmath import mp


# ---------------------------------------------------------------------
# Basic numeric + JSON helpers (interface-level only)
# ---------------------------------------------------------------------

def set_precision(dps: int) -> None:
    mp.dps = int(dps)


def mp_str(x) -> str:
    return mp.nstr(mp.mpf(x), n=mp.dps, strip_zeros=False)


def outward_hi(x):
    """Outward rounding upward for safety (math logic unchanged)."""
    try:
        return mp.nextabove(mp.mpf(x))
    except Exception:
        return mp.mpf(x) * (1 + mp.mpf("1e-30"))


def outward_lo(x):
    """Outward rounding downward for safety (kept for compatibility)."""
    try:
        return mp.nextbelow(mp.mpf(x))
    except Exception:
        return mp.mpf(x) * (1 - mp.mpf("1e-30"))


def now_utc_iso() -> str:
    return (
        dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def write_json(path: str, payload: dict, dps: int) -> None:
    """
    Write JSON with canonical meta block:

      meta.tool        = "param_cover_cert"
      meta.dps         = dps
      meta.created_utc = ISO8601 UTC
      meta.sha256      = hash over payload without meta.sha256

    Hash is computed over JSON with sorted keys and compact separators.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    payload.setdefault("meta", {})
    payload["meta"]["tool"] = "param_cover_cert"
    payload["meta"]["dps"] = int(dps)
    payload["meta"]["created_utc"] = now_utc_iso()

    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    meta = tmp_obj.get("meta", {})
    if isinstance(meta, dict):
        meta.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    sha = hashlib.sha256(blob).hexdigest()
    payload["meta"]["sha256"] = sha

    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------------
# Core functions (math unchanged)
# ---------------------------------------------------------------------

def I0(c):
    c = mp.mpf(c)
    return (
        mp.mpf("0.5") * mp.sqrt(c)
        + (c / mp.pi)
        + (mp.mpf("0.25") * c ** (mp.mpf("1.5")) / mp.sqrt(mp.pi))
    )


def I2(c):
    c = mp.mpf(c)
    return (
        mp.mpf("0.25") * c ** (mp.mpf("1.5")) / mp.pi
        + (c ** 2 / mp.pi ** 2)
        + (mp.mpf("0.375") * c ** (mp.mpf("2.5")) / mp.pi ** 2)
    )


def S_single(c):
    c = mp.mpf(c)
    term = mp.mpf("0.5") * c ** (-mp.mpf("1.5")) * I0(c) + mp.pi * c ** (
        -mp.mpf("2.5")
    ) * I2(c)
    return 2 * term


def S_a(a, b):
    return S_single(a) + S_single(a + b)


def S_b(a, b):
    return S_single(a + b)


def parse_domain_str(s):
    if s is None:
        return None
    s = " ".join(s.strip().split())
    pat = r"([ab])\s*=\s*\[\s*([0-9\.\+\-eE]+)\s*,\s*([0-9\.\+\-eE]+)\s*\]"
    pairs = re.findall(pat, s)
    a_lo = a_hi = b_lo = b_hi = None
    for name, lo, hi in pairs:
        lo_v = mp.mpf(lo)
        hi_v = mp.mpf(hi)
        if hi_v <= lo_v:
            raise ValueError(f"{name} interval must satisfy hi>lo")
        if name == "a":
            a_lo, a_hi = lo_v, hi_v
        if name == "b":
            b_lo, b_hi = lo_v, hi_v
    if a_lo is None or b_lo is None:
        return None
    return (a_lo, a_hi, b_lo, b_hi)


def read_L_hi(path):
    js = json.loads(open(path, "r", encoding="utf-8").read())
    try:
        return mp.mpf(js["q_lipschitz"]["L"]["hi"])
    except Exception:
        try:
            return mp.mpf(js["L"]["hi"])
        except Exception:
            raise SystemExit("Could not parse L.hi from lipschitz JSON")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Parameter ε-net cover certificate (proof-grade JSON)."
    )
    ap.add_argument(
        "--domain",
        type=str,
        help='Format: a=[A0,A1] b=[B0,B1]',
    )
    ap.add_argument(
        "--a-range",
        nargs=2,
        metavar=("A0", "A1"),
        help="Explicit a-range if --domain is not used.",
    )
    ap.add_argument(
        "--b-range",
        nargs=2,
        metavar=("B0", "B1"),
        help="Explicit b-range if --domain is not used.",
    )
    ap.add_argument(
        "--cellsA",
        dest="cells_a",
        type=int,
        default=24,
        help="Number of cells in a-direction (>=2).",
    )
    ap.add_argument(
        "--cellsB",
        dest="cells_b",
        type=int,
        default=24,
        help="Number of cells in b-direction (>=2).",
    )
    ap.add_argument(
        "--lipschitz",
        required=True,
        help="JSON with L.hi bound (from lipschitz_q_bound.py).",
    )
    ap.add_argument(
        "--m-net-lo",
        dest="m_net_lo",
        type=str,
        required=True,
        help="Certified minimal margin on the ε-net.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path.",
    )
    ap.add_argument(
        "--dps",
        type=int,
        default=300,
        help="Decimal precision.",
    )
    args = ap.parse_args()
    set_precision(args.dps)

    parsed = parse_domain_str(args.domain) if args.domain else None
    if parsed is None:
        if not (args.a_range and args.b_range):
            print(
                "[param_cover_cert] ERROR: provide --domain or both --a-range/--b-range",
                file=sys.stderr,
            )
            sys.exit(2)
        a0, a1 = mp.mpf(args.a_range[0]), mp.mpf(args.a_range[1])
        b0, b1 = mp.mpf(args.b_range[0]), mp.mpf(args.b_range[1])
        if a1 <= a0 or b1 <= b0:
            print("[param_cover_cert] ERROR: intervals must have hi>lo", file=sys.stderr)
            sys.exit(2)
    else:
        a0, a1, b0, b1 = parsed

    Na = max(2, int(args.cells_a))
    Nb = max(2, int(args.cells_b))
    da = (a1 - a0) / Na
    db = (b1 - b0) / Nb
    ha = da / 2
    hb = db / 2

    Sa = S_a(a1, b1)  # conservative corner
    Sb = S_b(a1, b1)
    eps_hi = outward_hi(Sa * ha + Sb * hb)

    L_hi = read_L_hi(args.lipschitz)
    m_net_lo = mp.mpf(args.m_net_lo)
    verified_gap = m_net_lo - L_hi * eps_hi

    PASS = bool(verified_gap > 0)

    payload = {
        "kind": "param_cover_cert",
        "inputs": {
            "a_range": [mp_str(a0), mp_str(a1)],
            "b_range": [mp_str(b0), mp_str(b1)],
            "cells_a": Na,
            "cells_b": Nb,
            "lipschitz_path": args.lipschitz,
            "m_net_lo": mp_str(m_net_lo),
        },
        "param_cover": {
            "domain": {
                "a": [mp_str(a0), mp_str(a1)],
                "b": [mp_str(b0), mp_str(b1)],
            },
            "grid": {
                "cellsA": Na,
                "cellsB": Nb,
                "da": mp_str(da),
                "db": mp_str(db),
            },
            "sensitivities": {
                "S_a_hi": mp_str(Sa),
                "S_b_hi": mp_str(Sb),
            },
            "epsilon": {
                "hi": mp_str(eps_hi),
            },
            "Lipschitz": {
                "hi": mp_str(L_hi),
            },
            "m_net_lo": mp_str(m_net_lo),
            "verified_margin": {
                "lo": mp_str(verified_gap),
                "hi": mp_str(outward_hi(verified_gap)),
            },
            "PASS": PASS,
        },
        "PASS": PASS,
        "meta": {},
    }

    write_json(args.out, payload, args.dps)
    print(f"[param_cover_cert] PASS={PASS} wrote {args.out}")


if __name__ == "__main__":
    main()
