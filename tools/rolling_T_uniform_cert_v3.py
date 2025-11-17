#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import sys
import json
import time
import hashlib
from typing import Any, Dict, List, Optional, Tuple

from mpmath import mp, mpf
from tqdm import tqdm


# ---------- time / JSON utils ----------

def nows() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: str, payload: dict) -> None:
    """
    Write JSON with deterministic sha256 in meta.sha256.

    Hash is computed over the payload with any existing meta.sha256 removed,
    using canonical JSON (sorted keys, compact separators).
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

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

    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


# -------- interval-safe parsing (no isinstance(..., mpi)) --------

def _looks_interval_str(s: str) -> bool:
    s = s.strip()
    return len(s) >= 5 and s[0] == "[" and s[-1] == "]" and "," in s


def _parse_interval_str(s: str) -> Tuple[mpf, mpf]:
    lo, hi = s.strip()[1:-1].split(",", 1)
    return mp.mpf(lo.strip()), mp.mpf(hi.strip())


def as_mpf_lo(x: Any) -> mpf:
    # strings: "[lo, hi]" or "number"
    if isinstance(x, str):
        if _looks_interval_str(x):
            lo, _ = _parse_interval_str(x)
            return lo
        return mp.mpf(x)
    # dict with 'lo'/'hi'
    if isinstance(x, dict) and ("lo" in x and "hi" in x):
        return mp.mpf(str(x["lo"]))
    # duck-typed mpi-like: has .a
    if hasattr(x, "a"):
        return mp.mpf(x.a)
    # plain number
    return mp.mpf(str(x))


def as_mpf_hi(x: Any) -> mpf:
    if isinstance(x, str):
        if _looks_interval_str(x):
            _, hi = _parse_interval_str(x)
            return hi
        return mp.mpf(x)
    if isinstance(x, dict) and ("lo" in x and "hi" in x):
        return mp.mpf(str(x["hi"]))
    if hasattr(x, "b"):
        return mp.mpf(x.b)
    return mp.mpf(str(x))


def nstr(x: mpf, ndp: int = 120) -> str:
    return mp.nstr(x, ndp)


# -------- core calculus (math unchanged) --------

def compute_delta_lo(T_lo: mpf, eps_eff_lo: mpf, grid_hi: mpf,
                     Cp: mpf, ap: mpf, Cg: mpf, ag: mpf) -> mpf:
    # tails decrease with T; worst case on [T_lo, T_hi] is at T_lo
    pt_hi = Cp / (T_lo ** ap)
    gt_hi = Cg / (T_lo ** ag)
    return eps_eff_lo - pt_hi - gt_hi - grid_hi


def adaptive_cert(T0: mpf, T1: mpf, eps_eff_lo: mpf, grid_hi: mpf,
                  Cp: mpf, ap: mpf, Cg: mpf, ag: mpf,
                  delta_target: mpf, mesh_initial: int, mesh_max: int
                  ) -> Tuple[bool, mpf, Optional[Dict[str, str]], int, int, mpf]:
    """
    Returns:
      PASS, global_min, witness_or_none, intervals, max_depth, argmin_T
    """
    xs = [T0 + (T1 - T0) * (i / mp.mpf(mesh_initial)) for i in range(mesh_initial + 1)]
    work: List[Tuple[mpf, mpf, int]] = [(xs[i], xs[i + 1], 0) for i in range(mesh_initial)]
    total = 0
    max_depth = 0
    global_min = mp.inf
    argmin_T: Optional[mpf] = None
    pbar = tqdm(total=mesh_max, desc="[rolling-T]", leave=False)

    while work and total < mesh_max:
        a, b, depth = work.pop()  # DFS
        d_lo = compute_delta_lo(a, eps_eff_lo, grid_hi, Cp, ap, Cg, ag)
        if d_lo < global_min:
            global_min = d_lo
            argmin_T = a
        total += 1
        if total <= mesh_max:
            pbar.update(1)
        max_depth = max(max_depth, depth)

        if d_lo >= delta_target:
            continue
        if total >= mesh_max:
            pbar.close()
            return (
                False,
                global_min,
                {"T_left": nstr(a), "T_right": nstr(b)},
                total,
                max_depth,
                (argmin_T or a),
            )
        mid = (a + b) / 2
        work.append((mid, b, depth + 1))
        work.append((a, mid, depth + 1))
        pbar.total = min(mesh_max, pbar.total + 1)

    pbar.close()
    if not work:
        # PASS case with no refinement
        return True, global_min, None, total, max_depth, (argmin_T or T0)
    a, b, _ = work[-1]
    return (
        False,
        global_min,
        {"T_left": nstr(a), "T_right": nstr(b)},
        total,
        max_depth,
        (argmin_T or a),
    )


# -------- main CLI (v2.1 normalized) --------

def main():
    ap = argparse.ArgumentParser(
        description="Adaptive rolling uniform certificate on [T0, T1] using analytic tail bounds."
    )
    ap.add_argument("--packet-dir", required=True, help="Root packet directory.")
    ap.add_argument("--T0", required=True, help="Left endpoint T0.")
    ap.add_argument("--T1", required=True, help="Right endpoint T1.")
    ap.add_argument("--delta-target", default="1e-12", help="Target margin delta.")
    ap.add_argument("--mesh-initial", type=int, default=128, help="Initial mesh intervals.")
    ap.add_argument("--mesh-max", type=int, default=131072, help="Max refinement intervals.")
    ap.add_argument("--dps", type=int, default=220, help="Decimal precision for mpmath.")
    ap.add_argument("--digits", type=int, default=120, help="Digits of precision for output strings.")
    ap.add_argument("--out", required=True, help="Output JSON path.")
    args = ap.parse_args()

    mp.dps = args.dps
    T0 = mp.mpf(args.T0)
    T1 = mp.mpf(args.T1)
    if not (T0 < T1):
        print("[rolling_T_uniform_cert_v3] ERROR: T0 must be < T1", file=sys.stderr)
        sys.exit(1)
    delta_target = mp.mpf(args.delta_target)

    ab_path = os.path.join(args.packet_dir, "analytic_bounds.json")
    if not os.path.exists(ab_path):
        print(
            f"[rolling_T_uniform_cert_v3] ERROR: missing {ab_path}. "
            f"Run analytic_tail_fit.py first.",
            file=sys.stderr,
        )
        sys.exit(1)
    ab = load_json(ab_path)
    b = ab.get("bounds", ab)  # accept flat or nested layouts

    try:
        eps_eff_lo = as_mpf_lo(b["eps_eff_lo"])
        grid_hi = as_mpf_hi(b["grid_error_hi"])

        pt = b["prime_tail"]
        gt = b["gamma_tail"]
        Cp = as_mpf_hi(pt["C"])      # error bound -> take HI
        apow = as_mpf_lo(pt["a"])    # exponent -> take LO (worst case)
        Cg = as_mpf_hi(gt["C"])
        agow = as_mpf_lo(gt["a"])
    except Exception as e:
        print(
            f"[rolling_T_uniform_cert_v3] ERROR parsing analytic_bounds.json: {e}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Respect a declared tail domain T0 if present (clamp up)
    for tail in (b.get("prime_tail", {}), b.get("gamma_tail", {})):
        try:
            t0_req = tail.get("T0", None)
            if t0_req is not None:
                t0_req = as_mpf_lo(t0_req)
                if T0 < t0_req:
                    T0 = t0_req
        except Exception:
            pass

    t0_wall = time.time()
    PASS, delta_min, witness, nint, depth, argmin_T = adaptive_cert(
        T0,
        T1,
        eps_eff_lo,
        grid_hi,
        Cp,
        apow,
        Cg,
        agow,
        delta_target,
        args.mesh_initial,
        args.mesh_max,
    )
    dt = time.time() - t0_wall

    # Always provide a witness (even in degenerate PASS case)
    if witness is None:
        try:
            pad = max(abs(argmin_T) * mp.mpf("1e-12"), mp.mpf("1.0"))
        except Exception:
            pad = mp.mpf("1.0")
        T_left = argmin_T - pad
        T_right = argmin_T + pad
        delta_at_T_star = compute_delta_lo(
            argmin_T, eps_eff_lo, grid_hi, Cp, apow, Cg, agow
        )
        witness = {
            "T_left": nstr(T_left, args.digits),
            "T_right": nstr(T_right, args.digits),
            "T_star": nstr(argmin_T, args.digits),
            "delta_at_T_star": nstr(delta_at_T_star, args.digits),
            "depth": 0,
            "mode": "argmin-degenerate",
        }

    payload = {
        "kind": "rolling_T_uniform",
        "inputs": {
            "packet_dir": args.packet_dir,
            "T0": nstr(T0, args.digits),
            "T1": nstr(T1, args.digits),
            "delta_target": nstr(delta_target, args.digits),
            "mesh_initial": int(args.mesh_initial),
            "mesh_max": int(args.mesh_max),
            "dps": int(args.dps),
            "digits": int(args.digits),
        },
        "bounds": {
            "eps_eff_lo": nstr(eps_eff_lo, args.digits),
            "grid_error_hi": nstr(grid_hi, args.digits),
            "prime_tail": {
                "C": nstr(Cp, args.digits),
                "a": nstr(apow, args.digits),
            },
            "gamma_tail": {
                "C": nstr(Cg, args.digits),
                "a": nstr(agow, args.digits),
            },
        },
        "mesh": {
            "intervals": int(nint),
            "max_depth": int(depth),
            "elapsed_sec": f"{dt:.6f}",
        },
        "result": {
            "PASS": bool(PASS),
            "delta_min": nstr(delta_min, args.digits),
            "witness": witness,
        },
        "meta": {
            "tool": "rolling_T_uniform_cert_v3",
            "dps": int(args.dps),
            "created_utc": nows(),
        },
    }

    write_json(args.out, payload)

    print(
        f"[rolling_T_uniform_cert_v3] "
        f"{'PASS' if PASS else 'FAIL'}  "
        f"delta_min={mp.nstr(delta_min, 40)}  "
        f"intervals={nint}  depth={depth}  time={round(dt, 2)}s"
    )
    if witness:
        print(
            "[rolling_T_uniform_cert_v3] witness interval: "
            f"{witness}"
        )


if __name__ == "__main__":
    main()
