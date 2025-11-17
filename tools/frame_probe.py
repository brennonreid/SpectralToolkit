#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
frame_probe.py — probe eigenvalues / condition number of Gaussian frames.

Purpose:
  Build an L2 Gram matrix for a finite Gaussian-notch frame over a real line
  grid, compute its eigenvalues, and summarize min/max eigenvalues and the
  resulting condition number and strict positivity.

CLI (v2.1 normalized for this tool):
  --dict        : dictionary type ("gaussian" only for now)
  --atoms       : number of atoms in the frame
  --sigma-min   : min sigma
  --sigma-max   : max sigma
  --k0-min      : min k0
  --k0-max      : max k0
  --A           : half-width of integration window
  --mgrid       : number of grid points on [-A, A]
  --dps         : decimal precision
  --threads     : number of threads for row-building
  --chunk       : reserved (unused)
  --tol         : override tolerance for strictly_positive (optional)
  --out         : JSON summary
  --csv         : CSV of eigenvalues

JSON (v2.1 normalized):
  {
    "kind": "frame_probe",
    "inputs": {
      "dict": "gaussian",
      "atoms": <int>,
      "sigma_min": "<string>",
      "sigma_max": "<string>",
      "k0_min": "<string>",
      "k0_max": "<string>",
      "A": "<string>",
      "mgrid": <int>,
      "dps": <int>,
      "threads": <int>,
      "tol": "<string>"
    },
    "results": {
      "min_eigenvalue": "<string>",
      "max_eigenvalue": "<string>",
      "condition_number": "<string>",
      "strictly_positive": true/false,
      "tolerance": "<string>"
    },
    "meta": {
      "tool": "frame_probe",
      "algo": "frame_probe/L2_gram",
      "dps": <int>,
      "created_utc": "<ISO8601 UTC>",
      "elapsed_sec": "<string>",
      "sha256": "<digest>"
    }
  }

All real-valued numeric fields are stored as strings; counts remain integers.
"""

import argparse
import csv
import json
import math
import time
from dataclasses import dataclass
from typing import List, Tuple

import os
import hashlib
import datetime as _dt
from datetime import timezone

from mpmath import mp

try:
    from tqdm import tqdm
    TQDM = True
except Exception:
    TQDM = False

from concurrent.futures import ThreadPoolExecutor, as_completed


# ---------------------------------------------------------------------
# Precision / numeric helpers (math unchanged)
# ---------------------------------------------------------------------

def set_mp(dps: int):
    mp.dps = dps
    mp.eps = mp.mpf(10) ** (-(dps - 3))


def tol_from_dps(dps: int, k: int = 6) -> mp.mpf:
    return mp.mpf(10) ** (-(dps - k))


def kahan_sum(vals: List[mp.mpf]) -> mp.mpf:
    s = mp.mpf("0")
    c = mp.mpf("0")
    for x in vals:
        y = x - c
        t = s + y
        c = (t - s) - y
        s = t
    return s


def nstr(x, ndp=20):
    return mp.nstr(x, ndp)


def to_real_scalar(x) -> mp.mpf:
    if isinstance(x, mp.mpf):
        return x
    if isinstance(x, str):
        return mp.mpf(x)
    if isinstance(x, mp.mpc):
        if mp.fabs(mp.im(x)) <= tol_from_dps(mp.dps, 8):
            return mp.re(x)
        raise TypeError(f"Eigenvalue has non-negligible imaginary part: {x}")
    if isinstance(x, mp.matrix):
        if x.rows == 1 and x.cols == 1:
            return to_real_scalar(x[0, 0])
        raise TypeError("Unexpected eigenvalue container shape.")
    return mp.mpf(x)


def mp_str(x, ndp=None) -> str:
    if ndp is None:
        return mp.nstr(mp.mpf(x), mp.dps, strip_zeros=False)
    return mp.nstr(mp.mpf(x), ndp, strip_zeros=False)


# ---------------------------------------------------------------------
# JSON / meta helpers
# ---------------------------------------------------------------------

def now_utc_iso() -> str:
    return (
        _dt.datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def sha256_payload(payload: dict) -> str:
    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    meta = tmp_obj.get("meta")
    if isinstance(meta, dict):
        meta.pop("sha256", None)
    blob = json.dumps(
        tmp_obj,
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    h = hashlib.sha256()
    h.update(blob)
    return h.hexdigest()


def write_json(path: str, payload: dict, dps: int, elapsed_sec: float) -> None:
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)

    if "meta" not in payload or not isinstance(payload["meta"], dict):
        payload["meta"] = {}

    meta = payload["meta"]
    meta["tool"] = "frame_probe"
    meta["algo"] = "frame_probe/L2_gram"
    meta["dps"] = int(dps)
    meta["created_utc"] = now_utc_iso()
    meta["elapsed_sec"] = f"{elapsed_sec:.6f}"

    meta["sha256"] = sha256_payload(payload)

    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as fj:
        json.dump(payload, fj, indent=2, ensure_ascii=False, sort_keys=True)
        fj.write("\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------------
# Frame construction and Gram matrix (math unchanged)
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class Atom:
    sigma: mp.mpf
    k0: mp.mpf


def make_atoms_gaussian(
    n_atoms: int,
    sigma_min: mp.mpf,
    sigma_max: mp.mpf,
    k0_min: mp.mpf,
    k0_max: mp.mpf,
) -> List[Atom]:
    m = int(mp.floor(mp.sqrt(n_atoms)))
    if m * m < n_atoms:
        m += 1
    ns = m
    nk = m
    sigmas = [
        sigma_min
        + (sigma_max - sigma_min) * mp.mpf(i) / mp.mpf(max(ns - 1, 1))
        for i in range(ns)
    ]
    k0s = [
        k0_min
        + (k0_max - k0_min) * mp.mpf(j) / mp.mpf(max(nk - 1, 1))
        for j in range(nk)
    ]
    atoms: List[Atom] = []
    for s in sigmas:
        for k in k0s:
            atoms.append(Atom(s, k))
            if len(atoms) == n_atoms:
                return atoms
    return atoms[:n_atoms]


def h_gauss_notch(x: mp.mpf, sigma: mp.mpf, k0: mp.mpf) -> mp.mpf:
    g = mp.e ** (-(x ** 2) / (sigma ** 2))
    notch = mp.mpf("1") - mp.e ** (-(x - k0) ** 2)
    return g * notch


def inner_prod_L2(ai: Atom, aj: Atom, xs: List[mp.mpf]) -> mp.mpf:
    vals = []
    for x in xs:
        vals.append(h_gauss_notch(x, ai.sigma, ai.k0) * h_gauss_notch(x, aj.sigma, aj.k0))
    s = kahan_sum(vals) - (vals[0] + vals[-1]) / 2
    A = mp.fabs(xs[-1])
    step = (2 * A) / (len(xs) - 1)
    return step * s


def build_grid(A: mp.mpf, mgrid: int) -> List[mp.mpf]:
    return [
        (-A + (2 * A) * mp.mpf(k) / mp.mpf(mgrid - 1))
        for k in range(mgrid)
    ]


def accumulate_row(i: int, atoms: List[Atom], xs: List[mp.mpf]) -> Tuple[int, List[mp.mpf]]:
    ai = atoms[i]
    row = []
    for j in range(i, len(atoms)):
        aj = atoms[j]
        gij = inner_prod_L2(ai, aj, xs)
        row.append(gij)
    return i, row


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Frame probe via L2 Gram (robust eigenspectrum)."
    )
    ap.add_argument("--dict", choices=["gaussian"], default="gaussian")
    ap.add_argument("--atoms", type=int, required=True)
    ap.add_argument("--sigma-min", type=str, required=True)
    ap.add_argument("--sigma-max", type=str, required=True)
    ap.add_argument("--k0-min", type=str, required=True)
    ap.add_argument("--k0-max", type=str, required=True)
    ap.add_argument("--A", type=str, default="50")
    ap.add_argument("--mgrid", type=int, default=4097)
    ap.add_argument("--dps", type=int, default=120)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--tol", type=str, default=None)
    ap.add_argument("--out", type=str, required=True)
    ap.add_argument("--csv", type=str, required=True)
    args = ap.parse_args()

    set_mp(args.dps)
    A = mp.mpf(args.A)
    mgrid = int(args.mgrid)
    tol = mp.mpf(args.tol) if args.tol is not None else tol_from_dps(args.dps, 6)

    sigma_min = mp.mpf(args.sigma_min)
    sigma_max = mp.mpf(args.sigma_max)
    k0_min = mp.mpf(args.k0_min)
    k0_max = mp.mpf(args.k0_max)

    n = int(args.atoms)
    if args.dict != "gaussian":
        raise SystemExit("Only --dict gaussian is currently supported.")

    atoms = make_atoms_gaussian(n, sigma_min, sigma_max, k0_min, k0_max)
    xs = build_grid(A, mgrid)
    t0 = time.time()

    G = mp.matrix(n, n)

    if args.threads <= 1:
        rng = range(n)
        iterator = tqdm(rng, desc="[frame-probe] rows", leave=False) if TQDM else rng
        for i in iterator:
            _, row = accumulate_row(i, atoms, xs)
            for k, gij in enumerate(row):
                j = i + k
                G[i, j] = gij
                G[j, i] = gij
    else:
        futures = []
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            for i in range(n):
                futures.append(ex.submit(accumulate_row, i, atoms, xs))
            iterator = (
                tqdm(
                    as_completed(futures),
                    total=len(futures),
                    desc="[frame-probe] rows",
                    leave=False,
                )
                if TQDM
                else as_completed(futures)
            )
            rows = [None] * n
            for fut in iterator:
                i, row = fut.result()
                rows[i] = row
        for i in range(n):
            row = rows[i]
            for k, gij in enumerate(row):
                j = i + k
                G[i, j] = gij
                G[j, i] = gij

    evals, _ = mp.eig(G)
    evals_real = [mp.mpf(str(e)) if not isinstance(e, mp.mpf) else e for e in evals]
    evals_real.sort()

    A_theta = evals_real[0]
    B_theta = evals_real[-1]
    elapsed = time.time() - t0

    # Write eigenvalue CSV
    with open(args.csv, "w", newline="", encoding="utf-8") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["idx", "eigenvalue"])
        for i, ev in enumerate(evals_real):
            w.writerow([i, nstr(ev, 30)])

    cond = B_theta / A_theta if A_theta != 0 else mp.inf
    strictly_positive = bool(A_theta > tol)

    payload = {
        "kind": "frame_probe",
        "inputs": {
            "dict": args.dict,
            "atoms": n,
            "sigma_min": nstr(sigma_min, 18),
            "sigma_max": nstr(sigma_max, 18),
            "k0_min": nstr(k0_min, 18),
            "k0_max": nstr(k0_max, 18),
            "A": nstr(A, 18),
            "mgrid": mgrid,
            "dps": int(args.dps),
            "threads": int(args.threads),
            "tol": nstr(tol, 20),
        },
        "results": {
            "min_eigenvalue": nstr(A_theta, 30),
            "max_eigenvalue": nstr(B_theta, 30),
            "condition_number": nstr(cond, 30),
            "strictly_positive": strictly_positive,
            "tolerance": nstr(tol, 20),
        },
        "meta": {},
    }

    write_json(args.out, payload, args.dps, elapsed)

    pos = "YES" if strictly_positive else "NO"
    print(
        f"[frame-probe] spectrum: min={nstr(A_theta,20)}  "
        f"max={nstr(B_theta,20)}  strictly_positive={pos}  "
        f"(tol={nstr(tol,10)})"
    )


if __name__ == "__main__":
    main()
