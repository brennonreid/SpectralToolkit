#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
subspace_psd_cholesky.py — Bochner Gram + (pivoted) Cholesky PSD cert

Purpose:
  PSD certification of subspaces via Bochner Gram matrix construction and
  (pivoted) Cholesky factorization.

CLI (v2.1 normalized):
  --basis        : basis family (currently "gaussian")
  --atoms        : number of atoms in the subspace
  --sigma-min    : minimum sigma
  --sigma-max    : maximum sigma
  --k0-min       : minimum k0
  --k0-max       : maximum k0
  --gridA        : half-width A for integration interval [-A, A]
  --mgrid        : number of trapezoid nodes (>= 2)
  --eta          : diagonal jitter added to Gram matrix (optional)
  --threads      : worker threads for Gram build
  --dps          : decimal precision for mpmath
  --out          : JSON output path
  --csv          : CSV output path (Gram entries)
  --progress     : show tqdm progress bar (optional)

JSON (v2.1 normalized):
  {
    "kind": "subspace_psd_cholesky",
    "inputs": {
      "basis_path": "<string>",      # here used as a basis identifier
      "atoms": <int>,
      "sigma_min": "<string>",
      "sigma_max": "<string>",
      "k0_min": "<string>",
      "k0_max": "<string>",
      "gridA": "<string>",
      "mgrid": <int>,
      "eta": "<string>",
      "threads": <int>
    },
    "result": {
      "chol_success": true/false,
      "min_diag_L": "<string or null>",
      "pivot_success": true/false,
      "min_pivot": "<string or null>",
      "rank": <int>,
      "psd_certified": true/false
    },
    "meta": {
      "tool": "subspace_psd_cholesky",
      "dps": <int>,
      "created_utc": "<ISO8601 UTC>",
      "elapsed_sec": "<string>",
      "sha256": "<digest>",
      ...  # additional diagnostic fields allowed
    }
  }

CSV:
  - Gram matrix entries with headers:
      i,j,h_ij
  - All headers are snake_case; h_ij is the Gram entry H[i,j].

All real-valued numeric fields in JSON are stored as strings. Integers remain
JSON integers; booleans remain JSON booleans.
"""

import argparse
import csv
import json
import math
import time
import os
import hashlib
from dataclasses import dataclass
from typing import List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime, timezone

from mpmath import mp
from tqdm import tqdm

# -------------------------
# Utilities
# -------------------------

def set_mp(dps: int):
    mp.dps = dps
    mp.eps = mp.mpf(10) ** (-(dps - 3))


def adaptive_tol_from_dps(dps: int, k: int = 8) -> mp.mpf:
    return mp.mpf(10) ** (-(dps - k))


def zero_tiny_imag(z, imag_tol):
    if isinstance(z, mp.mpf):
        return z
    zi = mp.im(z)
    if mp.fabs(zi) <= imag_tol:
        return mp.re(z)
    return z


def kahan_stream_trap(fn, a: mp.mpf, b: mp.mpf, mgrid: int) -> mp.mpf:
    """Trapezoid integral of fn over [a,b] with mgrid points, streaming Kahan sum."""
    if mgrid < 2:
        raise ValueError("mgrid must be >= 2")
    h = (b - a) / (mgrid - 1)
    s = mp.mpf("0")
    c = mp.mpf("0")
    # first point
    f0 = fn(a)
    # inner points
    x = a
    for k in range(1, mgrid - 1):
        x = a + h * k
        y = fn(x) - c
        t = s + y
        c = (t - s) - y
        s = t
    # last point
    fN = fn(b)
    s = s + (f0 + fN) * mp.mpf("0.5")
    return h * s


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

    # Clone without meta.sha256 for hashing.
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


# -------------------------
# Bochner kernel (real-symmetric form)
# -------------------------

@dataclass(frozen=True)
class Atom:
    sigma: mp.mpf
    k0: mp.mpf


def make_atoms_gaussian(n_atoms: int,
                        sigma_min: mp.mpf, sigma_max: mp.mpf,
                        k0_min: mp.mpf, k0_max: mp.mpf) -> List[Atom]:
    m = int(mp.floor(mp.sqrt(n_atoms)))
    if m * m < n_atoms:
        m += 1
    ns = m
    nk = m
    sigmas = [
        sigma_min + (sigma_max - sigma_min) * mp.mpf(i) / mp.mpf(max(ns - 1, 1))
        for i in range(ns)
    ]
    k0s = [
        k0_min + (k0_max - k0_min) * mp.mpf(j) / mp.mpf(max(nk - 1, 1))
        for j in range(nk)
    ]
    atoms = []
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


def bochner_gram_entry_real(atom_i: Atom, atom_j: Atom,
                            A: mp.mpf, mgrid: int) -> mp.mpf:
    # streaming trapezoid of hi(x)*hj(x) on [-A, A]
    def f(x: mp.mpf) -> mp.mpf:
        return h_gauss_notch(x, atom_i.sigma, atom_i.k0) * h_gauss_notch(
            x, atom_j.sigma, atom_j.k0
        )

    return kahan_stream_trap(f, -A, A, mgrid)


# -------------------------
# Pivoted Cholesky (PSD cert)
# -------------------------

def pivoted_cholesky_psd(H: mp.matrix, tol: mp.mpf) -> Tuple[bool, mp.mpf, int]:
    n = H.rows
    diag = [H[i, i] for i in range(n)]
    piv_order = list(range(n))
    L = [[mp.mpf("0") for _ in range(n)] for _ in range(n)]
    min_pivot = mp.inf
    rank = 0

    for k in range(n):
        p = max(range(k, n), key=lambda i: diag[i])
        if diag[p] < -tol:
            return False, mp.mpf(diag[p]), rank
        if diag[p] <= tol:
            break

        if p != k:
            diag[p], diag[k] = diag[k], diag[p]
            piv_order[p], piv_order[k] = piv_order[k], piv_order[p]
            for t in range(n):
                H[k, t], H[p, t] = H[p, t], H[k, t]
            for t in range(n):
                H[t, k], H[t, p] = H[t, p], H[t, k]
            L[k], L[p] = L[p], L[k]

        pivot = mp.sqrt(diag[k])
        L[k][k] = pivot
        min_pivot = min(min_pivot, pivot)
        rank += 1

        if k < n - 1:
            for i in range(k + 1, n):
                lij = H[i, k] / L[k][k] if L[k][k] != 0 else mp.mpf("0")
                L[i][k] = lij
                diag[i] = diag[i] - lij * lij
                for j in range(k + 1, i + 1):
                    H[i, j] = H[i, j] - lij * L[j][k]
                    H[j, i] = H[i, j]

    return True, min_pivot if min_pivot != mp.inf else mp.mpf("0"), rank


def try_plain_cholesky(H: mp.matrix) -> Tuple[bool, mp.mpf]:
    try:
        L = mp.cholesky(H)
        mind = min(L[i, i] for i in range(H.rows))
        return True, mind
    except Exception:
        return False, mp.mpf("0")


# -------------------------
# Main
# -------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Certified subspace PSD via Bochner Gram and pivoted Cholesky."
    )
    ap.add_argument("--basis", choices=["gaussian"], default="gaussian")
    ap.add_argument("--atoms", type=int, required=True)

    ap.add_argument("--sigma-min", type=str, required=True)
    ap.add_argument("--sigma-max", type=str, required=True)
    ap.add_argument("--k0-min", type=str, required=True)
    ap.add_argument("--k0-max", type=str, required=True)

    ap.add_argument("--gridA", type=str, default="50")  # log-domain half-width
    ap.add_argument("--mgrid", type=int, default=2049)  # integration nodes

    ap.add_argument("--dps", type=int, default=120)
    ap.add_argument("--eta", type=str, default="0")

    ap.add_argument("--threads", type=int, default=1, help="Worker threads for Gram build")
    ap.add_argument("--progress", action="store_true", help="Show tqdm progress bar")

    ap.add_argument("--out", type=str, required=True)
    ap.add_argument("--csv", type=str, required=True)

    args = ap.parse_args()
    set_mp(args.dps)

    sigma_min = mp.mpf(args.sigma_min)
    sigma_max = mp.mpf(args.sigma_max)
    k0_min = mp.mpf(args.k0_min)
    k0_max = mp.mpf(args.k0_max)
    A = mp.mpf(args.gridA)
    mgrid = int(args.mgrid)
    eta = mp.mpf(args.eta)

    imag_tol = adaptive_tol_from_dps(args.dps, k=6)

    n = int(args.atoms)
    atoms = make_atoms_gaussian(n, sigma_min, sigma_max, k0_min, k0_max)

    t0 = time.time()

    # Pre-allocate H
    H = mp.matrix(n, n)

    # Build list of upper-triangle index pairs
    pairs = [(i, j) for i in range(n) for j in range(i, n)]

    def task(pair):
        i, j = pair
        gij = bochner_gram_entry_real(atoms[i], atoms[j], A=A, mgrid=mgrid)
        gij = zero_tiny_imag(gij, imag_tol)
        if not isinstance(gij, mp.mpf):
            raise SystemExit(
                "Backend produced materially complex entry; real-symmetric PSD witness invalid."
            )
        return (i, j, gij)

    if args.threads and args.threads > 1:
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = [ex.submit(task, p) for p in pairs]
            it = as_completed(futs)
            if args.progress:
                it = tqdm(it, total=len(pairs), desc="Gram (upper)")
            for fut in it:
                i, j, gij = fut.result()
                H[i, j] = gij
                H[j, i] = gij
    else:
        it = pairs
        if args.progress:
            it = tqdm(it, total=len(pairs), desc="Gram (upper)")
        for i, j in it:
            gij = bochner_gram_entry_real(atoms[i], atoms[j], A=A, mgrid=mgrid)
            gij = zero_tiny_imag(gij, imag_tol)
            if not isinstance(gij, mp.mpf):
                raise SystemExit(
                    "Backend produced materially complex entry; real-symmetric PSD witness invalid."
                )
            H[i, j] = gij
            H[j, i] = gij

    if eta != 0:
        for i in range(n):
            H[i, i] = H[i, i] + eta

    chol_ok, min_diagL = try_plain_cholesky(H.copy())
    pivot_ok, min_pivot, rank = (False, mp.mpf("0"), 0)

    if not chol_ok:
        tol = adaptive_tol_from_dps(args.dps, k=6)
        pivot_ok, min_pivot, rank = pivoted_cholesky_psd(H.copy(), tol)

    elapsed = time.time() - t0
    cert_ok = chol_ok or pivot_ok

    # CSV (snake_case header)
    os.makedirs(os.path.dirname(os.path.abspath(args.csv)), exist_ok=True)
    with open(args.csv, "w", newline="", encoding="utf-8") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["i", "j", "h_ij"])
        for i in range(n):
            for j in range(n):
                w.writerow([i, j, mp.nstr(H[i, j], 20)])

    # JSON (v2.1-normalized shape)
    inputs_block = {
        "basis_path": str(args.basis),  # basis identifier (no external path used)
        "atoms": n,
        "sigma_min": mp_str(sigma_min),
        "sigma_max": mp_str(sigma_max),
        "k0_min": mp_str(k0_min),
        "k0_max": mp_str(k0_max),
        "gridA": mp_str(A),
        "mgrid": mgrid,
        "eta": mp_str(eta),
        "threads": int(args.threads),
    }

    result_block = {
        "chol_success": bool(chol_ok),
        "min_diag_L": mp_str(min_diagL) if chol_ok else None,
        "pivot_success": bool(pivot_ok),
        "min_pivot": mp_str(min_pivot) if pivot_ok else None,
        "rank": rank if pivot_ok else (n if chol_ok else 0),
        "psd_certified": bool(cert_ok),
    }

    meta_block = {
        "tool": "subspace_psd_cholesky",
        "dps": int(args.dps),
        "created_utc": now_utc_iso(),
        "elapsed_sec": "{:.6f}".format(elapsed),
        "algo": "subspace_psd_cholesky/pivoted",
        "n": n,
        "gridA": mp_str(A),
        "mgrid": mgrid,
        "eta": mp_str(eta),
        "threads": int(args.threads),
    }

    payload = {
        "kind": "subspace_psd_cholesky",
        "inputs": inputs_block,
        "result": result_block,
        "meta": meta_block,
    }

    write_json(args.out, payload)

    if cert_ok:
        if chol_ok:
            print(f"[subspace-psd] PASS (Cholesky). min_diag_L={min_diagL}")
        else:
            print(f"[subspace-psd] PASS (Pivoted Cholesky). min_pivot={min_pivot} rank={rank}/{n}")
    else:
        print("[subspace-psd] FAIL (not PSD under current precision/tolerances).")


if __name__ == "__main__":
    main()
