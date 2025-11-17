Here’s `param_search.py` normalized to your v2.1 spec (math / trial sequence unchanged; only CLI, naming, and JSON/CSV handling updated). 

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
param_search.py — grid search over (sigma, k0) with per-trial certificates.

Purpose:
  Sweep a rectangular grid in (sigma, k0), run the full local certification
  pipeline for each point, and log the resulting margins and gaps to CSV.

CLI (v2.1 normalized):
  --zeros        : path to zeros text file
  --N            : number of zeros to use
  --sigma-min    : minimum sigma in grid
  --sigma-max    : maximum sigma in grid
  --sigma-step   : sigma step size
  --k0-min       : minimum k0 in grid
  --k0-max       : maximum k0 in grid
  --k0-step      : k0 step size
  --T0           : base T0 for prime tail envelope (default 1e9)
  --dps          : decimal precision
  --outdir       : output directory root (per-trial subdirs are created)
  --workers      : reserved for future parallelization (currently unused)

Artifacts per trial (under outdir/trial_sXXX_kYYY):
  band_cert.json
  prime_tail_envelope.json
  prime_block_norm.json
  gamma_tails.json
  continuum_operator_rollup.json

CSV:
  outdir/param_search_results.csv with header:
    sigma,k0,lhs_total,epsilon_eff,gap,pass,band_margin,prime_cap,prime_tail,gamma_tail
"""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    tqdm = None


def frange(start, stop, step):
    x = float(start)
    stop = float(stop)
    step = float(step)
    # inclusive range with small epsilon
    while x <= stop + 1e-12:
        yield round(x, 12)
        x += step


def run(cmd, **kwargs):
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        **kwargs,
    )
    if p.returncode != 0:
        print("CMD FAIL:", " ".join(cmd))
        print(p.stdout)
        raise SystemExit(p.returncode)
    return p.stdout


def jload(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def main():
    ap = argparse.ArgumentParser(
        description="Grid search over (sigma, k0) with per-trial certificates."
    )
    ap.add_argument("--zeros", required=True)
    ap.add_argument("--N", type=int, default=20000)

    ap.add_argument("--sigma-min", type=float, required=True)
    ap.add_argument("--sigma-max", type=float, required=True)
    ap.add_argument("--sigma-step", type=float, required=True)

    ap.add_argument("--k0-min", type=float, required=True)
    ap.add_argument("--k0-max", type=float, required=True)
    ap.add_argument("--k0-step", type=float, required=True)

    ap.add_argument("--T0", type=str, default="1000000000")
    ap.add_argument("--dps", type=int, default=260)
    ap.add_argument("--outdir", default="PROOF_PACKET")
    ap.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Reserved for future parallel implementation.",
    )
    args = ap.parse_args()

    PY = sys.executable or "python"
    N = int(args.N)
    sig_a, sig_b, sig_s = args.sigma_min, args.sigma_max, args.sigma_step
    k0_a, k0_b, k0_s = args.k0_min, args.k0_max, args.k0_step
    ptT0 = str(args.T0)
    dps = int(args.dps)
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    # CSV log (normalized header)
    csv_path = outdir / "param_search_results.csv"
    csv_f = open(csv_path, "w", newline="", encoding="utf-8")
    writer = csv.writer(csv_f)
    writer.writerow(
        [
            "sigma",
            "k0",
            "lhs_total",
            "epsilon_eff",
            "gap",
            "pass",
            "band_margin",
            "prime_cap",
            "prime_tail",
            "gamma_tail",
        ]
    )

    best = None
    tried = 0

    # Count total trials for tqdm
    total_sigma = int(round((sig_b - sig_a) / sig_s)) + 1
    total_k0 = int(round((k0_b - k0_a) / k0_s)) + 1
    total = total_sigma * total_k0
    bar = tqdm(total=total, desc="param search") if (tqdm is not None) else None

    # Shared inputs paths (window + bands)
    win_path = Path("packs/rh/inputs/window.json")
    bands_path = Path("packs/rh/inputs/auto_bands.json")

    for sigma in frange(sig_a, sig_b, sig_s):
        for k0 in frange(k0_a, k0_b, k0_s):
            tried += 1
            tag = f"s{sigma:.6f}_k{k0:.6f}"
            run_dir = outdir / f"trial_{tag}"
            ensure_dir(run_dir)

            # 1) window (normalized window_gen.py CLI)
            run(
                [
                    PY,
                    "tools/window_gen.py",
                    "--mode",
                    "gauss",
                    "--sigma",
                    f"{sigma}",
                    "--k0",
                    f"{k0}",
                    "--dps",
                    f"{dps}",
                    "--out",
                    str(win_path),
                ]
            )

            # 2) bands (normalized bands_make.py CLI: --window-config)
            run(
                [
                    PY,
                    "tools/bands_make.py",
                    "--window-config",
                    str(win_path),
                    "--out",
                    str(bands_path),
                    "--critical-left",
                    "-0.50",
                    "--critical-right",
                    "0.50",
                    "--grid",
                    "6000",
                    "--digits",
                    "50",
                ]
            )

            # 3) band cert (normalized band_cert.py CLI: --window-config)
            bc_path = run_dir / "band_cert.json"
            run(
                [
                    PY,
                    "tools/band_cert.py",
                    "--window-config",
                    str(win_path),
                    "--bands",
                    str(bands_path),
                    "--out",
                    str(bc_path),
                    "--dps",
                    f"{dps}",
                    "--grid",
                    "6000",
                ]
            )
            bc = jload(bc_path)

            # 4) prime block (normalized prime_block_norm.py CLI)
            pb_path = run_dir / "prime_block_norm.json"
            run(
                [
                    PY,
                    "tools/prime_block_norm.py",
                    "--zeros",
                    args.zeros,
                    "--N",
                    f"{N}",
                    "--sigma",
                    f"{sigma}",
                    "--k0",
                    f"{k0}",
                    "--out",
                    str(pb_path),
                    "--dps",
                    f"{dps}",
                ]
            )

            # Flatten prime block to a simple cap file
            pb_flat = run_dir / "prime_block_norm.flat.json"
            with open(pb_path, "r", encoding="utf-8") as f:
                P = json.load(f)
            used = (
                (P.get("prime_block_norm", {}) or {}).get("used_operator_norm")
                or (P.get("prime_block_norm", {}) or {}).get("cap_total_hi")
                or (P.get("numbers", {}) or {}).get("cap_total_hi")
            )
            with open(pb_flat, "w", encoding="utf-8") as f_out:
                json.dump({"used_operator_norm": used}, f_out, indent=2)

            # 5) prime tail envelope (normalized prime_tail_envelope.py CLI)
            pt_path = run_dir / "prime_tail_envelope.json"
            run(
                [
                    PY,
                    "tools/prime_tail_envelope.py",
                    "--T0",
                    ptT0,
                    "--sigma",
                    f"{sigma}",
                    "--k0",
                    f"{k0}",
                    "--A-prime",
                    "1.2762",
                    "--K",
                    "3",
                    "--out",
                    str(pt_path),
                    "--dps",
                    f"{dps}",
                ]
            )
            with open(pt_path, "r", encoding="utf-8") as f:
                PT = json.load(f)
            prime_tail_val = (
                (PT.get("prime_tail", {}) or {}).get("env_T0_hi")
                or (PT.get("prime_tail", {}) or {}).get("norm")
                or (PT.get("numbers", {}) or {}).get("prime_tail_norm")
                or PT.get("prime_tail_norm")
            )

            # 6) gamma tails (T* = last zero used in this N)
            with open(args.zeros, "r", encoding="utf-8") as fz:
                lines = [ln for ln in fz if ln.strip()]
                if len(lines) < N:
                    raise SystemExit(
                        f"zeros file has only {len(lines)} lines; N={N}"
                    )
                Tstar = lines[N - 1].strip().split()[0]
            gt_path = run_dir / "gamma_tails.json"
            run(
                [
                    PY,
                    "tools/core_integral_prover.py",
                    "--T0",
                    Tstar,
                    "--window-config",
                    str(win_path),
                    "--out",
                    str(gt_path),
                    "--dps",
                    f"{dps}",
                ]
            )

            # 7) continuum rollup (normalized continuum_operator_rollup.py CLI)

            # write prime tail flat stub
            pt_flat = run_dir / "prime_tail_envelope.flat.json"
            with open(pt_flat, "w", encoding="utf-8") as f_out:
                json.dump({"prime_tail_norm": prime_tail_val}, f_out, indent=2)

            # band flat stub
            band_flat = run_dir / "band_cert.flat.json"
            band_margin_val = (
                (bc.get("numbers", {}) or {}).get("band_margin_lo")
                or ((bc.get("band_cert", {}) or {}).get("band_margin", {}) or {}).get(
                    "lo"
                )
                or bc.get("band_margin_lo")
                or (
                    bc.get("band_margin")
                    if isinstance(bc.get("band_margin"), str)
                    else None
                )
            )
            with open(band_flat, "w", encoding="utf-8") as f_out:
                json.dump(
                    {
                        "band_margin": band_margin_val,
                        "PASS": bc.get("PASS", False),
                    },
                    f_out,
                    indent=2,
                )

            # grid error: if missing, stub with zero
            grid_err = Path("PROOF_PACKET/grid_error_bound.json")
            if not grid_err.exists():
                with open(grid_err, "w", encoding="utf-8") as f_out:
                    json.dump({"grid_error_norm": "0"}, f_out, indent=2)

            roll_path = run_dir / "continuum_operator_rollup.json"
            run(
                [
                    PY,
                    "tools/continuum_operator_rollup.py",
                    "--band-cert",
                    str(band_flat),
                    "--gamma-tails",
                    str(gt_path),
                    "--prime-block",
                    str(pb_flat),
                    "--prime-tail",
                    str(pt_flat),
                    "--grid-error",
                    str(grid_err),
                    "--out",
                    str(roll_path),
                    "--dps",
                    f"{dps}",
                ]
            )

            R = jload(roll_path)
            lhs = float(R["numbers"]["lhs_total"])
            eps = float(R["numbers"]["epsilon_eff"])
            gap = eps - lhs
            passed = bool(R.get("PASS", False))
            gamma_tail_val = R["numbers"].get("gamma_tails", "0")

            writer.writerow(
                [
                    f"{sigma}",
                    f"{k0}",
                    f"{lhs}",
                    f"{eps}",
                    f"{gap}",
                    passed,
                    f"{band_margin_val}",
                    f"{used}",
                    f"{prime_tail_val}",
                    gamma_tail_val,
                ]
            )
            csv_f.flush()

            if best is None or gap > best["gap"]:
                best = {
                    "sigma": sigma,
                    "k0": k0,
                    "gap": gap,
                    "lhs": lhs,
                    "eps": eps,
                    "passed": passed,
                    "dir": str(run_dir),
                }

            if bar:
                bar.update(1)
            print(
                f"[trial {tag}] PASS={passed}  "
                f"gap={gap:.6e}  lhs={lhs:.6e}  eps={eps:.6e}"
            )

    csv_f.close()
    if bar:
        bar.close()

    if best:
        print(
            "\n[best-coarse] sigma={sigma:.6f}  k0={k0:.6f}  "
            "gap={gap:.6e}  PASS={passed}".format(**best)
        )
        print("See:", best["dir"])
        print(
            "\nTo lock numbers with N=90000, DPS=300, run:"
        )
        print(
            f"  SIGMA={best['sigma']} K0={best['k0']} "
            f"N=90000 DPS=300 ./regenerate_and_run_all_v2.sh"
        )
    else:
        print("No best found (unexpected). Check logs.")


if __name__ == "__main__":
    main()
```
