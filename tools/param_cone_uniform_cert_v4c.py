#!/usr/bin/env python3
# param_cone_uniform_cert_v4c.py
#
# Cone-shaped local search in (sigma, k0) around a candidate point.
#
# Purpose:
#   Sample a cone-shaped neighborhood around (sigma_mid, k0_mid), compute
#   band margins via window_gen + bands_make + band_cert, and compare
#   band_margin_lo to a fixed lhs_total. If any point in the cone has
#   band_margin_lo > lhs_total, the cone certificate PASSes.
#
# CLI (v2.1 normalized core):
#   --sigma-mid     : center sigma
#   --k0-mid        : center k0
#   --sigma-span    : total span in sigma (s_mid ± span/2)
#   --k0-width      : half-width in k0 (k0_mid ± width)
#   --sigma-steps   : number of sigma samples
#   --k0-steps      : number of k0 samples
#   --lhs-total     : fixed lhs_total to compare against (string/number)
#   --dps           : decimal precision (passed to band_cert)
#   --out           : output JSON (cone_uniform_cert)
#   --csv-dir       : directory for ok.csv / fail.csv
#
# Additional tool-specific controls (kept from v4c/v4d logic):
#   --grid, --digits,
#   --critical-left, --critical-right,
#   --inner-left, --inner-right,
#   --jobs, --executor,
#   --timeout-window, --timeout-bands, --timeout-cert,
#   --fname-frac-digits,
#   --reuse, --stop-on-first.
#
# JSON (v2.1 normalized):
#   kind  = "cone_uniform_cert"
#   inputs {
#     sigma_mid, k0_mid, sigma_span, k0_width,
#     sigma_steps, k0_steps, lhs_total
#   }
#   results {
#     PASS,
#     points_evaluated,
#     ok_points,
#     min_gap,
#     witness {
#       sigma,
#       k0,
#       band_margin_lo
#     },
#     stopped_early
#   }
#   PASS : boolean (same as results.PASS)
#   meta {
#     tool        = "param_cone_uniform_cert_v4c",
#     dps,
#     created_utc,
#     sha256
#   }

import os
import sys
import json
import csv
import time
import math
import subprocess
from pathlib import Path
from decimal import Decimal, getcontext
from concurrent.futures import (
    ThreadPoolExecutor,
    ProcessPoolExecutor,
    as_completed,
    CancelledError,
)
import hashlib

getcontext().prec = 200
try:
    from tqdm import tqdm
except Exception:
    tqdm = None

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
TOOLS = os.path.abspath(os.path.dirname(__file__))


# ---------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------

def utc_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_json(path: Path):
    with open(path, "r", encoding="utf-8-sig") as f:
        return json.load(f)


def to_dec(x, default=None):
    if x is None:
        return default
    try:
        return Decimal(str(x))
    except Exception:
        return default


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def newer(a: Path, b: Path) -> bool:
    """Is file a newer than file b?"""
    try:
        return a.stat().st_mtime > b.stat().st_mtime
    except FileNotFoundError:
        return False


def up_to_date(outp: Path, *inputs: Path) -> bool:
    if not outp.exists():
        return False
    try:
        omt = outp.stat().st_mtime
        for i in inputs:
            if i.exists() and i.stat().st_mtime > omt:
                return False
        return True
    except FileNotFoundError:
        return False


def short(s, maxlen=400):
    s = (s or "").strip()
    return s[-maxlen:] if len(s) > maxlen else s


def exe():
    return sys.executable


def tool(pathname: str) -> str:
    return os.path.join(TOOLS, pathname)


def mp_meta_sha(payload: dict, dps: int) -> dict:
    """
    Attach canonical meta {tool, dps, created_utc, sha256}.

    - Assumes payload['meta'] is a dict (created if missing).
    - sha256 is computed over the payload with meta.sha256 removed.
    """
    if "meta" not in payload or not isinstance(payload["meta"], dict):
        payload["meta"] = {}
    meta = payload["meta"]
    meta["tool"] = "param_cone_uniform_cert_v4c"
    meta["dps"] = int(dps)
    meta["created_utc"] = utc_iso()

    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    meta2 = tmp_obj.get("meta")
    if isinstance(meta2, dict):
        meta2.pop("sha256", None)

    blob = json.dumps(
        tmp_obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    digest = hashlib.sha256(blob).hexdigest()
    meta["sha256"] = digest
    return payload


def write_json(path: Path, payload: dict, dps: int) -> None:
    ensure_dir(path.parent)
    payload = mp_meta_sha(payload, dps)
    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------------
# Band margin extraction (unchanged)
# ---------------------------------------------------------------------

def read_band_margin_lo(obj):
    paths = [
        ("band_cert", "band_margin_lo"),
        ("band_cert", "band_margin", "lo"),
        ("numbers", "band_margin"),
        ("band", "margin_lo"),
        ("band_margin_lo",),
    ]
    for path in paths:
        o = obj
        ok = True
        for k in path:
            if not isinstance(o, dict) or k not in o:
                ok = False
                break
            o = o[k]
        if ok:
            d = to_dec(o)
            if d is not None:
                return d
    return None


# ---------------------------------------------------------------------
# Grid + worker helpers
# ---------------------------------------------------------------------

def safe_tag(x, max_frac=6):
    d = Decimal(str(x))
    if d == d.to_integral_value():
        s = format(d, "f")
    else:
        q = d.quantize(Decimal(10) ** -max_frac)
        s = format(q, "f")
        if "." in s:
            s = s.rstrip("0").rstrip(".")
    return s


def spiral_indices(nr, nc):
    """Yield (r,c) indices from center spiraling outward in a rectangle nr×nc."""
    cr = (nr - 1) / 2.0
    cc = (nc - 1) / 2.0
    pts = [(r, c) for r in range(nr) for c in range(nc)]
    pts.sort(
        key=lambda rc: (
            abs(rc[0] - cr) + abs(rc[1] - cc),
            (rc[0] - cr) ** 2 + (rc[1] - cc) ** 2,
        )
    )
    for rc in pts:
        yield rc


def run_cmd(name, cmd, timeout_sec):
    """Run a subprocess with a timeout; return (code, stdout, stderr).

    We keep tqdm output clean by not logging successful calls.
    Only nonzero exits, timeouts, or exceptions are reported, using
    tqdm.write when available so the progress bar is not mangled.
    """
    try:
        p = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=ROOT,
            timeout=timeout_sec,
        )
        if p.returncode != 0:
            msg = f"[cone:{name}] nonzero exit {p.returncode}"
            if tqdm is not None:
                tqdm.write(msg)
            else:
                print(msg, flush=True)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        msg = f"[cone:{name}] TIMEOUT after {timeout_sec}s"
        if tqdm is not None:
            tqdm.write(msg)
        else:
            print(msg, flush=True)
        return -9, "", msg
    except Exception as ex:
        msg = f"[cone:{name}] EXCEPTION: {ex}"
        if tqdm is not None:
            tqdm.write(msg)
        else:
            print(msg, flush=True)
        return -8, "", msg


def worker(args):
    """
    Run one (sigma, k0) point. Returns a dict with status and fields.

    Math/logic:
      - call window_gen.py, bands_make.py, band_cert.py
      - parse band_margin_lo
    """
    (
        work_root,
        grid,
        digits,
        dps,
        crit_l,
        crit_r,
        inner_l,
        inner_r,
        sigma,
        k0,
        t_win,
        t_bands,
        t_cert,
        frac_digits,
        reuse,
    ) = args
    work_root = Path(work_root)

    try:
        s_tag = safe_tag(sigma, frac_digits)
        k_tag = safe_tag(k0, frac_digits)
        cone_dir = work_root
        win_tmp = cone_dir / f"win_{s_tag}_{k_tag}.json"
        bands_tmp = cone_dir / f"bands_{s_tag}_{k_tag}.json"
        cert_tmp = cone_dir / f"band_cert_{s_tag}_{k_tag}.json"
        ensure_dir(win_tmp.parent)

        # window_gen (reuse if requested)
        if not (reuse and win_tmp.exists()):
            c1, o1, e1 = run_cmd(
                "window_gen",
                [
                    exe(),
                    tool("window_gen.py"),
                    "--mode",
                    "gauss",
                    "--sigma",
                    str(sigma),
                    "--k0",
                    str(k0),
                    "--out",
                    str(win_tmp),
                ],
                t_win,
            )
            if c1 != 0:
                return {
                    "ok": False,
                    "sigma": str(sigma),
                    "k0": str(k0),
                    "stage": "window_gen",
                    "code": c1,
                    "stderr": short(e1),
                }

        # bands_make (reuse if up-to-date wrt window)
        if not (reuse and up_to_date(bands_tmp, win_tmp)):
            c2, o2, e2 = run_cmd(
                "bands_make",
                [
                    exe(),
                    tool("bands_make.py"),
                    "--window-config",
                    str(win_tmp),
                    "--grid",
                    str(grid),
                    "--digits",
                    str(digits),
                    "--critical-left",
                    str(crit_l),
                    "--critical-right",
                    str(crit_r),
                    "--inner-left",
                    str(inner_l),
                    "--inner-right",
                    str(inner_r),
                    "--out",
                    str(bands_tmp),
                ],
                t_bands,
            )
            if c2 != 0:
                return {
                    "ok": False,
                    "sigma": str(sigma),
                    "k0": str(k0),
                    "stage": "bands_make",
                    "code": c2,
                    "stderr": short(e2),
                }

        # band_cert (reuse if up-to-date wrt window & bands)
        if not (reuse and up_to_date(cert_tmp, win_tmp, bands_tmp)):
            c3, o3, e3 = run_cmd(
                "band_cert",
                [
                    exe(),
                    tool("band_cert.py"),
                    "--bands",
                    str(bands_tmp),
                    "--window-config",
                    str(win_tmp),
                    "--out",
                    str(cert_tmp),
                    "--dps",
                    str(dps),
                ],
                t_cert,
            )
            if c3 != 0 or not cert_tmp.exists():
                return {
                    "ok": False,
                    "sigma": str(sigma),
                    "k0": str(k0),
                    "stage": "band_cert",
                    "code": c3,
                    "stderr": short(e3 or "no-output"),
                }

        cert_obj = load_json(cert_tmp)
        m_lo = read_band_margin_lo(cert_obj)
        if m_lo is None:
            return {
                "ok": False,
                "sigma": str(sigma),
                "k0": str(k0),
                "stage": "parse_margin",
                "code": 0,
                "stderr": "band_margin_lo missing",
            }

        return {
            "ok": True,
            "sigma": str(sigma),
            "k0": str(k0),
            "band_margin_lo": str(m_lo),
        }

    except Exception as ex:
        return {
            "ok": False,
            "sigma": str(sigma),
            "k0": str(k0),
            "stage": "exception",
            "code": 0,
            "stderr": short(str(ex)),
        }


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Cone-shaped local search in (sigma, k0) around a candidate point."
    )

    # v2.1 canonical arguments for this tool
    ap.add_argument("--sigma-mid", required=True, help="Center sigma.")
    ap.add_argument("--sigma-span", type=float, default=0.25,
                    help="Total sigma span; samples in [sigma_mid-span/2, sigma_mid+span/2].")
    ap.add_argument("--sigma-steps", type=int, default=25,
                    help="Number of sigma grid points.")

    ap.add_argument("--k0-mid", required=True, help="Center k0.")
    ap.add_argument("--k0-width", type=float, default=0.02,
                    help="Half-width in k0; samples in [k0_mid-width, k0_mid+width].")
    ap.add_argument("--k0-steps", type=int, default=25,
                    help="Number of k0 grid points.")

    ap.add_argument(
        "--lhs-total",
        required=True,
        help="Fixed lhs_total value to compare band_margin_lo against.",
    )

    ap.add_argument("--grid", type=int, default=6000, help="bands_make grid.")
    ap.add_argument("--digits", type=int, default=120, help="bands_make digits.")
    ap.add_argument("--dps", type=int, default=220, help="band_cert precision.")

    ap.add_argument("--critical-left", required=True, help="Critical line left endpoint.")
    ap.add_argument("--critical-right", required=True, help="Critical line right endpoint.")
    ap.add_argument("--inner-left", required=True, help="Inner line left endpoint.")
    ap.add_argument("--inner-right", required=True, help="Inner line right endpoint.")

    ap.add_argument(
        "--jobs",
        type=int,
        default=max(1, (os.cpu_count() or 4) - 1),
        help="Number of worker processes/threads.",
    )
    ap.add_argument(
        "--executor",
        choices=["thread", "process"],
        default="process",
        help="Executor type for parallelism.",
    )
    ap.add_argument("--csv-dir", required=True, help="Directory for ok.csv and fail.csv.")
    ap.add_argument("--out", required=True, help="Output cone_uniform_cert JSON path.")

    ap.add_argument("--timeout-window", type=int, default=600)
    ap.add_argument("--timeout-bands", type=int, default=1800)
    ap.add_argument("--timeout-cert", type=int, default=3600)
    ap.add_argument("--fname-frac-digits", type=int, default=6)

    ap.add_argument(
        "--reuse",
        action="store_true",
        help="Reuse stage outputs when present and up-to-date.",
    )
    ap.add_argument(
        "--stop-on-first",
        action="store_true",
        help="Stop sweep when first positive gap is found.",
    )

    args = ap.parse_args()

    # Tool sanity check
    needed = ["window_gen.py", "bands_make.py", "band_cert.py"]
    missing = [t for t in needed if not Path(tool(t)).exists()]
    if missing:
        print(f"[cone] ERROR: missing tools: {missing}", file=sys.stderr)
        sys.exit(2)

    csv_dir = Path(args.csv_dir)
    ensure_dir(csv_dir)

    fail_csv = csv_dir / "fail.csv"
    ok_csv = csv_dir / "ok.csv"

    # lhs_total is now an explicit CLI argument (canonical v2.1)
    lhs_total = to_dec(args.lhs_total)
    if lhs_total is None:
        print("[cone] ERROR: could not parse lhs_total", file=sys.stderr)
        sys.exit(2)

    sigma_mid = Decimal(str(args.sigma_mid))
    k0_mid = Decimal(str(args.k0_mid))

    s_lo = sigma_mid - Decimal(str(args.sigma_span)) / 2
    s_hi = sigma_mid + Decimal(str(args.sigma_span)) / 2
    k_lo = k0_mid - Decimal(str(args.k0_width))
    k_hi = k0_mid + Decimal(str(args.k0_width))

    # build grids
    if args.sigma_steps < 1 or args.k0_steps < 1:
        print("[cone] ERROR: sigma_steps and k0_steps must be >= 1", file=sys.stderr)
        sys.exit(2)

    sigmas = [
        s_lo + (s_hi - s_lo) * Decimal(i) / (args.sigma_steps - 1 or 1)
        for i in range(args.sigma_steps)
    ]
    k0s = [
        k_lo + (k_hi - k_lo) * Decimal(j) / (args.k0_steps - 1 or 1)
        for j in range(args.k0_steps)
    ]

    # spiral order (center-first)
    grid_points = []
    for (ri, ci) in spiral_indices(len(sigmas), len(k0s)):
        grid_points.append((str(sigmas[ri]), str(k0s[ci])))

    # CSV writers (v2.1 headers)
    f_fail = open(fail_csv, "w", newline="", encoding="utf-8")
    w_fail = csv.writer(f_fail)
    w_fail.writerow(["sigma", "k0", "stage", "code", "stderr_tail"])

    f_ok = open(ok_csv, "w", newline="", encoding="utf-8")
    w_ok = csv.writer(f_ok)
    w_ok.writerow(["sigma", "k0", "band_margin_lo", "lhs_total", "gap"])

    work_root = csv_dir / "cone"
    ensure_dir(work_root)

    task_args = [
        (
            str(work_root),
            args.grid,
            args.digits,
            args.dps,
            str(args.critical_left),
            str(args.critical_right),
            str(args.inner_left),
            str(args.inner_right),
            s,
            k,
            int(args.timeout_window),
            int(args.timeout_bands),
            int(args.timeout_cert),
            int(args.fname_frac_digits),
            bool(args.reuse),
        )
        for (s, k) in grid_points
    ]

    Executor = ThreadPoolExecutor if args.executor == "thread" else ProcessPoolExecutor
    successes = []
    total = len(task_args)

    # Print config first, then create tqdm so the bar appears after these lines
    print(f"[cone] ROOT={ROOT}")
    print(f"[cone] TOOLS={TOOLS}")
    print(
        f"[cone] jobs={args.jobs} executor={args.executor} points={total}",
        flush=True,
    )
    print(
        f"[cone] grid={args.grid} digits={args.digits} dps={args.dps}",
        flush=True,
    )
    print(
        f"[cone] timeouts: win={args.timeout_window}s "
        f"bands={args.timeout_bands}s cert={args.timeout_cert}s",
        flush=True,
    )

    pbar = tqdm(total=total, desc="cone grid") if tqdm else None

    cancelled = False
    try:
        with Executor(max_workers=args.jobs) as ex:
            futs = [ex.submit(worker, a) for a in task_args]
            for fu in as_completed(futs):
                try:
                    res = fu.result()
                except CancelledError:
                    break
                if pbar:
                    pbar.update(1)
                if res.get("ok"):
                    m_lo = to_dec(res["band_margin_lo"])
                    gap = (
                        m_lo - lhs_total
                        if (m_lo is not None and lhs_total is not None)
                        else None
                    )
                    if gap is not None and gap > 0:
                        successes.append((res["sigma"], res["k0"], m_lo, gap))
                        w_ok.writerow(
                            [res["sigma"], res["k0"], str(m_lo), str(lhs_total), str(gap)]
                        )
                        if args.stop_on_first:
                            cancelled = True
                            for f in futs:
                                f.cancel()
                            break
                else:
                    w_fail.writerow(
                        [
                            res.get("sigma", ""),
                            res.get("k0", ""),
                            res.get("stage", ""),
                            res.get("code", ""),
                            res.get("stderr", ""),
                        ]
                    )
    finally:
        if pbar:
            pbar.close()
        f_ok.close()
        f_fail.close()

    out = Path(args.out)

    if successes:
        # minimal positive gap
        s, k, m_lo, gap = min(successes, key=lambda t: t[3])
        results_block = {
            "PASS": True,
            "points_evaluated": total,
            "ok_points": len(successes),
            "min_gap": str(gap),
            "witness": {
                "sigma": str(s),
                "k0": str(k),
                "band_margin_lo": str(m_lo),
            },
            "stopped_early": bool(cancelled and args.stop_on_first),
        }
        PASS = True
    else:
        results_block = {
            "PASS": False,
            "points_evaluated": total,
            "ok_points": 0,
        }
        PASS = False

    payload = {
        "kind": "cone_uniform_cert",
        "inputs": {
            "sigma_mid": str(sigma_mid),
            "k0_mid": str(k0_mid),
            "sigma_span": str(args.sigma_span),
            "k0_width": str(args.k0_width),
            "sigma_steps": int(args.sigma_steps),
            "k0_steps": int(args.k0_steps),
            "lhs_total": str(lhs_total),
        },
        "results": results_block,
        "PASS": PASS,
        "meta": {},
    }

    write_json(out, payload, args.dps)  # or 'cert' / 'payload' if that’s the name in your script

    print(
        f"[cone] wrote {out}  PASS={payload['results']['PASS']}  "
        f"ok={payload['results'].get('ok_points', 0)}"
    )


if __name__ == "__main__":
    main()
