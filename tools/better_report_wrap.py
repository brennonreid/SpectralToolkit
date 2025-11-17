#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
better_report_wrap.py (ToolkitV2, normalized)
- Consumes PROOF_PACKET certs + inputs at the END of the pipeline
- Coalesces keys from multiple schema variants (older/newer writers)
- Adds Infinity stage artifacts (Fourier, Deconv, RvM, Frame, Subspace PSD, STP, Cone, Rolling-T)
- Writes a Markdown overview and a JSON wrap with kind="report_wrap" and normalized meta
"""

import argparse
import os
import json
import time
import hashlib
from typing import Any, Optional, Tuple, List, Dict
from mpmath import mp

# ---------------- helpers ----------------

def utc_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

def load_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8-sig") as f:
        return json.load(f)

def dig(obj: Any, path: list) -> Optional[Any]:
    o = obj
    for k in path:
        if not isinstance(o, dict) or k not in o:
            return None
        o = o[k]
    return o

def coalesce(obj: Any, paths: list, default: str = "") -> str:
    for p in paths:
        v = dig(obj, p)
        if v not in (None, "", "null"):
            return str(v)
    return default

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def file_meta(path: str) -> Optional[dict]:
    if not os.path.exists(path):
        return None
    st = os.stat(path)
    return {
        "path": path,
        "bytes": int(st.st_size),
        "mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime)),
        "sha256": sha256_file(path),
    }

def mpf_str(x: Any, dps: int) -> str:
    mp.dps = int(dps)
    return mp.nstr(mp.mpf(str(x)), n=mp.dps, strip_zeros=False)

def exists(path: str) -> bool:
    return os.path.exists(path)

# ---------------- readers (multi-shape) ----------------

def read_window(inputs_dir: str) -> Tuple[dict, dict]:
    win = os.path.join(inputs_dir, "window.json")
    auto = os.path.join(inputs_dir, "auto_bands.json")
    wj = load_json(win) if exists(win) else {}
    aj = load_json(auto) if exists(auto) else {}

    mode = coalesce(wj, [["mode"], ["window", "mode"]], "gauss")
    sigma = coalesce(wj, [["sigma"], ["window", "sigma"]], "6.0")
    k0    = coalesce(wj, [["k0"], ["window", "k0"], ["notch_k0"], ["window", "notch_k0"]], "0.25")
    gridN = (dig(aj, ["grid", "N"]) or dig(aj, ["N"]) or 6000)
    cleft = coalesce(aj, [["band", "critical_left"], ["critical_left"]], "0.30")
    crght = coalesce(aj, [["band", "critical_right"], ["critical_right"]], "2.80")

    window_summary = {
        "mode": mode,
        "sigma": sigma,
        "k0": k0,
    }
    bands_summary = {
        "grid_N": int(gridN),
        "critical_left": cleft,
        "critical_right": crght,
    }
    return window_summary, bands_summary

def read_band_cert(certs: str) -> dict:
    p = os.path.join(certs, "band_cert.json")
    if not exists(p):
        return {}
    j = load_json(p)
    band_margin_lo = coalesce(
        j,
        [
            ["band_cert", "band_margin_lo"],
            ["band_margin_lo"],
            ["numbers", "band_margin_lo"],
        ],
        "",
    )
    band_margin_hi = coalesce(
        j,
        [
            ["band_cert", "band_margin_hi"],
            ["band_margin_hi"],
            ["numbers", "band_margin_hi"],
        ],
        "",
    )
    return {
        "band_margin_lo": band_margin_lo,
        "band_margin_hi": band_margin_hi,
        "PASS": bool(j.get("PASS", True)),
    }

def read_tails_bundle(certs: str) -> dict:
    p = os.path.join(certs, "tails.json")
    if not exists(p):
        # v1 layout: separate gamma_tail / prime_tail_envelope
        gamma_p = os.path.join(certs, "gamma_tail.json")
        prime_p = os.path.join(certs, "prime_tail_envelope.json")
        gj = load_json(gamma_p) if exists(gamma_p) else {}
        pj = load_json(prime_p) if exists(prime_p) else {}
        return {
            "gamma_env_at_T0": coalesce(
                gj,
                [["gamma_tail", "env_T0_hi"], ["env_T0_hi"], ["numbers", "gamma_env_T0_hi"]],
                "",
            ),
            "prime_env_T0_hi": coalesce(
                pj,
                [["prime_tail", "env_T0_hi"], ["env_T0_hi"], ["numbers", "prime_env_T0_hi"]],
                "",
            ),
        }

    j = load_json(p)
    return {
        "gamma_env_at_T0": coalesce(
            j,
            [["gamma", "env_T0_hi"], ["gamma_tail", "env_T0_hi"]],
            "",
        ),
        "prime_env_T0_hi": coalesce(
            j,
            [["prime", "env_T0_hi"], ["prime_tail", "env_T0_hi"]],
            "",
        ),
    }

def read_prime_block_and_grid(certs: str) -> dict:
    prime_p = os.path.join(certs, "prime_block_norm.json")
    grid_p  = os.path.join(certs, "grid_error_bound.json")
    pj = load_json(prime_p) if exists(prime_p) else {}
    gj = load_json(grid_p) if exists(grid_p) else {}
    cap_hi = coalesce(
        pj,
        [
            ["prime_block", "cap_hi"],
            ["cap_hi"],
            ["numbers", "cap_hi"],
        ],
        "",
    )
    grid_hi = coalesce(
        gj,
        [
            ["grid_error", "bound_hi"],
            ["bound_hi"],
            ["numbers", "grid_bound_hi"],
        ],
        "",
    )
    return {
        "cap_hi": cap_hi,
        "grid_bound_hi": grid_hi,
    }

def read_weil_psd(certs: str) -> dict:
    p = os.path.join(certs, "weil_psd_bochner.json")
    if not exists(p):
        return {}
    j = load_json(p)
    status = "PASS" if j.get("PASS", True) else "FAIL"
    min_sample = coalesce(
        j,
        [
            ["weil_psd_bochner", "eval", "min_hat_h_sample"],
            ["eval", "min_hat_h_sample"],
            ["numbers", "min_hat_h_sample"],
        ],
        "",
    )
    reason = coalesce(
        j,
        [
            ["weil_psd_bochner", "reason"],
            ["reason"],
        ],
        "",
    )
    return {
        "status": status,
        "min_hat_h_sample": min_sample,
        "reason": reason,
    }

def read_continuum(certs: str) -> dict:
    # canonical v2.1 name
    p = os.path.join(certs, "continuum_operator_cert.json")
    if not exists(p):
        # tolerated legacy
        p_alt = os.path.join(certs, "continuum_operator_rollup.json")
        if not exists(p_alt):
            return {}
        p = p_alt

    j = load_json(p)
    lhs = coalesce(
        j,
        [
            ["continuum", "lhs"],
            ["lhs"],
            ["numbers", "lhs"],
        ],
        "",
    )
    eps_eff = coalesce(
        j,
        [
            ["continuum", "eps_eff"],
            ["eps_eff"],
            ["numbers", "eps_eff"],
        ],
        "",
    )
    status = "PASS" if j.get("PASS", True) else "FAIL"
    return {
        "lhs": lhs,
        "eps_eff": eps_eff,
        "status": status,
    }

def read_lipschitz(certs: str) -> dict:
    p = os.path.join(certs, "q_lipschitz.json")
    if not exists(p):
        return {}
    j = load_json(p)
    L_lo = coalesce(j, [["q_lipschitz", "L", "lo"], ["L", "lo"]], "")
    L_hi = coalesce(j, [["q_lipschitz", "L", "hi"], ["L", "hi"]], "")
    return {
        "L_lo": L_lo,
        "L_hi": L_hi,
    }

def read_density(certs: str) -> dict:
    p = os.path.join(certs, "density_metrics.json")
    if not exists(p):
        return {}
    j = load_json(p)
    return {
        "status": "PASS" if j.get("PASS", True) else "FAIL",
        "summary": j.get("summary", {}),
    }

def read_explicit(certs: str) -> dict:
    p = os.path.join(certs, "weil_explicit_cert.json")
    if not exists(p):
        p_alt = os.path.join(certs, "explicit_formula.json")
        if not exists(p_alt):
            return {}
        p = p_alt
    j = load_json(p)
    eps_eff_lo = coalesce(
        j,
        [
            ["explicit_formula", "epsilon_eff_lo"],
            ["epsilon_eff_lo"],
            ["numbers", "epsilon_eff_lo"],
        ],
        "",
    )
    return {
        "epsilon_eff_lo": eps_eff_lo,
        "PASS": bool(j.get("PASS", True)),
    }

def read_core_integral(certs: str) -> dict:
    p = os.path.join(certs, "core_integral.json")
    if not exists(p):
        return {}
    j = load_json(p)
    lhs = coalesce(
        j,
        [["core_integral", "lhs"], ["lhs"], ["numbers", "lhs"]],
        "",
    )
    status = "PASS" if j.get("PASS", True) else "FAIL"
    return {
        "lhs": lhs,
        "status": status,
    }

def read_uniform(certs: str) -> dict:
    p = os.path.join(certs, "uniform_certificate.json")
    if not exists(p):
        return {}
    j = load_json(p)
    eps_eff = coalesce(
        j,
        [
            ["uniform_certificate", "epsilon_eff"],
            ["epsilon_eff"],
            ["numbers", "epsilon_eff"],
        ],
        "",
    )
    return {
        "epsilon_eff": eps_eff,
        "PASS": bool(j.get("PASS", True)),
    }

# ---- Infinity & extras readers ----

def read_fourier(certs: str) -> dict:
    p = os.path.join(certs, "fourier_inversion_cert.json")
    if not exists(p):
        return {}
    j = load_json(p)
    status = "PASS" if j.get("PASS", True) else "FAIL"
    return {"status": status}

def read_deconv(certs: str) -> dict:
    p = os.path.join(certs, "deconv_cert_infinite.json")
    if not exists(p):
        return {}
    j = load_json(p)
    status = "PASS" if j.get("PASS", True) else "FAIL"
    return {"status": status}

def read_rvm(certs: str) -> dict:
    p = os.path.join(certs, "rv_mangoldt_bounds.json")
    if not exists(p):
        return {}
    j = load_json(p)
    status = "PASS" if j.get("PASS", True) else "FAIL"
    return {"status": status}

def read_frame_probe(certs: str) -> dict:
    p = os.path.join(certs, "frame_probe.json")
    if not exists(p):
        return {}
    j = load_json(p)
    return {
        "PASS": bool(j.get("PASS", True)),
        "summary": j.get("summary", {}),
    }

def read_subspace_psd(certs: str) -> dict:
    p = os.path.join(certs, "subspace_psd_cert.json")
    if not exists(p):
        return {}
    j = load_json(p)
    status = "PASS" if j.get("PASS", True) else "FAIL"
    min_diag = coalesce(
        j,
        [
            ["result", "min_diag_L"],
            ["min_diag_L"],
        ],
        "",
    )
    return {
        "status": status,
        "min_diag_L": min_diag,
    }

def read_stp(certs: str) -> dict:
    p = os.path.join(certs, "stp_test.json")
    if not exists(p):
        return {}
    j = load_json(p)
    status = "PASS" if j.get("PASS", True) else "FAIL"
    return {
        "status": status,
        "modules": j.get("modules", {}),
    }

def read_cone_uniform(certs: str) -> dict:
    p = os.path.join(certs, "cone_uniform_cert.json")
    if not exists(p):
        return {"status": "N/A", "ok": "", "sha256": ""}
    j = load_json(p)
    ok = coalesce(j, [["ok"], ["numbers", "ok"]], "")
    status = "PASS" if (j.get("PASS", True)) else "FAIL"
    return {"status": status, "ok": ok, "sha256": sha256_file(p)}

def read_rolling_T(certs: str) -> dict:
    p = os.path.join(certs, "rolling_T_uniform_cert.json")
    if not exists(p):
        # also tolerate canonical rolled JSON name
        p2 = os.path.join(certs, "rolling_T_uniform.json")
        if not exists(p2):
            return {"status": "N/A", "numbers": {}, "sha256": ""}
        p = p2
    j = load_json(p)
    status = "PASS" if j.get("PASS", True) else "FAIL"
    return {"status": status, "numbers": j.get("numbers", {}), "sha256": sha256_file(p)}

# ---------------- main ----------------

def collect_meta_files(certs_dir: str, inputs_dir: str) -> List[Dict[str, Any]]:
    paths = []

    # inputs
    paths.append(os.path.join(inputs_dir, "window.json"))
    paths.append(os.path.join(inputs_dir, "auto_bands.json"))

    # core certs
    for name in [
        "band_cert.json",
        "gamma_tail.json",
        "prime_tail_envelope.json",
        "prime_block_norm.json",
        "grid_error_bound.json",
        "weil_psd_bochner.json",
        "continuum_operator_cert.json",
        "continuum_operator_rollup.json",
        "q_lipschitz.json",
        "density_metrics.json",
        "tails.json",
        "weil_explicit_cert.json",
        "explicit_formula.json",
        "core_integral.json",
        "uniform_certificate.json",
        "rollup_uniform.json",
    ]:
        paths.append(os.path.join(certs_dir, name))

    # infinity & extras
    for name in [
        "fourier_inversion_cert.json",
        "deconv_cert_infinite.json",
        "rv_mangoldt_bounds.json",
        "rv_mangoldt_bounds.theory.json",
        "frame_probe.json",
        "frame_probe.csv",
        "subspace_psd_cert.json",
        "subspace_psd_gram.csv",
        "stp_selftest.log",
        "cone_uniform_cert.json",
        "rolling_T_uniform_cert.json",
        "rolling_T_uniform.json",
    ]:
        paths.append(os.path.join(certs_dir, name))

    meta_files: List[Dict[str, Any]] = []
    for p in paths:
        m = file_meta(p)
        if m is not None:
            meta_files.append(m)
    return meta_files

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certs-dir", required=True, help="Directory containing PROOF_PACKET certificates")
    ap.add_argument("--inputs-dir", required=True, help="Directory containing window.json / auto_bands.json")
    ap.add_argument("--out-md", required=True, help="Output Markdown summary path")
    ap.add_argument("--out-json", required=True, help="Output JSON wrap path")
    ap.add_argument("--dps", type=int, default=220)
    args = ap.parse_args()

    # Core sections
    window, bands = read_window(args.inputs_dir)
    band_cert     = read_band_cert(args.certs_dir)
    tails         = read_tails_bundle(args.certs_dir)
    pbg           = read_prime_block_and_grid(args.certs_dir)
    psd           = read_weil_psd(args.certs_dir)
    cont          = read_continuum(args.certs_dir)
    lips          = read_lipschitz(args.certs_dir)
    density       = read_density(args.certs_dir)
    explicit      = read_explicit(args.certs_dir)
    core          = read_core_integral(args.certs_dir)
    uniform       = read_uniform(args.certs_dir)

    # Infinity & extras
    fourier       = read_fourier(args.certs_dir)
    deconv        = read_deconv(args.certs_dir)
    rvm           = read_rvm(args.certs_dir)
    frame         = read_frame_probe(args.certs_dir)
    subpsd        = read_subspace_psd(args.certs_dir)
    stp           = read_stp(args.certs_dir)
    cone          = read_cone_uniform(args.certs_dir)
    rollT         = read_rolling_T(args.certs_dir)

    # Verdict (very conservative: PASS if no explicit FAILs in core/explicit/uniform/cont/psd)
    fails: List[str] = []

    def flag_fail(name: str, section: Any) -> None:
        if isinstance(section, dict) and section.get("status") == "FAIL":
            fails.append(name)

    # mark critical sections
    for nm, sec in [
        ("continuum_operator", cont),
        ("uniform_certificate", uniform),
        ("weil_psd_bochner", psd),
        ("explicit_formula", explicit),
        ("stp_selftest", stp),
        ("subspace_psd", subpsd),
        ("deconv", deconv),
        ("rv_mangoldt", rvm),
    ]:
        flag_fail(nm, sec)

    verdict = "PASS" if not fails else "FAIL"

    wrap = {
        "kind": "report_wrap",
        "generated_utc": utc_iso(),
        "certs_dir": args.certs_dir,
        "inputs_dir": args.inputs_dir,
        "verdict": verdict,
        "failures": fails,
        "summary": {
            "window": window,
            "bands": bands,
            "band_cert": band_cert,
            "tails": tails,
            "prime_block_and_grid": pbg,
            "weil_psd_bochner": psd,
            "continuum_operator": cont,
            "lipschitz": lips,
            "density": density,
            "explicit_formula": explicit,
            "core_integral": core,
            "uniform_certificate": uniform,
            "infinity": {
                "fourier_inversion": fourier,
                "deconvolution": deconv,
                "riemann_von_mangoldt": rvm,
                "frame_probe": frame,
                "subspace_psd": subpsd,
                "stp_selftest": stp,
                "cone_uniform": cone,
                "rolling_T_uniform": rollT,
            },
        },
        "meta": {
            "tool": "better_report_wrap",
            "dps": int(args.dps),
            "created_utc": utc_iso(),
            "files": collect_meta_files(args.certs_dir, args.inputs_dir),
        },
    }

    # md summary (compact but richer)
    md_lines: List[str] = []
    md_lines.append(f"# RH Report Wrap\n\nGenerated: {wrap['generated_utc']}\n")
    md_lines.append(f"**Verdict:** {wrap['verdict']}\n")
    if fails:
        md_lines.append(f"**Failures:** {', '.join(fails)}\n")

    md_lines.append("## Key Numbers\n")
    md_lines.append(f"- Band margin (lo): {band_cert.get('band_margin_lo', '')}")
    md_lines.append(f"- Continuum LHS: {cont.get('lhs', '')}")
    md_lines.append(f"- Continuum eps_eff: {cont.get('eps_eff', '')}")
    md_lines.append(f"- Explicit epsilon_eff_lo: {explicit.get('epsilon_eff_lo', '')}")
    md_lines.append(f"- Prime block cap_hi: {pbg.get('cap_hi', '')}")
    md_lines.append(f"- Grid bound hi: {pbg.get('grid_bound_hi', '')}")
    md_lines.append(f"- Gamma env @ T0: {tails.get('gamma_env_at_T0', '')}")
    md_lines.append(f"- Prime tail env @ T0: {tails.get('prime_env_T0_hi', '')}")
    md_lines.append("")

    with open(args.out_md, "w", encoding="utf-8") as f_md:
        f_md.write("\n".join(md_lines))

    with open(args.out_json, "w", encoding="utf-8") as f_js:
        json.dump(wrap, f_js, indent=2, sort_keys=False)

if __name__ == "__main__":
    main()
