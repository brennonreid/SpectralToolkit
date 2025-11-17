#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bands_make.py — exact bands/grids writer (proof-grade I/O, no logic change).

Purpose:
  Create exact bands and uniform grids for the spectral window, with
  proof-grade I/O:
    - exact decimal nodes using Python's Decimal
    - deterministic JSON layout and sha256
    - explicit linkage to the source window config

CLI (normalized v2.1):
  --window-config    : path to window.json (kind="window")
  --out              : output bands JSON (e.g., packs/rh/inputs/auto_bands.json)
  --dps / --digits   : decimal digits for endpoints/steps/nodes (default 80)
  --grid             : uniform grid size per band (>= 2, default 6000)
  --critical-left    : left endpoint of critical band (decimal string)
  --critical-right   : right endpoint of critical band (decimal string)
  --inner-left       : optional inner band left
  --inner-right      : optional inner band right
  --outer-left       : optional outer band left
  --outer-right      : optional outer band right

JSON (v2.1-style, non-certificate artifact):
  {
    "kind": "bands",
    "version": "1.2",
    "source_window": {
      "window_config_path": "...",
      "mode": "gauss" | ...
    },
    "bands": [
      { "label": "critical", "left": "<string>", "right": "<string>" },
      ...
    ],
    "named_grids": {
      "critical": {
        "left": "<string>",
        "right": "<string>",
        "grid": <int>,
        "h": "<string>",
        "nodes": ["<string>", ...]
      },
      ...
    },
    "critical_left": "<string>",
    "critical_right": "<string>",
    "grid_N": <int>,
    "bands_count": <int>,
    "meta": {
      "tool": "bands_make",
      "dps": <int>,
      "grid": <int>,
      "python": "...",
      "os": "...",
      "cpu": "...",
      "workers_detected": <int>,
      "runtime_sec": "<string>",
      "created_utc": "<ISO8601 UTC>",
      "sha256": "<digest>"
    }
  }

All band endpoints, steps, and nodes are decimal strings. Integers such as
grid size and bands_count remain JSON integers.
"""

import argparse
import json
import hashlib
import os
import sys
import time
import platform
import multiprocessing
from pathlib import Path
from decimal import Decimal as D, getcontext
from datetime import datetime, timezone


def now_utc_iso() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def sha256_canonical(obj) -> str:
    """
    Compute sha256 over a canonical JSON serialization (sorted keys, no
    pretty-printing). Used for meta.sha256.
    """
    blob = json.dumps(
        obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def write_json(path, payload):
    """
    Write JSON with newline-terminated UTF-8 and sha256 in meta.sha256.
    The hash is computed over the payload with any existing meta.sha256
    field removed.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    # Ensure meta exists
    meta = payload.get("meta")
    if not isinstance(meta, dict):
        meta = {}
        payload["meta"] = meta

    # Compute hash on a copy without meta.sha256
    tmp_obj = json.loads(json.dumps(payload, ensure_ascii=False))
    if "meta" in tmp_obj and isinstance(tmp_obj["meta"], dict):
        tmp_obj["meta"].pop("sha256", None)
    digest = sha256_canonical(tmp_obj)
    meta["sha256"] = digest

    tmp = str(p) + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, p)


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def dec(s):
    if isinstance(s, str):
        return D(s)
    return D(str(s))


def str_dec(x, ndigits):
    q = D(10) ** (-ndigits)
    return str(x.quantize(q)) if ndigits >= 0 else str(x)


def exact_step(left, right, grid):
    if grid < 2:
        raise ValueError("grid must be >= 2")
    return (right - left) / D(grid - 1)


def build_nodes(left, right, grid, ndigits):
    h = exact_step(left, right, grid)
    nodes = [str_dec(left + h * D(i), ndigits) for i in range(grid)]
    left_check = left + h * D(grid - 1)
    if left_check != right:
        raise RuntimeError(
            f"grid mismatch: left + h*(n-1) != right (got {left_check} vs {right})"
        )
    return str_dec(h, ndigits), nodes


def capture_meta(dps_val, grid):
    return {
        "python": platform.python_version(),
        "os": platform.platform(),
        "cpu": platform.processor(),
        "workers_detected": multiprocessing.cpu_count(),
        "dps": dps_val,
        "grid": grid,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Create exact bands and grids from window.json (proof-grade IO)."
    )
    ap.add_argument(
        "--window-config",
        required=True,
        help="Path to window.json (kind='window').",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON (e.g., packs/rh/inputs/auto_bands.json).",
    )
    # Canonical flag is --dps, with --digits as a backwards-compatible alias
    ap.add_argument(
        "--dps",
        "--digits",
        dest="dps",
        type=int,
        default=80,
        help="Decimal digits to store for endpoints/steps/nodes. Default 80.",
    )
    ap.add_argument(
        "--grid",
        type=int,
        default=6000,
        help="Uniform grid size per band. Default 6000.",
    )
    ap.add_argument(
        "--critical-left",
        type=str,
        required=True,
        help="Left endpoint of critical band (decimal string).",
    )
    ap.add_argument(
        "--critical-right",
        type=str,
        required=True,
        help="Right endpoint of critical band (decimal string).",
    )
    ap.add_argument(
        "--inner-left",
        type=str,
        default=None,
        help="Optional inner band left.",
    )
    ap.add_argument(
        "--inner-right",
        type=str,
        default=None,
        help="Optional inner band right.",
    )
    ap.add_argument(
        "--outer-left",
        type=str,
        default=None,
        help="Optional outer band left.",
    )
    ap.add_argument(
        "--outer-right",
        type=str,
        default=None,
        help="Optional outer band right.",
    )

    args = ap.parse_args()

    dps_val = int(args.dps)
    grid = int(args.grid)
    getcontext().prec = max(100, dps_val + 20)
    t0 = time.time()

    # Sanity on window artifact; used only for traceability.
    win = read_json(args.window_config)
    if "window" in win and isinstance(win["window"], dict):
        win = {
            "kind": "window",
            "mode": win["window"].get("mode"),
            "params": {k: v for k, v in win["window"].items() if k not in ("mode",)},
        }
    if win.get("kind") != "window":
        raise RuntimeError(
            "window_config does not look like a SpectralToolkit window artifact (expected kind='window')."
        )

    named_specs = []
    for label, Ls, Rs in [
        ("critical", args.critical_left, args.critical_right),
        ("inner", args.inner_left, args.inner_right),
        ("outer", args.outer_left, args.outer_right),
    ]:
        if Ls is None or Rs is None:
            continue
        L = dec(Ls)
        R = dec(Rs)
        if not (L < R):
            raise ValueError(f"{label} band requires left < right, got {L} >= {R}")
        named_specs.append((label, L, R))

    if not named_specs:
        raise RuntimeError(
            "no bands were specified; at least one band (e.g., critical) is required"
        )

    grids = {}
    flat_bands = []
    for label, L, R in named_specs:
        h_str, nodes = build_nodes(L, R, grid, dps_val)
        grids[label] = {
            "left": str_dec(L, dps_val),
            "right": str_dec(R, dps_val),
            "grid": grid,
            "h": h_str,
            "nodes": nodes,
        }
        flat_bands.append(
            {
                "label": label,
                "left": str_dec(L, dps_val),
                "right": str_dec(R, dps_val),
            }
        )

    critical_left = next(
        (b["left"] for b in flat_bands if b["label"] == "critical"), None
    )
    critical_right = next(
        (b["right"] for b in flat_bands if b["label"] == "critical"), None
    )

    meta = capture_meta(dps_val, grid)
    meta["tool"] = "bands_make"
    meta["created_utc"] = now_utc_iso()
    meta["runtime_sec"] = str(round(time.time() - t0, 6))

    payload = {
        "kind": "bands",
        "version": "1.2",
        "source_window": {
        "window_config_path": args.window_config,
        "mode": win.get("mode", None),
        },
        "bands": flat_bands,
        "named_grids": grids,
        "critical_left": critical_left,
        "critical_right": critical_right,
        "grid_N": grid,
        "bands_count": len(flat_bands),
        "meta": meta,
    }

    write_json(args.out, payload)
    dt = time.time() - t0
    print(
        f"[OK] wrote {args.out} sha256={payload['meta']['sha256'][:16]}... in {dt:.2f}s"
    )
    print(
        f"[info] bands: {', '.join(b['label'] for b in flat_bands)} "
        f"(grid={grid}, dps={dps_val})"
    )


if __name__ == "__main__":
    main()
