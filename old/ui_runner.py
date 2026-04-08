# ui_runner.py
# Minimal local web UI runner for heatsonar + testbash orchestration.
#
# STRICT POLICY:
# - If required executables/scripts are missing, hard-fail.

import os
import time
import json
import shutil
import signal
import subprocess
from pathlib import Path
from typing import Dict, Any, Generator, Optional
from threading import Lock

from fastapi import FastAPI
from fastapi.responses import (
    HTMLResponse,
    StreamingResponse,
    JSONResponse,
    PlainTextResponse,
    Response,
)

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "ui_runs"
RUNS.mkdir(exist_ok=True)

REQUIRED = [
    ROOT / "heatsonar.exe",
    ROOT / "testbash.sh",
]

app = FastAPI()

# Track the currently running process per run_id (either heatsonar or testbash)
RUN_PROCS: Dict[str, subprocess.Popen] = {}
RUN_LOCK = Lock()


def hard_fail_missing():
    missing = [str(p) for p in REQUIRED if not p.exists()]
    if missing:
        raise RuntimeError("Missing required files: " + ", ".join(missing))


def sh(cmd, cwd: Path, env: Optional[Dict[str, str]] = None) -> subprocess.Popen:
    # Use bash so your existing .sh works in MSYS2.
    # Start a new process group on Windows so we can stop the whole tree.
    creationflags = 0
    if os.name == "nt":
        creationflags = subprocess.CREATE_NEW_PROCESS_GROUP

    return subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
        creationflags=creationflags,
        env=env,
    )


def stream_process(p: subprocess.Popen) -> Generator[str, None, int]:
    assert p.stdout is not None
    for line in p.stdout:
        yield line
    return p.wait()


def new_run_dir() -> Path:
    ts = time.strftime("%Y-%m-%d_%H-%M-%S")
    rd = RUNS / ts
    rd.mkdir(parents=True, exist_ok=False)
    return rd


def copy_artifacts(src_dir: Path, run_dir: Path):
    """
    Copy artifacts created in ROOT into this run's artifacts folder.
    heatsonar/testbash still run in ROOT; we copy outputs after completion.
    """
    out = run_dir / "artifacts"
    out.mkdir(exist_ok=True)

    # Always include UI config for this run (viz requires it)
    cfg = run_dir / "config.json"
    if cfg.exists():
        shutil.copy2(cfg, out / "config.json")

    # Core singletons
    candidates = [
        # heatsonar outputs (forced paths)
        src_dir / "heatmap.csv",
        src_dir / "seeds_run.txt",

        # sniper/testbash outputs
        src_dir / "sniper_runs" / "refined_Ts.txt",
        src_dir / "sniper_runs" / "refined_raw.txt",

        # your script currently writes finalts.txt (NOT tight_pairs.txt)
        src_dir / "sniper_runs" / "finalts.txt",

        # keep backward compatibility if you later rename back
        src_dir / "sniper_runs" / "tight_pairs.txt",
    ]

    for p in candidates:
        if p.exists():
            shutil.copy2(p, out / p.name)

    # Copy ALL per-seed step traces
    steps_dir = src_dir / "sniper_runs"
    if steps_dir.exists():
        for p in steps_dir.glob("run_*.steps.txt"):
            if p.is_file():
                shutil.copy2(p, out / p.name)


def set_status(run_dir: Path, status: str):
    (run_dir / "status.txt").write_text(status, encoding="utf-8")


def kill_process_tree_windows(pid: int):
    # Most reliable way on Windows to kill the whole subtree.
    subprocess.run(
        ["taskkill", "/F", "/T", "/PID", str(pid)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )


INDEX_HTML = r"""
<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>CToolkitV1 Scan UI</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 16px; }
    .row { display: flex; gap: 12px; flex-wrap: wrap; align-items: end; }
    label { display: block; font-size: 12px; margin-bottom: 4px; }
    input { width: 160px; padding: 6px; }
    button { padding: 8px 14px; }
    pre { background: #111; color: #ddd; padding: 12px; height: 420px; overflow: auto; }
    .box { border: 1px solid #ccc; padding: 12px; border-radius: 8px; }
    .muted { color: #666; font-size: 12px; }
    .pill { display: inline-block; padding: 4px 8px; border: 1px solid #ccc; border-radius: 999px; font-size: 12px; }
  </style>
</head>
<body>
  <h2>Heatsonar + Sniper Full Scan</h2>

  <div class="box">
    <div class="row">
      <div><label>int_start</label><input id="int_start" value="195"></div>
      <div><label>int_end</label><input id="int_end" value="204"></div>
      <div><label>T_step</label><input id="T_step" value="0.01"></div>
      <div><label>top</label><input id="top" value="500"></div>
      <div><label>dps (heatsonar)</label><input id="dps" value="80"></div>
    </div>

    <div class="row" style="margin-top: 12px;">
      <div><label>depth (sniper)</label><input id="depth" value="1000"></div>
      <div><label>max_step (sniper)</label><input id="max_step" value="500"></div>
      <div><label>pair_thresh</label><input id="pair_thresh" value="0.05"></div>

      <button onclick="startRun()">Run</button>
      <button onclick="stopRun()">Stop</button>

      <div id="run_id" style="padding: 0 0 2px 10px;"></div>
    </div>

    <div class="muted" style="margin-top: 8px;">
      <span class="pill">heatsonar writes heatmap.csv + seeds.txt</span>
      <span class="pill">testbash runs nsniper per seed and writes run_*.steps.txt</span>
      <span class="pill">Stop kills the active process tree</span>
    </div>
  </div>

  <h3>Live Log</h3>
  <pre id="log"></pre>

  <h3>Artifacts</h3>
  <div id="artifacts"></div>

<script>
let CURRENT_RUN_ID = "";
let ES = null;

async function startRun() {
  document.getElementById("log").textContent = "";
  document.getElementById("artifacts").innerHTML = "";
  document.getElementById("run_id").textContent = "";

  const payload = {
    int_start: document.getElementById("int_start").value,
    int_end: document.getElementById("int_end").value,
    T_step: document.getElementById("T_step").value,
    top: document.getElementById("top").value,
    dps: document.getElementById("dps").value,
    depth: document.getElementById("depth").value,
    max_step: document.getElementById("max_step").value,
    pair_thresh: document.getElementById("pair_thresh").value
  };

  const resp = await fetch("/start", {
    method: "POST",
    headers: {"Content-Type":"application/json"},
    body: JSON.stringify(payload)
  });

  if (!resp.ok) {
    const txt = await resp.text();
    document.getElementById("log").textContent = "ERROR: " + txt;
    return;
  }

  const data = await resp.json();
  CURRENT_RUN_ID = data.run_id;
  document.getElementById("run_id").textContent = "Run: " + data.run_id;

  if (ES) { ES.close(); ES = null; }
  ES = new EventSource("/stream?run_id=" + encodeURIComponent(data.run_id));

  ES.onmessage = (ev) => {
    const log = document.getElementById("log");
    log.textContent += ev.data + "\n";
    log.scrollTop = log.scrollHeight;
  };

  ES.addEventListener("done", async () => {
    if (ES) { ES.close(); ES = null; }

    const a = await fetch("/artifacts?run_id=" + encodeURIComponent(CURRENT_RUN_ID));
    const j = await a.json();

    const div = document.getElementById("artifacts");
    div.innerHTML = "";

    for (const item of j.items) {
      const p = document.createElement("div");
      const link = document.createElement("a");
      link.href = item.url;
      link.textContent = item.name;
      link.target = "_blank";
      p.appendChild(link);
      div.appendChild(p);
    }
  });
}

async function stopRun() {
  if (!CURRENT_RUN_ID) return;
  await fetch("/stop", {
    method: "POST",
    headers: {"Content-Type":"application/json"},
    body: JSON.stringify({run_id: CURRENT_RUN_ID})
  });
}
</script>
</body>
</html>
"""


# Your proposed viz page as a raw template; we substitute __RUN_ID__ at runtime.
VIZ_HTML_TEMPLATE = r"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Heatmap Viz</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 16px; }
    canvas { border: 1px solid #ccc; border-radius: 8px; }
    .row { display: flex; gap: 12px; flex-wrap: wrap; align-items: center; margin-bottom: 10px; }
    .pill { display: inline-block; padding: 4px 8px; border: 1px solid #ccc; border-radius: 999px; font-size: 12px; }
    .muted { color: #666; font-size: 12px; }
    a { text-decoration: none; }
    .small { font-size: 12px; color: #333; margin-top: 8px; }
  </style>
</head>
<body>
  <div class="row">
    <h2 style="margin:0;">Heatmap + Work Curves</h2>
    <span class="pill">Run: __RUN_ID__</span>
    <a class="pill" href="/">Back</a>
    <a class="pill" href="/artifacts?run_id=__RUN_ID__" target="_blank">Artifacts JSON</a>
  </div>

  <div class="row muted">
    <span class="pill">Top: heatmap log10(|Z|) curve with integer ticks + minima overlays</span>
    <span class="pill">Bottom: per-run sniper work curves from ALL run_*.steps.txt found in artifacts</span>
  </div>

  <canvas id="cv" width="1200" height="720"></canvas>
  <div id="debug" class="small"></div>

<script>
const RUN_ID = "__RUN_ID__";

function clamp01(x) {
  if (x < 0) return 0;
  if (x > 1) return 1;
  return x;
}

async function tryFetchText(url) {
  try {
    const r = await fetch(url);
    if (!r.ok) return null;
    return await r.text();
  } catch (e) {
    return null;
  }
}

function parseCSV(text) {
  const lines = text.split(/\r?\n/).filter(l => l.trim().length > 0);
  if (lines.length < 2) return null;
  const header = lines[0].split(",");
  const rows = [];
  for (let i=1; i<lines.length; i++) {
    const parts = lines[i].split(",");
    if (parts.length < header.length) continue;
    const obj = {};
    for (let k=0; k<header.length; k++) obj[header[k]] = parts[k];
    rows.push(obj);
  }
  return rows;
}

function parseTListTextKeepStrings(text) {
  const out = [];
  const lines = text.split(/\r?\n/);
  for (const l of lines) {
    const s = l.trim();
    if (!s) continue;
    if (s.startsWith("#")) continue;
    const tok = s.split(/\s+/)[0];
    out.push(tok);
  }
  return out;
}

function parseTListAsFloat(text) {
  const out = [];
  const lines = text.split(/\r?\n/);
  for (const l of lines) {
    const s = l.trim();
    if (!s) continue;
    if (s.startsWith("#")) continue;
    const tok = s.split(/\s+/)[0];
    const v = parseFloat(tok);
    if (!isNaN(v)) out.push(v);
  }
  return out;
}

function lowerBound(arr, x) {
  let lo = 0, hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (arr[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

// UPDATED: parse your new report format.
// We support lines like:
// "FRAC: 17" etc by maintaining current block state.
function parseStepsNewReport(text) {
  const lines = text.split(/\r?\n/);
  const pts = [];

  let frac = null;
  let elog = null;

  function flush() {
    if (frac !== null && elog !== null && isFinite(frac) && isFinite(elog)) {
      pts.push({x: frac, y: elog});
    }
    frac = null;
    elog = null;
  }

  for (const raw of lines) {
    const s = raw.trim();
    if (!s) continue;

    if (s.startsWith("----------------------------------------------------------------")) {
      flush();
      continue;
    }

    let m;

    m = s.match(/^FRAC:\s*([0-9]+)/i);
    if (m) { frac = parseInt(m[1], 10); continue; }

    // Prefer explicit E_LOG10 first
    m = s.match(/^E_LOG10:\s*\[\s*([-0-9.+eE]+)/i);
    if (m) { elog = parseFloat(m[1]); continue; }

    // Back-compat: "log10E=..."
    m = s.match(/log10E\s*=\s*([-0-9.+eE]+)/i);
    if (m) { elog = parseFloat(m[1]); continue; }

    // "E~1e-65"
    m = s.match(/E\s*~\s*1e([+-]?[0-9]+)/i);
    if (m) { elog = parseInt(m[1], 10); continue; }
  }

  flush();

  pts.sort((a,b) => a.x - b.x);
  const out = [];
  let lastX = null;
  for (const p of pts) {
    if (lastX === null || p.x !== lastX) { out.push(p); lastX = p.x; }
    else out[out.length - 1] = p;
  }
  return out;
}

function drawAxes(ctx, x0, y0, w, h) {
  ctx.strokeStyle = "#000000";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(x0, y0);
  ctx.lineTo(x0, y0 + h);
  ctx.lineTo(x0 + w, y0 + h);
  ctx.stroke();
}

function drawHGrid(ctx, x0, y0, w, h, n) {
  ctx.strokeStyle = "rgba(0,0,0,0.08)";
  ctx.lineWidth = 1;
  for (let i=0; i<=n; i++) {
    const yy = y0 + Math.floor((i / n) * h);
    ctx.beginPath();
    ctx.moveTo(x0, yy);
    ctx.lineTo(x0 + w, yy);
    ctx.stroke();
  }
}

async function main() {
  const base = "/dl/" + encodeURIComponent(RUN_ID) + "/";

  const cfgText = await tryFetchText(base + "config.json");
  if (!cfgText) {
    alert("Missing config.json in artifacts for run " + RUN_ID);
    return;
  }
  const cfg = JSON.parse(cfgText);

  const heatText = await tryFetchText(base + "heatmap.csv");
  if (!heatText) { alert("Missing heatmap.csv for run " + RUN_ID); return; }
  const heatRows = parseCSV(heatText);
  if (!heatRows || heatRows.length < 3) { alert("heatmap.csv parse failed."); return; }

  const Ts = [];
  const logEs = [];
  for (const r of heatRows) {
    const t = parseFloat(r["T"]);
    const le = parseFloat(r["log10E"]);
    if (!isNaN(t) && !isNaN(le)) { Ts.push(t); logEs.push(le); }
  }
  if (Ts.length < 3) { alert("heatmap.csv Ts too small."); return; }

  let minLE = +Infinity, maxLE = -Infinity;
  for (const v of logEs) { if (v < minLE) minLE = v; if (v > maxLE) maxLE = v; }

  const seedsText = await tryFetchText(base + "seeds.txt");
  const refinedText = await tryFetchText(base + "refined_Ts.txt");

  const seedStrings = seedsText ? parseTListTextKeepStrings(seedsText) : [];
  const seedFloats  = seedsText ? parseTListAsFloat(seedsText) : [];
  const refinedFloats = refinedText ? parseTListAsFloat(refinedText) : [];

  function interpLogE(t) {
    if (t <= Ts[0]) return logEs[0];
    if (t >= Ts[Ts.length - 1]) return logEs[logEs.length - 1];
    const i = lowerBound(Ts, t);
    const i0 = Math.max(0, i - 1);
    const i1 = Math.min(Ts.length - 1, i);
    const x0 = Ts[i0], x1 = Ts[i1];
    const y0 = logEs[i0], y1 = logEs[i1];
    if (x1 === x0) return y0;
    const a = (t - x0) / (x1 - x0);
    return y0 + a * (y1 - y0);
  }

  const cv = document.getElementById("cv");
  const ctx = cv.getContext("2d");
  const W = cv.width, H = cv.height;

  const padL = 80, padR = 24, padT = 26, padB = 60;
  const gap = 36;
  const topH = 300;
  const botH = H - padT - padB - topH - gap;
  const topY0 = padT;
  const botY0 = padT + topH + gap;
  const plotW = W - padL - padR;

  const tMin = Ts[0];
  const tMax = Ts[Ts.length - 1];

  function xOfT(t) {
    const u = (t - tMin) / (tMax - tMin);
    return padL + clamp01(u) * plotW;
  }
  function yTopOfLogE(le) {
    const u = (le - minLE) / (maxLE - minLE);
    return topY0 + (1 - clamp01(u)) * topH;
  }
  function intensityOfLogE(le) {
    const u = (le - minLE) / (maxLE - minLE);
    return 1 - clamp01(u);
  }

  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, W, H);

  // TOP
  drawHGrid(ctx, padL, topY0, plotW, topH, 6);

  const n = Ts.length;
  for (let i=0; i<n; i++) {
    const le = logEs[i];
    const inten = intensityOfLogE(le);
    const gg = Math.floor(inten * 255);
    const x0 = xOfT(Ts[i]);
    const x1 = (i + 1 < n) ? xOfT(Ts[i + 1]) : (padL + plotW);
    const w = Math.max(1, Math.floor(x1 - x0));
    ctx.fillStyle = "rgba(" + gg + "," + gg + "," + gg + ",0.55)";
    ctx.fillRect(Math.floor(x0), topY0, w, topH);
  }

  drawAxes(ctx, padL, topY0, plotW, topH);

  function drawTopCurve(style, width) {
    ctx.strokeStyle = style;
    ctx.lineWidth = width;
    ctx.beginPath();
    for (let i=0; i<n; i++) {
      const x = xOfT(Ts[i]);
      const y = yTopOfLogE(logEs[i]);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }
  drawTopCurve("rgba(255,255,255,0.90)", 5);
  drawTopCurve("rgba(0,80,255,0.95)", 2);

  const intStart = parseInt(cfg.int_start, 10);
  const intEnd = parseInt(cfg.int_end, 10);
  if (isFinite(intStart) && isFinite(intEnd) && intEnd >= intStart && (intEnd - intStart) <= 2000) {
    ctx.strokeStyle = "rgba(0,0,0,0.25)";
    ctx.fillStyle = "rgba(0,0,0,0.80)";
    ctx.lineWidth = 1;
    ctx.font = "12px Arial";

    for (let k=intStart; k<=intEnd; k++) {
      const x = xOfT(k);
      ctx.beginPath();
      ctx.moveTo(x, topY0 + topH);
      ctx.lineTo(x, topY0 + topH + 8);
      ctx.stroke();
      if ((intEnd - intStart) <= 40 || (k - intStart) % 2 === 0) {
        ctx.fillText(String(k), x - 10, topY0 + topH + 24);
      }
    }
  }

  function drawTopPoints(ts, strokeColor, fillColor, radius) {
    ctx.strokeStyle = strokeColor;
    ctx.fillStyle = fillColor;
    ctx.lineWidth = 2;
    for (const t of ts) {
      if (isNaN(t)) continue;
      const x = xOfT(t);
      const le = interpLogE(t);
      const y = yTopOfLogE(le);
      ctx.beginPath();
      ctx.arc(x, y, radius, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }
  }
  drawTopPoints(seedFloats, "#cc9900", "rgba(255,204,0,0.60)", 5);
  drawTopPoints(refinedFloats, "#00aa44", "rgba(0,255,102,0.50)", 6);

  ctx.fillStyle = "#000000";
  ctx.font = "16px Arial";
  ctx.fillText("Top: Heatmap log10(|Z|) vs T (integer ticks + minima overlays)", padL, topY0 - 6);

  // BOTTOM
  drawHGrid(ctx, padL, botY0, plotW, botH, 6);
  drawAxes(ctx, padL, botY0, plotW, botH);

  const allCols = [];
  for (let i = 0; i < seedStrings.length; i++) {
    const s = seedStrings[i];
    const tCol = parseFloat(s);
    if (isNaN(tCol)) continue;

    const fn = "run_" + s + ".steps.txt";
    const txt = await tryFetchText(base + fn);
    if (!txt) continue;

    const pts = parseStepsNewReport(txt);
    if (pts.length < 1) continue;

    allCols.push({ seed: s, t: tCol, pts: pts });
  }

  if (allCols.length === 0) {
    ctx.fillStyle = "#000000";
    ctx.font = "14px Arial";
    ctx.fillText("Bottom: no per-seed steps found. (Did you copy run_*.steps.txt into artifacts?)", padL, botY0 + 20);
    return;
  }

  let botMinE = +Infinity, botMaxE = -Infinity;
  for (const c of allCols) {
    for (const p of c.pts) {
      if (p.y < botMinE) botMinE = p.y;
      if (p.y > botMaxE) botMaxE = p.y;
    }
  }
  if (!isFinite(botMinE) || !isFinite(botMaxE) || botMaxE === botMinE) {
    botMinE = -1000; botMaxE = 0;
  } else {
    const yPad = 0.08 * (botMaxE - botMinE);
    botMinE -= yPad;
    botMaxE += yPad;
  }

  function xBotT(t) {
    const u = (t - tMin) / (tMax - tMin);
    return padL + clamp01(u) * plotW;
  }
  function yBotE(le) {
    const u = (le - botMinE) / (botMaxE - botMinE);
    return botY0 + (1 - clamp01(u)) * botH;
  }

  ctx.fillStyle = "#000000";
  ctx.font = "16px Arial";
  ctx.fillText("Bottom: Sniper lattice samples (x=T, y=E_LOG10 per step)", padL, botY0 - 6);

  ctx.font = "12px Arial";
  ctx.fillText("log10(|E|)", 10, botY0 + 12);

  if (isFinite(intStart) && isFinite(intEnd) && intEnd >= intStart && (intEnd - intStart) <= 2000) {
    ctx.strokeStyle = "rgba(0,0,0,0.25)";
    ctx.fillStyle = "rgba(0,0,0,0.80)";
    ctx.lineWidth = 1;
    ctx.font = "12px Arial";

    for (let k=intStart; k<=intEnd; k++) {
      const x = xBotT(k);
      ctx.beginPath();
      ctx.moveTo(x, botY0 + botH);
      ctx.lineTo(x, botY0 + botH + 8);
      ctx.stroke();
      if ((intEnd - intStart) <= 40 || (k - intStart) % 2 === 0) {
        ctx.fillText(String(k), x - 10, botY0 + botH + 24);
      }
    }
  }

  const palette = [
    "rgba(0, 170, 68, 0.90)",
    "rgba(0, 80, 255, 0.90)",
    "rgba(204, 153, 0, 0.90)",
    "rgba(180, 0, 180, 0.90)",
    "rgba(0, 160, 160, 0.90)",
    "rgba(220, 60, 0, 0.90)"
  ];

  for (let i = 0; i < allCols.length; i++) {
    const c = allCols[i];
    const color = palette[i % palette.length];
    const x = xBotT(c.t);

    let cMin = +Infinity, cMax = -Infinity;
    for (const p of c.pts) { if (p.y < cMin) cMin = p.y; if (p.y > cMax) cMax = p.y; }
    if (isFinite(cMin) && isFinite(cMax)) {
      ctx.strokeStyle = "rgba(0,0,0,0.10)";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(x, yBotE(cMin));
      ctx.lineTo(x, yBotE(cMax));
      ctx.stroke();
    }

    ctx.fillStyle = color;
    ctx.strokeStyle = "rgba(0,0,0,0.30)";
    ctx.lineWidth = 1;

    for (let j = 0; j < c.pts.length; j++) {
      const y = yBotE(c.pts[j].y);
      ctx.beginPath();
      ctx.arc(x, y, 3, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }

    const last = c.pts[c.pts.length - 1];
    if (last) {
      ctx.fillStyle = color;
      ctx.strokeStyle = "rgba(0,0,0,0.35)";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(x, yBotE(last.y), 6, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }
  }

  ctx.fillStyle = "#000000";
  ctx.font = "12px Arial";
  ctx.fillText(
    "seeds=" + seedStrings.length +
    "  step_files_loaded=" + allCols.length +
    "  depth_cfg=" + String(cfg.depth) +
    "  heat_samples=" + String(n),
    padL, H - 20
  );
}

main();
</script>
</body>
</html>
"""


def viz_html(run_id: str) -> str:
    return VIZ_HTML_TEMPLATE.replace("__RUN_ID__", run_id)


@app.get("/", response_class=HTMLResponse)
def index():
    hard_fail_missing()
    return INDEX_HTML


@app.get("/viz", response_class=HTMLResponse)
def viz(run_id: str):
    hard_fail_missing()
    run_dir = RUNS / run_id
    if not run_dir.exists():
        return HTMLResponse(f"<pre>Run not found: {run_id}</pre>", status_code=404)
    return viz_html(run_id)


@app.post("/start")
async def start(payload: Dict[str, Any]):
    hard_fail_missing()

    run_dir = new_run_dir()
    (run_dir / "log.txt").write_text("", encoding="utf-8")

    cfg = {
        "int_start": str(payload.get("int_start", "")),
        "int_end": str(payload.get("int_end", "")),
        "T_step": str(payload.get("T_step", "")),
        "top": str(payload.get("top", "")),
        "dps": str(payload.get("dps", "")),
        "depth": str(payload.get("depth", "")),
        "max_step": str(payload.get("max_step", "")),
        "pair_thresh": str(payload.get("pair_thresh", "")),
    }
    (run_dir / "config.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    set_status(run_dir, "pending")

    return JSONResponse({"run_id": run_dir.name})


@app.post("/stop")
async def stop(payload: Dict[str, Any]):
    hard_fail_missing()
    run_id = str(payload.get("run_id", ""))
    if not run_id:
        return JSONResponse({"error": "missing run_id"}, status_code=400)

    run_dir = RUNS / run_id
    if run_dir.exists():
        set_status(run_dir, "stopping")

    with RUN_LOCK:
        p = RUN_PROCS.get(run_id)

    if p is None:
        if run_dir.exists():
            set_status(run_dir, "not_running")
        return JSONResponse({"status": "not_running"})

    try:
        if os.name == "nt":
            try:
                p.send_signal(signal.CTRL_BREAK_EVENT)
                time.sleep(0.2)
            except Exception:
                pass
            kill_process_tree_windows(p.pid)
        else:
            p.terminate()
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=500)

    with RUN_LOCK:
        RUN_PROCS.pop(run_id, None)

    if run_dir.exists():
        set_status(run_dir, "stopped")

    return JSONResponse({"status": "stopped"})


@app.get("/stream")
async def stream(run_id: str):
    hard_fail_missing()
    run_dir = RUNS / run_id
    if not run_dir.exists():
        return StreamingResponse(iter(["Run not found\n"]), media_type="text/plain")

    cfg = json.loads((run_dir / "config.json").read_text(encoding="utf-8"))
    log_path = run_dir / "log.txt"

    def emit(line: str):
        with log_path.open("a", encoding="utf-8") as f:
            f.write(line + "\n")

    def event_gen():
        set_status(run_dir, "running")

        env = os.environ.copy()
        env["DEPTH"] = cfg["depth"]
        env["MAX_STEP"] = cfg["max_step"]
        env["PAIR_THRESH"] = cfg["pair_thresh"]

        heatmap_file = ROOT / "heatmap.csv"
        seeds_file = ROOT / "seeds.txt"

        heat_cmd = [
            str(ROOT / "heatsonar.exe"),
            "--mode", "range",
            "--int_start", cfg["int_start"],
            "--int_end", cfg["int_end"],
            "--T_step", cfg["T_step"],
            "--top", cfg["top"],
            "--dps", cfg["dps"],
            "--heatmap", str(heatmap_file),
            "--seeds", str(seeds_file),
        ]

        bash_cmd = [
            "bash",
            str(ROOT / "testbash.sh"),
            cfg["depth"],
            cfg["max_step"],
            cfg["pair_thresh"],
        ]

        try:
            emit("[UI] heatsonar: " + " ".join(heat_cmd))
            yield f"data: [UI] heatsonar: {' '.join(heat_cmd)}\n\n"

            p1 = sh(heat_cmd, ROOT, env=env)
            with RUN_LOCK:
                RUN_PROCS[run_id] = p1

            for line in stream_process(p1):
                line = line.rstrip("\n")
                emit(line)
                yield f"data: {line}\n\n"

            rc1 = p1.poll()
            if rc1 is None:
                rc1 = p1.wait()

            if rc1 != 0:
                msg = f"[UI] heatsonar failed rc={rc1}"
                emit(msg)
                yield f"data: {msg}\n\n"
                set_status(run_dir, "failed_heatsonar")
                with RUN_LOCK:
                    RUN_PROCS.pop(run_id, None)
                yield "event: done\ndata: done\n\n"
                return

            emit("[UI] testbash: " + " ".join(bash_cmd))
            yield f"data: [UI] testbash: {' '.join(bash_cmd)}\n\n"

            p2 = sh(bash_cmd, ROOT, env=env)
            with RUN_LOCK:
                RUN_PROCS[run_id] = p2

            for line in stream_process(p2):
                line = line.rstrip("\n")
                emit(line)
                yield f"data: {line}\n\n"

            rc2 = p2.poll()
            if rc2 is None:
                rc2 = p2.wait()

            if rc2 != 0:
                msg = f"[UI] testbash failed rc={rc2}"
                emit(msg)
                yield f"data: {msg}\n\n"
                set_status(run_dir, "failed_testbash")
            else:
                emit("[UI] done")
                yield "data: [UI] done\n\n"
                copy_artifacts(ROOT, run_dir)
                set_status(run_dir, "done")

        finally:
            with RUN_LOCK:
                RUN_PROCS.pop(run_id, None)

        yield "event: done\ndata: done\n\n"

    return StreamingResponse(event_gen(), media_type="text/event-stream")


@app.get("/artifacts")
async def artifacts(run_id: str):
    run_dir = RUNS / run_id
    art = run_dir / "artifacts"
    items = []

    if run_dir.exists():
        items.append({"name": "VIEW_HEATMAP.html", "url": f"/viz?run_id={run_id}"})

    if art.exists():
        for p in sorted(art.iterdir()):
            items.append({"name": p.name, "url": f"/dl/{run_id}/{p.name}"})
    return JSONResponse({"items": items})


@app.get("/dl/{run_id}/{name}")
async def dl(run_id: str, name: str):
    run_dir = RUNS / run_id
    p = run_dir / "artifacts" / name
    if not p.exists():
        return JSONResponse({"error": "not found"}, status_code=404)

    # Serve html as HTML so it renders
    if name.lower().endswith(".html"):
        return HTMLResponse(p.read_text(encoding="utf-8"))

    # Everything else as plain text
    return PlainTextResponse(p.read_text(encoding="utf-8"))


if __name__ == "__main__":
    hard_fail_missing()
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
