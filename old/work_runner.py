# work_runner.py
# Orchestration + filesystem layout + strict missing-file policy.
# FastAPI layer calls into this module.

import os
import time
import json
import shutil
import signal
import subprocess
from pathlib import Path
from typing import Dict, Any, Generator, Optional
from threading import Lock

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "ui_runs"
RUNS.mkdir(exist_ok=True)

REQUIRED = [
    ROOT / "heatsonar.exe",
    ROOT / "nsniperv6.exe",
    ROOT / "run_pipeline.sh",
]

RUN_PROCS: Dict[str, subprocess.Popen] = {}
RUN_LOCK = Lock()


def hard_fail_missing():
    missing = [str(p) for p in REQUIRED if not p.exists()]
    if missing:
        raise RuntimeError("Missing required files: " + ", ".join(missing))


def new_run_dir() -> Path:
    ts = time.strftime("%Y-%m-%d_%H-%M-%S")
    rd = RUNS / ts
    rd.mkdir(parents=True, exist_ok=False)
    return rd


def set_status(run_dir: Path, status: str):
    (run_dir / "status.txt").write_text(status, encoding="utf-8")


def append_log(run_dir: Path, line: str):
    log_path = run_dir / "log.txt"
    with log_path.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def kill_process_tree_windows(pid: int):
    subprocess.run(
        ["taskkill", "/F", "/T", "/PID", str(pid)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )


def sh(cmd, cwd: Path, env: Optional[Dict[str, str]] = None) -> subprocess.Popen:
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
        yield line.rstrip("\n")
    return p.wait()


def start_run(payload: Dict[str, Any]) -> str:
    """
    Creates run dir + config.json + empty log.
    The actual pipeline is executed by run_pipeline_stream(run_id).
    """
    hard_fail_missing()
    run_dir = new_run_dir()

    # Store everything the UI may need later (viz + artifact selection).
    cfg = {
        # heatsonar range inputs
        "int_start": str(payload.get("int_start", "")),
        "int_end": str(payload.get("int_end", "")),
        "T_step": str(payload.get("T_step", "")),
        "dps": str(payload.get("dps", "")),
        "pretty_digits": str(payload.get("pretty_digits", "")),
        "max_step_heat": str(payload.get("max_step_heat", "")),
        "target_energy": str(payload.get("target_energy", "")),
        "order": str(payload.get("order", "")),
        "csv_hits": str(payload.get("csv_hits", "")),

        # outputs and pipeline layout
        "heatmap_path": str(payload.get("heatmap_path", "heatmap.csv")),
        "seeds_path": str(payload.get("seeds_path", "seeds.txt")),
        "outdir": str(payload.get("outdir", "sniper_runs")),

        # nsniperv6
        "sniper_depth": str(payload.get("sniper_depth", "")),
        "sniper_max_step": str(payload.get("sniper_max_step", "")),
        "sniper_extra": str(payload.get("sniper_extra", "")),

        # pipeline flags
        "skip_heatsonar": str(payload.get("skip_heatsonar", "")),
        "quiet": str(payload.get("quiet", "")),
    }

    (run_dir / "config.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    (run_dir / "log.txt").write_text("", encoding="utf-8")
    set_status(run_dir, "pending")
    return run_dir.name


def stop_run(run_id: str) -> str:
    """
    Stops active process tree for this run_id.
    Returns one of: not_running | stopped | stop_error
    """
    hard_fail_missing()
    run_dir = RUNS / run_id
    if run_dir.exists():
        set_status(run_dir, "stopping")

    with RUN_LOCK:
        p = RUN_PROCS.get(run_id)

    if p is None:
        if run_dir.exists():
            set_status(run_dir, "not_running")
        return "not_running"

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
    except Exception:
        if run_dir.exists():
            set_status(run_dir, "stop_error")
        return "stop_error"

    with RUN_LOCK:
        RUN_PROCS.pop(run_id, None)

    if run_dir.exists():
        set_status(run_dir, "stopped")
    return "stopped"


def copy_artifacts(src_dir: Path, run_dir: Path, cfg: Optional[Dict[str, Any]] = None):
    """
    Copy outputs created in ROOT into run_dir/artifacts.

    New model:
      - Per-seed logs are outdir/run_*.stdout.txt (stdout is the artifact)
      - Combined log is outdir/all_runs.stdout.txt
      - heatmap + seeds are configured by cfg["heatmap_path"], cfg["seeds_path"]
    """
    out = run_dir / "artifacts"
    out.mkdir(exist_ok=True)

    # Always copy config.json (viz requires it)
    cfg_path = run_dir / "config.json"
    if cfg_path.exists():
        shutil.copy2(cfg_path, out / "config.json")

    heatmap_name = "heatmap.csv"
    seeds_name = "seeds.txt"
    outdir_name = "sniper_runs"

    if cfg:
        hm = str(cfg.get("heatmap_path", "")).strip()
        sd = str(cfg.get("seeds_path", "")).strip()
        od = str(cfg.get("outdir", "")).strip()
        if hm:
            heatmap_name = Path(hm).name
        if sd:
            seeds_name = Path(sd).name
        if od:
            outdir_name = Path(od).name

    # Copy heatmap + seeds if present at ROOT
    for p in [src_dir / heatmap_name, src_dir / seeds_name]:
        if p.exists():
            shutil.copy2(p, out / p.name)

    # Copy sniper stdout artifacts
    sniper_dir = src_dir / outdir_name
    if sniper_dir.exists():
        all_log = sniper_dir / "all_runs.stdout.txt"
        if all_log.exists():
            shutil.copy2(all_log, out / all_log.name)

        for p in sorted(sniper_dir.glob("run_*.stdout.txt")):
            shutil.copy2(p, out / p.name)


def run_pipeline_stream(run_id: str) -> Generator[str, None, None]:
    """
    Executes:
      - bash run_pipeline.sh ... which runs:
          1) heatsonar -> heatmap + seeds
          2) nsniperv6 per seed -> outdir/run_*.stdout.txt + outdir/all_runs.stdout.txt

    Yields plain log lines (FastAPI wraps into SSE).
    Always terminates (even on failure), so UI sees done event.
    """
    hard_fail_missing()
    run_dir = RUNS / run_id
    if not run_dir.exists():
        yield f"[UI] ERROR: run not found: {run_id}"
        return

    cfg = json.loads((run_dir / "config.json").read_text(encoding="utf-8"))
    set_status(run_dir, "running")

    env = os.environ.copy()

    def cfgv(key: str, default: str) -> str:
        v = str(cfg.get(key, "")).strip()
        return v if v else default

    # Build command for your pipeline script
    bash_cmd = [
        "bash",
        str(ROOT / "run_pipeline.sh"),
        "--outdir", cfgv("outdir", "sniper_runs"),
        "--heatmap", cfgv("heatmap_path", "heatmap.csv"),
        "--seeds", cfgv("seeds_path", "seeds.txt"),
        "--int_start", cfgv("int_start", "10"),
        "--int_end", cfgv("int_end", "100"),
        "--T_step", cfgv("T_step", "0.001"),
        "--dps", cfgv("dps", "80"),
        "--pretty_digits", cfgv("pretty_digits", "25"),
        "--max_step_heat", cfgv("max_step_heat", "200"),
        "--order", cfgv("order", "T"),
        "--sniper_depth", cfgv("sniper_depth", "25"),
        "--sniper_max_step", cfgv("sniper_max_step", "100"),
    ]

    # Optional flags (only add if set)
    target_energy = str(cfg.get("target_energy", "")).strip()
    if target_energy:
        bash_cmd += ["--target_energy", target_energy]

    csv_hits = str(cfg.get("csv_hits", "")).strip()
    if csv_hits:
        bash_cmd += ["--csv", csv_hits]

    sniper_extra = str(cfg.get("sniper_extra", "")).strip()
    if sniper_extra:
        bash_cmd += ["--sniper_extra", sniper_extra]

    skip_heatsonar = str(cfg.get("skip_heatsonar", "")).strip()
    if skip_heatsonar in ("1", "true", "True", "yes", "YES", "on", "ON"):
        bash_cmd += ["--skip_heatsonar"]

    quiet = str(cfg.get("quiet", "")).strip()
    if quiet in ("1", "true", "True", "yes", "YES", "on", "ON"):
        bash_cmd += ["--quiet"]

    def log(line: str):
        append_log(run_dir, line)

    try:
        line = "[UI] pipeline: " + " ".join(bash_cmd)
        log(line)
        yield line

        p = sh(bash_cmd, ROOT, env=env)
        with RUN_LOCK:
            RUN_PROCS[run_id] = p

        for ln in stream_process(p):
            log(ln)
            yield ln

        rc = p.poll()
        if rc is None:
            rc = p.wait()

        with RUN_LOCK:
            RUN_PROCS.pop(run_id, None)

        if rc != 0:
            msg = f"[UI] pipeline failed rc={rc}"
            log(msg)
            yield msg
            set_status(run_dir, "failed_pipeline")
            return

        msg = "[UI] done"
        log(msg)
        yield msg

        copy_artifacts(ROOT, run_dir, cfg)
        set_status(run_dir, "done")

    finally:
        with RUN_LOCK:
            RUN_PROCS.pop(run_id, None)
