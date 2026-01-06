# ui_app.py
# FastAPI layer only. No process/orchestration logic lives here.

from pathlib import Path
from typing import Dict, Any

from fastapi import FastAPI
from fastapi.responses import HTMLResponse, StreamingResponse, JSONResponse, PlainTextResponse

import work_runner

ROOT = Path(__file__).resolve().parent
TEMPLATES = ROOT / "templates"

app = FastAPI()


def load_template(name: str) -> str:
    p = TEMPLATES / name
    if not p.exists():
        raise RuntimeError(f"Missing template: {p}")
    return p.read_text(encoding="utf-8")


@app.get("/", response_class=HTMLResponse)
def index():
    work_runner.hard_fail_missing()
    return load_template("index.html")


@app.get("/viz", response_class=HTMLResponse)
def viz(run_id: str):
    work_runner.hard_fail_missing()
    run_dir = work_runner.RUNS / run_id
    if not run_dir.exists():
        return HTMLResponse(f"<pre>Run not found: {run_id}</pre>", status_code=404)

    html = load_template("viz.html").replace("__RUN_ID__", run_id)
    return HTMLResponse(html)


@app.post("/start")
async def start(payload: Dict[str, Any]):
    work_runner.hard_fail_missing()
    run_id = work_runner.start_run(payload)
    return JSONResponse({"run_id": run_id})


@app.post("/stop")
async def stop(payload: Dict[str, Any]):
    work_runner.hard_fail_missing()
    run_id = str(payload.get("run_id", "")).strip()
    if not run_id:
        return JSONResponse({"error": "missing run_id"}, status_code=400)

    status = work_runner.stop_run(run_id)
    if status == "stop_error":
        return JSONResponse({"error": "stop_error"}, status_code=500)

    return JSONResponse({"status": status})


@app.get("/stream")
async def stream(run_id: str):
    work_runner.hard_fail_missing()

    run_dir = work_runner.RUNS / run_id
    if not run_dir.exists():
        return StreamingResponse(iter(["Run not found\n"]), media_type="text/plain")

    def event_gen():
        # work_runner is responsible for setting status/logging.
        for line in work_runner.run_pipeline_stream(run_id):
            # SSE requires each message to be prefixed with "data: "
            yield f"data: {line}\n\n"
        yield "event: done\ndata: done\n\n"

    return StreamingResponse(event_gen(), media_type="text/event-stream")


@app.get("/artifacts")
async def artifacts(run_id: str):
    run_dir = work_runner.RUNS / run_id
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
    run_dir = work_runner.RUNS / run_id
    p = run_dir / "artifacts" / name
    if not p.exists():
        return JSONResponse({"error": "not found"}, status_code=404)

    return PlainTextResponse(p.read_text(encoding="utf-8"))


if __name__ == "__main__":
    work_runner.hard_fail_missing()
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
