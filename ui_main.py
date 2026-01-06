# ui_main.py
# Entry point for running the UI server.

import work_runner
from ui_app import app

if __name__ == "__main__":
    work_runner.hard_fail_missing()
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
