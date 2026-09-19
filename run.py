#!/usr/bin/env python3
"""Launch the explorer: python run.py [--data-root DIR] [--port 8766] [--host 127.0.0.1]"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=None, help="workspace root containing pyepi_results/, default_300cells/, pyEpiAneufinder/")
    ap.add_argument("--port", type=int, default=int(os.environ.get("PYEPI_SCIENTIFIC_PORT", 8766)))
    ap.add_argument("--host", default=os.environ.get("PYEPI_SCIENTIFIC_HOST", "127.0.0.1"))
    ap.add_argument("--reload", action="store_true")
    args = ap.parse_args()
    if args.data_root:
        os.environ["PYEPI_SCIENTIFIC_ROOT"] = os.path.abspath(args.data_root)
    import uvicorn
    uvicorn.run("server.app:app", host=args.host, port=args.port, reload=args.reload, log_level="info")


if __name__ == "__main__":
    main()
