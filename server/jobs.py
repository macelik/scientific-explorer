"""Bounded background execution: an interactive pool (permutation tests) and a
lower-priority background pool (cohort indexing). Jobs are deduplicated by key
and results cached (in memory + JSON on disk for seeded permutation tests)."""
from __future__ import annotations

import json
import multiprocessing
import os
import threading
import time
import uuid
from concurrent.futures import Future, ProcessPoolExecutor
from typing import Any, Callable, Dict, Optional


class JobManager:
    def __init__(self, cache_dir: str, interactive_workers: int = 4, background_workers: int = 4):
        self.cache_dir = cache_dir
        # Fresh processes must not inherit the server socket or its signal handlers.
        context = multiprocessing.get_context("spawn")
        self.interactive = ProcessPoolExecutor(max_workers=interactive_workers, mp_context=context)
        self.background = ProcessPoolExecutor(max_workers=background_workers, mp_context=context)
        self.jobs: Dict[str, dict] = {}
        self.by_key: Dict[str, str] = {}
        self.lock = threading.Lock()
        self.result_cache: Dict[str, Any] = {}
        self._disk_cache_path = os.path.join(cache_dir, "permutation_tests.json")
        self._load_disk_cache()

    def _load_disk_cache(self) -> None:
        if os.path.exists(self._disk_cache_path):
            try:
                with open(self._disk_cache_path) as f:
                    self.result_cache.update(json.load(f))
            except Exception:
                pass

    def _persist(self, key: str, result: Any) -> None:
        if not key.startswith("test:"):
            return
        try:
            with self.lock:
                data = {}
                if os.path.exists(self._disk_cache_path):
                    with open(self._disk_cache_path) as f:
                        data = json.load(f)
                data[key] = result
                tmp = self._disk_cache_path + ".tmp"
                with open(tmp, "w") as f:
                    json.dump(data, f)
                os.replace(tmp, self._disk_cache_path)
        except Exception:
            pass

    def cached(self, key: str) -> Optional[Any]:
        return self.result_cache.get(key)

    def submit(self, key: str, fn: Callable, *args, kind: str = "interactive", meta: Optional[dict] = None) -> dict:
        with self.lock:
            if key in self.result_cache:
                return {"id": None, "status": "done", "result": self.result_cache[key], "key": key, "cached": True}
            if key in self.by_key:
                jid = self.by_key[key]
                return self._view(jid)
            jid = uuid.uuid4().hex[:12]
            pool = self.interactive if kind == "interactive" else self.background
            job = {"id": jid, "key": key, "status": "running", "result": None, "error": None,
                   "submitted": time.time(), "finished": None, "meta": meta or {}, "progress": None}
            self.jobs[jid] = job
            self.by_key[key] = jid
        fut: Future = pool.submit(fn, *args)

        def _done(f: Future, jid=jid, key=key):
            with self.lock:
                job = self.jobs[jid]
                job["finished"] = time.time()
                try:
                    res = f.result()
                    job["result"] = res
                    job["status"] = "done"
                    self.result_cache[key] = res
                except Exception as e:  # noqa: BLE001
                    job["error"] = f"{type(e).__name__}: {e}"
                    job["status"] = "error"
                self.by_key.pop(key, None)
            if job["status"] == "done":
                self._persist(key, job["result"])

        fut.add_done_callback(_done)
        return self._view(jid)

    def _view(self, jid: str) -> dict:
        j = self.jobs[jid]
        return {"id": jid, "key": j["key"], "status": j["status"], "result": j["result"], "error": j["error"],
                "meta": j["meta"], "elapsed": (j["finished"] or time.time()) - j["submitted"], "cached": False}

    def get(self, jid: str) -> Optional[dict]:
        with self.lock:
            if jid not in self.jobs:
                return None
            return self._view(jid)

    def shutdown(self) -> None:
        self.interactive.shutdown(wait=False, cancel_futures=True)
        self.background.shutdown(wait=False, cancel_futures=True)
