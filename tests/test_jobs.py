"""Regression: computation workers must not keep the server socket alive."""
import os
from pathlib import Path
import socket
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server.jobs import JobManager


def owns_socket(socket_id):
    for fd in Path('/proc/self/fd').iterdir():
        try:
            if os.readlink(fd) == socket_id:
                return True
        except FileNotFoundError:
            pass
    return False


def test_computation_worker_does_not_inherit_server_socket(tmp_path):
    with socket.socket() as listener:
        socket_id = os.readlink(f'/proc/self/fd/{listener.fileno()}')
        jobs = JobManager(str(tmp_path), interactive_workers=1, background_workers=1)
        try:
            job = jobs.submit('socket-inheritance', owns_socket, socket_id)
            until = time.monotonic()+15
            while job['status'] == 'running' and time.monotonic()<until:
                time.sleep(.05)
                job = jobs.get(job['id'])
            assert job['status'] == 'done',job
            assert job['result'] is False
        finally:
            jobs.shutdown()
