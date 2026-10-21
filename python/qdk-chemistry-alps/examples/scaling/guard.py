# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/scaling/guard.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 bab8396d12e81b48abf6b276a9457775d712e0a1818555a61a1f8627270e9372.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Run one benchmark with a sampled RAM watchdog (needs psutil).

Example: python guard.py -- python run.py --case hubbard ...
Run resource-intensive benchmarks sequentially. RSS and macOS physical
footprint include the child process tree (conservatively double-counting
shared pages). This sampled watchdog is not an OS-enforced allocation limit.
"""

import argparse
import math
import os
import signal
import subprocess
import sys
import ctypes
import time

import psutil

class _MacUsage(ctypes.Structure):
    # rusage_info_v0, from <sys/resource.h>. Footprint includes compressed
    # memory, which macOS RSS omits; guard both on this platform.
    _fields_ = [("uuid", ctypes.c_uint8 * 16)] + [
        (name, ctypes.c_uint64) for name in (
            "user", "system", "idle", "interrupt", "pageins", "wired",
            "resident", "footprint", "start", "exit")]


def footprint_reader():
    if sys.platform != "darwin":
        return lambda pid: 0
    lib = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
    lib.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(_MacUsage)]
    lib.proc_pid_rusage.restype = ctypes.c_int

    def read(pid):
        usage = _MacUsage()
        if lib.proc_pid_rusage(pid, 0, ctypes.byref(usage)) != 0:
            error = ctypes.get_errno()
            if error == 3:  # ESRCH: child exited between samples.
                return 0
            raise OSError(error, os.strerror(error))
        return usage.footprint
    return read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-gb", type=float, default=12)
    parser.add_argument("--min-available-gb", type=float, default=8)
    parser.add_argument("--max-seconds", type=float, default=None)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if (not command or not math.isfinite(args.max_gb) or args.max_gb <= 0
            or not math.isfinite(args.min_available_gb) or args.min_available_gb < 0):
        parser.error("provide a command after -- and valid nonnegative memory limits")
    if args.max_seconds is not None and (not math.isfinite(args.max_seconds) or args.max_seconds <= 0):
        parser.error("max-seconds must be positive")
    footprint = footprint_reader()
    started = time.monotonic()
    child = subprocess.Popen(command, start_new_session=True)
    process = psutil.Process(child.pid)
    print(f"Watching PID {child.pid}: RSS/footprint cap {args.max_gb:g} GB, "
          f"system available floor {args.min_available_gb:g} GB", flush=True)
    peak = peak_footprint = 0
    try:
        while child.poll() is None:
            try:
                descendants = process.children(recursive=True)
            except psutil.NoSuchProcess:
                break
            rss = physical = 0
            for member in [process, *descendants]:
                try:
                    rss += member.memory_info().rss
                    physical += footprint(member.pid)
                except psutil.NoSuchProcess:
                    pass  # A short-lived compiler/helper exited between samples.
            peak = max(peak, rss)
            peak_footprint = max(peak_footprint, physical)
            available = psutil.virtual_memory().available
            if max(rss, physical) > args.max_gb * 1e9 or available < args.min_available_gb * 1e9:
                print(f"MEMORY GUARD STOP: RSS {rss / 1e9:.2f} GB, "
                      f"footprint {physical / 1e9:.2f} GB, available {available / 1e9:.2f} GB", flush=True)
                raise SystemExit(75)
            if args.max_seconds is not None and time.monotonic() - started > args.max_seconds:
                print(f"TIME GUARD STOP: {args.max_seconds:g} s limit; "
                      f"watched peak RSS {peak / 1e9:.2f} GB", flush=True)
                raise SystemExit(124)
            time.sleep(0.25)
        print(f"Watched peak RSS: {peak / 1e9:.2f} GB; "
              f"footprint: {peak_footprint / 1e9:.2f} GB", flush=True)
        raise SystemExit(child.wait())
    finally:
        # A failed or interrupted watchdog must not leave its benchmark running.
        if child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()


if __name__ == "__main__":
    main()
