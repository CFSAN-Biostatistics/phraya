#!/usr/bin/env python3
"""Run a command, reporting wall-clock seconds and peak child RSS (MB).

Portable replacement for `/usr/bin/time -v`, which is not installed everywhere
(e.g. minimal WSL/containers). Peak RSS comes from getrusage(RUSAGE_CHILDREN).
ru_maxrss is in kilobytes on Linux.

Usage: time_run.py <cmd> [args...]
Prints one line to stdout: WALL_SECONDS<TAB>PEAK_RSS_MB<TAB>EXIT_CODE
The child's own stdout/stderr are inherited (pass through).
"""
import resource
import subprocess
import sys
import time


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: time_run.py <cmd> [args...]", file=sys.stderr)
        return 2
    cmd = sys.argv[1:]
    start = time.perf_counter()
    proc = subprocess.run(cmd)
    wall = time.perf_counter() - start
    # RUSAGE_CHILDREN accumulates over all reaped children; for a single run this
    # is the peak RSS of the one child. ru_maxrss is KB on Linux.
    peak_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak_mb = peak_kb / 1024.0
    print(f"{wall:.3f}\t{peak_mb:.1f}\t{proc.returncode}")
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main())
