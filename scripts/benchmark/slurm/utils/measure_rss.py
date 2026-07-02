#!/usr/bin/env python3
"""
Run a command and report its peak RSS (VmPeak from /proc/PID/status).
Writes a time_verbose.txt compatible line:
    Maximum resident set size (kbytes): N

Usage: measure_rss.py <output_file> -- <cmd> [args...]
"""
import os
import subprocess
import sys
import time


def peak_rss_kb(pid):
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except FileNotFoundError:
        pass
    return 0


def main():
    if len(sys.argv) < 4 or sys.argv[2] != "--":
        print("Usage: measure_rss.py <output_file> -- <cmd> [args...]", file=sys.stderr)
        sys.exit(1)

    out_file = sys.argv[1]
    cmd = sys.argv[3:]

    proc = subprocess.Popen(cmd)
    max_rss = 0
    while proc.poll() is None:
        rss = peak_rss_kb(proc.pid)
        if rss > max_rss:
            max_rss = rss
        time.sleep(0.2)

    # Final read after process ends
    returncode = proc.returncode

    with open(out_file, "w") as f:
        f.write(f"\tMaximum resident set size (kbytes): {max_rss}\n")
        f.write(f"\tExit status: {returncode}\n")

    sys.exit(returncode)


if __name__ == "__main__":
    main()
