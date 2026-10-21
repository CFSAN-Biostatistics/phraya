#!/usr/bin/env python3
"""
Run a command and report the peak RSS of the command and its descendants.
Writes a time_verbose.txt compatible line:
    Maximum resident set size (kbytes): N

Peak RSS comes from getrusage(RUSAGE_CHILDREN).ru_maxrss after the child is reaped: the
kernel's high-water mark for the largest waited-for descendant. The previous implementation
sampled VmRSS of the direct child every 0.2 s, but the direct child is the `bash -c` wrapper
every caller uses, so it reported ~3 MB (peak_rss_gb=0.003) for every aligner.

Usage: measure_rss.py <output_file> -- <cmd> [args...]
"""
import resource
import subprocess
import sys


def main():
    if len(sys.argv) < 4 or sys.argv[2] != "--":
        print("Usage: measure_rss.py <output_file> -- <cmd> [args...]", file=sys.stderr)
        sys.exit(1)

    out_file = sys.argv[1]
    cmd = sys.argv[3:]

    proc = subprocess.Popen(cmd)
    returncode = proc.wait()
    # Linux reports ru_maxrss in KB. This process spawns exactly one child tree, so
    # RUSAGE_CHILDREN is that tree's maximum.
    max_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    with open(out_file, "w") as f:
        f.write(f"\tMaximum resident set size (kbytes): {max_rss}\n")
        f.write(f"\tExit status: {returncode}\n")

    sys.exit(returncode)


if __name__ == "__main__":
    main()
