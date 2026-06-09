#!/usr/bin/env python3
"""
Parse /usr/bin/time -v output to JSON.

Extracts:
- wall_time_s: Elapsed wall clock time (minutes:seconds → seconds)
- peak_rss_gb: Maximum resident set size (kbytes → GB)
"""
import re
import sys
import json


def parse_time_output(text):
    """Parse /usr/bin/time -v output."""
    wall_match = re.search(r'Elapsed \(wall clock\) time.*?(\d+):(\d+\.\d+)', text)
    rss_match = re.search(r'Maximum resident set size \(kbytes\): (\d+)', text)

    if not wall_match:
        return {"error": "Could not parse wall time from /usr/bin/time output"}
    if not rss_match:
        return {"error": "Could not parse RSS from /usr/bin/time output"}

    # Wall time: MM:SS.ss → seconds
    minutes = int(wall_match.group(1))
    seconds = float(wall_match.group(2))
    wall_s = minutes * 60 + seconds

    # RSS: kbytes → GB
    rss_kb = int(rss_match.group(1))
    rss_gb = rss_kb / 1_048_576

    return {
        "wall_time_s": round(wall_s, 2),
        "peak_rss_gb": round(rss_gb, 3)
    }


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: parse_time.py <timing.txt>", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        text = f.read()

    result = parse_time_output(text)
    print(json.dumps(result, indent=2))
