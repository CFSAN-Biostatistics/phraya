#!/usr/bin/env python3
"""Normalize a phraya variant TSV so it is a stable correctness oracle.

The `all_alleles` column (field 3) renders a HashMap<u8,u32>, whose iteration
order is randomized per process — so `C:3,T:2` and `T:2,C:3` denote the same
multiset but differ textually. Sorting the comma-separated tokens in that column
(and then the caller sorts lines) removes that noise, leaving a digest that
changes only if the actual variant data changes.

Reads TSV on stdin, writes normalized TSV on stdout.
"""
import sys


def main() -> None:
    for line in sys.stdin:
        fields = line.rstrip("\n").split("\t")
        if len(fields) > 2:
            fields[2] = ",".join(sorted(fields[2].split(",")))
        print("\t".join(fields))


if __name__ == "__main__":
    main()
