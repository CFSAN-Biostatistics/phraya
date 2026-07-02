#!/usr/bin/env python3
"""
Count aligned reads from a .phraya.queries file.

.phraya.queries format: zstd-compressed msgpack of HashMap<String, Vec<(u32, f64)>>
A query appears in the map iff it placed at least one alignment with score ≥ 0.95.
Outputs count to stdout.

Usage: count_phraya_aligned.py <file.phraya.queries>
"""
import sys
import zstandard
import msgpack


def main():
    if len(sys.argv) != 2:
        print("Usage: count_phraya_aligned.py <file.phraya.queries>", file=sys.stderr)
        sys.exit(1)

    path = sys.argv[1]
    try:
        with open(path, "rb") as f:
            compressed = f.read()
        dctx = zstandard.ZstdDecompressor()
        raw = dctx.decompress(compressed)
        data = msgpack.unpackb(raw, raw=False)
        print(len(data))
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
