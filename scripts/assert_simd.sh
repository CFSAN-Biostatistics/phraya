#!/usr/bin/env bash
# Anti-fraud gate: prove the WFA diagonal fill actually compiles to SIMD.
#
# Disassembles the `fill_simd` symbol in phraya-align and asserts the function
# body contains target-native vector instructions (NEON on aarch64, SSE/AVX on
# x86_64). A scalar delegation — like the fake NEON/SSE paths this replaced —
# contains none, and fails this check.
#
# Builds NATIVELY for the host arch, so the x86 CI job verifies SSE and the
# aarch64 CI job verifies NEON, with no cross-compilation required.
set -euo pipefail

cd "$(dirname "$0")/.."

ARCH="$(uname -m)"
case "$ARCH" in
    aarch64 | arm64)
        # NEON: signed-min / compare-equal / 4x32-bit & 16x8-bit vector forms.
        PATTERN='smin|umin|cmeq|\.4s|\.16b'
        LABEL="NEON (aarch64)"
        ;;
    x86_64 | amd64)
        # SSE/AVX: packed signed min, packed compare-equal, vector registers.
        PATTERN='pminsd|pminsw|pcmpeq|xmm|ymm'
        LABEL="SSE/AVX (x86_64)"
        ;;
    *)
        echo "assert_simd: unsupported host arch '$ARCH'" >&2
        exit 2
        ;;
esac

OBJDUMP="${OBJDUMP:-llvm-objdump}"
if ! command -v "$OBJDUMP" >/dev/null 2>&1; then
    echo "assert_simd: '$OBJDUMP' not found (set OBJDUMP=... or install llvm)" >&2
    exit 2
fi

echo "assert_simd: building phraya-align (release, symbols kept)…"
# -C target-cpu=native: required by simd-minimizers' ensure_simd guard in
#   release, and lets `wide` lower to the widest native vectors (AVX2 / NEON).
# strip=false so the symbol survives; opt so wide lowers to single instructions.
# lto=false so the rlib carries real machine code (thin-LTO rlibs are
# bitcode-only, deferring codegen and leaving nothing to disassemble).
export RUSTFLAGS="${RUSTFLAGS:-} -C target-cpu=native"
cargo build -p phraya-align --release \
    --config 'profile.release.strip=false' \
    --config 'profile.release.lto=false' >&2

RLIB="$(find target/release/deps -name 'libphraya_align-*.rlib' -print -quit)"
if [[ -z "${RLIB:-}" ]]; then
    echo "assert_simd: could not locate libphraya_align rlib" >&2
    exit 2
fi
echo "assert_simd: disassembling $RLIB"

# Isolate the fill_simd function body (from its label to the next label).
BODY="$("$OBJDUMP" -d --demangle "$RLIB" 2>/dev/null | awk '
    /<.*fill_simd.*>:/ { cap = 1; next }
    cap && /<.*>:/      { cap = 0 }
    cap                 { print }
')"

if [[ -z "$BODY" ]]; then
    echo "assert_simd: FAIL — no fill_simd symbol found in disassembly" >&2
    exit 1
fi

COUNT="$(printf '%s\n' "$BODY" | grep -Eci "$PATTERN" || true)"
if [[ "$COUNT" -eq 0 ]]; then
    echo "assert_simd: FAIL — fill_simd contains no $LABEL vector instructions." >&2
    echo "             Expected one of: $PATTERN" >&2
    echo "             This is what a scalar fake looks like. Disassembly excerpt:" >&2
    printf '%s\n' "$BODY" | head -25 >&2
    exit 1
fi

echo "assert_simd: PASS — fill_simd contains $COUNT $LABEL vector instruction(s)."
