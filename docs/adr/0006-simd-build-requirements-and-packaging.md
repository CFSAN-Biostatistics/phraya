# 6. SIMD build requirements and packaging (`target-cpu` / `ensure_simd`)

- **Status**: Proposed
- **Date**: 2026-07-01

## Context

Phraya has two independent SIMD stories, and they have different build requirements:

1. **Match extension** (WFA/Myers inner loop, [ADR-0002](0002-simd-match-extension.md)) uses
   the SSE2 (x86-64) and NEON (AArch64) *baselines*, selected at compile time. These are
   mandatory on their architectures, so this path needs **no special build flags** and always
   compiles.
2. **k-mer sketching** (`simd-minimizers`) wants AVX2 on x86-64 / NEON on AArch64 for full
   throughput (~2× on sketching). Its transitive dependency `ensure_simd` **hard-errors the
   build** unless AVX2/native is enabled, *or* its `scalar` feature is activated:

   ```
   error: … uses AVX2 (on x64) or NEON (on aarch64) SIMD instructions …
   To get the expected performance, compile/install using e.g.:
   RUSTFLAGS="-C target-cpu=native" cargo …
   Alternatively, silence this error by activating the `scalar` feature.
   ```

This surfaced during the `perf/align-hot-path` work: a plain `cargo build --release` **fails to
compile** with no `RUSTFLAGS` set. (Debug/test builds happened to succeed only because an
`ensure_simd` artifact compiled earlier with the flag was cached in `target/`.) So the
requirement is currently implicit — it works on machines/CI that already export the flag and
breaks silently on those that don't.

Current packaging (see `README.md`):
- Prebuilt releases ship a **native** binary (`-C target-cpu=x86-64-v3`, AVX2) and a **portable**
  binary (`-C target-feature=+sse4.2`, broad compatibility).
- The Docker image is built on the **SSE4.2 baseline** for portability.

The portable *distribution* path therefore exists, but building it still has to satisfy
`ensure_simd` — the SSE4.2 flags do **not** enable AVX2, so a portable/Docker build must either
raise the floor to `x86-64-v3` (losing the "runs anywhere" property) or activate the `scalar`
feature (slower sketching). That interaction has not been written down or deliberately chosen.

## Decision

*Proposed — not yet decided.* This ADR records the constraint and the open question so it is not
rediscovered by accident. The options to weigh:

1. **Require and document `RUSTFLAGS` for all from-source builds** (`-C target-cpu=native` for
   local/HPC, `-C target-cpu=x86-64-v3` for reproducible/distributable). Make CI and the
   Dockerfile set it explicitly so the build never depends on ambient environment.
2. **Activate the `ensure_simd`/`simd-minimizers` `scalar` feature for the portable + Docker
   build** so it compiles on the SSE4.2 baseline without AVX2, accepting slower sketching there.
3. **Expose a Phraya cargo feature** that selects (1) vs (2), wired into the release matrix so
   "native" and "portable" are explicit build modes rather than an implicit `RUSTFLAGS` contract.

## Consequences / open questions to think through

- **Docker reproducibility**: `-C target-cpu=native` compiles for the *builder's* CPU, so a
  `native` image is not reproducible and may `SIGILL` on older hardware. A distributable image
  must pin a concrete level (`x86-64-v3`) or use the `scalar` feature — not `native`.
- **Runtime SIGILL risk**: an AVX2 binary run on a pre-2013 CPU crashes with an illegal
  instruction. The portable tier exists to avoid this; the build story must keep it genuinely
  AVX2-free (hence the `scalar`-feature question).
- **CI fragility**: without an explicit flag/feature in CI config, the build depends on ambient
  `RUSTFLAGS` and can break on a clean runner. This is the concrete breakage that prompted the ADR.
- **Sketching speed on the portable tier**: the `scalar` path is ~2× slower at sketching; that is
  the cost of the broadest-compatibility image and should be a conscious choice.
- Whatever is chosen should be captured in the Dockerfile, CI config, and `README.md` so the
  requirement is explicit rather than an implicit environment contract.

## Alternatives considered

- **Leave it implicit** (status quo): rejected as the thing to *keep* — it produces silent,
  environment-dependent build failures. Documenting the constraint is the minimum.
- **Vendor-patch `ensure_simd` to drop the guard**: rejected — fights the upstream contract and
  would mask genuinely slow scalar sketching instead of making the tradeoff explicit.
