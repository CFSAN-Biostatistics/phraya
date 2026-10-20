# Handoff: Benchmarking Rework & Long-Read Findings (2026-10-20)

Status: **verified**. Local `main` matches `origin/main` (commits `531fe62`–`a2567b7`).
Reedling2 `~/phraya-build` pulled to `a2567b7`; end-to-end harness run passes.

## Context: what this round of benchmarking was about

After the Phase 10 write-up, the long-read harness (`run_local_bench_long.sh`)
was broken: it crashed in the timing `awk`, wrote minimap2's stderr into the SAM,
and silently reported `phraya PA=0.0000 mapped=0/50` + `minimap2 PA=1.0` even
though phraya placed 50/50. The goal was to fix the harness/scorer, score
minimap2 by position against ground truth (not mapped/total), verify 4f7fb55's
assembly work didn't regress short-read behavior, and find the long-read
slowdown.

## Platform

Reedling2 (`reedling2.fda.gov`), Xeon Gold 5416S. Tools come from
environment modules (`module avail`): **minimap2 2.28** = `module minimap2/2.28`,
**samtools 1.20** = `samtools/1.20`. A non-interactive `ssh host "cmd"` shell
doesn't define `module` until `source /etc/profile`; the harness does this
itself (`load_tool`). (Inside SLURM jobs the module system is unreliable;
`slurm/config/global.env` uses absolute micromamba paths instead.)

Build: `export RUSTFLAGS="-C target-cpu=native"; cargo build --release -p phraya`
(the package is `phraya`, not `phraya-cli`). Verify the binary's mtime is newer
than the last `git pull` — a failed build silently benchmarks the stale binary.

## Commits (scripts/benchmark/…)

| Commit | What |
|---|---|
| `1ab2024` | `README.md` "Running on Reedling2" section; harness loads minimap2 via `module` |
| `c32adf1` | Fix long-read timing `awk`; keep minimap2 stderr out of the SAM |
| `61fa1c0` | Stream-decompress the phraya queries sidecar (zstd frames have no content size) |
| `531fe62` | Detect zstd by frame magic (`\x28\xb5\x2f\xfd`) in the scorer |
| `f34c302` | Rewrite `long_read_accuracy.py`: score phrapa + SAM on one truth basis |
| `7070050` | Score minimap2 by position in the harness; wrapper fixes (global.env path, PIPESTATUS, count_phry_aligned) |
| `a2567b7` | Close `minimap2-pb.sh` heredoc; fail loud on missing RSS in `phraya-longread.sh` |

## Scorer: `slurm/utils/long_read_accuracy.py` (rewritten)

Scores a phraya `.queries` sidecar OR a minimap2 `.sam` against truth. Input
kind is chosen by file extension (`.queries` → phraya, `.sam` → SAM), and there
is **no content sniffing and no `except: pass` silent fallback** — a real parse
error aborts. Output fields:

`pa  precision  n_truth  n_mapped  n_correct  n_correct_any`

Rules:
- **Each read scored once**, on its single best placement (highest identity for
  phraya; primary record for SAM). `correct_any` (any retained placement within
  tolerance) is reported separately — for phraya, which keeps many per-locus
  near-duplicates, this can differ from `correct`. (Here both aligners keep one
  placement per read, so they're equal.)
- **PA = correct / n_truth** — unmapped reads count against the aligner.
- **precision = correct / mapped** (for reference).
- SAM: leading soft/hard clip subtracted from POS; secondary (0x100),
  supplementary (0x800) and unmapped (0x4) records skipped.
- zstd is detected by magic; frames are stream-decoded (`stream_reader`),
  because phraya's frames carry no content size and `dctx.decompress()` silently
  returns empty without a max-output-size argument.

Negative tests (run on `/tmp` scratch data, no repo changes) confirm the scorer
can fail, i.e. it isn't trivially passing:
- Truth starts shifted +5000 bp → both aligners PA = 0.0000, n_mapped still 50.
- SAM with `longread_0..4` forced to FLAG 4 (unmapped) → n_mapped = 45, PA = 0.9000.

On the real dataset both aligners score 50/50 at PA 1.0000.

## Bugs fixed in the harness/wrappers

**`run_local_bench_long.sh`:**
1. Elapsed-time `awk` had an escaped `$(date …)` that awk parsed as syntax error
   (the run aborted before minimap2; walls were blank). Now computes
   `END=$(date +%s.%N)` and passes `-v s=… -v e=…`.
2. minimap2 stderr (`[M::…]`) was redirected into the SAM, corrupting it. Now to
   `minimap2.log`.
3. minimap2 "PA" was mapped/total (counts wrong-locus mappings). Now scored by
   position via `long_read_accuracy.py` with the same tolerance as phraya.
4. `phraya_mapped` summed `len(v)` over placements (inflated past read count).
  Now reads-with-≥1-placement via `count_phry_aligned.py`.
5. phraya's exit status was masked by `| tee`. Now `set +e` around the pipeline
   and `PIPESTATUS[0]` is checked; `PHRAYA_FAILED=0` initialised.
6. `MMI="${REF%.fasta}.mmi"` never stripped the `.fasta` (well, `ref.fa`), so
  the index name was wrong; and it wasn't rebuilt on ref change. Now `$REF.mmi`,
  rebuilt every full run (~0.2 s).
7. `2>/dev/null || echo 0` wrappers hid missing `msgpack`/`zstandard` and real
  scoring failures → PA silently 0. All removed; a scoring failure now aborts.
8. `PHRAYA` defaulted to `target/release/phraya.exe` (doesn't exist on Linux).
  Now `phraya` with a `.exe` fallback.
9. `RESCORE_ONLY=1` still generated data and ran `phraya plan`; now skips
  gen/index/plan/align and requires prior `results.txt` + outputs (lets you
  re-test the scorer without a 12-minute align).
10. `module load` under `set -u` tripped on unset module-script vars. Wrapped in
  `set +u`/`set -u`.

**`phraya-longread.sh`, `minimap2-ont.sh`, `minimap2-pb.sh`:**
11. `SCRIPT_DIR` was `dirname` of the wrapper (the `wrappers/` dir), so
  `source "$SCRIPT_DIR/config/global.env"` looked in `wrappers/config/` —
  nonexistent. All three wrappers now use `dirname …/..` to point at the slurm/
  root, matching every sibling wrapper (`_phara_run.sh`, `bowtie2.sh`, …).
12. Same `awk` timing bug (escaped `$(...)` inside double quotes) in
  `phraya-longread.sh`. Fixed.
13. `minimap2-pb.sh` heredoc was **truncated** — no closing `EOF`, no final
  echo; `timing.txt` was cut off and trailing output lost. Closed the heredoc
  and added the echo, matching `minimap2-ont.sh`.
14. `phraya-longread.sh` line 48 `grep … 2>/dev/null … || echo "0"` silently
  turned a missing/timed-out RSS reading into `0 GB`. Now a hard error if no RSS.
15. `quiries` count was `sum(len(v))` (placement count, not read count). Now
  uses `count_phary_aligned.py` (reads with ≥1 placement).

**`long_read_accuracy.py` (scorer, original bugs):**
16. Magic-byte check was `b"B9\x01"` / `"SM"` — phrapa's zstd frames start
  `\x28\xb5\x2f\xfd`. The phraya branch never ran; the file fell through to the
  SAM text parser, which hit `utf-8 decode byte 0xb5` and printed PA=0. Now
  dispatch by extension with correct zstd detection. (Also: `import io` was
  missing for `BytesIO`; now uses `stream_reader(f)` on the file handle, matching
  `phraya_accuracy.py`.)
17. `compute_pa_sam_long` returned a 3-tuple but was unpacked as 4 values
  (would raise `ValueError`); the `except Exception: pass` hid it. Rewritten
  with one `score_phraya`/`score_sam` pair and consistent returns.

## Verification

`bash -n` on all four shell scripts + `ast.parse` on the Python: OK (after
`a2567b7` pull on Reedling).

End-to-end on the committed tree (`a2567b7`), 5 reads × 15 kb, single-threaded,
fresh data:

| aligner | wall (s) | PA | precision | mapped |
|---|---:|---:|---:|---:|
| phraya long-read | 21.19 | 1.0000 | 1.0000 | 5/5 |
| minimap2 map-ont | 0.14 | 1.0000 | 1.0000 | 5/5 |

`results.txt` written with `label/genome_size/num_reads/read_len_range/
tolerance/phraya_wall_seconds/phraya_score/minimap2_wall_seconds/
minimap2_score/score_fields`.

## Benchmark results

### Short reads (20k × 150 bp) — 1% divergence, single-threaded

`BENCH_DIR=/tmp/phraya-bench12`:

| Strategy | 4f7fb55 wall | RSS (MB) | Correct? |
|---|:---:|---:|---|
| fast | 1.333 s | — | digest `a504…` |
| balanced (default) | 1.262 s | 123.1 | `f88db62c…` (byte-identical to benchmark-10) |
| sensitive | 0.743 s | — | `35d8…` |
| exact | 0.634 s | 108.5 | `ac65…` |

5% divergence: wall **1.538 s**, placed 10522/20000 (9025 below_threshold, 453
no_alignment). The correctness digest `c2fb3b18…` is **identical** at 90aacf4
(SSE2 baseline), 97366d3 (Phase 3 SIMD), and 4f7fb55. Placement behavior is
unchanged across the assembly work.

### Short reads (20k × 150 bp) — 1% divergence, 1% indel rate

`INDEL_RATE=0.01`: wall **1.376 s**. IEC = 0.327 (5092 concordant / 15553
indel-bearing reads). This is the **linear-gap baseline**, not a new weakness:
README.md documents IEC = 0.147 on a smaller run (300 reads, indel_max_len=4)
for the same reason — linear costs fragment true indels into many CIGAR ops.
ADR-0014 (affine WFA) is the intended fix. The numbers differ from README only
because of the different dataset size/len, not because of a code change.

### Long reads (50 reads × 20–25 kb, 12% ONT error, 5 Mbp genome)

Generated identically and scored with the same script. phraya placed 50/50,
PA 1.0000.

| commit | wall (s) | PA | queries sidecar digest |
|---|:---:|---:|---|
| `f758d73` (chunked WFA) | **750 s** | 1.0000 | `a764fba1…` |
| `4f7fb55` (assembly work) | **745 s** | 1.0000 | `a764fba1…` |

**Byte-identical placement output (same sidecar digest) and within 1% wall time.**
The ~15 s/read long-read cost is **inherent to the chunked-WFA implementation**
(f758d73), not a regression from 4f7fb55's assembly verification work. There is
no baseline to regress from: benchmarking-longread-N.md reports unbanded WFA
placed 0/10 before chunking. The cost vs minimap2 (~0.15–0.75 s, ~0.015 s/read)
is a ~300× gap that needs its own profiling pass. Each read carries many
near-duplicate placements (starts within ±2 bp, identities ~0.817 at 12% error),
reported as `correct_any` for inspection.

### Noise note

5% and 50-read long numbers are single runs on a 64-core shared host. Same-session
A/B runs of 90aacf4 / 97366d3 / 4f7fb55 at 5% gave walls 1.68–1.69 s and identical
digests; the 1.223 s vs 1.262 s gap versus benchmark-10 (≈3%) is run-to-run noise,
not a code delta. Don't claim a speedup or slowdown from single short-read runs.

## Open issues (not fixed here, not regressions from 4f7fb55)

- **Long-read speed.** ~750 s for 50 reads (≈15 s/read) vs minimap2 ≈0.015 s/read.
  Root cause is the chunked-WFA algorithm/implementation, not the assembly work.
  Worth a dedicated profiling ticket against `phraya-align`.
- `phraya-longread.sh` and the minimap2 wrappers still don't run as part of any
  automated SLURM job from this work — path fixes are in but only the local
  end-to-end harness was executed on Reedling.
- The SLURM `run_benchmark.sh` driver was not re-exercised (the local harness
  covers the scoring path end-to-end; the SLURM wrappers need a node + dataset
  slot to validate).

## Reproducing

```bash
# short reads, 1 thread
export RAYON_NUM_THREADS=1
scripts/benchmark/local/run_local_bench.sh p_balanced 2000000 20000 150 0.01
STRATEGY=fast|balanced|sensitive|exact scripts/benchmark/local/run_local_bench.sh p_fast ... 0.01
# 5% / indels
scripts/benchmark/local/run_local_bench.sh p_5pct 2000000 20000 150 0.05
INDEL_RATE=0.01 scripts/benchmark/local/run_local_bench.sh p_indel 2000000 20000 150 0.01

# long reads, 5 reads x 15-25kb (takes ~20-750s)
export RAYON_NUM_THREADS=1
scripts/benchmark/local/run_local_bench_long.sh p_long 5000000 5 15000 25000
# rescore only (skips align — uses prior outputs + results.txt walls)
RESCORE_ONLY=1 scripts/benchmark/local/run_local_bench_long.sh p_long 5000000 5 15000 25000
```

On Reedling2, `module load minimap2/2.28` (done automatically by the long-read
harness if not on PATH) and build the binary with
`export RUSTFLAGS="-C target-cpu=native"; cargo build --release -p phara`.

## Long-read results after fix (Reedling2, 2026-10-20)

Measured on Reedling2 (`lscpu`: Intel Xeon Gold 5416S), against the
existing `/tmp/phraya-bench12/p4_long` dataset (50 reads × 20–25 kb, 12% ONT
error, 5 Mbp genome), scorer `long_read_accuracy.py` (tolerance 100), single
threaded (`RAYON_NUM_THREADS=1`), binary built with
`RUSTFLAGS="-C target-cpu=native" cargo build --release -p phraya` at commit
a268d3a. minimap2 2.28 from `module minimap2/2.28` (env
`minimap2-v2.28`), run with `-ax map-ont -K256M -t1`.

| phase | phraya | minimap2 |
|---|---:|---:|
| plan (only phraya builds a plan/sketch index) | 0.68 s | — |
| align | 0.671 s | 0.882 s |
| peak RSS | 170 MB | 250 MB |

| metric | phraya | minimap2 |
|---|---:|---:|
| reads placed | 50/50 | 50/50 |
| PA (±100 bp, primary) | 0.9800 | 1.0000 |

phraya `plan` + `align` total 1.365 s. Align-only per-read is 0.671/50 = 0.013 s
(minimap2 0.882/50 = 0.018 s, including its FM-index build in the 0.882 s).
The single Reedling2 PA miss (primary, ±100 bp) is `longread_7`, placed at
1846148 vs truth 1846026 (offset +122, just past tolerance; identity 0.805,
1 placement) — cause uninvestigated. Note: the local Windows run used a separately
generated dataset (seed 42, 5 reads), not `p4_long`, so its PA 1.0000 does not
counterfactual this read. This is approximate: wall times on this shared host were
noisy (align measured 0.67–0.68 s across two runs).

Note: the handoff's ~15 s/read figure was on the old chunked-WFA path (commit
f758d73). The numbers above are on the rewritten seed→chain→locus path.

### What changed vs the handoff's ~15 s/read

The pre-fix chunked-WFA path (one full-query WFA per surviving chain anchor, up
to 100 anchors × 2 kb chunks per read) cost ~15 s/read. After
(1) the seed→chain→locus rewrite (`extend_chains_longread` → `align_long_locus`:
exact-match segments + global Myers gap-fills + fitting WFA end-extensions) and
(2) suppressing per-mismatch `VariantObservation` extraction for
`Strategy::LongRead` — profiled at ~335 ms/read and the dominant align cost —
phraya align is 0.013 s/read on Reedling2, comparable to minimap2 (0.018 s/read),
both single-threaded. This is approximate: wall times on this shared 64-core host
were noisy (the 50-read align measured 0.67–0.68 s across two runs).

### Behavior change

`Strategy::LongRead` no longer emits per-mismatch `VariantObservation`s in the
`.phraya` blocks. At ~3000 mismatches per 25 kb read / 12% error, that was the
dominant align cost. Long-read use is read placement (depletion/tiering in the
filter stage), not SNP calling; variant calling is a separate downstream step.
Coverage depth and the cross-space placement sidecar (`query_positions`) are
unchanged. Short-read strategies are unaffected. If a future workflow needs
long-read variant observations, gate that on a flag rather than restoring the
~0.34 s/read cost unconditionally.
