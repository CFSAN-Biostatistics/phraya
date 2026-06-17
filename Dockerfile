FROM rust:1.75-slim AS builder

WORKDIR /build

# Cache dependencies layer
COPY Cargo.toml Cargo.lock ./
COPY phraya-core/Cargo.toml phraya-core/
COPY phraya-align/Cargo.toml phraya-align/
COPY phraya-io/Cargo.toml phraya-io/
COPY phraya-filter/Cargo.toml phraya-filter/
COPY phraya-cli/Cargo.toml phraya-cli/

# Stub source files to warm the dependency cache
RUN for crate in phraya-core phraya-align phraya-io phraya-filter phraya-cli; do \
      mkdir -p ${crate}/src && echo "fn main() {}" > ${crate}/src/lib.rs; \
    done && \
    echo "fn main() {}" > phraya-cli/src/main.rs && \
    cargo build --release --bin phraya 2>/dev/null || true && \
    for crate in phraya-core phraya-align phraya-io phraya-filter phraya-cli; do \
      rm -f ${crate}/src/lib.rs; \
    done && \
    rm -f phraya-cli/src/main.rs

# Build the real binary (portable SIMD — runtime CPU unknown in container)
COPY phraya-core/src phraya-core/src
COPY phraya-align/src phraya-align/src
COPY phraya-io/src phraya-io/src
COPY phraya-filter/src phraya-filter/src
COPY phraya-cli/src phraya-cli/src

RUN cargo build --release --bin phraya

# ── Minimal runtime image ─────────────────────────────────────────────────
FROM gcr.io/distroless/cc-debian12

COPY --from=builder /build/target/release/phraya /usr/local/bin/phraya

ENTRYPOINT ["/usr/local/bin/phraya"]
