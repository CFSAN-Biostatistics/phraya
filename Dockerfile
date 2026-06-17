# Build stage
FROM rust:1.75-slim AS builder

WORKDIR /build

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy workspace manifests first for layer caching
COPY Cargo.toml Cargo.lock ./
COPY phraya-core/Cargo.toml phraya-core/
COPY phraya-align/Cargo.toml phraya-align/
COPY phraya-io/Cargo.toml phraya-io/
COPY phraya-filter/Cargo.toml phraya-filter/
COPY phraya-cli/Cargo.toml phraya-cli/

# Create stub lib/main files so cargo can resolve the workspace
RUN mkdir -p phraya-core/src phraya-align/src phraya-io/src phraya-filter/src phraya-cli/src \
    && echo "fn main() {}" > phraya-cli/src/main.rs \
    && for crate in phraya-core phraya-align phraya-io phraya-filter; do \
         echo "" > $crate/src/lib.rs; \
       done

# Pre-fetch dependencies (portable SIMD — SSE4.2 baseline, no native CPU features)
RUN RUSTFLAGS="-C target-feature=+sse4.2" cargo build --release --bin phraya 2>/dev/null || true

# Copy full source and build for real
COPY . .

# Touch to force rebuild of non-stub sources
RUN touch phraya-core/src/lib.rs phraya-align/src/lib.rs phraya-io/src/lib.rs \
         phraya-filter/src/lib.rs phraya-cli/src/main.rs

# Portable SSE4.2 baseline build — runtime CPU is unknown in containers
RUN RUSTFLAGS="-C target-feature=+sse4.2" cargo build --release --bin phraya

# Runtime stage — minimal distroless image
FROM gcr.io/distroless/cc-debian12

COPY --from=builder /build/target/release/phraya /usr/local/bin/phraya

ENTRYPOINT ["/usr/local/bin/phraya"]
