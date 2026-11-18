# syntax=docker/dockerfile:1

ARG RUST_VERSION=1.99
ARG DEBIAN_RELEASE=bookworm

# ---- build -----------------------------------------------------------------
FROM rust:${RUST_VERSION}-slim-${DEBIAN_RELEASE} AS build
WORKDIR /src

# Fetch and compile dependencies in their own layer so source edits stay fast.
COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo 'fn main() {}' > src/main.rs && touch src/lib.rs \
 && cargo build --release --locked \
 && rm -rf src

COPY src ./src
COPY data ./data
COPY tests ./tests
RUN touch src/main.rs src/lib.rs && cargo build --release --locked

# ---- test (opt-in: docker build --target test .) ---------------------------
FROM build AS test
RUN rustup component add rustfmt clippy \
 && cargo fmt --check \
 && cargo clippy --all-targets --locked -- -D warnings \
 && cargo test --release --locked

# ---- runtime ---------------------------------------------------------------
FROM debian:${DEBIAN_RELEASE}-slim AS runtime
RUN useradd --system --uid 10001 --no-create-home --shell /usr/sbin/nologin app
COPY --from=build /src/target/release/acetyl-scan /usr/local/bin/acetyl-scan
USER app
WORKDIR /data
ENTRYPOINT ["acetyl-scan"]
CMD ["--help"]
