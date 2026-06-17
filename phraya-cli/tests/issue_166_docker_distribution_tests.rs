/// Acceptance tests for Issue #166: Docker distribution channel
/// Tests verify the existence and correctness of Dockerfile, .dockerignore,
/// release.yml (multi-arch buildx), and README Docker section.
/// All tests must fail (RED) until the implementation is complete.
use std::fs;
use std::path::PathBuf;

/// Resolve the workspace root from the phraya-cli manifest directory.
fn workspace_root() -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest_dir)
        .parent()
        .expect("phraya-cli has a parent workspace root")
        .to_path_buf()
}

// ============================================================================
// SECTION A: Dockerfile existence and multi-stage structure
// ============================================================================

/// Test: Dockerfile exists at the workspace root
#[test]
fn issue_166_dockerfile_exists() {
    let dockerfile = workspace_root().join("Dockerfile");
    assert!(
        dockerfile.exists(),
        "Dockerfile must exist at workspace root (got path: {})",
        dockerfile.display()
    );
}

/// Test: Dockerfile uses rust:1.75-slim as the builder stage base image
#[test]
fn issue_166_dockerfile_builder_stage_uses_rust_1_75_slim() {
    let dockerfile = workspace_root().join("Dockerfile");
    let content = fs::read_to_string(&dockerfile)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerfile.display()));

    assert!(
        content.contains("FROM rust:1.75-slim"),
        "Dockerfile must use 'FROM rust:1.75-slim' as the builder base image; got:\n{}",
        content
    );
}

/// Test: Dockerfile uses gcr.io/distroless/cc-debian12 as the runtime stage
#[test]
fn issue_166_dockerfile_runtime_stage_uses_distroless() {
    let dockerfile = workspace_root().join("Dockerfile");
    let content = fs::read_to_string(&dockerfile)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerfile.display()));

    assert!(
        content.contains("FROM gcr.io/distroless/cc-debian12"),
        "Dockerfile must use 'FROM gcr.io/distroless/cc-debian12' as the runtime image; got:\n{}",
        content
    );
}

/// Test: Dockerfile is multi-stage (has at least two FROM instructions)
#[test]
fn issue_166_dockerfile_is_multistage() {
    let dockerfile = workspace_root().join("Dockerfile");
    let content = fs::read_to_string(&dockerfile)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerfile.display()));

    let from_count = content
        .lines()
        .filter(|line| line.trim_start().starts_with("FROM "))
        .count();

    assert!(
        from_count >= 2,
        "Dockerfile must have at least 2 FROM instructions (multi-stage build), found {}",
        from_count
    );
}

// ============================================================================
// SECTION B: Portable SSE4.2 baseline build (no target-cpu=native)
// ============================================================================

/// Test: Dockerfile does not use -C target-cpu=native (must use portable SSE4.2 baseline)
#[test]
fn issue_166_dockerfile_no_target_cpu_native() {
    let dockerfile = workspace_root().join("Dockerfile");
    let content = fs::read_to_string(&dockerfile)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerfile.display()));

    assert!(
        !content.contains("target-cpu=native"),
        "Dockerfile must NOT use '-C target-cpu=native'; \
         Docker images require a portable baseline (SSE4.2). Got:\n{}",
        content
    );
}

/// Test: Dockerfile sets RUSTFLAGS for SSE4.2 as the SIMD baseline
#[test]
fn issue_166_dockerfile_uses_sse42_baseline() {
    let dockerfile = workspace_root().join("Dockerfile");
    let content = fs::read_to_string(&dockerfile)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerfile.display()));

    assert!(
        content.contains("sse4.2") || content.contains("target-feature=+sse4.2"),
        "Dockerfile must set RUSTFLAGS with SSE4.2 as the portable SIMD baseline; got:\n{}",
        content
    );
}

// ============================================================================
// SECTION C: .dockerignore content
// ============================================================================

/// Test: .dockerignore exists at the workspace root
#[test]
fn issue_166_dockerignore_exists() {
    let dockerignore = workspace_root().join(".dockerignore");
    assert!(
        dockerignore.exists(),
        ".dockerignore must exist at workspace root (got path: {})",
        dockerignore.display()
    );
}

/// Test: .dockerignore excludes target/ (build artifacts)
#[test]
fn issue_166_dockerignore_excludes_target_dir() {
    let dockerignore = workspace_root().join(".dockerignore");
    let content = fs::read_to_string(&dockerignore)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerignore.display()));

    assert!(
        content.lines().any(|line| {
            let l = line.trim();
            l == "target/" || l == "target" || l == "/target/" || l == "/target"
        }),
        ".dockerignore must exclude 'target/' to keep image build context small; got:\n{}",
        content
    );
}

/// Test: .dockerignore excludes .git/
#[test]
fn issue_166_dockerignore_excludes_git_dir() {
    let dockerignore = workspace_root().join(".dockerignore");
    let content = fs::read_to_string(&dockerignore)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerignore.display()));

    assert!(
        content.lines().any(|line| {
            let l = line.trim();
            l == ".git/" || l == ".git" || l == "/.git/" || l == "/.git"
        }),
        ".dockerignore must exclude '.git/' to keep image build context small; got:\n{}",
        content
    );
}

/// Test: .dockerignore excludes *.md files
#[test]
fn issue_166_dockerignore_excludes_markdown_files() {
    let dockerignore = workspace_root().join(".dockerignore");
    let content = fs::read_to_string(&dockerignore)
        .unwrap_or_else(|_| panic!("cannot read {}", dockerignore.display()));

    assert!(
        content.lines().any(|line| {
            let l = line.trim();
            l == "*.md" || l == "**/*.md"
        }),
        ".dockerignore must exclude '*.md' files; got:\n{}",
        content
    );
}

// ============================================================================
// SECTION D: release.yml workflow with docker buildx
// ============================================================================

/// Test: release.yml workflow file exists
#[test]
fn issue_166_release_workflow_exists() {
    let release_yml = workspace_root()
        .join(".github")
        .join("workflows")
        .join("release.yml");
    assert!(
        release_yml.exists(),
        ".github/workflows/release.yml must exist (got path: {})",
        release_yml.display()
    );
}

/// Test: release.yml uses docker/setup-buildx-action for multi-arch builds
#[test]
fn issue_166_release_workflow_uses_buildx() {
    let release_yml = workspace_root()
        .join(".github")
        .join("workflows")
        .join("release.yml");
    let content = fs::read_to_string(&release_yml)
        .unwrap_or_else(|_| panic!("cannot read {}", release_yml.display()));

    assert!(
        content.contains("docker/setup-buildx-action"),
        "release.yml must use docker/setup-buildx-action for multi-arch builds; got:\n{}",
        content
    );
}

/// Test: release.yml targets both amd64 and arm64 platforms
#[test]
fn issue_166_release_workflow_targets_multiarch() {
    let release_yml = workspace_root()
        .join(".github")
        .join("workflows")
        .join("release.yml");
    let content = fs::read_to_string(&release_yml)
        .unwrap_or_else(|_| panic!("cannot read {}", release_yml.display()));

    assert!(
        content.contains("linux/amd64") && content.contains("linux/arm64"),
        "release.yml must target both linux/amd64 and linux/arm64 platforms; got:\n{}",
        content
    );
}

/// Test: release.yml pushes to ghcr.io/cfsan-biostatistics/phraya:latest
#[test]
fn issue_166_release_workflow_pushes_latest_tag() {
    let release_yml = workspace_root()
        .join(".github")
        .join("workflows")
        .join("release.yml");
    let content = fs::read_to_string(&release_yml)
        .unwrap_or_else(|_| panic!("cannot read {}", release_yml.display()));

    assert!(
        content.contains("ghcr.io/cfsan-biostatistics/phraya") && content.contains(":latest"),
        "release.yml must push to ghcr.io/cfsan-biostatistics/phraya:latest; got:\n{}",
        content
    );
}

/// Test: release.yml pushes versioned tags (not only :latest)
#[test]
fn issue_166_release_workflow_pushes_versioned_tags() {
    let release_yml = workspace_root()
        .join(".github")
        .join("workflows")
        .join("release.yml");
    let content = fs::read_to_string(&release_yml)
        .unwrap_or_else(|_| panic!("cannot read {}", release_yml.display()));

    // Versioned tag: a ref to the git tag version (e.g. ${{ github.ref_name }} or semver pattern)
    let has_versioned_tag = content.contains("ref_name")
        || content.contains("github.event.release.tag_name")
        || content.contains("semver")
        || (content.contains("tags:") && content.contains("ghcr.io/cfsan-biostatistics/phraya"));

    assert!(
        has_versioned_tag,
        "release.yml must push versioned tags (e.g. using github.ref_name or semver); got:\n{}",
        content
    );
}

/// Test: release.yml logs in to GHCR (GitHub Container Registry)
#[test]
fn issue_166_release_workflow_logs_in_to_ghcr() {
    let release_yml = workspace_root()
        .join(".github")
        .join("workflows")
        .join("release.yml");
    let content = fs::read_to_string(&release_yml)
        .unwrap_or_else(|_| panic!("cannot read {}", release_yml.display()));

    assert!(
        content.contains("docker/login-action") && content.contains("ghcr.io"),
        "release.yml must use docker/login-action to log in to ghcr.io; got:\n{}",
        content
    );
}

// ============================================================================
// SECTION E: README Docker section
// ============================================================================

/// Test: README.md contains a Docker section heading
#[test]
fn issue_166_readme_has_docker_section() {
    let readme = workspace_root().join("README.md");
    let content = fs::read_to_string(&readme)
        .unwrap_or_else(|_| panic!("cannot read {}", readme.display()));

    assert!(
        content.contains("## Docker") || content.contains("# Docker"),
        "README.md must have a '## Docker' section; got no Docker heading in README"
    );
}

/// Test: README Docker section references ghcr.io/cfsan-biostatistics/phraya image
#[test]
fn issue_166_readme_docker_section_references_image() {
    let readme = workspace_root().join("README.md");
    let content = fs::read_to_string(&readme)
        .unwrap_or_else(|_| panic!("cannot read {}", readme.display()));

    assert!(
        content.contains("ghcr.io/cfsan-biostatistics/phraya"),
        "README.md Docker section must reference 'ghcr.io/cfsan-biostatistics/phraya'; got no such reference"
    );
}

/// Test: README Docker section includes a docker pull or docker run example
#[test]
fn issue_166_readme_docker_section_has_usage_example() {
    let readme = workspace_root().join("README.md");
    let content = fs::read_to_string(&readme)
        .unwrap_or_else(|_| panic!("cannot read {}", readme.display()));

    assert!(
        content.contains("docker pull") || content.contains("docker run"),
        "README.md Docker section must include a 'docker pull' or 'docker run' usage example"
    );
}
