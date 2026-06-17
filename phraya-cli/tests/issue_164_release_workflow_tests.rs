/// Acceptance tests for Issue #164: GitHub Releases: Linux portable build (MVP)
///
/// Validates that .github/workflows/release.yml exists and contains the required
/// configuration: tag trigger, x86_64-linux-gnu target, SSE4.2 baseline RUSTFLAGS,
/// strip step, portable tarball with correct naming, and GitHub Release upload.
/// Also validates that README contains prebuilt binary download instructions.
use std::path::PathBuf;

fn workspace_root() -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    std::path::Path::new(manifest_dir)
        .parent()
        .expect("phraya-cli has a parent directory")
        .to_path_buf()
}

fn release_workflow_path() -> PathBuf {
    workspace_root().join(".github").join("workflows").join("release.yml")
}

fn read_release_workflow() -> String {
    let path = release_workflow_path();
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!(".github/workflows/release.yml must exist and be readable: {}", e))
}

// ============================================================================
// File existence
// ============================================================================

#[test]
fn issue_164_release_workflow_file_exists() {
    let path = release_workflow_path();
    assert!(
        path.exists(),
        ".github/workflows/release.yml must exist; got path: {}",
        path.display()
    );
}

// ============================================================================
// Trigger: push on version tags
// ============================================================================

#[test]
fn issue_164_workflow_triggers_on_push_tags() {
    let content = read_release_workflow();
    assert!(
        content.contains("push:"),
        "workflow must have a 'push:' trigger; got:\n{}", content
    );
    assert!(
        content.contains("tags:"),
        "workflow push trigger must include 'tags:'; got:\n{}", content
    );
}

#[test]
fn issue_164_workflow_tag_filter_matches_version_tags() {
    let content = read_release_workflow();
    // The tag filter must start with 'v' (e.g. 'v*' or 'v**' or "'v*'")
    assert!(
        content.contains("'v*'") || content.contains("\"v*\"") || content.contains("- v*") || content.contains("- 'v*'"),
        "tags filter must match version tags starting with 'v' (e.g. 'v*'); got:\n{}", content
    );
}

// ============================================================================
// Build: x86_64-unknown-linux-gnu target
// ============================================================================

#[test]
fn issue_164_workflow_builds_x86_64_linux_gnu_target() {
    let content = read_release_workflow();
    assert!(
        content.contains("x86_64-unknown-linux-gnu"),
        "workflow must reference target x86_64-unknown-linux-gnu; got:\n{}", content
    );
}

// ============================================================================
// Build: SSE4.2 portable SIMD baseline via RUSTFLAGS
// ============================================================================

#[test]
fn issue_164_workflow_uses_sse42_rustflag() {
    let content = read_release_workflow();
    assert!(
        content.contains("target-feature=+sse4.2"),
        "workflow must set RUSTFLAGS with target-feature=+sse4.2 for portable SSE4.2 baseline; got:\n{}", content
    );
}

// ============================================================================
// Strip step
// ============================================================================

#[test]
fn issue_164_workflow_strips_binary() {
    let content = read_release_workflow();
    assert!(
        content.contains("strip"),
        "workflow must include a 'strip' step to reduce binary size; got:\n{}", content
    );
}

// ============================================================================
// Tarball naming: phraya-{version}-x86_64-linux-gnu-portable.tar.gz
// ============================================================================

#[test]
fn issue_164_workflow_creates_portable_tarball() {
    let content = read_release_workflow();
    assert!(
        content.contains(".tar.gz"),
        "workflow must create a .tar.gz tarball; got:\n{}", content
    );
}

#[test]
fn issue_164_workflow_tarball_name_contains_phraya_prefix() {
    let content = read_release_workflow();
    assert!(
        content.contains("phraya-"),
        "tarball name must start with 'phraya-'; got:\n{}", content
    );
}

#[test]
fn issue_164_workflow_tarball_name_contains_portable_suffix() {
    let content = read_release_workflow();
    assert!(
        content.contains("x86_64-linux-gnu-portable"),
        "tarball name must contain 'x86_64-linux-gnu-portable' per the naming spec; got:\n{}", content
    );
}

// ============================================================================
// GitHub Release upload
// ============================================================================

#[test]
fn issue_164_workflow_uploads_asset_to_github_release() {
    let content = read_release_workflow();
    // Accept any standard upload mechanism: softprops/action-gh-release,
    // actions/upload-release-asset, or gh CLI
    let has_upload = content.contains("softprops/action-gh-release")
        || content.contains("upload-release-asset")
        || content.contains("gh release upload")
        || content.contains("gh release create");
    assert!(
        has_upload,
        "workflow must upload the tarball to a GitHub Release using \
         softprops/action-gh-release, upload-release-asset, or gh CLI; got:\n{}", content
    );
}

// ============================================================================
// README: prebuilt binary download instructions
// ============================================================================

#[test]
fn issue_164_readme_mentions_prebuilt_binaries() {
    let readme_path = workspace_root().join("README.md");
    let content = std::fs::read_to_string(&readme_path)
        .expect("README.md must exist and be readable");
    let content_lower = content.to_lowercase();
    assert!(
        content_lower.contains("prebuilt") || content_lower.contains("pre-built"),
        "README must mention prebuilt binaries; got README without 'prebuilt' or 'pre-built'"
    );
}

#[test]
fn issue_164_readme_has_download_instructions() {
    let readme_path = workspace_root().join("README.md");
    let content = std::fs::read_to_string(&readme_path)
        .expect("README.md must exist and be readable");
    let content_lower = content.to_lowercase();
    // Must have some kind of download/install section
    assert!(
        content_lower.contains("download") || content_lower.contains("install"),
        "README must include download or install instructions for prebuilt binaries"
    );
}

#[test]
fn issue_164_readme_references_release_tarball_pattern() {
    let readme_path = workspace_root().join("README.md");
    let content = std::fs::read_to_string(&readme_path)
        .expect("README.md must exist and be readable");
    assert!(
        content.contains("x86_64-linux-gnu-portable") || content.contains("phraya-") && content.contains("linux"),
        "README must reference the portable Linux tarball (x86_64-linux-gnu-portable or similar); \
         found no matching pattern in README"
    );
}
