/// Acceptance tests for Issue #169: Automated packaging and release orchestration
///
/// These tests verify the static presence and structure of release infrastructure:
/// - .github/workflows/release.yml (all release channels, 5 platform binaries, Docker, crates.io)
/// - README.md "For Maintainers" section
/// - CHANGELOG.md with release notes template
///
/// All tests read files from the repository root (one directory above phraya-cli/).
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    // CARGO_MANIFEST_DIR is phraya-cli/; repo root is one level up
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest.parent().unwrap().to_path_buf()
}

fn release_yml_path() -> PathBuf {
    repo_root().join(".github/workflows/release.yml")
}

fn read_release_yml() -> String {
    std::fs::read_to_string(release_yml_path())
        .expect("release.yml must exist at .github/workflows/release.yml")
}

// ============================================================================
// SECTION A: release.yml existence and trigger
// ============================================================================

/// release.yml must exist in .github/workflows/
#[test]
fn issue_169_release_yml_exists() {
    assert!(
        release_yml_path().exists(),
        ".github/workflows/release.yml must exist to automate releases"
    );
}

/// release.yml must trigger on version tag pushes (v* pattern)
#[test]
fn issue_169_release_yml_triggers_on_version_tags() {
    let content = read_release_yml();
    // Accept both 'v*' and 'v[0-9]*' tag patterns
    let has_vtag_pattern = content.contains("'v*'")
        || content.contains("\"v*\"")
        || content.contains("v[0-9]")
        || content.contains("'v[0-9]");
    assert!(
        has_vtag_pattern,
        "release.yml must trigger on version tag pushes (e.g. 'v*'); \
         found no recognizable tag pattern in:\n{}",
        &content[..content.len().min(500)]
    );
}

// ============================================================================
// SECTION B: GitHub Releases — 5 platform binaries
// ============================================================================

/// release.yml must build a linux-x86_64 binary
#[test]
fn issue_169_release_yml_targets_linux_x86_64() {
    let content = read_release_yml();
    let has_target = content.contains("x86_64-unknown-linux")
        || content.contains("x86_64-linux")
        || (content.contains("linux") && content.contains("x86_64"));
    assert!(
        has_target,
        "release.yml must include a linux x86_64 build target for one of the 5 release binaries"
    );
}

/// release.yml must build a linux-aarch64 binary
#[test]
fn issue_169_release_yml_targets_linux_aarch64() {
    let content = read_release_yml();
    let has_target = content.contains("aarch64-unknown-linux")
        || content.contains("aarch64-linux")
        || (content.contains("linux") && content.contains("aarch64"));
    assert!(
        has_target,
        "release.yml must include a linux aarch64 build target for one of the 5 release binaries"
    );
}

/// release.yml must build a macOS binary (x86_64 or aarch64)
#[test]
fn issue_169_release_yml_targets_macos() {
    let content = read_release_yml();
    let has_target = content.contains("apple-darwin")
        || content.contains("macos-latest")
        || content.contains("macos-14")
        || content.contains("macos-13")
        || content.contains("macos-12");
    assert!(
        has_target,
        "release.yml must include at least one macOS build target for the 5 release binaries"
    );
}

/// release.yml must build a Windows x86_64 binary
#[test]
fn issue_169_release_yml_targets_windows_x86_64() {
    let content = read_release_yml();
    let has_target = content.contains("x86_64-pc-windows")
        || content.contains("windows-latest")
        || content.contains("windows-2022")
        || content.contains("windows-2019");
    assert!(
        has_target,
        "release.yml must include a Windows x86_64 build target for one of the 5 release binaries"
    );
}

/// release.yml must upload binaries to GitHub Releases
#[test]
fn issue_169_release_yml_uploads_github_release_binaries() {
    let content = read_release_yml();
    let has_upload = content.contains("softprops/action-gh-release")
        || content.contains("actions/upload-release-asset")
        || content.contains("gh release create")
        || content.contains("gh release upload");
    assert!(
        has_upload,
        "release.yml must upload binaries to GitHub Releases using \
         softprops/action-gh-release, actions/upload-release-asset, or gh CLI"
    );
}

// ============================================================================
// SECTION C: Docker multi-arch — ghcr.io with :latest and versioned tags
// ============================================================================

/// release.yml must have a Docker build step
#[test]
fn issue_169_release_yml_has_docker_build() {
    let content = read_release_yml();
    let has_docker = content.contains("docker/build-push-action")
        || content.contains("docker buildx build")
        || (content.contains("docker") && content.contains("buildx"));
    assert!(
        has_docker,
        "release.yml must build Docker images using docker buildx or docker/build-push-action"
    );
}

/// release.yml must push Docker images to ghcr.io
#[test]
fn issue_169_release_yml_pushes_to_ghcr() {
    let content = read_release_yml();
    assert!(
        content.contains("ghcr.io"),
        "release.yml must push Docker images to ghcr.io (GitHub Container Registry)"
    );
}

/// release.yml must tag the Docker image with :latest
#[test]
fn issue_169_release_yml_docker_has_latest_tag() {
    let content = read_release_yml();
    assert!(
        content.contains(":latest") || content.contains("latest"),
        "release.yml must tag Docker image with :latest"
    );
}

/// release.yml must tag the Docker image with a versioned tag (not just :latest)
#[test]
fn issue_169_release_yml_docker_has_versioned_tag() {
    let content = read_release_yml();
    // Common patterns: ${{ github.ref_name }}, tags from docker/metadata-action, type=semver
    let has_versioned = content.contains("github.ref_name")
        || content.contains("type=semver")
        || content.contains("docker/metadata-action")
        || content.contains("ref_name");
    assert!(
        has_versioned,
        "release.yml must tag Docker image with a versioned tag (e.g. github.ref_name or semver from metadata-action)"
    );
}

/// release.yml must target multiple architectures for the Docker image
#[test]
fn issue_169_release_yml_docker_multiarch() {
    let content = read_release_yml();
    // Multi-arch Docker requires platforms specification
    let has_multiarch = content.contains("linux/amd64")
        || content.contains("linux/arm64")
        || content.contains("platforms:");
    assert!(
        has_multiarch,
        "release.yml must build multi-arch Docker image; expected 'platforms:' with linux/amd64 and/or linux/arm64"
    );
}

// ============================================================================
// SECTION D: crates.io publish
// ============================================================================

/// release.yml must publish to crates.io via cargo publish
#[test]
fn issue_169_release_yml_triggers_cargo_publish() {
    let content = read_release_yml();
    assert!(
        content.contains("cargo publish"),
        "release.yml must run 'cargo publish' to push to crates.io"
    );
}

/// cargo publish step must use CARGO_REGISTRY_TOKEN or CARGO_TOKEN secret
#[test]
fn issue_169_release_yml_cargo_publish_uses_registry_token() {
    let content = read_release_yml();
    let has_token = content.contains("CARGO_REGISTRY_TOKEN")
        || content.contains("CARGO_TOKEN")
        || content.contains("secrets.CARGO");
    assert!(
        has_token,
        "release.yml cargo publish step must authenticate via CARGO_REGISTRY_TOKEN or CARGO_TOKEN secret"
    );
}

// ============================================================================
// SECTION E: README "For Maintainers" section
// ============================================================================

/// README.md must contain a "For Maintainers" section
#[test]
fn issue_169_readme_has_for_maintainers_section() {
    let readme_path = repo_root().join("README.md");
    let content = std::fs::read_to_string(&readme_path)
        .expect("README.md must exist");
    assert!(
        content.contains("For Maintainers"),
        "README.md must contain a 'For Maintainers' section describing the release process"
    );
}

/// README.md "For Maintainers" section must mention the release tag process
#[test]
fn issue_169_readme_maintainers_section_describes_release_tag() {
    let readme_path = repo_root().join("README.md");
    let content = std::fs::read_to_string(&readme_path).expect("README.md must exist");

    assert!(
        content.contains("For Maintainers"),
        "README.md must have 'For Maintainers' section"
    );

    // Find the section and verify it contains tagging instructions
    let section_start = content
        .find("For Maintainers")
        .expect("'For Maintainers' heading must be present");
    let section_content = &content[section_start..];

    let mentions_tag = section_content.contains("git tag")
        || section_content.contains("v0.")
        || section_content.contains("v1.")
        || section_content.contains("tag v");
    assert!(
        mentions_tag,
        "README.md 'For Maintainers' section must describe how to create a release tag \
         (e.g. 'git tag v0.x.x' or similar)"
    );
}

// ============================================================================
// SECTION F: CHANGELOG.md with release notes template
// ============================================================================

/// CHANGELOG.md must exist at the repository root
#[test]
fn issue_169_changelog_exists() {
    let changelog_path = repo_root().join("CHANGELOG.md");
    assert!(
        changelog_path.exists(),
        "CHANGELOG.md must exist at the repository root"
    );
}

/// CHANGELOG.md must contain an [Unreleased] section or version header
#[test]
fn issue_169_changelog_has_unreleased_section() {
    let changelog_path = repo_root().join("CHANGELOG.md");
    let content = std::fs::read_to_string(&changelog_path)
        .expect("CHANGELOG.md must exist");
    // Keep a Changelog format: ## [Unreleased] is the standard
    assert!(
        content.contains("[Unreleased]") || content.contains("## Unreleased"),
        "CHANGELOG.md must contain an '## [Unreleased]' section following Keep a Changelog convention"
    );
}

/// CHANGELOG.md must have version section structure (Added/Changed/Fixed categories)
#[test]
fn issue_169_changelog_has_release_notes_template() {
    let changelog_path = repo_root().join("CHANGELOG.md");
    let content = std::fs::read_to_string(&changelog_path)
        .expect("CHANGELOG.md must exist");
    // Keep a Changelog defines standard subsections
    let has_categories = content.contains("### Added")
        || content.contains("### Changed")
        || content.contains("### Fixed")
        || content.contains("### Removed");
    assert!(
        has_categories,
        "CHANGELOG.md must contain release notes template with standard categories \
         (### Added, ### Changed, ### Fixed, ### Removed)"
    );
}
