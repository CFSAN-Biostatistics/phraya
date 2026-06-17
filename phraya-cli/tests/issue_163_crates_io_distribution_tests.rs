/// RED acceptance tests for issue #163: crates.io distribution channel
///
/// Verifies:
///   - README.md has an Installation section with the correct `cargo install` commands
///   - README.md documents the performance difference between portable and SIMD-native builds
///   - Cargo.toml workspace metadata includes `description` (required for cargo publish)
///   - phraya-cli crate resolves a description for cargo publish

use std::path::PathBuf;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("phraya-cli must have a parent directory")
        .to_path_buf()
}

fn readme_contents() -> String {
    let path = workspace_root().join("README.md");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("could not read README.md: {e}"))
}

fn workspace_cargo_toml_contents() -> String {
    let path = workspace_root().join("Cargo.toml");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("could not read Cargo.toml: {e}"))
}

fn cli_cargo_toml_contents() -> String {
    let path = workspace_root().join("phraya-cli").join("Cargo.toml");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("could not read phraya-cli/Cargo.toml: {e}"))
}

// ============================================================================
// README — Installation section
// ============================================================================

/// issue #163: README must have a dedicated Installation heading
#[test]
fn issue_163_readme_has_installation_section() {
    let readme = readme_contents();
    assert!(
        readme.contains("## Installation"),
        "README.md must contain a '## Installation' section; found none.\n\
         Add a section that documents how to install phraya via cargo install."
    );
}

/// issue #163: README must document the standard (portable) cargo install command using the git URL
#[test]
fn issue_163_readme_documents_standard_git_install_command() {
    let readme = readme_contents();
    let expected = "cargo install --git https://github.com/CFSAN-Biostatistics/phraya --locked phraya-cli";
    assert!(
        readme.contains(expected),
        "README.md must contain the standard install command:\n  {expected}\n\
         This lets users install without waiting for a crates.io publish."
    );
}

/// issue #163: README must document the SIMD-optimised install command with RUSTFLAGS
#[test]
fn issue_163_readme_documents_simd_optimized_install_command() {
    let readme = readme_contents();
    // The README must show RUSTFLAGS="-C target-cpu=native" applied to a cargo install invocation
    let has_rustflags_install = readme.contains("RUSTFLAGS=\"-C target-cpu=native\" cargo install")
        || readme.contains("RUSTFLAGS='-C target-cpu=native' cargo install");
    assert!(
        has_rustflags_install,
        "README.md must document the SIMD-optimised install command, e.g.:\n\
         RUSTFLAGS=\"-C target-cpu=native\" cargo install --git ... phraya-cli\n\
         Found only a build-time variant or none at all."
    );
}

/// issue #163: README must explain the performance cost of portable (non-native) builds
///
/// Acceptance criterion: "Explains performance difference: portable builds ~40-60% speed of native SIMD"
#[test]
fn issue_163_readme_explains_portable_build_performance_penalty() {
    let readme = readme_contents();
    // The spec says ~40-60%; we accept either the exact range or the individual boundaries
    let mentions_range = readme.contains("40-60%")
        || readme.contains("40–60%")
        || readme.contains("40 to 60%")
        || (readme.contains("40%") && readme.contains("60%"));
    assert!(
        mentions_range,
        "README.md must explain that portable builds run at ~40-60% the speed of native SIMD builds.\n\
         The 'Installation' section should set user expectations about the performance trade-off."
    );
}

// ============================================================================
// Cargo metadata — required for `cargo publish`
// ============================================================================

/// issue #163: [workspace.package] must supply a description so crates can inherit it
///
/// `cargo publish` rejects crates without a description. The workspace convention
/// already centralises version/license/repository in [workspace.package]; description
/// should follow the same pattern.
#[test]
fn issue_163_workspace_package_has_description_field() {
    let toml = workspace_cargo_toml_contents();

    // Find the [workspace.package] section and confirm description appears inside it
    let workspace_pkg_section = toml
        .find("[workspace.package]")
        .expect("root Cargo.toml must have a [workspace.package] section");

    // Grab everything after [workspace.package] up to the next section header
    let after_section = &toml[workspace_pkg_section..];
    let section_body = match after_section.find("\n[") {
        Some(end) => &after_section[..end],
        None => after_section,
    };

    assert!(
        section_body.contains("description"),
        "[workspace.package] must include a 'description' field.\n\
         cargo publish requires a non-empty description for every published crate.\n\
         Example: description = \"General-purpose pairwise sequence aligner for bacterial genomics\""
    );
}

/// issue #163: phraya-cli must resolve a description for cargo publish
///
/// Either the crate sets `description.workspace = true` (inheriting from [workspace.package])
/// or it supplies its own `description = "..."`.  Both satisfy the acceptance criterion.
#[test]
fn issue_163_phraya_cli_resolves_description_for_publish() {
    let toml = cli_cargo_toml_contents();
    let has_workspace_inherit = toml.contains("description.workspace")
        || toml.contains("description = { workspace = true }");
    let has_inline_description = {
        // A bare `description = "..."` line (not workspace = true)
        toml.lines().any(|line| {
            let trimmed = line.trim();
            trimmed.starts_with("description") && trimmed.contains('=') && !trimmed.contains("workspace")
        })
    };
    assert!(
        has_workspace_inherit || has_inline_description,
        "phraya-cli/Cargo.toml must have a description for cargo publish.\n\
         Either add:\n\
           description.workspace = true\n\
         (after adding description to [workspace.package]), or supply a direct:\n\
           description = \"...\""
    );
}
