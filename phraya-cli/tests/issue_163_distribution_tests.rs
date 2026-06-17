/// Acceptance tests for issue #163: crates.io distribution channel
///
/// Verifies that:
/// - README has an Installation section with standard and SIMD-optimized cargo install commands
/// - README explains the portable vs native SIMD performance difference
/// - phraya-cli/Cargo.toml has a [[bin]] entry naming the binary "phraya"
/// - Workspace Cargo.toml has the metadata fields required for crates.io publication

use std::path::Path;

fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().to_path_buf()
}

fn cli_manifest_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn read_readme() -> String {
    std::fs::read_to_string(workspace_root().join("README.md")).expect("README.md not found")
}

fn read_cli_cargo_toml() -> String {
    std::fs::read_to_string(cli_manifest_dir().join("Cargo.toml"))
        .expect("phraya-cli/Cargo.toml not found")
}

fn read_workspace_cargo_toml() -> String {
    std::fs::read_to_string(workspace_root().join("Cargo.toml"))
        .expect("workspace Cargo.toml not found")
}

#[test]
fn issue_163_readme_has_installation_section() {
    let readme = read_readme();
    assert!(
        readme.contains("## Installation"),
        "README must have an '## Installation' section"
    );
}

#[test]
fn issue_163_readme_documents_standard_cargo_install() {
    let readme = read_readme();
    assert!(
        readme.contains("cargo install --git https://github.com/CFSAN-Biostatistics/phraya"),
        "README must document: cargo install --git https://github.com/CFSAN-Biostatistics/phraya"
    );
    assert!(
        readme.contains("--locked phraya"),
        "README must include --locked phraya in the install command"
    );
}

#[test]
fn issue_163_readme_documents_simd_optimized_install() {
    let readme = read_readme();
    assert!(
        readme.contains("RUSTFLAGS") && readme.contains("target-cpu=native"),
        "README must document SIMD-optimized build with RUSTFLAGS=\"-C target-cpu=native\""
    );
}

#[test]
fn issue_163_readme_explains_performance_difference() {
    let readme = read_readme();
    let has_pct = readme.contains("40") && readme.contains("60");
    let has_speed = readme.contains("speed") || readme.contains("slower") || readme.contains("faster");
    assert!(
        has_pct && has_speed,
        "README must explain portable vs native SIMD performance difference (40-60% speed)"
    );
}

#[test]
fn issue_163_cli_cargo_toml_has_bin_named_phraya() {
    let toml = read_cli_cargo_toml();
    assert!(
        toml.contains("[[bin]]"),
        "phraya-cli/Cargo.toml must have a [[bin]] section"
    );
    assert!(
        toml.contains("name = \"phraya\""),
        "phraya-cli/Cargo.toml [[bin]] must set name = \"phraya\""
    );
}

#[test]
fn issue_163_cli_cargo_toml_has_description() {
    let toml = read_cli_cargo_toml();
    assert!(
        toml.contains("description"),
        "phraya-cli/Cargo.toml must have a description field for crates.io"
    );
}

#[test]
fn issue_163_workspace_has_keywords() {
    let toml = read_workspace_cargo_toml();
    assert!(
        toml.contains("keywords"),
        "workspace Cargo.toml must have keywords for crates.io"
    );
}

#[test]
fn issue_163_workspace_has_categories() {
    let toml = read_workspace_cargo_toml();
    assert!(
        toml.contains("categories"),
        "workspace Cargo.toml must have categories for crates.io"
    );
}

#[test]
fn issue_163_workspace_has_readme() {
    let toml = read_workspace_cargo_toml();
    assert!(
        toml.contains("readme"),
        "workspace Cargo.toml must have readme for crates.io"
    );
}

#[test]
fn issue_163_workspace_has_homepage() {
    let toml = read_workspace_cargo_toml();
    assert!(
        toml.contains("homepage"),
        "workspace Cargo.toml must have homepage for crates.io"
    );
}

#[test]
fn issue_163_workspace_has_documentation() {
    let toml = read_workspace_cargo_toml();
    assert!(
        toml.contains("documentation"),
        "workspace Cargo.toml must have documentation for crates.io"
    );
}

#[test]
fn issue_163_readme_installation_before_philosophy() {
    let readme = read_readme();
    let install_pos = readme.find("## Installation");
    let philosophy_pos = readme.find("## Philosophy");
    assert!(install_pos.is_some(), "README must have ## Installation section");
    assert!(philosophy_pos.is_some(), "README must have ## Philosophy section");
    assert!(
        install_pos.unwrap() < philosophy_pos.unwrap(),
        "## Installation must appear before ## Philosophy in README"
    );
}
