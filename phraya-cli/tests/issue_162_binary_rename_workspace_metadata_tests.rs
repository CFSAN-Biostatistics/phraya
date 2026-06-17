/// Issue #162: Foundation - Binary rename and workspace metadata
///
/// Acceptance criteria:
/// - phraya-cli/Cargo.toml has [[bin]] name = "phraya" entry
/// - Workspace Cargo.toml [workspace.package] includes: description, keywords,
///   categories, readme, homepage, documentation
/// - cargo build --release produces target/release/phraya binary
/// - Binary runs: target/release/phraya --version succeeds
use std::fs;
use std::path::PathBuf;

fn cli_manifest_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml")
}

fn workspace_cargo_toml() -> PathBuf {
    // CARGO_MANIFEST_DIR is phraya-cli/; parent is workspace root
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("phraya-cli must have a parent directory (workspace root)")
        .join("Cargo.toml")
}

/// issue #162: phraya-cli/Cargo.toml must declare [[bin]] with name = "phraya"
///
/// Without this entry, cargo defaults the binary name to the package name
/// ("phraya-cli"), which fails the distribution requirement.
#[test]
fn issue_162_cli_cargo_toml_has_bin_section() {
    let contents = fs::read_to_string(cli_manifest_path())
        .expect("phraya-cli/Cargo.toml must exist");
    assert!(
        contents.contains("[[bin]]"),
        "phraya-cli/Cargo.toml must have a [[bin]] section to explicitly name the binary"
    );
}

/// issue #162: the [[bin]] section in phraya-cli/Cargo.toml must set name = "phraya"
#[test]
fn issue_162_cli_bin_name_is_phraya() {
    let contents = fs::read_to_string(cli_manifest_path())
        .expect("phraya-cli/Cargo.toml must exist");
    // Verify the bin name is "phraya", not "phraya-cli" or anything else.
    // We check the literal TOML value; the default package name "phraya-cli" would
    // never produce this string in a [[bin]] section without an explicit entry.
    assert!(
        contents.contains("name = \"phraya\""),
        "phraya-cli/Cargo.toml [[bin]] must set `name = \"phraya\"`, \
         but the file contains:\n{}",
        contents
    );
}

/// issue #162: workspace [workspace.package] must include a non-empty description
#[test]
fn issue_162_workspace_package_has_description() {
    let contents = fs::read_to_string(workspace_cargo_toml())
        .expect("workspace Cargo.toml must exist");
    assert!(
        contents.contains("description"),
        "[workspace.package] must include a 'description' field for crates.io publication"
    );
    // Verify it's under [workspace.package], not just anywhere in the file
    let after_workspace_package = contents
        .split("[workspace.package]")
        .nth(1)
        .expect("[workspace.package] section must exist");
    // The description must appear before the next section header
    let before_next_section = after_workspace_package
        .split("\n[")
        .next()
        .unwrap_or(after_workspace_package);
    assert!(
        before_next_section.contains("description"),
        "'description' must be inside [workspace.package], not elsewhere in the file"
    );
}

/// issue #162: workspace [workspace.package] must include a keywords array
#[test]
fn issue_162_workspace_package_has_keywords() {
    let contents = fs::read_to_string(workspace_cargo_toml())
        .expect("workspace Cargo.toml must exist");
    let after_workspace_package = contents
        .split("[workspace.package]")
        .nth(1)
        .expect("[workspace.package] section must exist");
    let before_next_section = after_workspace_package
        .split("\n[")
        .next()
        .unwrap_or(after_workspace_package);
    assert!(
        before_next_section.contains("keywords"),
        "[workspace.package] must include a 'keywords' field for crates.io discoverability"
    );
}

/// issue #162: workspace [workspace.package] must include a categories array
#[test]
fn issue_162_workspace_package_has_categories() {
    let contents = fs::read_to_string(workspace_cargo_toml())
        .expect("workspace Cargo.toml must exist");
    let after_workspace_package = contents
        .split("[workspace.package]")
        .nth(1)
        .expect("[workspace.package] section must exist");
    let before_next_section = after_workspace_package
        .split("\n[")
        .next()
        .unwrap_or(after_workspace_package);
    assert!(
        before_next_section.contains("categories"),
        "[workspace.package] must include a 'categories' field for crates.io classification"
    );
}

/// issue #162: workspace [workspace.package] must include a readme field
#[test]
fn issue_162_workspace_package_has_readme() {
    let contents = fs::read_to_string(workspace_cargo_toml())
        .expect("workspace Cargo.toml must exist");
    let after_workspace_package = contents
        .split("[workspace.package]")
        .nth(1)
        .expect("[workspace.package] section must exist");
    let before_next_section = after_workspace_package
        .split("\n[")
        .next()
        .unwrap_or(after_workspace_package);
    assert!(
        before_next_section.contains("readme"),
        "[workspace.package] must include a 'readme' field pointing to README.md"
    );
}

/// issue #162: workspace [workspace.package] must include a homepage field
#[test]
fn issue_162_workspace_package_has_homepage() {
    let contents = fs::read_to_string(workspace_cargo_toml())
        .expect("workspace Cargo.toml must exist");
    let after_workspace_package = contents
        .split("[workspace.package]")
        .nth(1)
        .expect("[workspace.package] section must exist");
    let before_next_section = after_workspace_package
        .split("\n[")
        .next()
        .unwrap_or(after_workspace_package);
    assert!(
        before_next_section.contains("homepage"),
        "[workspace.package] must include a 'homepage' field for crates.io publication"
    );
}

/// issue #162: workspace [workspace.package] must include a documentation field
#[test]
fn issue_162_workspace_package_has_documentation() {
    let contents = fs::read_to_string(workspace_cargo_toml())
        .expect("workspace Cargo.toml must exist");
    let after_workspace_package = contents
        .split("[workspace.package]")
        .nth(1)
        .expect("[workspace.package] section must exist");
    let before_next_section = after_workspace_package
        .split("\n[")
        .next()
        .unwrap_or(after_workspace_package);
    assert!(
        before_next_section.contains("documentation"),
        "[workspace.package] must include a 'documentation' field for crates.io publication"
    );
}

/// issue #162: the binary must be invocable as "phraya" (not "phraya-cli")
///
/// This test runs `cargo run --bin phraya -- --version` and verifies it
/// succeeds. It fails if [[bin]] name is absent or set to anything other
/// than "phraya".
#[test]
fn issue_162_binary_invocable_as_phraya_with_version_flag() {
    let manifest_path = cli_manifest_path();
    let output = std::process::Command::new("cargo")
        .args([
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--bin",
            "phraya",
            "--",
            "--version",
        ])
        .output()
        .expect("cargo run must be executable");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        output.status.success(),
        "`cargo run --bin phraya -- --version` must exit 0.\nstdout: {}\nstderr: {}",
        stdout,
        stderr
    );

    // --version output must contain the binary name "phraya" and a version string
    assert!(
        stdout.contains("phraya"),
        "`phraya --version` stdout must contain the name 'phraya', got: {}",
        stdout
    );
}
