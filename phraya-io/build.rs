use std::fs;
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let manifest_dir =
        PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR not set"));
    let fixtures_dir = manifest_dir.join("tests/fixtures");

    // Generate compressed versions
    if fixtures_dir.join("sample.fasta").exists() {
        let fasta_content = fs::read(fixtures_dir.join("sample.fasta"))
            .expect("Failed to read sample.fasta");

        // Generate gzipped versions
        {
            use flate2::Compression;
            use flate2::write::GzEncoder;

            let gz_path = fixtures_dir.join("sample.fasta.gz");
            let file = fs::File::create(&gz_path).expect("Failed to create sample.fasta.gz");
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder
                .write_all(&fasta_content)
                .expect("Failed to write gzipped FASTA");
        }

        // Generate bzip2 versions
        {
            use bzip2::Compression;
            use bzip2::write::BzEncoder;

            let bz2_path = fixtures_dir.join("sample.fasta.bz2");
            let file = fs::File::create(&bz2_path).expect("Failed to create sample.fasta.bz2");
            let mut encoder = BzEncoder::new(file, Compression::default());
            encoder
                .write_all(&fasta_content)
                .expect("Failed to write bzip2 FASTA");
        }
    }

    if fixtures_dir.join("sample.fastq").exists() {
        let fastq_content = fs::read(fixtures_dir.join("sample.fastq"))
            .expect("Failed to read sample.fastq");

        // Generate gzipped versions
        {
            use flate2::Compression;
            use flate2::write::GzEncoder;

            let gz_path = fixtures_dir.join("sample.fastq.gz");
            let file = fs::File::create(&gz_path).expect("Failed to create sample.fastq.gz");
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder
                .write_all(&fastq_content)
                .expect("Failed to write gzipped FASTQ");
        }

        // Generate bzip2 versions
        {
            use bzip2::Compression;
            use bzip2::write::BzEncoder;

            let bz2_path = fixtures_dir.join("sample.fastq.bz2");
            let file = fs::File::create(&bz2_path).expect("Failed to create sample.fastq.bz2");
            let mut encoder = BzEncoder::new(file, Compression::default());
            encoder
                .write_all(&fastq_content)
                .expect("Failed to write bzip2 FASTQ");
        }
    }
}
