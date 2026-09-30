class Phraya < Formula
  desc "General-purpose pairwise sequence aligner for bacterial genomics"
  homepage "https://github.com/CFSAN-Biostatistics/phraya"
  version "0.1.0"

  on_linux do
    if Hardware::CPU.intel?
      url "https://github.com/CFSAN-Biostatistics/phraya/releases/download/v#{version}/phraya-v#{version}-x86_64-linux-gnu-portable.tar.gz"
      sha256 "TBD"  # Update after release
    elsif Hardware::CPU.arm?
      url "https://github.com/CFSAN-Biostatistics/phraya/releases/download/v#{version}/phraya-v#{version}-aarch64-linux-gnu.tar.gz"
      sha256 "TBD"  # Update after release
    end
  end

  on_macos do
    if Hardware::CPU.intel?
      url "https://github.com/CFSAN-Biostatistics/phraya/releases/download/v#{version}/phraya-v#{version}-x86_64-darwin.tar.gz"
      sha256 "TBD"  # Update after release
    elsif Hardware::CPU.arm?
      url "https://github.com/CFSAN-Biostatistics/phraya/releases/download/v#{version}/phraya-v#{version}-aarch64-darwin.tar.gz"
      sha256 "TBD"  # Update after release
    end
  end

  def install
    bin.install "phraya"
  end

  def test
    system "#{bin}/phraya", "--version"
  end
end
