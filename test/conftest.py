import pytest
import os
import hashlib


@pytest.fixture
def root_dir():
    return os.path.dirname(__file__)


@pytest.fixture
def test_dir():
    return os.path.join(os.path.dirname(__file__), "test_files")


@pytest.fixture
def mock_subprocess_run(mocker):
    return mocker.patch("subprocess.run")


@pytest.fixture
def get_md5():
    # Calculate the MD5 hash of the file
    def get_md5_file(file_path):
        md5_hash = hashlib.md5(open(file_path, "rb").read()).hexdigest()
        return md5_hash

    return get_md5_file


@pytest.fixture
def test_files(test_dir):
    return {
        "bam": f"{test_dir}/bams/NA12878.bam",
        "bai": f"{test_dir}/bams/NA12878.bai",
        "unmapped_bam": f"{test_dir}/bams/NA12878_unmapped.bam",
        "dedup_bam": f"{test_dir}/bams/NA12878_dedup.bam",
        "shifted_dedup_bam": f"{test_dir}/bams/NA12878_shifted_dedup.bam",
        "vcf": f"{test_dir}/vcfs/NA12878.vcf",
        "vcf_stats": f"{test_dir}/vcfs/NA12878.vcf.stats",
        "shifted_vcf": f"{test_dir}/vcfs/NA12878_shifted.vcf",
        "shifted_vcf_stats": f"{test_dir}/vcfs/NA12878_shifted.vcf.stats",
        "coverage_csv": f"{test_dir}/vis/coverage.csv",
    }
