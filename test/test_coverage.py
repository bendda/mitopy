from mitopy.coverage import do_coverage
import pytest
import warnings


def test_do_coverage(test_files, tmp_path, get_md5):
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    cov = do_coverage(
        test_files["dedup_bam"], test_files["shifted_dedup_bam"], out_dir=tmp_path
    )

    # Check main outputs
    assert get_md5(cov["coverage_csv"]) == "3c4073b9073e14ae6b10f745a5c96c84"
