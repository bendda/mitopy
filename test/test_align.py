from mitopy.align import do_align
import pytest
import pysam


# shifted and unshifted mode
@pytest.mark.parametrize(
    "shifted_mode, expected_md5_bam",
    [
        (False, "454b6e5c2f58c6254429664f32f5dbb1"),
        (True, "4cf9524a77c439987a764ad4a3120ad2"),
    ],
)
def test_do_align(test_files, tmp_path, get_md5, shifted_mode, expected_md5_bam):
    aln = do_align(test_files["unmapped_bam"], out_dir=tmp_path, shifted=shifted_mode)

    # convert to sam (to ignore header)
    sam = f"{tmp_path}/test.sam"
    pysam.view(aln["dedup_sorted_bam"], "-o", sam, catch_stdout=False)
    # Check main output
    assert get_md5(sam) == expected_md5_bam
