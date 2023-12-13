from mitopy.call import do_call
import pytest
import subprocess


# shifted and unshifted mode
@pytest.mark.parametrize(
    "shifted_mode, expected_md5_vcf,  expected_md5_stats",
    [
        (False, "8e851dd866f99c1ecaef75d25885ff30", "c76d2113af8dda3f76d616598a236e59"),
        (True, "c1b4d20814f7cc4eeeafd0d0f71afea7", "e7d6c05baa9dede12aa3041af9f73258"),
    ],
)
def test_do_call(
    test_files, tmp_path, get_md5, shifted_mode, expected_md5_vcf, expected_md5_stats
):
    input_bam = (
        test_files["shifted_dedup_bam"] if shifted_mode else test_files["dedup_bam"]
    )
    call = do_call(input_bam, out_dir=tmp_path, shifted=shifted_mode)

    # ignore header
    no_header_vcf = f"{tmp_path}/no_header.vcf"
    subprocess.run(f"egrep -v '^#' {call['raw_vcf']} > {no_header_vcf}", shell=True)

    # Check main output
    assert get_md5(no_header_vcf) == expected_md5_vcf
    assert get_md5(call["raw_vcf_stats"]) == expected_md5_stats
