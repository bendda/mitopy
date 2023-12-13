from mitopy.postprocess import do_postprocess
import pytest
import logging
import subprocess


@pytest.mark.parametrize(
    "autosomal_coverage, contamination_filter, expected_md5_vcf",
    [
        (0, True, "00fab9cfd45e2dd6d4eec1c911888ffb"),  # activate contamination filter
        (30, False, "00fab9cfd45e2dd6d4eec1c911888ffb"),  # activate numt filter
        (0, False, "00fab9cfd45e2dd6d4eec1c911888ffb"),  # deactivate both filters
    ],
)
def test_do_postprocess(
    test_files,
    tmp_path,
    get_md5,
    caplog,
    autosomal_coverage,
    contamination_filter,
    expected_md5_vcf,
):
    caplog.set_level(logging.INFO)

    postprocess = do_postprocess(
        test_files["vcf"],
        test_files["vcf_stats"],
        autosomal_coverage=autosomal_coverage,
        contamination_filter=contamination_filter,
        out_dir=tmp_path,
    )

    no_header_vcf = f"{tmp_path}/no_header.vcf"
    subprocess.run(
        f"egrep -v '^#' {postprocess['postprocessed_vcf']} > {no_header_vcf}",
        shell=True,
    )
    # Check main output
    assert get_md5(no_header_vcf) == expected_md5_vcf

    if contamination_filter:
        assert (
            "Estimating sample contamination level using Haplocheck..." in caplog.text
        )

    if autosomal_coverage != 0:
        assert "Filtering Numts..." in caplog.text
