from mitopy.haplogroup import do_identify_haplogroup
import pytest
import logging


@pytest.mark.parametrize(
    "mt_ref, expected_md5",
    [
        ("rcrs", "a3612e09d3bc0ec1abb936bdbfc1a38a"),
        ("rsrs", "6b69f48e25b87855820146685190b516"),
    ],
)
def test_do_identify_haplogroup(
    test_files, tmp_path, get_md5, caplog, mt_ref, expected_md5
):
    caplog.set_level(logging.INFO)

    haplo_report = do_identify_haplogroup(test_files["vcf"], mt_ref, out_dir=tmp_path)

    # Check main outputs
    assert get_md5(haplo_report["haplogroups"]) == expected_md5

    if mt_ref == "rsrs":
        assert "Downloading phylotree for RSRS reference..." in caplog.text
