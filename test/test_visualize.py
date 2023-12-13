from mitopy.visualize import do_visualize
import pytest
import logging


@pytest.mark.parametrize(
    "include_coverage, expected_md5",
    [
        (True, "ce55e661f267ee916b983344ab1ecd15"),
        (False, "256f3ab0b5198924576359c525e9eaa1"),
    ],
)
def test_do_visualize(
    test_files, tmp_path, get_md5, caplog, include_coverage, expected_md5
):
    caplog.set_level(logging.INFO)
    if include_coverage:
        vis = do_visualize(
            test_files["vcf"],
            test_files["coverage_csv"],
            out_dir=tmp_path,
            save_as_png=True,
        )
    else:
        vis = do_visualize(test_files["vcf"], out_dir=tmp_path, save_as_png=True)

    # Check main outputs
    assert get_md5(vis["vis_png"]) == expected_md5

    if include_coverage:
        assert "Plotting per base coverage..." in caplog.text
