from mitopy.preprocess import do_preprocess


def test_do_preprocess(test_files, tmp_path, get_md5):
    prep = do_preprocess(test_files["bam"], test_files["bai"], out_dir=tmp_path)

    # Check main output
    assert get_md5(prep["unmapped_bam"]) == "aab17066bb85e07ff89733f055847f79"
