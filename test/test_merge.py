from mitopy.merge import do_merge


def test_do_merge(test_files, tmp_path, get_md5):
    merge = do_merge(test_files["vcf"], test_files["shifted_vcf"], out_dir=tmp_path)

    # Check main output
    assert get_md5(merge["merged_vcf"]) == "79f41b8d71dcf941e70667eb44f9b3ad"
    assert get_md5(merge["merged_vcf_stats"]) == "6eee922d3d362325d9bdaab81af2d1cc"
