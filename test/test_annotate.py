from mitopy.annotate import do_annotate
import pytest


def test_do_annotate(test_files, tmp_path, get_md5):
    ann = do_annotate(test_files["vcf"], out_dir=tmp_path)

    # Check main outputs
    assert get_md5(ann["annotated_vcf"]) == "5e65270bdc4281598675c20740ba5e4e"
    assert get_md5(ann["annot_csv"]) == "3090fe9b532780ccfd29da3e15aa0099"
