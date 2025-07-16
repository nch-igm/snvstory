import os
from igm_churchill_ancestry.utilities.utilities import get_extension, filter_extension


def test_get_extension(tmp_path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text("")
    assert get_extension(str(vcf)) == ".vcf"


def test_filter_extension(tmp_path):
    paths = [str(tmp_path / "a.vcf"), str(tmp_path / "b.txt"), str(tmp_path / "c.g.vcf")]
    result = filter_extension(paths)
    assert str(tmp_path / "a.vcf") in result
    assert str(tmp_path / "c.g.vcf") in result
    assert str(tmp_path / "b.txt") not in result
