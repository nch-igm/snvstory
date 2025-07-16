import json
from igm_churchill_ancestry.utilities.vcf2sparse import vcf_to_json


def test_vcf_to_json_basic(tmp_path):
    attr_dir = tmp_path / "attrs"
    attr_dir.mkdir()
    variants = {
        "1_123_A_T": 0,
        "1_456_G_A": 0,
    }
    with open(attr_dir / "vars.json", "w") as f:
        json.dump(variants, f)

    parsed_vcf = [
        "1\t123\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
        "1\t456\t.\tG\tA\t.\tPASS\t.\tGT\t1/1",
    ]

    result = vcf_to_json(parsed_vcf, str(attr_dir), None)
    assert result["1_123_A_T"] == 1
    assert result["1_456_G_A"] == 2
