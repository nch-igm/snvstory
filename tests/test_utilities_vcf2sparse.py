import json
import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")
from scipy import sparse

from igm_churchill_ancestry.utilities.vcf2sparse import vcf_to_json, json_to_sparse_matrix


def test_vcf_to_json_and_sparse(tmp_path):
    attr_dir = tmp_path / "attr"
    attr_dir.mkdir()
    variant_file = attr_dir / "variants.json"
    with variant_file.open("w") as fh:
        json.dump({"1_1000_A_G": 0}, fh)
    parsed_vcf = ["1\t1000\t.\tA\tG\t.\t.\t.\tGT\t0/1"]
    g_container = vcf_to_json(parsed_vcf, str(attr_dir), None)
    assert g_container["1_1000_A_G"] == 1
    s_matrix = json_to_sparse_matrix(g_container, ["1_1000_A_G"])
    assert sparse.issparse(s_matrix)
    assert s_matrix.shape == (1, 1)
    assert s_matrix.toarray()[0, 0] == 1

