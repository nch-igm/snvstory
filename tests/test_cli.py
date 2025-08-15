import os
import sys
import types

boto3_stub = types.ModuleType("boto3")
boto3_stub.client = lambda *a, **k: types.SimpleNamespace(download_file=lambda *a, **k: None,
                                                         upload_file=lambda *a, **k: None)
sys.modules.setdefault("boto3", boto3_stub)

awscli_module = types.ModuleType("awscli")
clidriver_module = types.ModuleType("awscli.clidriver")
def _create_clidriver():
    class _Driver:
        def main(self, cmd):
            return 0
    return _Driver()
clidriver_module.create_clidriver = _create_clidriver
sys.modules.setdefault("awscli", awscli_module)
sys.modules.setdefault("awscli.clidriver", clidriver_module)

numpy_stub = types.ModuleType("numpy")
class _Random:
    def seed(self, x):
        return None
numpy_stub.random = _Random()
sys.modules.setdefault("numpy", numpy_stub)

vcf2sparse_stub = types.ModuleType("igm_churchill_ancestry.utilities.vcf2sparse")
vcf2sparse_stub.vcf_to_json = lambda *a, **k: None
vcf2sparse_stub.load_snp_order = lambda *a, **k: None
vcf2sparse_stub.json_to_sparse_matrix = lambda *a, **k: None
sys.modules.setdefault("igm_churchill_ancestry.utilities.vcf2sparse", vcf2sparse_stub)

plot_stub = types.ModuleType("igm_churchill_ancestry.utilities.plot_ancestry")
plot_stub.plot_parser = lambda *a, **k: None
sys.modules.setdefault("igm_churchill_ancestry.utilities.plot_ancestry", plot_stub)

umap_stub = types.ModuleType("igm_churchill_ancestry.utilities.plot_umap")
umap_stub.plot_umap_parser = lambda *a, **k: None
sys.modules.setdefault("igm_churchill_ancestry.utilities.plot_umap", umap_stub)

ancestry_stub = types.ModuleType("igm_churchill_ancestry.pipelines.ancestry_prediction")
ancestry_stub.run_ancestry_pipeline = lambda *a, **k: None
sys.modules.setdefault("igm_churchill_ancestry.pipelines.ancestry_prediction", ancestry_stub)

from igm_churchill_ancestry.cli import setup_workspace


def test_setup_workspace(tmp_path, monkeypatch):
    monkeypatch.setenv("TMP_DIR", str(tmp_path))
    monkeypatch.setenv("AWS_BATCH_JOB_ID", "testjob")
    wrk_dir = setup_workspace()
    assert wrk_dir == str(tmp_path / "testjob")
    assert os.path.isdir(wrk_dir)

