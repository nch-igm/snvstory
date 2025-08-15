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

from igm_churchill_ancestry.utilities.flex import flex_input, flex_output


def test_flex_input_local_file(tmp_path):
    src = tmp_path / "file.txt"
    src.write_text("data")
    result = flex_input(str(src))
    assert result == str(src)


def test_flex_input_force_copy(tmp_path):
    src = tmp_path / "file.txt"
    src.write_text("data")
    out_dir = tmp_path / "dest"
    out_dir.mkdir()
    result = flex_input(str(src), str(out_dir), force_copy=True)
    expected = out_dir / "file.txt"
    assert result == str(expected)
    assert expected.read_text() == "data"


def test_flex_output_local(tmp_path):
    src = tmp_path / "file.txt"
    src.write_text("data")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    result = flex_output(str(src), str(out_dir))
    expected = out_dir / "file.txt"
    assert result == str(expected)
    assert expected.exists()

