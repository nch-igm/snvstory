from igm_churchill_ancestry.utilities.flex import flex_input


def test_flex_input_local_return(tmp_path):
    f = tmp_path / "file.txt"
    f.write_text("data")
    assert flex_input(str(f)) == str(f)


def test_flex_input_force_copy(tmp_path):
    f = tmp_path / "file.txt"
    f.write_text("data")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    result = flex_input(str(f), out_dir=str(out_dir), force_copy=True)
    assert result == str(out_dir / "file.txt")
    assert (out_dir / "file.txt").exists()
