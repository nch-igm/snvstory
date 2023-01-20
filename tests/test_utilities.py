import os
import uuid

import pytest

from igm_churchill_cnvloh.utilities.utilities import check_resource_folder, compute_sex, make_index
from igm_churchill_cnvloh.pipelines.variables import variables

def create_tmp_file(name):
    file_name = f"{name}"
    with open(file_name, 'w') as fn:
        fn.write("temp file data")
    return file_name

def test_check_resource_folder():
    """Test that normal pipeline returns True if a local file exists."""
    fn1 = create_tmp_file('targets.preprocessed.interval_list')
    fn2 = create_tmp_file('biallelic.vcf')
    fn3 = create_tmp_file('female.normals.pon.hdf5')
    fn4 = create_tmp_file('male.normals.pon.hdf5')
    fn5 = create_tmp_file('gene_table.txt')
    output = check_resource_folder(path='.')
    assert output == True
    os.remove(fn1)
    os.remove(fn2)
    os.remove(fn3)
    os.remove(fn4)
    os.remove(fn5)

def test_variables_return_true():
    PLT_PATH = variables.PLT_PATH
    assert os.path.isfile(PLT_PATH) is True

def test_make_index_return_true():
    path_input_0 = 'test.bam'
    path_input_1 = 'test.sam'
    path_output =  'test.bai'
    with pytest.raises(RuntimeError):
        make_index(path_input_0, path_output)
    with pytest.raises(ValueError):
        make_index(path_input_1, path_output)
    
