# import os
# import uuid

# import pytest

# from igm_churchill_cnvloh.pipelines.tumor import tumor_pipeline
# from igm_churchill_cnvloh.pipelines.normal import normal_pipeline
# from igm_churchill_cnvloh.pipelines.tumor_normal import tumor_normal_pipeline


# def create_tmp_file():
#     file_name = f"{uuid.uuid4().hex}.txt"
#     with open(file_name, 'w') as fn:
#         fn.write("temp file data")
#     return file_name

# def test_normal_pipeline_return_true():
#     """Test that normal pipeline returns True if a local file exists."""
#     fn = create_tmp_file()
#     output = normal_pipeline(fn)
#     assert output == True
#     os.remove(fn)

# def test_tumor_pipeline_return_true():
#     """Test that normal pipeline returns True if a local file exists."""
#     fn = create_tmp_file()
#     output = tumor_pipeline(fn)
#     assert output == True
#     os.remove(fn)

# def test_tumor_normal_pipeline_return_true():
#     """Test that normal pipeline returns True if a local file exists."""
#     fn1 = create_tmp_file()
#     fn2 = create_tmp_file()
#     output = tumor_normal_pipeline(fn1, fn2)
#     assert output == True
#     os.remove(fn1)
#     os.remove(fn2)