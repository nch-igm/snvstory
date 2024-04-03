import sys
import numpy as np
import pandas as pd
from unittest import TestCase, main
from unittest.mock import patch


# Import the functions to be tested
sys.path.append(str('../'))
from Feature_importance import get_chrom_sizes 

class FeatureImportanceTestCase(TestCase):
    def setUp(self):
        self.resource_dir = "resource_dir/"
        self.genome_ver = "hg38"
        self.output_path = "output_test/"
        self.shap_values = [np.array([[1, 2, 5], [3, 4, 6]]), np.array([[5, 6, 9], [7, 8, 10]])]
        self.sample_list = ["sample1", "sample2"]
        self.var_ids = ["var1", "var2", "var3"]
        self.label_dict = {0:'afr', 1:'amr', 2:'asj', 3:'eas', 4:'eur', 5:'sas'}

    def test_get_chrom_sizes(self):
        # Mock the chrom_sizes file
        with patch("pandas.read_csv") as mock_read_csv:
            mock_read_csv.return_value = pd.DataFrame(
                {"chrom": ["chr1", "chr2", "chr3"], "start": [0, 0, 1], "end": [100, 200, 300]}
            )
            df_chrom_sizes = get_chrom_sizes(f"path/to/chrom_sizes.txt")

        # Assert the returned DataFrame
        self.assertIsInstance(df_chrom_sizes, pd.DataFrame)
        self.assertEqual(df_chrom_sizes.shape, (3, 5))
        self.assertEqual(df_chrom_sizes["chrom"].tolist(), ["chr1", "chr2", "chr3"])
        self.assertEqual(df_chrom_sizes["start"].tolist(), [0, 0, 1])
        self.assertEqual(df_chrom_sizes["end"].tolist(), [100, 200, 300])
        self.assertEqual(df_chrom_sizes["width"].tolist(), [100, 200, 299])
        self.assertEqual(df_chrom_sizes["color"].tolist(), ["#D3D3D3", "#D3D3D3", "#D3D3D3"])


if __name__ == "__main__":
    main()
