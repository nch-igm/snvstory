import sys
import numpy as np
import pandas as pd
from unittest import TestCase, main
from unittest.mock import patch

sys.path.append(str('../'))
from Feature_importance_collection import calculate_feature_importance, agg_shap_genome_region, get_shap_abs_mean, get_shap_abs

#shap_values, variant_ids, sample_names, var_convert, label_dict
class FeatureImportanceCollectionTestCase(TestCase):
    def setUp(self):
        self.shap_values = [np.array([[1, -2, 5], [3, 4, 6]]), np.array([[5, -6, 9], [7, 8, 10]])]
        self.sample_list = ["sample1", "sample2"]
        self.var_ids = ["var1", "var2", "var3"]
        self.var_convert = {"var1": "gene1", "var2": "gene2", "var3": "gene3"}
        self.label_dict = {0:'afr', 1:'amr'}

    def test_agg_shap_genome_region(self):
        shap_values, gene_ids = agg_shap_genome_region(self.shap_values, self.var_ids, self.sample_list, self.var_convert, self.label_dict)

        # Assert the returned values
        self.assertIsInstance(shap_values, list)
        self.assertEqual(len(shap_values), 2)
        self.assertEqual(shap_values[0].shape, (2, 3))
        self.assertEqual(gene_ids, ['gene1', 'gene2', 'gene3'])

    def test_get_shap_abs(self):
        shap_abs = get_shap_abs(self.shap_values)

        # Assert the returned values
        self.assertIsInstance(shap_abs, np.ndarray)
        self.assertEqual(shap_abs.shape, (2, 3))
        self.assertEqual(shap_abs.tolist(), [[2.0, 3.0, 5.5], [6.0, 7.0, 9.5]])


    def test_get_shap_abs_mean(self):
        shap_abs = get_shap_abs(self.shap_values)
        shap_abs_mean = get_shap_abs_mean(shap_abs)

        # Assert the returned values
        self.assertIsInstance(shap_abs_mean, np.ndarray)
        self.assertEqual(shap_abs_mean.shape, (3,))
        self.assertEqual(shap_abs_mean.tolist(), [4.0, 5.0, 7.5])
    


if __name__ == "__main__":
    main()