import sys
from unittest import TestCase, main

sys.path.append(str('../'))
from Parse_vcf import get_sample_names, get_file_handle

class ParseVcfTestCase(TestCase):
    def setUp(self):
        self.path_input = "input_test.vcf"
        self.return_sample_names = False

    def test_get_sample_names(self):
        
        indices, sample_names = get_sample_names(self.path_input)

        # Assert the returned values
        self.assertIsInstance(indices, list)
        self.assertEqual(indices, [9, 10, 11])
        self.assertIsInstance(sample_names, list)
        self.assertEqual(sample_names, ["sample1", "sample2", "sample3"])


    def test_get_sample_names_invalid_file_extension(self):
        # Set an invalid file path
        path_input = "path/to/input.txt"

        # Call the function and assert that it exits with an error message
        with self.assertRaises(SystemExit):
            get_sample_names(path_input)


    def test_get_file_handle(self):
        o, gz_file = get_file_handle(self.path_input)

        # Assert the returned values
        self.assertIsInstance(o, list)
        self.assertEqual(o, ["##fileformat=VCFv4.2\n", "##FILTER=<ID=PASS,Description=\"All filters passed\">\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n"])
        self.assertFalse(gz_file)


    def test_get_file_handle_invalid_file(self):
        # Set an invalid file path
        path = "/path/to/invalid/file.txt"

        # Call the function and assert that it exits with an error message
        with self.assertRaises(SystemExit):
            get_file_handle(path)

if __name__ == "__main__":
    main()