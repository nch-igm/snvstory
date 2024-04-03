import sys, os
import gzip
import logging
import numpy as np
from scipy.sparse import csr_matrix, vstack
from vcf_to_sparse import vcf_to_json, json_to_sparse_matrix


def get_sample_names(path_input):
    """
    Retrieves the sample names from a VCF file.

    Args:
        path_input (str): The path to the VCF file.
        return_sample_names (bool, optional): Whether to return the sample names as a list. 
            Defaults to False.

    Returns:
        list or None: The sample names as a list if `return_sample_names` is True, 
            otherwise None.
    """

    if path_input.endswith(tuple(['.vcf', '.vcf.gz'])):
        logging.debug("Checking vcf sample composition")
        o, gz_file = get_file_handle(path_input)
        for line in o:
            if gz_file:
                line = line.decode('utf-8').strip()
            if line.startswith('#CHROM'):
                line = line.strip()
                break

        aline = np.asarray(line.split('\t'))
        index = np.where(aline == 'FORMAT')[0][0]
        indices = np.arange(index + 1, len(aline)).tolist()

        return (indices, list(aline[indices]))

def get_file_handle(path):
    """
    Runtime check if file exists. Opens vcf(.gz) file and reads lines.

    Args:
        path (str): The path to the file.

    Returns:
        tuple: A tuple containing the list of lines from the file and a boolean indicating if the file is gzipped.

    Raises:
        Exception: If the file cannot be opened or does not exist.

    """
    gz_file = False
    try:
        if path.endswith('.gz'):
            with gzip.open(path, 'r') as fin:
                o = fin.readlines()
            gz_file = True
            return o, gz_file
        else:
            with open(path, 'r') as fin:
                o = fin.readlines()
            return o, gz_file
    except Exception as e:
        print(f'Unable to open file or file does not exist. {path}. {e}')
        return


def parse_vcf(o, gz_file, sample_index):
    """
    Parse a VCF file and extract the genotype information for a specific sample.

    Args:
        o (iterable): An iterable object containing the lines of the VCF file.
        gz_file (bool): A boolean indicating whether the VCF file is gzipped.
        sample_index (int): The index of the sample to extract the genotype information from.

    Returns:
        list: A list of parsed VCF lines, where each line contains the genotype information for the specified sample.

    Raises:
        Exception: If there is an error parsing the sample at the specified position.

    """
    parsed_vcf = []
    for line in o:
        if gz_file:
            line = line.decode('utf-8').strip()
        if not line.startswith('#'):
            line = line.strip()
            # assume the sample has been converted to the numeric column
            try:
                sline = line.split('\t')
                sample_gt = sline[sample_index]
                nline = '\t'.join(sline[:9] + [sample_gt])  # I assume that the samples start after the FORMAT field which is the 8th field by zero index.
                parsed_vcf.append(nline)
            except Exception as e:
                print(f'Failed to parse sample at position {sample_index} from line {line}. {e}')
                return
    return parsed_vcf


def process_vcf(vcf_path, ordered_snps, variant_container, locus_converter_path):
    """
    Process a VCF file and convert it into a sparse matrix.

    Args:
        vcf_path (str): The path to the VCF file.
        ordered_snps (str): The path to the file containing ordered SNPs.
        variant_container (object): The variant container object.
        locus_converter (object): The locus converter path.

    Returns:
        tuple: A tuple containing the sparse matrix of shape (n_samples, n_snps) and the list of sample names.
    """
    sample_indices, sample_names = get_sample_names(vcf_path)
    o, gz_file = get_file_handle(vcf_path)

    s_matrix_list = []
    for sample_index in sample_indices:
        parsed_vcf = parse_vcf(o, gz_file, sample_index)
        t_vcf_json = vcf_to_json(parsed_vcf, variant_container, locus_converter_path)
        s_matrix = json_to_sparse_matrix(t_vcf_json, ordered_snps)
        s_matrix_list.append(s_matrix)

    s_matrix = vstack(s_matrix_list)

    return s_matrix, sample_names

    

    

        





