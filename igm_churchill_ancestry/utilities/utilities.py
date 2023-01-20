from igm_churchill_ancestry.pipelines.variables import variables
import gzip
import numpy as np
import logging
import string
import os
import sys
log = logging.getLogger(__name__)
np.random.seed(0)  # guaranteed to have at least 10000 unique filecodes


def check_resources(var):
    """
    Each model expects certain auxillary files necessary to complete prediction
    This function checks to see whether those files are accessible based on the generated
    variable class.
    """
    for x in var.MATRIX_ATT_DIRS + var.MODEL_DIRS:
        if os.path.isdir(x) is False:
            print(f"Failed to find the resource directory {x}")
            sys.exit()


def filter_extension(path):
    """
    filter for files that only end with
    approved extension and return
    these files.

    args
    ----
    path - list of paths (str)

    returns
    -------
    ext - returns the matched extension (str)
          or signals that the path is a directory

    """
    path = filter(lambda x: x.endswith(variables.EXTENSIONS), path)
    return list(path)


def get_extension(path):
    """
    check for files that only end with
    approved extension and return
    the extension.

    args
    ----
    path - path to vcf type  /path/to/sample.vcf (str)

    returns
    -------
    ext - returns the matched extension (str)
          or signals that the path is a directory

    """
    extensions = [x for x in variables.EXTENSIONS if path.endswith(x)]
    if len(extensions) == 1:
        return extensions[0]
    else:
        log.debug(f"{path} matches multiple or no extensions {variables.EXTENSIONS}. Assuming directory")
        return 'dir'


def random_file_code():
    rnd_str = ''.join(np.random.choice(list(string.ascii_uppercase) + list(string.ascii_lowercase), size=5).tolist())
    return rnd_str


def genotype_dictionary(var):
    """
    A function to build the genotype dictionary that
    converts the raw vcf version of a genotype into a
    numeric representation

    args
    ----
    var - variant e.g. 0/0, 0|1

    returns
    -------
    int - int representation of genotype
    """
    if var[1] == '/':
        split_gt = var.strip().split("/")
    else:
        split_gt = var.strip().split("|")
    if split_gt[0] == '0' and split_gt[1] == '0':
        return 0
    elif split_gt[0] == '.' and split_gt[1] == '.':
        return 0
    elif split_gt[0] != split_gt[1]:
        return 1
    elif split_gt[0] == split_gt[1]:
        return 2


def get_file_handle(path):
    """
    A method to determine gzip status of file
    and load the file into memory

    args
    ----
    path - path to vcf

    returns
    -------
    gz_file - boolean represents whether file is gzipped
    o - vcf file read into memory and broken into lines
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


def is_vcf_multisample(path_input, return_sample_names=False):
    """
    Runtime check for vcf type e.g. multisample or single

    args
    ----
    path_input - path to vcf/gvcf file
    return_sample_names - bool to return sample names

    returns
    -------
    Boolean - True if it is multisample and false if it is not
    sample_names - optionally will return a list of the sample name(s)
    """
    if path_input.endswith(variables.EXTENSIONS):
        logging.debug("Checking vcf sample composition")
        o, gz_file = get_file_handle(path_input)
        for num, line in enumerate(o):
            if gz_file:
                line = line.decode('utf-8').strip()
            if line.startswith('#CHROM'):
                line = line.strip()
                break
        aline = np.asarray(line.split('\t'))
        # we assume samples always come after the FORMAT field
        index = np.where(aline == 'FORMAT')[0][0]
        # calculate the number of samples
        n_samples = len(aline[index + 1:])
        if n_samples > 1:
            if return_sample_names is False:
                return True
            else:
                return True, list(aline[index + 1:])
        else:
            if return_sample_names is False:
                return False
            else:
                return False, list(aline[index + 1:])
