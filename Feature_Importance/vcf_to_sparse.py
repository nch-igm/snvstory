import sys
import json
import warnings
import pandas as pd
from scipy import sparse

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


def vcf_to_json(parsed_vcf, variant_container, locus_converter_json_path):
    """
    Transforms genotypes from VCF into numeric representation and
    fills a JSON containing ancestry informative alleles with those
    numeric values.

    Returns a JSON filled with the VCF data
    """
    
    if locus_converter_json_path is not None:
        try:
            with open(locus_converter_json_path, 'r') as myfile:
                locus_converter = json.load(myfile)
        except Exception as e:
            print(f"Failed to load locus converter error: {e}")
            sys.exit(1)
    else:
        locus_converter = None
    gt_sum = 0
    for k in parsed_vcf:
        values = k.split("\t")
        chrom = values[0]
        pos = values[1]
        ref = values[3]
        alt = values[4]
        try:
            locus_id = chrom + "_" + pos + "_" + ref + "_" + alt
        except Exception as e:
            print(f"VCF is malformed, not able to generate locus id. {e} ")
            return
        # convert the locus if we are using locus vcf
        if locus_converter is not None and locus_id in locus_converter:
            locus_id = locus_converter[locus_id]
        if locus_id in variant_container:
            try:
                genotype = values[9][:3:]
            except Exception as e:
                print(f"VCF is malformed, not able to extract genotype. {e} ")
                return
            try:
                genotype_int = genotype_dictionary(genotype)
                gt_sum += genotype_int
            except Exception as e:
                print(f"Unknown genotype: {genotype}. {e}")
                return
            variant_container[locus_id] = genotype_int
    if gt_sum > 0:
        return variant_container
    else:
        warnings.warn("No genotypes present. Check to make sure VCF has variants")
        return variant_container


def json_to_sparse_matrix(g_container, o_snps):
    """
    Convert the raw JSON data into a sparse row matrix
    that is patelable for XGBoost.
    """
    try:
        df = pd.DataFrame([g_container])
    except Exception as e:
        print(f'DataFrame unable to be built, something wrong with json length: {len(g_container)}. {e}')
        return
    try:
        df = df[o_snps]
    except Exception as e:
        print(f'Check ordered length of SNPs. {e}')
        return
    # coerce to numeric
    try:
        vals = df.values.astype(int).flatten()
    except Exception as e:
        print(f'Unable to coerce to numeric. {e}')
        return
    try:
        s_matrix = sparse.csr_matrix(vals)
    except Exception as e:
        print(f'Cannot make matrix sparse. {e}')
        return
    return s_matrix