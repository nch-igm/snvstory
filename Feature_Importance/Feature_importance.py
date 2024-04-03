import sys, os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import matplotlib.patches as mpatches
import seaborn as sns

import shap
import xgboost as xgb
from scipy import sparse
import json
import argparse

#from shap_calculations import validate_shap_values, calculate_feature_importance
from Parse_vcf import process_vcf
from Feature_importance_collection import calculate_feature_importance, agg_shap_genome_region, get_shap_abs_mean, get_shap_abs, bar_plot_max_abs, stacked_bar_plot, IdeogramPlotter


def load_files(resource_dir, genome_ver):
    """
    Load annotation files and return necessary data.

    Args:
        resource_dir (str): The directory path where the annotation files are located.
        genome_ver (str): The version of the genome.

    Returns:
        tuple: A tuple containing the following elements:
            - clf (xgb.Booster): The XGBoost model loaded from the specified path.
            - variant_container (dict): A dictionary containing variant information loaded from a JSON file.
            - o_snps (list): A list of ordered SNPs loaded from a text file.
            - var_convert_cyto (dict): A dictionary containing variant-to-cytoband conversion information loaded from a JSON file.
            - var_convert_genes (dict): A dictionary containing variant-to-gene conversion information loaded from a JSON file.
            - locus_converter_path (str or None): The path to the locus converter file, or None if the genome version is not 'hg38'.
            - label_dict (dict): A dictionary mapping label indices to label names.

    Raises:
        SystemExit: If there is an error loading any of the files.
    """

    print('Loading annotation files...')

    gnomad_cont_model = f'{resource_dir}/continental/machine_learning_models/xgb_model_7.50E-49.bin'
    clf = xgb.Booster({'nthread': 1})
    try:
        clf.load_model(gnomad_cont_model)
    except Exception as e:
        print(f'Cannot load model, check path: {gnomad_cont_model}. {e}')
        sys.exit(1)

    gnomad_features = f'{resource_dir}/continental/matrix_attributes/continental_features.json'
    try:
        with open(gnomad_features, 'r') as json_file:
            variant_container = json.load(json_file)
    except Exception as e:
        print(f'Error reading file: {gnomad_features}. {e}')
        sys.exit(1)

    ordered_snps = f'{resource_dir}/continental/matrix_attributes/77402_ordered_snps.txt'
    o_snps = []
    try:
        with open(ordered_snps) as f:
            for line in f:
                line = line.strip()
                o_snps.append(line)
    except Exception as e:
        print(f'File path probably incorrect: {ordered_snps}. {e}')
        sys.exit(1)
 
    cyto_conversion = f'{resource_dir}/feature_importance/var_cyto_conversion.json'
    try:
        with open(cyto_conversion, 'r') as json_file:
            var_convert_cyto = json.load(json_file)
    except Exception as e:
        print(f'Error reading file: {e}')
        sys.exit(1)

    gene_conversion = f'{resource_dir}/feature_importance/var_genes_conversion.json'
    try:
        with open(gene_conversion, 'r') as json_file:
            var_convert_genes = json.load(json_file)
    except Exception as e:
        print(f'Error reading file: {gene_conversion}. {e}')
        sys.exit(1)
    
    if genome_ver == 'hg38':
        locus_converter_path = f'{resource_dir}/genome_ver_converters/hg38tob37/hg38_liftoverAIMs.json'
    else:
        locus_converter_path = None

    label_dict = {0:'afr', 1:'amr', 2:'asj', 3:'eas', 4:'eur', 5:'sas'}

    return (clf, variant_container, o_snps, var_convert_cyto, var_convert_genes, locus_converter_path, label_dict)
    

def get_chrom_sizes(chrom_sizes):
    """
    Reads a file containing chromosome sizes and returns a DataFrame.

    Parameters:
    chrom_sizes (str): The file path of the chromosome sizes file.

    Returns:
    pandas.DataFrame: A DataFrame containing the chromosome sizes, start and end positions, width, and color.

    """
    try:
        df_chrom_sizes = pd.read_csv(chrom_sizes, delim_whitespace=True)
    except Exception as e:
        print(f'Error reading file: {chrom_sizes}. {e}')
        sys.exit(1)

    df_chrom_sizes.columns = ['chrom', 'end']
    df_chrom_sizes['start'] = 0
    df_chrom_sizes['end'] = df_chrom_sizes['end'].astype(int)
    df_chrom_sizes = df_chrom_sizes[['chrom', 'start', 'end']]
    df_chrom_sizes['width'] = df_chrom_sizes['end'] - df_chrom_sizes['start']  
    df_chrom_sizes['color'] = '#D3D3D3'

    return df_chrom_sizes


def write_npz(shap_values, sample_list, var_ids, label_dict, output_path):
    """
    Save aggregated SHAP values to a compressed numpy array (.npz) file.

    Args:
        shap_values (list): List of SHAP values for each feature.
        sample_list (list): List of sample names.
        var_ids (list): List of variable IDs.
        label_dict (dict): Dictionary mapping variable IDs to labels.
        output_path (str): Path to save the output file.

    Returns:
        None

    Raises:
        Exception: If there is an error saving the compressed numpy array.
    """

    dim_labels = np.array([
        list(label_dict.values()),
        sample_list,
        var_ids
    ], dtype=object)
    try:
        np.savez_compressed(output_path, array=np.stack(shap_values, axis=0), labels=dim_labels)
        print(f'Saved aggregated shap values to {output_path}')
    except Exception as e:
        print(f'Error saving compressed numpy array: {output_path}. {e}')
        sys.exit(1)
    

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate feature importance aggregated to genes and cytolocations. Returns two .npz files with shap values. Optionally create summary plots.')
    parser.add_argument('vcf', type=str, help='Path to input VCF file.')
    parser.add_argument('output', type=str, help='Output folder to write results. Will exit if folder exists.')
    parser.add_argument('resource_dir', type=str, help='Path to resource directory.')
    parser.add_argument('genome_ver', type=str, choices=['hg19', 'hg38'], help='Genome version (hg19 or hg38) of input vcf.')
    parser.add_argument('-b', '--bar-plot', action='store_true', help='Flag to indicate whether to create mean(|SHAP val|) bar plot. All samples in the VCF are aggregated together.')
    parser.add_argument('-s', '--stacked-bar-plot', action='store_true', help='Flag to indicate whether to create stacked bar plot. All samples in the VCF are aggregated together.')
    parser.add_argument('-i', '--ideogram-plot', action='store_true', help='Flag to indicate whether to create mean(|SHAP val|) bar plot. This will create a separate plot for each sample in the VCF.')
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Check if the output path exists
    if os.path.exists(args.output):
        print(f"Error: The path '{args.output}' already exists.")
        sys.exit(1)
    os.makedirs(args.output, exist_ok=False)
    
    # Load in the model and anotation files
    clf, variant_container, ordered_snps, var_convert_cyto, var_convert_genes, locus_converter_path, label_dict = load_files(args.resource_dir, args.genome_ver)

    print('Converting VCF to sparse...')
    s_matrix, sample_list = process_vcf(args.vcf, ordered_snps, variant_container, locus_converter_path)
    shap_values = calculate_feature_importance(s_matrix, clf, ordered_snps, sample_list)
    shap_cyto, cyto_ids = agg_shap_genome_region(shap_values, ordered_snps, sample_list, var_convert_cyto, label_dict)
    shap_gene, gene_ids = agg_shap_genome_region(shap_values, ordered_snps, sample_list, var_convert_genes, label_dict)

    # Write out aggregated shap values
    write_npz(shap_cyto, sample_list, cyto_ids, label_dict, f'{args.output}/shap_cyto.npz')
    write_npz(shap_gene, sample_list, gene_ids, label_dict, f'{args.output}/shap_gene.npz')
    
    if args.bar_plot:
        shap_genes_abs = get_shap_abs(shap_gene)
        shap_genes_mean = get_shap_abs_mean(shap_genes_abs)
        bar_plot_max_abs(shap_genes_mean, gene_ids, 20, f'{args.output}/Top_20_genes.png')
    
    if args.stacked_bar_plot:
        shap_genes_abs = get_shap_abs(shap_gene)
        shape_genes_mean = get_shap_abs_mean(shap_genes_abs)
        stacked_bar_plot(shap_genes_abs, shap_genes_mean, gene_ids, 20, label_dict, f'{args.output}/Top_20_genes_stacked_label.png')
    
    if args.ideogram_plot:
        chromosome_list = ['chr%s' % i for i in list(range(1, 23)) + ['X', 'Y']]
        df_chrom_sizes = get_chrom_sizes(f'{args.resource_dir}/feature_importance/Hg19_chrom_sizes.txt')

        ideogram = IdeogramPlotter(chromosome_list, df_chrom_sizes)
        for sample_ind, sample in enumerate(sample_list):
            shap_indiv = []
            # Get shap values for each sample separately
            for i in range(len(shap_cyto)):
                shap_indiv.append(shap_cyto[i][sample_ind])
            df_indiv = pd.DataFrame(shap_indiv).T
            df_indiv.columns = label_dict.values()
            df_indiv.index = cyto_ids
            ideogram.plot_ideo(df_indiv, f'{args.output}/{sample}_ideogram.png')

if __name__ == '__main__':
    main()


