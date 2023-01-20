from igm_churchill_ancestry.utilities.parsing import parse_vcf, parse_multisample_vcf_sample, parse_multisample_vcf
from igm_churchill_ancestry.utilities.vcf2sparse import vcf_to_json, load_snp_order, json_to_sparse_matrix
from igm_churchill_ancestry.utilities.utilities import get_file_handle
from igm_churchill_ancestry.utilities.plot_ancestry import plot_parser
from igm_churchill_ancestry.utilities.plot_umap import plot_umap_parser
# from igm_churchill_ancestry.utilities.sgdp_knn import sgdp_knn
import pandas as pd
import pickle
import xgboost as xgb
import numpy as np
import scipy
from scipy import sparse
import glob
import umap
import os

## Sklearn falsely warns about unpickling an estimator from a version other than what was used to build the model.
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)


def predict_ancestry(s_matrix, ml_dir, n_classes, model_type='xgb'):
    """
    Predicts ancestry at different geographic resolution


    args
    ----
    s_matrix - sparse matrix generated from the vcf
    data_dir - location of the data directory
    ml_dir - location of the model directroy
    n_classes - number of classes
    model_type - str value can be c for continental, s for
                 subcontinental, or k for 1000genomes

    returns
    -------
    yprob - a numpy array of the probabilities for each ancestry
    ylabel - the numeric label that maps to the ancestral label

    """
    if model_type == 'xgb':
        # Prep the model
        model_path = glob.glob(ml_dir + '/*.bin')[0]
        bst = xgb.Booster({'nthread': 1})  # static
        try:
            bst.load_model(model_path)
        except Exception as e:
            print(f'Cannot load model, check path: {model_path}. {e}')
            return
        # Sparse and numpy matrix needs to be converted
        if isinstance(s_matrix, xgb.core.DMatrix):
            dtest = s_matrix
        else:
            try:
                dtest = xgb.DMatrix(s_matrix, label=np.zeros(s_matrix.shape[0]))
            except Exception as e:
                print(f'Cannot coerce data. Check sparse matrix: {e}')
                return

        yprob = bst.predict(dtest).reshape(s_matrix.shape[0], n_classes)
        ylabel = np.argmax(yprob, axis=1)[0]

    elif model_type == 'svm':
        # Prep the model
        model_path = glob.glob(ml_dir + '/*.p')[0]
        try:
            with open(model_path, 'rb') as fin:
                clf = pickle.load(fin)
        except Exception as e:
            print(f'Cannot load model, check path: {model_path}. {e}')
            return
        # Sparse and numpy matrix needs to be converted
        if isinstance(s_matrix, scipy.sparse.csr.csr_matrix):
            dtest = s_matrix
        else:
            print(f'Cannot coerce data. Check sparse matrix: {s_matrix}')
            return

        try:
            yprob = clf.predict_proba(dtest).reshape(s_matrix.shape[0], n_classes)
        except ValueError:
            yprob = clf.predict_proba(dtest.A).reshape(s_matrix.shape[0], n_classes)
        ylabel = np.argmax(yprob, axis=1)[0]

    return (yprob, ylabel)


def run_ancestry_pipeline(vcf_path, multi_sample_status, sample, sample_position, var, outdir, genome_ver, mode, ofn):

    # Read VCF-type file into memory
    o, gz_file = get_file_handle(vcf_path)
    print(f'Gzipped status: {gz_file}')

    # Single sample analysis
    if multi_sample_status is False:
        DATA_TO_PLOT = {}
        p_vcf = parse_vcf(o, gz_file)
        print('Getting model predictions:')
        for att_dir, ml_dir, n_classes, m_type in var.R_DIRS:
            print(m_type)
            # Make the vcf ml compatible
            t = m_type.split('_')[0]
            t_vcf_json = vcf_to_json(parsed_vcf=p_vcf, attribute_dir=att_dir, locus_converter_json_path=var.JSON_CONVERTS[genome_ver][mode][t])
            o_snps = load_snp_order(attribute_dir=att_dir)
            s_matrix = json_to_sparse_matrix(t_vcf_json, o_snps)
            # Ancestry prediction
            if 'gnomAD' in m_type:
                yprob, ylabel = predict_ancestry(s_matrix, ml_dir=ml_dir, n_classes=n_classes)
            else:
                yprob, ylabel = predict_ancestry(s_matrix, ml_dir=ml_dir, n_classes=n_classes, model_type='svm')
            # UMAP plotting
            if 'continental' in m_type:
                plot_umap_parser(s_matrix, ml_dir=ml_dir, att_dir=att_dir, sample_name=sample, outdir=outdir, m_type=m_type)

            DATA_TO_PLOT[m_type] = yprob.flatten().tolist()
        df_data = pd.DataFrame([DATA_TO_PLOT], sample)

    # Multisample anaylsis
    else:
        DATA_TO_PLOT_LIST = []
        if sample_position == 'all':
            samples, sample_names = parse_multisample_vcf_sample(o, gz_file, sample_position)
            for s in samples:
                DATA_TO_PLOT = {}
                p_mvcf = parse_multisample_vcf(o, gz_file, s)
                for att_dir, ml_dir, n_classes, m_type in var.R_DIRS:
                    print(m_type)
                    # Make the vcf ml compatible
                    t = m_type.split('_')[0]
                    t_vcf_json = vcf_to_json(parsed_vcf=p_mvcf, attribute_dir=att_dir, locus_converter_json_path=var.JSON_CONVERTS[genome_ver][mode][t])
                    o_snps = load_snp_order(attribute_dir=att_dir)
                    s_matrix = json_to_sparse_matrix(t_vcf_json, o_snps)
                    if 'gnomAD' in m_type:
                        yprob, ylabel = predict_ancestry(s_matrix, ml_dir=ml_dir, n_classes=n_classes)
                    else:
                        yprob, ylabel = predict_ancestry(s_matrix, ml_dir=ml_dir, n_classes=n_classes, model_type='svm')
                    DATA_TO_PLOT[m_type] = yprob.flatten().tolist()
                DATA_TO_PLOT_LIST.append(DATA_TO_PLOT)
            df_data = pd.DataFrame(DATA_TO_PLOT_LIST, sample_names)
        else:
            # User chosen sample(s); if handpicked should be comma separated.
            sample_position_s = sample_position.split(',')
            for s in sample_position_s:
                DATA_TO_PLOT = {}
                sample, sample_names = parse_multisample_vcf_sample(o, gz_file, s)
                p_mvcf = parse_multisample_vcf(o, gz_file, sample)
                for att_dir, ml_dir, n_classes, m_type in var.R_DIRS:
                    t = m_type.split('_')[0]
                    # Make the vcf ml compatible
                    t_vcf_json = vcf_to_json(parsed_vcf=p_mvcf, attribute_dir=att_dir, locus_converter_json_path=var.JSON_CONVERTS[genome_ver][mode][t])
                    o_snps = load_snp_order(attribute_dir=att_dir)
                    s_matrix = json_to_sparse_matrix(t_vcf_json, o_snps)
                    if 'gnomAD' in m_type:
                        yprob, ylabel = predict_ancestry(s_matrix, ml_dir=ml_dir, n_classes=n_classes)
                    else:
                        yprob, ylabel = predict_ancestry(s_matrix, ml_dir=ml_dir, n_classes=n_classes, model_type='svm')
                    DATA_TO_PLOT[m_type] = yprob.flatten().tolist()
                DATA_TO_PLOT_LIST.append(DATA_TO_PLOT)
            df_data = pd.DataFrame(DATA_TO_PLOT_LIST, sample_names)

    # begin the plotting and figure writing
    if not df_data.empty:
        # Placeholder plotting
        plot_parser(df_data, var, outdir, ofn)
    else:
        return
