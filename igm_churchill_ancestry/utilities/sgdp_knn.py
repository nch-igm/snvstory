import sys, os
from igm_churchill_ancestry.pipelines.variables import variables
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree as KDTree
import json
import pickle
from scipy import sparse

''' Plot sample with simons data, output table with top nearest neighbors '''

def create_distance_table(emb, att_dir):
    npx = emb.to_numpy()
    # This will create the index
    k = KDTree(npx)
    # Get all pairwise distances to create adjacency matrix
    (dists, idxs) = k.query(npx, npx.shape[0])

    with open(os.path.join(att_dir,'sgdp_country.json'), 'r') as f:
        country_dict = json.load(f)
    with open(os.path.join(att_dir,'sgdp_country.json'), 'r') as f:
        sample_dict = json.load(f)

    # last entry is target sample
    countries = np.vectorize(country_dict.get)(idxs[-1])
    samples = np.vectorize(sample_dict.get)(idxs[-1])
    distances = dists[-1]
    inds = idxs[-1]

    kdist_df = pd.DataFrame({'idxs' : inds, 
                            'Sample' : samples, 
                            'Country' : countries, 
                            'Distance' : distances})
    return(kdist_df)

# Broken
# The way umap is saved is not compatible? There's some issue with reading in umap file
def sgdp_knn(s_matrix, ml_dir, attribute_dir, outdir, sample):
    try:
        with open(os.path.join(ml_dir, "sgdp_svd_fit.pkl"), 'rb') as f:
            svd_fit = pickle.load(f)
        # Does not work
        with open(os.path.join(ml_dir, "sgdp_umap_fit.pkl"), 'rb') as f:
            umap_fit = pickle.load(f)
        model_embedding = pd.read_csv(os.path.join(ml_dir, "sgdp_embedding.txt"), index_col=0)
    except Exception as e:
        print(f'Cannot load models, check path. {e}')
        return
    # Not sure why I did this...
    if isinstance(s_matrix, sparse.csr.csr_matrix):
        sparse_matrix = s_matrix
    else:
        print(f'Cannot coerce data. Check sparse matrix: {s_matrix}')
        return

    # Transform data
    xsvd = svd_fit.transform(sparse_matrix)
    input_embedding = umap_fit.transform(xsvd)

    np_concat = np.concatenate([model_embedding, input_embedding])
    kdist_df = create_distance_table(np_concat, attribute_dir)
    kdist_df.to_csv(os.path.join(outdir, sample + "_knn.csv")) 
