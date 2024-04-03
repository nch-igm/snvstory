import sys, os

import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import matplotlib.patches as mpatches
import seaborn as sns

def validate_shap_values(shap_values, gnomad_vars, sample_list):
    """
    Validates the shap_values list for a given set of gnomad_vars and sample_list.

    Args:
        shap_values (list): A list of length 6 containing the shapley values.
        gnomad_vars (list): A list of variables.
        sample_list (list): A list of samples.

    Raises:
        ValueError: If shap_values is not a list of length 6 or if any element in shap_values is not a numpy array of shape (len(sample_list), len(gnomad_vars)).

    Returns:
        None
    """

    # Check if shap_values is a list of length 6
    if not isinstance(shap_values, list) or len(shap_values) != 6:
        raise ValueError("shap_values is not a list of length 6")
    
    # Check if each element in shap_values is a numpy array of shape (len(sample_list), len(gnomad_vars))
    for i, arr in enumerate(shap_values):
        if not isinstance(arr, np.ndarray) or arr.shape != (len(sample_list), len(gnomad_vars)):
            raise ValueError(f"Error: shap_values[{i}] is not a numpy array of shape ({len(sample_list)}, {len(gnomad_vars)})")
    
    print(f"SHAP values are valid array of shape ({len(sample_list)}, {len(gnomad_vars)}).")


def calculate_feature_importance(s_matrix, clf, ordered_snps, sample_list):
    """
    Calculate the feature importance using SHAP values.

    Parameters:
    s_matrix (numpy.ndarray): The feature matrix.
    clf (object): The trained classifier model.
    ordered_snps (list): The ordered list of SNPs.
    sample_list (list): The list of samples.

    Returns:
    numpy.ndarray: The SHAP values representing the feature importance.
    """

    explainer = shap.TreeExplainer(clf)
    shap_values = explainer.shap_values(s_matrix)
    validate_shap_values(shap_values, ordered_snps, sample_list)

    return shap_values


def agg_shap_genome_region(shap_values, variant_ids, sample_names, var_convert, label_dict):
    """
    Aggregate SHAP values across genomic region by the mean.

    Args:
        shap_values (list of numpy arrays): List of SHAP values for each SNV feature.
        variant_ids (list): List of SNV feature IDs.
        sample_names (list): List of sample names.
        var_convert (dict): Dictionary mapping SNVs to their respective genomic region (e.g., genes or cytolocations).
        label_dict (dict): Dictionary mapping label indices to label names.

    Returns:
        list of numpy arrays: List of aggregated SHAP values across genomic region with list length of 6 for each label and numpy array
        with shape (n_samples, n_genomic_regions).
        list: List of genomic region IDs corresponding to numpy array order.
    """

    df_list = []
    for i in range(len(shap_values)):
        sub_np = shap_values[i]
        df = pd.DataFrame(sub_np)
        df.index = [f'{x}_{label_dict[i]}' for x in sample_names]
        df_list.append(df)

    # Aggregate SHAP values across genomic region by the mean
    df_shap = pd.concat(df_list)
    df_shap.columns = variant_ids
    df_shap.columns = df_shap.columns.map(var_convert)
    df_shap = df_shap.groupby(by=df_shap.columns, axis=1).mean()

    # Convert DataFrame to np array
    n_labels = len(label_dict)
    n_samples = int(len(df_shap)/n_labels)
    shap_values = []
    for i in range(n_labels):
        i1 = i * n_samples
        i2 = (i + 1) * n_samples
        npa = df_shap.iloc[i1:i2, :].to_numpy()
        shap_values.append(npa)

    return (shap_values, df_shap.columns.to_list())


def get_shap_abs(shap_values):
    """
    Calculate the mean absolute SHAP values for each feature across all samples.

    Parameters:
    shap_values (list): A list of SHAP values for each label.

    Returns:
    numpy.ndarray: The absolute SHAP values for each feature, with shape (6, n_genomic_regions).
    """

    shap_abs = []

    for i in range(len(shap_values)):
        np_abs = np.absolute(shap_values[i])
        # If multiple samples, take the mean of the absolute SHAP values
        shap_m1 = np_abs.mean(axis=0)
        shap_abs.append(shap_m1)

    shap_abs = np.stack(shap_abs, axis=0)

    return shap_abs


def get_shap_abs_mean(shap_abs):
    """
    Calculate the mean absolute SHAP values for each feature across labels.

    Parameters:
    shap_abs (numpy.ndarray): A 2D array of SHAP values for each label. Each row represents a label, and each column represents a feature.

    Returns:
    numpy.ndarray: A 1D array containing the mean absolute SHAP values for each feature with length n_genomic_regions.
    """

    shap_abs_mean = shap_abs.mean(axis=0)

    return shap_abs_mean


def bar_plot_max_abs(shap_abs_mean, feat_ids, n_indices, output):
    """
    Plot bar chart of absolute(SHAP values).

    Parameters:
    shap_vals (numpy.ndarray): Array of SHAP values.
    feat_ids (list): List of feature IDs.
    n_indices (int): Number of top features to plot.
    output (str): Output file path to save the plot.

    Returns:
    None
    """

    feat_indices = np.argsort(shap_abs_mean)[::-1][:n_indices]
    shap_abs_mean[feat_indices]

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(y=[feat_ids[i] for i in feat_indices], x=shap_abs_mean[feat_indices], orient='h', color='steelblue')
    ax.set_xlabel('Mean(|SHAP value|)', fontsize=14)
    ax.set_ylabel('Gene(s)', fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    plt.title('Feature Importance', fontsize=14)
    plt.tight_layout()
    plt.savefig(output, dpi=300)

    print(f'Feature importance bar plot saved to {output}')


def stacked_bar_plot(shap_abs, shap_abs_mean, feat_ids, n_indices, label_dict, output):
    """
    Plot a stacked bar chart of SHAP values across genomic regions.

    Args:
        shap_abs (list of numpy arrays): List of absolute SHAP values for each SNV feature.
        shap_abs_mean (numpy array): Mean absolute SHAP values for each SNV feature.
        feat_ids (list): List of feature IDs.
        n_indices (int): Number of top features to plot.
        label_dict (dict): Dictionary mapping label indices to label names.
        output (str): Output file path to save the plot.

    Returns:
        None
    """
    color_lookup = {
        'afr':'#B847A3',
        'amr':'#FBDF6C',
        'asj':'#2EDB7E',
        'eas':'#ED592A',
        'eur':'#2FA4DC',
        'sas':'#DC2E31'
        }  

    df_shap_abs = pd.DataFrame(shap_abs).T
    df_shap_abs.columns = label_dict.values()
    df_shap_abs.index = feat_ids

    feat_indices = np.argsort(shap_abs_mean)[::-1][:n_indices]
    df_plot = df_shap_abs.iloc[feat_indices]
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color_lookup[label]) for label in df_shap_abs.columns]

    plt.figure()
    df_plot.plot(kind='barh', stacked=True, figsize=(10,6), color=[color_lookup[label] for label in df_plot.columns]).invert_yaxis()
    plt.xlabel('|SHAP value|', fontsize=14)
    plt.ylabel('Gene(s)', fontsize=14)
    plt.tick_params(axis='both', labelsize=12)
    plt.title('Feature Importance', fontsize=14)
    plt.legend(legend_handles, df_shap_abs.columns)
    plt.tight_layout()
    plt.savefig(output, dpi=300)

    print(f'Stacked bar plot saved to {output}')


class IdeogramPlotter:
    """
    A class for plotting genome features.

    Parameters
    ----------
    chromosome_list : list
        List of chromosomes to be plotted.
    df_chrom_sizes : pandas.DataFrame
        DataFrame containing chromosome sizes.
    chrom_height : float, optional
        Height of each chromosome, by default 1.
    chrom_spacing : float, optional
        Spacing between chromosomes, by default 1.
    gene_height : float, optional
        Height of each gene, by default 0.4.
    gene_padding : float, optional
        Padding between genes, by default 0.1.
    figsize : tuple, optional
        Figure size, by default (6, 8).

    Methods
    -------
    chromosome_collections(df, y_positions, height, **kwargs)
        Yields BrokenBarHCollection of features that can be added to an Axes object.
    chromosome_base_tracker(chromosome_list)
        Calculates the base positions for SHAP values and chromosome visualization.
    plot_genome(df_positive_shap_values, df_negative_shap_values)
        Plots the genome features.
    """

    def __init__(self, chromosome_list, df_chrom_sizes, chrom_height=1, chrom_spacing=1, gene_height=0.4, gene_padding=0.1, shap_height=0.5,
                 color_lookup = {'afr':'#B847A3','amr':'#FBDF6C','asj':'#2EDB7E','eas':'#ED592A','eur':'#2FA4DC','sas':'#DC2E31'}, figsize=(12, 16)):
        self.chrom_height = chrom_height
        self.chrom_spacing = chrom_spacing
        self.gene_height = gene_height
        self.gene_padding = gene_padding
        self.shap_height = shap_height
        self.figsize = figsize
        self.color_lookup = color_lookup
        self.chromosome_list = chromosome_list
        self.df_chrom_sizes = df_chrom_sizes
    
    def get_shap_vals(self, df_shap):
        """
        Get the positive and negative SHAP values.

        Parameters:
        - df_shap (DataFrame): The DataFrame containing SHAP values.

        Returns:
        - df_shap_ideo_pos (DataFrame): The DataFrame containing positive SHAP values.
        - df_shap_ideo_neg (DataFrame): The DataFrame containing negative SHAP values.
        """
        df_shap['max_name'] = df_shap.iloc[:, 0:6].apply(lambda x: x.abs().idxmax(), axis=1)
        df_shap['max_val'] = df_shap.iloc[:, 0:6].agg({lambda x: max(x, key=abs)}, axis=1)

        df_shap_ideo = df_shap[['max_name', 'max_val']].copy()
        df_shap_ideo['chrom'] = df_shap_ideo.index.to_series().str.split('_').str[0].str.split('-').str[1]
        df_shap_ideo['start'] = df_shap_ideo.index.to_series().str.split('_').str[1].astype(int)
        df_shap_ideo['end'] = df_shap_ideo.index.to_series().str.split('_').str[2].astype(int)
        df_shap_ideo['color'] = df_shap_ideo['max_name'].map(self.color_lookup) 
        df_shap_ideo = df_shap_ideo[df_shap_ideo['max_val'] != 0]

        # Separate positive and negative shap values
        df_shap_ideo_pos = df_shap_ideo[df_shap_ideo['max_val'] > 0]
        df_shap_ideo_neg = df_shap_ideo[df_shap_ideo['max_val'] < 0]

        return df_shap_ideo_pos, df_shap_ideo_neg
        
    # Referenced from https://gist.github.com/daler/c98fc410282d7570efc3#file-ideograms-py
    def chromosome_collections(self, height, df, y_positions, **kwargs):
        """
        Yields BrokenBarHCollection of features that can be added to an Axes
        object.

        Parameters
        ----------
        df : pandas.DataFrame
            Must at least have columns ['chrom', 'start', 'end', 'color']. If no
            column 'width', it will be calculated from start/end.
        y_positions : dict
            Keys are chromosomes, values are y-value at which to anchor the
            BrokenBarHCollection
        height : float
            Height of each BrokenBarHCollection

        Additional kwargs are passed to BrokenBarHCollection
        """

        del_width = False
        if 'width' not in df.columns:
            del_width = True
            df['width'] = df['end'] - df['start']
        for chrom, group in df.groupby('chrom'):
            yrange = (y_positions[chrom], height)
            xranges = group[['start', 'width']].values
            yield BrokenBarHCollection(
                xranges, yrange, facecolors=group['color'], **kwargs)
        if del_width:
            del df['width']

    # Referenced from https://gist.github.com/daler/c98fc410282d7570efc3#file-ideograms-py
    def chromosome_base_tracker(self, chromosome_list):
        """
        Calculates the base positions for SHAP values and chromosome visualization.

        Args:
            chromosome_list (list): A list of chromosome names.

        Returns:
            tuple: A tuple containing three dictionaries:
                - chrom_centers (dict): A dictionary mapping chromosome names to the center positions for chromosome visualization.
                - chrom_ybase (dict): A dictionary mapping chromosome names to the base positions for chromosome visualization.
                - pos_shap_ybase (dict): A dictionary mapping chromosome names to the base positions for positive SHAP values.
                - neg_shap_ybase (dict): A dictionary mapping chromosome names to the base positions for negative SHAP values.
                
        """
        ybase = 0
        chrom_ybase = {}
        chrom_centers = {}

        # Iterate in reverse so that items in the beginning of `chromosome_list` will
        # appear at the top of the plot
        for chrom in chromosome_list[::-1]:
            chrom_ybase[chrom] = ybase
            chrom_centers[chrom] = ybase + self.chrom_height / 2.
            ybase += self.chrom_height + self.chrom_spacing

        # Offset positive and negative SHAP values from chromosome
        ybase = 0
        pos_shap_ybase = {}
        neg_shap_ybase = {}

        for chrom in chromosome_list[::-1]:
            pos_shap_ybase[chrom] = ybase + .5
            neg_shap_ybase[chrom] = ybase
            ybase += self.chrom_height + self.chrom_spacing

        return (chrom_centers, chrom_ybase, pos_shap_ybase, neg_shap_ybase)
    
    def plot_ideo(self, df_shap, output):
        """
        Plots the genome features.

        This method creates a figure and axes, adds the chromosome background,
        adds the gene collections, tweaks the axes, adds a legend, and shows the plot.

        Parameters:
        - df_shap: DataFrame
            The DataFrame containing the SHAP values for the genome features.

        Returns:
        None
        """

        # Create figure and axes
        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(111)

        chrom_centers, chrom_ybase, pos_shap_ybase, neg_shap_ybase = self.chromosome_base_tracker(self.chromosome_list)

        # Add grey chromosome background
        for collection in self.chromosome_collections(self.chrom_height, self.df_chrom_sizes, chrom_ybase):
            ax.add_collection(collection)

        # Add SHAP values
        df_shap_ideo_pos, df_shap_ideo_neg = self.get_shap_vals(df_shap)
        for collection in self.chromosome_collections(self.shap_height, df_shap_ideo_pos, pos_shap_ybase):
            ax.add_collection(collection)
        for collection in self.chromosome_collections(self.shap_height, df_shap_ideo_neg, neg_shap_ybase):
            ax.add_collection(collection)

        # Axes tweaking
        ax.set_yticks([chrom_centers[i] for i in self.chromosome_list])
        ax.set_yticklabels(self.chromosome_list)
        ax.axis('tight')

        afr_patch = mpatches.Patch(color=self.color_lookup['afr'], label='afr')
        amr_patch = mpatches.Patch(color=self.color_lookup['amr'], label='amr')
        asj_patch = mpatches.Patch(color=self.color_lookup['asj'], label='asj')
        eas_patch = mpatches.Patch(color=self.color_lookup['eas'], label='eas')
        eur_patch = mpatches.Patch(color=self.color_lookup['eur'], label='eur')
        sas_patch = mpatches.Patch(color=self.color_lookup['sas'], label='sas')
        plt.legend(handles=[afr_patch, amr_patch, asj_patch, eas_patch, eur_patch, sas_patch], loc='lower right')

        plt.savefig(output, dpi=300)

        print(f'Ideogram plot saved to {output}')
