import sys, os
from igm_churchill_ancestry.pipelines.variables import variables
import pandas as pd
import numpy as np
import pickle
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
import umap
import glob

from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, CDSView, GroupFilter, HoverTool, Legend

'''
Plot UMAP of sample fitted by all three models
'''


def load_pca(ml_dir):
    pca_path = glob.glob(ml_dir + f'/*svd*.pkl')[0]
    try:
        with open(pca_path, 'rb') as fin:
            pca_transformer = pickle.load(fin)
        return pca_transformer
    except Exception as e:
        print(f'Cannot load pickled PCA: {e}')
        return


def load_umap(ml_dir):
    umap_path = glob.glob(ml_dir + f'/*umap*.pkl')[0]
    try:
        with open(umap_path, 'rb') as fin:
            umap_transformer = pickle.load(fin)
        #umap_transformer = joblib.load(umap_path)
        return umap_transformer
    except Exception as e:
        print(f'Cannot load pickled UMAP: {e}')
        return

def load_plot_attr(att_dir):
    df_plot_path = glob.glob(att_dir + f'/*umap_attributes.csv')[0]
    try:
        df_plot_attr = pd.read_csv(df_plot_path, index_col=0)
        return df_plot_attr
    except Exception as e:
        print(f'Cannot load umap plot attributes: {e}')
        return

def transform_input(s_matrix, pca, umap):
    xsvd = pca.transform(s_matrix)
    embedding = umap.transform(xsvd)
    return(embedding)


def bokeh_gnomad(embedding, plot_attr, sample_name, outdir):
    # Hover labels
    source = ColumnDataSource(plot_attr)

    # Create a new plot with a title and axis labels
    p = figure(title="Gnomad UMAP", x_axis_label='umap 1', y_axis_label='umap 2', width=800, height=600)
    p.add_layout(Legend(), 'right')

    # Add the scatter plot
    for continent in plot_attr['Continent'].unique():
        sub_plot = CDSView(source=source, filters=[GroupFilter(column_name='Continent', group=continent)])
        p.circle(x='x', y='y',
            color='color_code',
            size=4,
            alpha=0.6, 
            legend_label=continent,
            source=source,
            view=sub_plot,
            name='reference_samples')
    
    # Add sample
    p.circle(x=embedding[0, 0], y=embedding[0, 1],
            color='black',
            size=6,
            alpha=1,
            name='sample')

    p.add_tools(HoverTool(names=['reference_samples'],
            tooltips=[('subcont', '@Subcontinent')]),
        HoverTool(names=['sample'],
            tooltips=[('User Input Sample Name', sample_name[0])]))
    p.legend.title = "Continental Labels"

    output_file(os.path.join(outdir, f'{sample_name[0]}_gnomAD_umap.html'))
    save(p)


def bokeh_1kgp(embedding, plot_attr, sample_name, outdir):
    # Hover labels
    source = ColumnDataSource(plot_attr)

    # Create a new plot with a title and axis labels
    p = figure(title="1kGP UMAP", x_axis_label='umap 1', y_axis_label='umap 2', width=800, height=600)
    p.add_layout(Legend(), 'right')

    # Add the scatter plot
    for continent in plot_attr['Continent'].unique():
        sub_plot = CDSView(source=source, filters=[GroupFilter(column_name='Continent', group=continent)])
        p.circle(x='x', y='y',
            color='color_code',
            size=4,
            alpha=0.6, 
            legend_label=continent,
            source=source,
            view=sub_plot,
            name='reference_samples')
    
    # Add sample
    p.circle(x=embedding[0, 0], y=embedding[0, 1],
            color='black',
            size=6,
            alpha=1,
            name='sample')

    p.add_tools(HoverTool(names=['reference_samples'],
            tooltips=[('population', '@{Population name}')]),
        HoverTool(names=['sample'],
            tooltips=[('User Input Sample Name', sample_name[0])]))
    p.legend.title = "Continental Labels"

    output_file(os.path.join(outdir, f'{sample_name[0]}_1kGP_umap.html'))
    save(p)


def bokeh_sgdp(embedding, plot_attr, sample_name, outdir):
    # Hover labels
    source = ColumnDataSource(plot_attr)

    # Create a new plot with a title and axis labels
    p = figure(title="SGDP UMAP", x_axis_label='umap 1', y_axis_label='umap 2', width=800, height=600)
    p.add_layout(Legend(), 'right')

    # Add the scatter plot
    for continent in plot_attr['Region'].unique():
        sub_plot = CDSView(source=source, filters=[GroupFilter(column_name='Region', group=continent)])
        p.circle(x='x', y='y',
            color='color_code',
            size=6,
            alpha=0.8, 
            legend_label=continent,
            source=source,
            view=sub_plot,
            name='reference_samples')
    
    # Add sample
    p.circle(x=embedding[0, 0], y=embedding[0, 1],
            color='black',
            size=6,
            alpha=1,
            name='sample')

    p.add_tools(HoverTool(names=['reference_samples'],
            tooltips=[('country', '@Country'),
                ('population', '@{Population ID}')]),
        HoverTool(names=['sample'],
            tooltips=[('User Input Sample Name', sample_name[0])]))
    p.legend.title = "Continental Labels"

    output_file(os.path.join(outdir, f'{sample_name[0]}_SGDP_umap.html'))
    save(p)




def plot_umap_parser(s_matrix, ml_dir, att_dir, sample_name, outdir, m_type):
    pca = load_pca(ml_dir)
    umap = load_umap(ml_dir)
    embedding = transform_input(s_matrix, pca, umap)
    plot_attr = load_plot_attr(att_dir)

    if m_type == 'gnomAD_continental':
        bokeh_gnomad(embedding, plot_attr, sample_name, outdir)
    elif m_type == '1kGP_continental':
        bokeh_1kgp(embedding, plot_attr, sample_name, outdir)
    elif m_type == 'SGDP_continental':
        bokeh_sgdp(embedding, plot_attr, sample_name, outdir)

