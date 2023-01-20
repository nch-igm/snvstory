import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval


'''
Plots all lolipop and donut plot results.
'''


r1 = 1.3  # Wedge plot radius
color_dict = {'purple':'#B847A3', 'yellow':'#FBDF6C', 'orange':'#ED592A', 'grey':'#646464', 'blue':'#2FA4DC', 'pink':'#FFC1C1', 'red':'#DC2E31', 'green':'#2EDB7E'}


# Normalize the sub continental model predictions using the continental probabilities
def get_score_normalizers(data, var):
    normalize_sub_continent = np.asarray(data['gnomAD_continental'])
    normalize_1kg_continent = np.asarray(data['1kGP_continental'])
    # Normalizing the subcontinent data
    data['gnomAD_eur'] = np.asarray(data['gnomAD_eur']) * normalize_sub_continent[4] # Indicies are the positions of the sub populations in the coninental array
    data['gnomAD_eas'] = np.asarray(data['gnomAD_eas']) * normalize_sub_continent[3]
    data['1kGP_amr'] = np.asarray(data['1kGP_amr']) * normalize_1kg_continent[2]
    data['1kGP_afr'] = np.asarray(data['1kGP_afr']) * normalize_1kg_continent[4]
    data['1kGP_eur'] = np.asarray(data['1kGP_eur']) * normalize_1kg_continent[0]
    data['1kGP_sas'] = np.asarray(data['1kGP_sas']) * normalize_1kg_continent[3]
    data['1kGP_eas'] = np.asarray(data['1kGP_eas']) * normalize_1kg_continent[1]
    # Get the label associated with the max prob and generate an independent col
    top_hits = []
    for c in data.index:
        top_i = np.argsort(data[c])[-2:]
        top_val = np.asarray(data[c])[top_i]
        top_labels = [var.LABS_CONVERTER[c][x][0] for x in top_i]
        lab_val = (top_labels, top_val.tolist())
        top_hits.append(lab_val)
    return data, top_hits


def write_out_normed_df(data_list, t_list):
    normed_df = pd.concat(data_list, axis=1).T
    normed_df['total'] = t_list
    return normed_df


def transform_df(p):
    df = pd.read_csv(p, index_col=0)
    df = df.applymap(literal_eval)
    df = df.applymap(np.asarray)


# Plot all three models in separate donut plots
def get_values_for_donut_separate(data, model_lab, var, ax):
    labs = [x[0] for x in list(var.LABS_CONVERTER[model_lab].values())]
    lab_zip = dict(zip(labs, data[model_lab]))
    df = pd.DataFrame(data=lab_zip.values(), columns=[model_lab], index=lab_zip.keys())
    if model_lab == '1kGP_continental':
        df['colors'] = [color_dict['blue'], color_dict['orange'], color_dict['yellow'], color_dict['red'], color_dict['purple']]
    elif model_lab == 'gnomAD_continental':
        df['colors'] = [color_dict['purple'], color_dict['yellow'], color_dict['green'], color_dict['orange'], color_dict['blue'], color_dict['red']]
    else:
        df['colors'] = [color_dict['purple'], color_dict['yellow'], color_dict['orange'], color_dict['grey'], color_dict['blue'], color_dict['pink'], color_dict['red']]
    wedges, _ = ax.pie(df[model_lab], radius=r1, colors=df['colors'], wedgeprops=dict(width=0.5), startangle=-40)
    kw = dict(arrowprops=dict(arrowstyle="-"), zorder=1, va="center")
    for i, p in enumerate(wedges):
        if df.iloc[i][model_lab] >= 0.05:
            ang = (p.theta2 - p.theta1) / 2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))
            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = "angle,angleA=0,angleB={}".format(ang)
            kw["arrowprops"].update({"connectionstyle": connectionstyle})
            ax.annotate(df.index[i], xy=(x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
                        horizontalalignment=horizontalalignment, **kw)


def hex2rbg(hex_color):
    h = matplotlib.colors.to_rgba(hex_color)
    return h


def shaded_area(vals, alt_order):
    x_above = np.where(vals >= 0.05)[0]
    vals = vals[x_above]
    if x_above.size != 0:
        if alt_order.size != 0:
            x_above = alt_order[x_above]
        shade_coords = []
        for x in x_above:
            y0 = x - .5
            y1 = x + .5
            shade_coords.append((y0, y1))
        return shade_coords, vals, x_above
    else:
        return None, None, None


def plot_predictions_separate(data, var, sample_name: str, outdir: str):
    fig = plt.figure(figsize=(11,15))
    gs = plt.GridSpec(7,6)
    with plt.style.context('classic'):
        # Add lolipop plots
        for plot_ind, idx in enumerate(data.index):
            if plot_ind%2 == 0:
                ax = fig.add_subplot(gs[plot_ind//2, 0:3])
            else:
                ax = fig.add_subplot(gs[plot_ind//2, 3:])

            vals = np.asarray(data.loc[idx])
            anc_order = np.asarray([x[1] for x in list(var.LABS_CONVERTER[idx].values())])
            for i, x in enumerate(vals):
                ax.plot((0, x), (anc_order[i], anc_order[i]), color='grey', linewidth=2, zorder=1)

            labels = [var.ABBR[x[0]][0] for x in list(var.LABS_CONVERTER[idx].values())]
            colors = [hex2rbg(var.ABBR[x[0]][1]) for x in list(var.LABS_CONVERTER[idx].values())]
            ax.scatter(vals, anc_order, c=colors, linewidths=0.5, edgecolors='k', s=100, zorder=10)
            shaded_coords, top_v, top_y = shaded_area(vals, anc_order)
            if shaded_coords is not None:
                for e, shd in enumerate(shaded_coords):
                    ax.axhspan(shd[0], shd[1], facecolor='grey', alpha=0.25, zorder=-1)
                    ax.annotate(s=f"{top_v[e]:.2f}", xy=(top_v[e] + 0.05, shd[0] + 0.25))
                    
            ax.set_yticks(ticks=anc_order)
            ax.set_yticklabels(labels=labels, rotation=0, fontsize=11)
            ax.set_title(var.TITLES[idx])
            ax.set_xlim(-0.05, 1.25, auto=True)

        # Add donut plots
        ax = fig.add_subplot(gs[5, 0:2])
        get_values_for_donut_separate(data, 'gnomAD_continental', var, ax)
        ax.set_title('gnomAD\n')
        ax = fig.add_subplot(gs[5, 2:4])
        get_values_for_donut_separate(data, '1kGP_continental', var, ax)
        ax.set_title('1kGP\n')
        ax = fig.add_subplot(gs[5, 4:6])
        get_values_for_donut_separate(data, 'SGDP_continental', var, ax)
        ax.set_title('SGDP\n')
        
        fig.suptitle(sample_name, fontsize=16)
        fig.patch.set_facecolor('#EDEDED')
        plt.subplots_adjust(left=None, bottom=0, right=None, top=1, wspace=None, hspace=None)
        gs.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
        
    plt.savefig(os.path.join(outdir, f'{sample_name}.pdf'), bbox_inches='tight')


def plot_parser(df, var, outdir, ofn):
    data_list = []
    t_list = []
    for row in df.iterrows():
        sample_name = row[0]
        data = row[1]
        data, top_hits = get_score_normalizers(data, var)
        data_list.append(data)
        t_list.append(top_hits)
        plot_predictions_separate(data, var, sample_name, outdir)
    normed_df = write_out_normed_df(data_list, t_list)
    normed_df.to_csv(os.path.join(outdir, ofn))
