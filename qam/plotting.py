# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 16:15:51 2017

@author: Raluca Sandu
"""
from collections import OrderedDict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

matplotlib.use('Agg')

np.seterr(divide='ignore', invalid='ignore')
cmap = sns.color_palette("colorblind")  # colorblind friendly palette


def plot_histogram_surface_distances(pat_name, lesion_id, output_file, distance_map, title, print_case_details=True):
    """

    Plots and saves the surface distances (traffic-light color schemes) between tumor and ablation.
    :param pat_name: Patient Name
    :param lesion_id: Lesion 1, 2, 3...etc
    :param output_file: Name of the img saved to the disk
    :param distance_map: Array containing all the surface distances computed.
    :param title: Title of the plot.
    :param print_case_details: True/False to plot an extended title on the image.
    :return: The percentages of non-ablated, insufficiently ablated and completely ablated tumor surface.
    """
    fontsize = 20

    fig, ax = plt.subplots(figsize=(12, 10))
    if len(distance_map) == 0:
        min_val = -15
        max_val = 15
    else:
        min_val = int(np.floor(min(distance_map)))
        max_val = int(np.ceil(max(distance_map)))

    bins = np.arange(min_val, max_val + 1.5, 1)
    col_height, bins, patches = ax.hist(distance_map, ec='black', align='mid', bins=bins)

    voxels_nonablated = []
    voxels_insuffablated = []
    voxels_ablated = []

    for b, p, col_val in zip(bins, patches, col_height):
        if b < 0:
            voxels_nonablated.append(col_val)
        elif 0 <= b < 5:
            voxels_insuffablated.append(col_val)
        elif b >= 5:
            voxels_ablated.append(col_val)

    num_voxels = np.sum(voxels_ablated) + np.sum(voxels_insuffablated) + np.sum(voxels_nonablated)

    # %% calculate the total percentage of surface for ablated, non-ablated, insufficiently ablated
    voxels_nonablated = np.asarray(voxels_nonablated)
    voxels_insuffablated = np.asarray(voxels_insuffablated)
    voxels_ablated = np.asarray(voxels_ablated)

    sum_perc_nonablated = ((voxels_nonablated / num_voxels) * 100).sum()
    sum_perc_insuffablated = ((voxels_insuffablated / num_voxels) * 100).sum()
    sum_perc_ablated = ((voxels_ablated / num_voxels) * 100).sum()
    # %% iterate through the bins to change the colors of the patches bases on the range [mm]
    for b, p, col_val in zip(bins, patches, col_height):
        if b < 0 and col_val > 0:
            plt.setp(p, 'facecolor', cmap[3],
                     label='Ablation Margin ' + r'$x < 0$' + 'mm :' + " %.2f" % sum_perc_nonablated + '%')
        elif 0 <= b < 5 and col_val > 0:
            plt.setp(p, 'facecolor', cmap[8],
                     label='Ablation Margin ' + r'$0 \leq x < 5$' + 'mm: ' + "%.2f" % sum_perc_insuffablated + '%')
        elif b >= 5 and col_val > 0:
            plt.setp(p, 'facecolor', cmap[2],
                     label='Ablation Margin ' + r'$x \geq 5$' + 'mm: ' + " %.2f" % sum_perc_ablated + '%')
    # %% edit the axes limits and labels
    plt.xlabel('Surface-to-Surface Exact Euclidean Distances (mm)', fontsize=fontsize, color='black')
    plt.tick_params(labelsize=fontsize, color='black')
    ax.tick_params(colors='black', labelsize=fontsize)
    ax.set_xlim([-15, 15])
    # edit the y-ticks: change to percentage of surface
    yticks, locs = plt.yticks()
    percent = (yticks / num_voxels) * 100
    percentage_surface_rounded = np.round(percent)
    yticks_percent = [str(x) + '%' for x in percentage_surface_rounded]
    new_yticks = (percentage_surface_rounded * yticks) / percent
    new_yticks[0] = 0
    plt.yticks(new_yticks, yticks_percent)

    plt.ylabel('Tumor Surface Covered (%)', fontsize=fontsize, color='black')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=fontsize, loc='upper left')
    plt.xticks(fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.grid(False)
    # %% save the fig to disk as png and eps
    if print_case_details:
        plt.title(title + '. Case ' + str(pat_name) + '. Lesion ' + str(lesion_id), fontsize=fontsize)
    else:
        plt.title(title, fontsize=fontsize)

    ax.set_rasterized(True)
    plt.savefig(output_file + '.png', dpi=600, bbox_inches='tight')
    plt.savefig(output_file + '.svg', dpi=600)
    plt.savefig(output_file + '.eps', dpi=600)

    plt.close()

    return sum_perc_nonablated, sum_perc_insuffablated, sum_perc_ablated
