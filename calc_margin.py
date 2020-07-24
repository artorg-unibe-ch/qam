# -*- coding: utf-8 -*-
"""
Created on Mon July 20 2020

@author: Raluca Sandu, Iwan Paolucci
"""

import argparse
import os
import sys
from datetime import date

import numpy as np
import pandas as pd

import scripts.plot_ablation_margin_hist as pm
from surface_distances.margin import compute_distances
from utils.niftireader import load_image

np.set_printoptions(suppress=True, precision=4)
today = date.today()


def is_running_in_snakemake():
    try:
        snakemake
        return True
    except NameError:
        return False


def get_args():
    if is_running_in_snakemake():
        # noinspection PyUnresolvedReferences
        args = {
            'tumor': snakemake.input['tumor'],
            'ablation': snakemake.input['ablation'],
            'liver': snakemake.input['liver'],
            'patient_id': snakemake.params['patient_id'],
            'lesion_id': snakemake.params['lesion_id'],
            'OUTPUT': snakemake.output[0]
        }
    else:
        ap = argparse.ArgumentParser()
        ap.add_argument("-t", "--tumor", required=True, help="path to the tumor segmentation")
        ap.add_argument("-a", "--ablation", required=True, help="path to the ablation segmentation")
        ap.add_argument("-l", "--liver", required=False, help="path to the liver segmentation")
        ap.add_argument("-p", "--patient-id", required=False, help="patient id from study")
        ap.add_argument("-i", "--lesion-id", required=False, help="lesion id")
        ap.add_argument("-d", "--ablation-date", required=False, help="ablation date from study")
        ap.add_argument("-o", "--OUTPUT", required=False, help="output filename (csv)")
        ap.add_argument("-r", "--DIR", required=False, help="directory path to write plots and csv")
        args = vars(ap.parse_args())
    return args


if __name__ == '__main__':
    # get the args
    args = get_args()
    tumor_file = args['tumor']
    ablation_file = args['ablation']
    liver_file = args['liver']
    patient_id = args['patient_id']
    lesion_id = args['lesion_id']
    ablation_date = args['ablation_date']
    output_file = args['OUTPUT']
    dir_path = args['DIR']

    # check whether the input has been provided
    if patient_id is None:
        patient_id = 'Test'
    if lesion_id is None:
        lesion_id = 1
    if ablation_date is None:
        ablation_date = today.strftime("%d-%m-%Y")
    if output_file is None:
        csv_dir = 'surface_distances_csv'
        if not os.path.exists(csv_dir):
            os.mkdir(csv_dir)
        output_file = os.path.join(csv_dir, str(patient_id) + '_' + str(lesion_id) + '_surface_distances.xlsx')
    if dir_path is None:
        dir_path = 'figures'
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)

    tumor, tumor_np = load_image(tumor_file)
    # check if there is actually a segmentation in the file
    has_tumor_segmented = np.sum(tumor_np.astype(np.uint8)) > 0
    if has_tumor_segmented is False:
        print('No tumor segmentation mask found in the file provided...program exiting')
        sys.exit()
    ablation, ablation_np = load_image(ablation_file)
    has_ablation_segmented = np.sum(ablation_np.astype(np.uint8)) > 0

    if has_ablation_segmented is False:
        print('No ablation segmentation mask found in the file provided...program exiting')
        sys.exit()

    if liver_file is not None:
        # load the image file
        liver, liver_np = load_image(liver_file)
        # check whether there is a liver segmentation or not.
        has_liver_segmented = np.sum(liver_np.astype(np.uint8)) > 0
    else:
        has_liver_segmented = False

    # extract the spacing from the ablation file
    pixdim = ablation.header['pixdim']
    spacing = (pixdim[1], pixdim[2], pixdim[3])
    # compute the surface distances based on tumor and ablation segmentations
    surface_distance = compute_distances(mask_gt=tumor_np, mask_pred=ablation_np,
                                         exclusion_zone=liver_np if has_liver_segmented else None,
                                         spacing_mm=spacing, connectivity=1, crop=True)

    if surface_distance['distances_gt_to_pred'].size > 0:
        # if surface distances returned are not empty
        non_ablated, insuffiecient_ablated, completely_ablated =\
            pm.plot_histogram_surface_distances(patient_id, lesion_id, dir_path,
                                                distance_map=surface_distance['distances_gt_to_pred'],
                                                num_voxels=len(surface_distance['distances_gt_to_pred']),
                                                title='Quantitative Ablation Margin',
                                                ablation_date=ablation_date)
        df = pd.DataFrame(data={
            'Patient': [patient_id] * len(surface_distance['distances_gt_to_pred']),
            'Lesion': [lesion_id] * len(surface_distance['distances_gt_to_pred']),
            'Distances': surface_distance['distances_gt_to_pred']})

        writer = pd.ExcelWriter(output_file)
        df.to_excel(writer, sheet_name='surface_distances', index=False, float_format='%.4f')
        max_distance = np.nanmax(surface_distance['distances_gt_to_pred'])
        min_distance = np.nanmin(surface_distance['distances_gt_to_pred'])
        q25_distance = np.nanquantile(surface_distance['distances_gt_to_pred'], 0.25)
        median_distance = np.nanmedian(surface_distance['distances_gt_to_pred'])
        q75_distance = np.nanquantile(surface_distance['distances_gt_to_pred'], 0.75)

        patient_data = {'Patient': patient_id,
                        'Lesion': lesion_id,
                        'max_distance': max_distance,
                        'min_distance': min_distance,
                        'q25_distance': q25_distance,
                        'median_distance': median_distance,
                        'q75_distance': q75_distance,
                        'x_less_than_0mm': non_ablated,
                        'x_equal_greater_than_0m': insuffiecient_ablated,
                        'x_equal_greater_than_5m': completely_ablated}

        df_percentages_coverage = pd.DataFrame([patient_data])
        df_percentages_coverage.to_excel(writer, sheet_name='percentages_coverage', index=False, float_format='%.4f')
        writer.save()

    else:
        print('No surface distance computed for patient ' + str(patient_id) + ' lesion ' + str(
            lesion_id) + '. Lesion could be completely within the subcapsular exclusion zone. '
                         'Please visualize your segmentation files for more insight.')
