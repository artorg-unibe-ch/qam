# -*- coding: utf-8 -*-
"""
Created on Mon July 20 2020

@author: Iwan Paolucci
"""

import argparse

import numpy as np
import pandas as pd
import nibabel as nib
import scripts.plot_ablation_margin_hist as pm
from surface_distances.margin import compute_distances
from utils.niftireader import load_image

np.set_printoptions(suppress=True, precision=4)


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
        ap.add_argument("-i", "--lesion-id", required=True, help="lesion id")
        ap.add_argument("-p", "--patient-id", required=True, help="patient id from study")
        ap.add_argument("-o", "--OUTPUT", required=False, help="output file (csv)")
        ap.add_argument("-r", "--DIR", required=False, help="rootdir to write plots and csv")
        args = vars(ap.parse_args())
    return args


if __name__ == '__main__':
    args = get_args()
    tumor_file = args['tumor']
    ablation_file = args['ablation']
    liver_file = args['liver']
    patient_id = args['patient_id']
    lesion_id = args['lesion_id']
    output_file = args['OUTPUT']
    rootdir = args['DIR']

    tumor, tumor_np = load_image(tumor_file)
    ablation, ablation_np = load_image(ablation_file)
    # todo: add liver flag
    liver, liver_np = load_image(liver_file)

    # check whether there is a liver segmentation or not.
    has_liver_segmented = np.sum(liver_np.astype(np.uint8)) > 0

    pixdim = liver.header['pixdim']
    spacing = (pixdim[1], pixdim[2], pixdim[3])
    surface_distance = compute_distances(tumor_np, ablation_np, exclusion_zone=liver_np if has_liver_segmented else None,
                                         spacing_mm=spacing, connectivity=1, crop=False)
                                         
    # save the contour of the tumor and the distance map of the ablation for visualization purposes
    borders_gt = nib.Nifti1Image(surface_distance["borders_gt"].astype(np.uint8), tumor.affine)
    nib.save(borders_gt, 'borders_gt.nii.gz')
    distmap_pred = nib.Nifti1Image(surface_distance["distmap_pred"], ablation.affine)
    nib.save(distmap_pred, 'distmap_pred.nii.gz')

    #TODO: check that the surface distances are not empty
    pm.plot_histogram_surface_distances(patient_id, lesion_id, rootdir, distance_map=surface_distance['distances_gt_to_pred'],
                                        num_voxels=len(surface_distance['distances_gt_to_pred']),
                                        title='Quantitative Ablation Margin', ablation_date="20200101")


    df = pd.DataFrame(data={
        'Patient': [patient_id] * len(surface_distance['distances_gt_to_pred']),
        'Lesion': [lesion_id] * len(surface_distance['distances_gt_to_pred']),
        'Distances': surface_distance['distances_gt_to_pred']})
    df.to_csv(output_file, index=False)
