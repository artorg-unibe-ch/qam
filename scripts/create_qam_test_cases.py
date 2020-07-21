# -*- coding: utf-8 -*-
"""
@author: Raluca Sandu
"""

import os

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

from utils.niftireader import load_image

TEST_DATA_PATH = 'test_data'


def save_test_case(data_tumor, data_ablation, data_liver, path, sample, name):
    """

    :param data_tumor:
    :param data_ablation:
    :param data_liver:
    :param path:
    :param sample:
    :param name:
    :return:
    """

    if not os.path.exists(os.path.join(path, sample)):
        os.mkdir(os.path.join(path, sample))
    if not os.path.exists(os.path.join(path, sample, name)):
        os.mkdir(os.path.join(path, sample, name))
    # create_paths
    data_tumor_path = os.path.join(path, sample, name, '{0}_L{1}_Tumor.nii.gz'.format(sample, name))
    data_ablation_path = os.path.join(path, sample, name, '{0}_L{1}_Ablation.nii.gz'.format(sample, name))
    data_liver_path = os.path.join(path, sample, name, '{0}_L{1}_Liver.nii.gz'.format(sample, name))

    affine = np.eye(4)
    tumor_ni = nib.Nifti1Image(data_tumor.astype(np.uint8), affine=affine)
    ablation_ni = nib.Nifti1Image(data_ablation.astype(np.uint8), affine=affine)
    liver_ni = nib.Nifti1Image(data_liver.astype(np.uint8), affine=affine)

    nib.save(tumor_ni, data_tumor_path)
    nib.save(ablation_ni, data_ablation_path)
    nib.save(liver_ni, data_liver_path)


def create_sphere(shape, radius, position):
    """

    :param shape:
    :param radius:
    :param position:
    :return: array

    Examples:
            # >>> arr = create_sphere((256, 256, 256), 10, (127, 127, 127))

    """
    # assume shape and position are both a 3-tuple of int or float
    # the units are pixels / voxels (px for short)
    # radius is a int or float in px
    semisizes = (radius,) * 3
    grid = [slice(-x0, dim - x0) for x0, dim in zip(position, shape)]
    position = np.ogrid[grid]
    # calculate the distance of all points from `position` center
    # scaled by the radius
    arr = np.zeros(shape, dtype=float)
    for x_i, semisize in zip(position, semisizes):
        # this can be generalized for exponent != 2
        # in which case `(x_i / semisize)`
        # would become `np.abs(x_i / semisize)`
        arr += (x_i / semisize) ** 2

    # the inner part of the create_sphere will have distance below 1
    return arr <= 1.0


if __name__ == '__main__':

    if not os.path.exists(TEST_DATA_PATH):
        os.mkdir(TEST_DATA_PATH)

    #%% Non-Subcapsular Cases
    # perfect overlap , radius_ablation=radius_tumor=10)
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    # replace everything with False for liver
    data_liver[data_liver == True] = False
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T01', 'perfect_overlap_0mm_margin')

    # centered, radius_ablation=10 mm, radius_tumor=5mm, QAM=5
    data_tumor = create_sphere((256, 256, 256), 5, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    # replace everything with False for liver
    data_liver[data_liver == True] = False
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T01', 'perfect_overlap_5mm_margin')

    # centered, radius_ablation=5mm, radius_tumor=10, QAM=-5
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 5, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    # replace everything with False for liver
    data_liver[data_liver == True] = False
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T01', 'perfect_overlap_-5mm_margin')

    # centered, radius_ablation=15mm, radius_tumor=5, QAM=10
    data_tumor = create_sphere((256, 256, 256), 5, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 15, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 15, (0, 0, 0))
    # replace everything with False for liver
    data_liver[data_liver == True] = False
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T01', 'perfect_overlap_10mm_margin')


    # no overlap, same radius
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 10, (15, 15, 15))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    # replace everything with False for liver
    data_liver[data_liver == True] = False
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T01', 'perfect_no_overlap_margin')

    # shifted xy 5mm

    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 15, (5, 5, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    # replace everything with False for liver
    data_liver[data_liver == True] = False
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T01', 'perfect_shifted_5mm_xy_margin')

    # shifted z 5mm
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 15, (0, 0, 5))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    # replace everything with False for liver
    data_liver[data_liver == True] = False
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T01', 'perfect_shifted_5mm_z_margin')

    #%% Subcapsular Cases
    # perfect overlap , radius_ablation=radius_tumor=10)
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T02', 'perfect_overlap_0mm_margin')

    # centered, radius_ablation=10 mm, radius_tumor=5mm, QAM=5
    data_tumor = create_sphere((256, 256, 256), 5, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T02', 'perfect_overlap_5mm_margin')

    # centered, radius_ablation=5mm, radius_tumor=10, QAM=-5
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 5, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T02', 'perfect_overlap_-5mm_margin')

    # centered, radius_ablation=15mm, radius_tumor=5, QAM=10
    data_tumor = create_sphere((256, 256, 256), 5, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 15, (0, 0, 0))
    data_liver = create_sphere((256, 256, 256), 15, (0, 0, 0))
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T02', 'perfect_overlap_10mm_margin')

    # no overlap, same radius
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 10, (15, 15, 15))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T02', 'perfect_no_overlap_margin')

    # shifted xy 5mm

    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 15, (5, 5, 0))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T02', 'perfect_shifted_5mm_xy_margin')

    # shifted z 5mm
    data_tumor = create_sphere((256, 256, 256), 10, (0, 0, 0))
    data_ablation = create_sphere((256, 256, 256), 15, (0, 0, 5))
    data_liver = create_sphere((256, 256, 256), 10, (0, 0, 0))
    save_test_case(data_tumor, data_ablation, data_liver, TEST_DATA_PATH, 'T02', 'perfect_shifted_5mm_z_margin')