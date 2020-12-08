# -*- coding: utf-8 -*-
"""
Created on Mon July 20 2020

@author: Iwan Paolucci
"""


import numpy as np
import nibabel as nib


def image_to_np(image):
    image_np = np.asanyarray(image.dataobj)
    if np.min(image_np) != np.max(image_np):
        image_np = image_np == np.max(image_np)
    return image_np


def load_image(file):
    image = nib.load(file)
    image = nib.as_closest_canonical(image)
    image_np = image_to_np(image)

    return image, image_np
