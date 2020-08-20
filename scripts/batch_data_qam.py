# -*- coding: utf-8 -*-
"""
@author: Raluca Sandu
"""
import argparse
import os
import time

import pandas as pd

if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--rootdir", required=True,
                    help="path to the patient folder with Radiomics XLSX to be processed")
    args = vars(ap.parse_args())

    if args["rootdir"] is not None:
        print("Path to folder with Radiomics Files for each patient and each lesion: ", args["rootdir"])

    result = []  # list to store all df per lesion.

    for subdir, dirs, files in os.walk(args["rootdir"]):
        ablation_path = None
        tumor_path = None
        liver_path = None
        patient_id = None
        lesion_id = None
        for file in sorted(files):
            if file.endswith('.gz'):
                patient_id = file[0:3]
                if "Ablation" in file and "Source" not in file:
                    ablation_path = os.path.join(subdir, file)
                    lesion_id = file.split("_")[1]
                if "Tumor" in file and "Source" not in file:
                    tumor_path = os.path.join(subdir, file)
                if "Liver" in file:
                    liver_path = os.path.join(subdir, file)
        if ablation_path is not None:
            frames = {'Patient_ID': patient_id,
                      'Lesion_ID': lesion_id,
                      'ablation_path': ablation_path,
                      'tumor_path': tumor_path,
                      'liver_path': liver_path
                      }
            result.append(frames)

file_paths = pd.DataFrame(result)
filepath_excel = 'Input_Clinical_Maverric_FilePaths.xlsx'
writer = pd.ExcelWriter(filepath_excel)
file_paths.to_excel(writer, sheet_name='radiomics', index=False, float_format='%.4f')
writer.save()
