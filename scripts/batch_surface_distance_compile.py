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

    frames = []  # list to store all df per lesion.

    for subdir, dirs, files in os.walk(args["rootdir"]):
        for file in sorted(files):
            if file.endswith('.xlsx'):
                # check file extension is xlsx
                excel_input_file_per_lesion = os.path.join(subdir, file)
                df_single_lesion = pd.read_excel(excel_input_file_per_lesion, sheet_name='percentages_coverage')
                df_single_lesion.rename(columns={'Lesion': 'Lesion_ID', 'Patient': 'Patient_ID'}, inplace=True)
                frames.append(df_single_lesion)


result = pd.concat(frames, ignore_index=True)
timestr = time.strftime("%H%M%S-%Y%m%d")
filepath_excel = 'Surface_Distances_Pooled_' + timestr + '_.xlsx'
writer = pd.ExcelWriter(filepath_excel)
result.to_excel(writer, sheet_name='radiomics', index=False, float_format='%.4f')
writer.save()
