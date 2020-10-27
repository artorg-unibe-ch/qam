import sys
sys.path.insert(0, "..")
import os
import unittest
from qam import margin
from utils import niftireader

OUTPUT_FILE = 'data/_output/Grouped.csv'

def _get_file_name(case_id, lesion_id, type):
    return os.path.join('test_data', case_id, lesion_id,
        '{0}_L{1}_{2}.nii.gz'.format(case_id, lesion_id, type))


class TestMargins(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)

    def test_01_no_overlap_10mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '01_no_overlap_10mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '01_no_overlap_10mm_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertEqual(record["min_distance"], -20.0)
        self.assertEqual(record["max_distance"], -10.0)

    def test_02_perfect_overlap_0mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '02_perfect_overlap_0mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '02_perfect_overlap_0mm_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]

        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertAlmostEqual(record["min_distance"], 0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 0, delta=0.01)

    def test_03_perfect_overlap_10mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '03_perfect_overlap_10mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '03_perfect_overlap_10mm_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertAlmostEqual(record["min_distance"], 9.05, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 10.0, delta=0.01)

    def test_04_perfect_overlap_5mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '04_perfect_overlap_5mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '04_perfect_overlap_5mm_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertAlmostEqual(record["min_distance"], 4.12, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 5.0, delta=0.01)

    def test_05_perfect_overlap_n5mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '05_perfect_overlap_-5mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '05_perfect_overlap_-5mm_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 978.0)
        self.assertAlmostEqual(record["min_distance"], -5.19, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], -4.12, delta=0.01)

    def test_06_perfect_shifted_5mm_xy_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '06_perfect_shifted_5mm_xy_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '06_perfect_shifted_5mm_xy_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertAlmostEqual(record["min_distance"], -2.23, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 7.0, delta=0.01)

    def test_07_perfect_shifted_5mm_x_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '07_perfect_shifted_5mm_x_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '07_perfect_shifted_5mm_x_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertAlmostEqual(record["min_distance"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 9.05, delta=0.01)

    def test_08_perfect_shifted_5mm_y_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '08_perfect_shifted_5mm_y_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '08_perfect_shifted_5mm_y_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertAlmostEqual(record["min_distance"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 9.05, delta=0.01)

    def test_09_perfect_shifted_5mm_z_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '09_perfect_shifted_5mm_z_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '09_perfect_shifted_5mm_z_margin', 'Ablation'))
        distances = margin.compute_distances(tumor_np, ablation_np, None, 1)
        df = margin.summarize('T01', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T01") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 222.0)
        self.assertAlmostEqual(record["min_distance"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 9.05, delta=0.01)


class TestSubcapsularLesions(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)

    def test_01_0mm_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '01_0mm_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '01_0mm_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '01_0mm_margin_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L01', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L01')].iloc[0]
        self.assertEqual(record["nr_voxels"], 97.0)
        self.assertAlmostEqual(record["min_distance"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 0.0, delta=0.01)

    def test_02_5mm_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '02_5mm_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '02_5mm_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '02_5mm_margin_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L02', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L02')].iloc[0]
        self.assertEqual(record["nr_voxels"], 97.0)
        self.assertAlmostEqual(record["min_distance"], 2.23, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 5.0, delta=0.01)

    def test_03_n5mm_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '03_-5mm_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '03_-5mm_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '03_-5mm_margin_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L03', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L03')].iloc[0]
        self.assertEqual(record["nr_voxels"], 317.0)
        self.assertAlmostEqual(record["min_distance"], -5.09, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], -2.23, delta=0.01)

    def test_04_2mm_shifted_tumor_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '04_2mm_shifted_tumor_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '04_2mm_shifted_tumor_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '04_2mm_shifted_tumor_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L04', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L04')].iloc[0]
        self.assertEqual(record["nr_voxels"], 149.0)
        self.assertAlmostEqual(record["min_distance"], 1.41, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 5.0, delta=0.01)

    def test_05_4mm_shifted_tumor_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '05_4mm_shifted_tumor_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '05_4mm_shifted_tumor_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '05_4mm_shifted_tumor_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L05', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L05')].iloc[0]
        self.assertEqual(record["nr_voxels"], 193.0)
        self.assertAlmostEqual(record["min_distance"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 5.0, delta=0.01)

    def test_06_shifted_ablation_5mm_xy_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '06_shifted_ablation_5mm_xy_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '06_shifted_ablation_5mm_xy_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '06_shifted_ablation_5mm_xy_margin_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L06', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L06')].iloc[0]
        self.assertEqual(record["nr_voxels"], 97.0)
        self.assertAlmostEqual(record["min_distance"], -3.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 4.47, delta=0.01)

    def test_07_shifted_5mm_x_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '07_shifted_5mm_x_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '07_shifted_5mm_x_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '07_shifted_5mm_x_margin_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L07', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L07')].iloc[0]
        self.assertEqual(record["nr_voxels"], 97.0)
        self.assertAlmostEqual(record["min_distance"], -1.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 4.12, delta=0.01)

    def test_08_shifted_5mm_y_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '08_shifted_ablation_5mm_y_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '08_shifted_ablation_5mm_y_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '08_shifted_ablation_5mm_y_margin_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L08', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L08')].iloc[0]
        self.assertEqual(record["nr_voxels"], 97.0)
        self.assertAlmostEqual(record["min_distance"], -2.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 5.74, delta=0.01)

    def test_09_shifted_ablation_5mm_z_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '09_shifted_ablation_5mm_z_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '09_shifted_ablation_5mm_z_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '09_shifted_ablation_5mm_z_margin_subcapsular', 'Liver'))

        distances = margin.compute_distances(tumor_np, ablation_np, liver_np, 1)
        df = margin.summarize('T02', 'L09', distances)
        print(df)

        record = df.loc[(df['Patient'] == "T02") & (df['Lesion'] == 'L09')].iloc[0]
        self.assertEqual(record["nr_voxels"], 97.0)
        self.assertAlmostEqual(record["min_distance"], -2.0, delta=0.01)
        self.assertAlmostEqual(record["max_distance"], 5.74, delta=0.01)

if __name__ == '__main__':
    unittest.main()
