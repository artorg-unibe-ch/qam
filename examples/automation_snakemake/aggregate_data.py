import pandas as pd

if __name__ == '__main__':
    input_files = snakemake.input
    output_file = snakemake.output[0]
    print('input files', input_files)
    df = pd.DataFrame(columns=['Patient', 'Lesion', 'max_distance', 'min_distance',
        'q25_distance', 'median_distance', 'q75_distance', 'x_less_than_0mm',
        'x_equal_greater_than_0m', 'x_equal_greater_than_5m'])

    data_frames = []

    for input_file in input_files:
        df = pd.read_excel(input_file, sheet_name='percentages_coverage')
        data_frames.append(df)

    data_aggregated = pd.concat(data_frames)
    data_aggregated.to_excel(output_file, sheet_name='percentages_coverage', index=False)
