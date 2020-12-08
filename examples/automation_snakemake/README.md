# Automation using Snakemake

This is a portion of a Snakefile that we used to automate QAM computations for a larger clinical study.

    snakemake --profile profiles/local --directory ../../tests/test_data _output/Aggregated.xlsx

The folder structure has to be as follows (same as in the test data used in the repository):

```
+-- PatientID
|   +-- LesionNr
    |   +-- [PatientID]_L[LesionNR]_Ablation.nii.gz
    |   +-- [PatientID]_L[LesionNR]_Liver.nii.gz
    |   +-- [PatientID]_L[LesionNR]_Tumor.nii.gz
```

For more information about Snakemake refer to their website https://snakemake.readthedocs.io.
