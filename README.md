# Ligand-based 3D Virtual Screening
## Lig3DLens
Bill Tatsis, Matt Seddon, Dan Mason, Dan O'Donovan, Gio Cincilla, Azedine Zoufir and Nath Brown

Lig3DLens performs the following tasks:
1. Prepares a commercial compound library to be used for a VS campaign. This task
involves: i) compound standardisation and ii) optional filtering of compounds outside a predefined range of physicochemical properties.
2. Generates conformers for all of the compounds in the commercial library and calculates their 3D similarity (shape & electrostatics) to a reference compound.
3. Finally, it can cluster the highest scoring hits and select a set number of representative compounds that can be ordered and tested.


## Installation

```
python -m pip install -r requirements.txt .
```
... for development
```
python -m pip install -r requirements.txt -r dev-requirements -e .
```

To run the test suite we use [`tox`](https://tox.wiki/en/4.24.1/) which can be
called with
```
tox
```

## Running a ligand-based 3D VS campaign

1. Prepare a chemical library for a 3D VS campaign
> **Note**
> In order to keep track of the library cmpds the input file should have a column containing the text "ID"

**With physicochemical filtering (default):**
```
lig3dlens-prepare --in input_SD_file --filter physchem_yaml_file --out output_SD_file
```

**Without filtering (standardization only):**
```
lig3dlens-prepare --in input_SD_file --no-filter --out output_SD_file
```
Use `--no-filter` when filtering has already been performed or when no filtering is desired. This significantly improves performance by skipping descriptor calculations.

2. Generates 3D conformers for both the library and reference compounds and scores the library compounds using a 3D shape & electrostatics similarity function to the reference molecule
> **Note**
> In order to keep track of the library cmpds the input file should have a column containing the text "ID"

```
lig3dlens-align --ref input_reference_molecule_file --lib input_library_file_name --conf num_conformers --out output_SD_file
```

> **New feature**
> Lig3DLens can work with input files that contain the 3D conformations for the reference compound and/or the compound library

```
lig3dlens-align --ref input_3D-reference_molecule_file --lib input_3D-library_file_name  --out output_SD_file
```

3. Clusters the highest scoring molecules and selects a representative (diverse) set of compounds. The user can input the number of clusters (`num_clusters`), the fingerprint type (`fingerprint_type`) and its dimension (`fingerprint_dimension`) used for the clustering.
```
lig3dlens-cluster –-in input_SD_file –-clusters num_clusters –-out output_file -–dim fingerprint_dimension -–fp_type fingerprint_type
```

### Example
To run `lig3dlens` with the input files (reference molecule and compound library) included in this repository;

... preparing the compound library for virtual screening (with filtering)
```shell
lig3dlens-prepare --in tests/test_data/Enamine_hts_collection_202303_first500_VS_results.sdf \
    --filter lig3dlens/physchem_properties.yaml \
    --out curated_compounds.sdf
```

... or without filtering (standardization only, faster)
```shell
lig3dlens-prepare --in tests/test_data/Enamine_hts_collection_202303_first500_VS_results.sdf \
    --no-filter \
    --out standardized_compounds.sdf
```

... generate the 3D conformers for both reference and library compounds, score the compounds in the library set using a 3D shape and electrostatics similarity
```shell
lig3dlens-align --ref input_reference_molecule_file --lib input_library_file_name --conf num_conformers --out output_SD_file
```

... cluster the results from the 3D virtual screening campaing
```shell
lig3dlens-cluster –-in input_SD_file –-clusters num_clusters –-out output_file -–dim fingerprint_dimension -–fp_type fingerprint_type
```


## Development

Use the Makefile commands to help tidy the codebase periodically. The following will reformat the code according to PEP8, and logically sort the imported modules:
```
make tidy
```

## Tests
Run pytest in lig3dlens directory
```
pytest tests
```

## Requirements

Note: The whole VS workflow was tested in a linux (ubuntu) environment and this environment variable had to be set:
Tell MKL (used by NumPy) to use the GNU OpenMP runtime instead of the Intel OpenMP runtime by setting the following environment variable:
```
export MKL_THREADING_LAYER=GNU
```
The open source quantum chemistry package *Psi4* is required to for the QM calculation of the partial charges.
More information about installing Psi4 in different CPU architectures (arm64 included) is provided in  [Psi4's](https://psicode.org/installs/v182/) website

## Future improvements
- Compound library preparation:
    * Apply a set of structural filters (for example [REOS](https://www.nature.com/articles/nrd1063) or [PAINS](https://pubs.acs.org/doi/10.1021/jm901137j)) - either remove or flag compounds.
    * Provide more autonomy to the drug designer when setting the physicochemical properties filters.

- Compound selection:
    * Multi-parameter selection of compounds using a score function that includes the 3D score, 2D similarity to the reference compound, and the physchem properties. The aim is to get an even distribution between highly scored cmpds and other properties.
    * Select an optimal number of clusters instead of a predefined one (e.g. using Silhouette or affinity propagation methods). Alternatively, using another method for maximum score-diversity selection problem (e.g. Score Erosion algorithm).
    * Provide tools to analyse the chemical diversity of the final selection compound set.
