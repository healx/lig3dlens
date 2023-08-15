# Ligand-based 3D Virtual Screening
## Lig3DLens

Lig3DLens performs the following tasks:
1. Prepares a commercial compound library to be used for a VS campaign
2. Generated conformers for the compounds in the commercial library and calculated the 3D similarity (shape & electrostatics) of these compounds against a reference compound.
3. Finally, it can cluster the highest scoring hits and select a set number of representative compounds that can be ordered and tested.


## Installation

```
conda env create -f environment.yml
conda activate lig3dlens
```

The above are also executable via the makefile command:
```
make create-env
```

## Running a ligand-based 3D VS campaign

1. Prepare a library of compounds for virtual screening
```
python prep_cmpds_library.py --in input_SD_file --filter physchem_yaml_file --out output_SD_file
```

2. Generated 3D conformers for both the library and refernce compounds and score the library compounds 
using a 3D shape & electrostatics similarity function
```
python main.py --ref input_reference_molecule_file --lib input_library_file_name --conf num_conformers --out output_SD_file
```

3. Cluster the highest scoring molecules and select a representative set
```
python kmeans_clustering.py –in input_SD_file –clusters num_clusters –out output_file –dim fingerprint_dimension –fp_type fingerprint_type
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

## Future improvements
- Compound library preparation: 
    * apply a set of structural filters (for example [REOS](https://www.nature.com/articles/nrd1063) or [PAINS](https://pubs.acs.org/doi/10.1021/jm901137j)) - either remove or flag compounds.
    * Provide more autonomy to drug designer when setting the physicochemical properties filters.

- Compound selection: 
    * Multi-parameter selection of compounds using a score function that includes the 3D score, 2D similarity to the reference compound, and the physchem properties. The aim is to get an even distribution between highly scored cmpds and other properties.
    * Select an optimal number of clusters instead of a predefined one (e.g. using Silhouette). Alternatively, using another method for maximum score-diversity selection problem (e.g. Score Erosion algorithm).


