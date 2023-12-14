# ----------------------------------------------------------------------------
# Authors: Bill Tatsis, Matt Seddon, Dan Mason, and Dan O'Donovan
# Company: Healx
# Date: August 17, 2023
#
# Description: This script standardises a chemical library and filters compounds
#              based on user defined physchem criteria
# ----------------------------------------------------------------------------

import gc

import click
import datamol as dm
import pandas as pd
import yaml
from loguru import logger
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def _preprocess(i: int, row: pd.Series, mol_column: str = "ROMol") -> pd.Series:
    """
    Standardise/"Sanitise" a chemical library

    Parameters
    ----------
    i: int
        used to iterate through the pandas df rows
    row: pd.Series
        row from pandas df
    mol_column: str
        pandas df column that holds the RDKit Mol objects

    Returns
    -------
    row
        a pd.series with the standardised smiles
    """

    with dm.without_rdkit_log():
        mol = dm.fix_mol(row[mol_column])
        mol = dm.sanitize_mol(mol, sanifix=True, charge_neutral=False)
        mol = dm.standardize_mol(
            mol,
            disconnect_metals=True,
            normalize=True,
            reionize=True,
            uncharge=False,
            stereo=True,
        )
        mol = dm.keep_largest_fragment(mol)

        row["smiles_sdt"] = dm.standardize_smiles(dm.to_smiles(mol))

        return row


def _calculate_descriptors(
    i: int, row: pd.Series, smiles_column: str = "smiles_sdt"
) -> pd.Series:
    """
    Calculate physicochemical properties of all the cmpds in a dataframe

    Parameters
    ----------
    i: int
        used to iterate throught pandas df rows
    row: pd.Series
        row from pandas df
    smiles_column: str
        pandas df SMILES column

    Returns
    -------
    row
        a pd.series with the standardised smiles
    """

    mol = Chem.MolFromSmiles(row[smiles_column])

    if mol is None:
        raise ValueError("Invalid SMILES:",  row[smiles_column])

    row["mw"] = dm.descriptors.mw(mol)
    row["hba"] = dm.descriptors.n_lipinski_hba(mol)
    row["hbd"] = dm.descriptors.n_lipinski_hbd(mol)
    row["rot_bonds"] = dm.descriptors.n_rotatable_bonds(mol)
    row["arom_rings"] = dm.descriptors.n_aromatic_rings(mol)
    row["fsp3"] = dm.descriptors.fsp3(mol)
    row["tpsa"] = dm.descriptors.tpsa(mol)
    row["clogp"] = dm.descriptors.clogp(mol)
    row["qed"] = dm.descriptors.qed(mol)
    row["heavy_atms"] = dm.descriptors.n_heavy_atoms(mol)
    row["stereo_cntrs"] = dm.descriptors.n_stereo_centers(mol)

    return row


@click.command(name="prepare")
@click.option(
    "--in", "input_cmpd_lib", type=str, required=True, help="Input compound library"
)
@click.option(
    "--filter",
    "input_physchem_props",
    type=str,
    required=True,
    help="yaml file with physicochemical properties filters",
)
@click.option(
    "--out",
    "output_file",
    type=str,
    required=True,
    help="Output SD file with Mols & ID columns",
)
def main(input_cmpd_lib, input_physchem_props, output_file):
    logger.info(
        "Initialising the compound library preparation workflow",
        colorize=True,
        format="<green>{time}</green> <level>{message}</level>",
    )

    logger.info(f"Loading cmpds from {input_cmpd_lib} to a Pandas dataframe")

    # Load input file to a pandas df
    mols_lib = dm.read_sdf(input_cmpd_lib, as_df=True, mol_column="ROMol")
    logger.debug(f"There are {mols_lib.shape[0]} cmpds in {input_cmpd_lib}")

    # -----------------------------------------------------------------
    # Standardising/Sanitising compound library                      #
    # -----------------------------------------------------------------
    logger.info("Standardising compound library")
    mols_lib_sdt = dm.parallelized(
        _preprocess,
        mols_lib.iterrows(),
        arg_type="args",
        progress=True,
        total=len(mols_lib),
        tqdm_kwargs={"desc": "Standardising compound library"},
    )

    mols_lib_sdt = pd.DataFrame([r for r in mols_lib_sdt if r is not None])

    # Drop the initial 'ROMol' column to reduce the allocated memory
    mols_lib_sdt.drop(columns=["ROMol"], inplace=True)
    mols_lib_sdt.reset_index(inplace=True, drop=True)

    logger.debug(f"There are {mols_lib_sdt.shape[0]} cmpds after standardization")

    # Free up memory by deleting the previous dataframe
    mols_lib = None
    gc.collect()

    logger.info(f"Parsing physchem filters from {input_physchem_props}")
    with open(input_physchem_props) as f:
        physchem_properties = yaml.safe_load(f)

    # -----------------------------------------------------------------
    # Calculating the physicochemical properties of the compounds
    # -----------------------------------------------------------------
    logger.info("Calculating physchem properties for the library cmpds")
    mols_lib_descs = dm.parallelized(
        _calculate_descriptors,
        mols_lib_sdt.iterrows(),
        arg_type="args",
        progress=True,
        total=len(mols_lib_sdt),
        tqdm_kwargs={"desc": "Calculating physchem properties"},
    )

    mols_lib_descs = pd.DataFrame(mols_lib_descs)

    # Free up memory by deleting the previous dataframe
    mols_lib_sdt = None
    gc.collect()

    # -----------------------------------------------------------------
    # Filtering cmpds using physicochemical properties criteria
    # -----------------------------------------------------------------

    mols_lib_filt = mols_lib_descs[
        (
            mols_lib_descs["mw"].between(
                physchem_properties["MIN_MOLWT"], physchem_properties["MAX_MOLWT"]
            )
        )
        & (
            mols_lib_descs["clogp"].between(
                physchem_properties["MIN_CLOGP"], physchem_properties["MAX_CLOGP"]
            )
        )
        & (
            mols_lib_descs["hbd"].between(
                physchem_properties["MIN_HBD"], physchem_properties["MAX_HBD"]
            )
        )
        & (
            mols_lib_descs["hba"].between(
                physchem_properties["MIN_HBA"], physchem_properties["MAX_HBA"]
            )
        )
        & (
            mols_lib_descs["tpsa"].between(
                physchem_properties["MIN_TPSA"], physchem_properties["MAX_TPSA"]
            )
        )
        & (
            mols_lib_descs["rot_bonds"].between(
                physchem_properties["MIN_ROT_BONDS"],
                physchem_properties["MAX_ROT_BONDS"],
            )
        )
        & (
            mols_lib_descs["arom_rings"].between(
                physchem_properties["MIN_AROM_RINGS"],
                physchem_properties["MAX_AROM_RINGS"],
            )
        )
        & (
            mols_lib_descs["stereo_cntrs"].between(
                physchem_properties["MIN_STEREO_CNTRS"],
                physchem_properties["MAX_STEREO_CNTRS"],
            )
        )
    ].copy()

    logger.debug(
        f"There are {mols_lib_filt.shape[0]} cmpds after applying physchem filters"
    )

    # Adding now the ROMol objects
    mols_lib_filt["ROMol"] = mols_lib_filt.smiles_sdt.apply(Chem.MolFromSmiles)

    # Get columns in dataframe that might contain IDs
    columns_with_id = mols_lib_filt.columns[
        mols_lib_filt.columns.str.contains("ID", case=False)
    ].to_list()
    columns_to_store = columns_with_id + ["smiles_sdt", "ROMol"]

    # Write the pandas df to a SD file
    dm.to_sdf(mols_lib_filt[columns_to_store], output_file)
    logger.success(f"Curated compound library stored in {output_file}")


if __name__ == "__main__":
    main()
