import multiprocessing
from pathlib import Path
from typing import Optional, Dict, Generator, Tuple

from loguru import logger
from rdkit import Chem
from tqdm import tqdm

from lig3dlens import parsers
from lig3dlens.align3D_score import score_alignment
from lig3dlens.gen_conformers import generate_conformers, map_conf_gen


def save_best_mol(
    output_file: Chem.SDWriter, best_mol: Chem.Mol, mol_id: str, alignment_score
) -> None:
    """
    Save the best molecule conformer with scores to the output file.

    Parameters
    ----------
    output_file : Chem.SDWriter
        The output file writer.
    best_mol : Chem.Mol
        The molecule with the best conformer.
    mol_id : str
        The ID of the molecule.
    alignment_score : object
        The alignment score object containing various scores.
    """
    best_mol.SetProp("CPD_ID", mol_id)
    best_mol.SetProp("SC_RDKit_score", f"{alignment_score.rdkit_score:.4f}")
    best_mol.SetProp("ESPSim", f"{alignment_score.esp_score:.4f}")
    best_mol.SetProp("ShapeSim", f"{alignment_score.shape_score:.4f}")
    output_file.write(best_mol)


def align_n_score_mols(
    molecules: Generator[Tuple[str, Chem.Mol], None, None],
    href_mol: Chem.Mol,
    output_file: Optional[Chem.SDWriter],
) -> None:
    """
    Process and score molecules, saving the best conformer to the output file.

    Parameters
    ----------
    molecules : Generator
        A generator of molecule ID and molecule tuples.
    href_mol : Chem.Mol
        The reference molecule for alignment.
    output_file : Optional[Chem.SDWriter]
        The output file writer, if specified.
    """
    for mol_id, mol in molecules:
        if mol is not None:
            alignment_score = score_alignment(mol, href_mol)
            save_best_mol(
                output_file, alignment_score.best_mol, mol_id, alignment_score
            )


def run_alignment(
    search_config: parsers.SearchConfig,
    num_conformers: int,
    output_file: Optional[str] = None,
) -> None:
    """
    Run alignment virtual screening pipeline.

    Parameters
    ----------
    search_config : parsers.SearchConfig
        A SearchConfig with the reference molecule and the molecules to be aligned.
    num_conformers : int
        Number of conformers to generate.
    output_file : Optional[str]
        File path to store the best overlays, if specified.

    Returns
    -------
    None, writes output to file.
    """
    if search_config.lig_3D_flag:
        href_mol = search_config.ref_mol

        logger.info(
            "Using the 3D conformations (of the reference and library cmpds) from the provided SD files"
        )
    else:
        href_mol = generate_conformers(
            search_config.ref_mol, num_conformers=1, optimize=True
        )
        output_path = Path(output_file).resolve().parent if output_file else Path(".")
        Chem.SDWriter(f"{output_path}/ref_mol-3D.sdf").write(href_mol)
        logger.info(
            "3D conformation of the reference compound has been stored in current working dir."
        )

    output_writer = Chem.SDWriter(str(output_file)) if output_file else None

    if search_config.lib_3D_flag:
        align_n_score_mols(
            tqdm(search_config.mols_to_align.items(), desc="3D alignment & scoring"),
            href_mol,
            output_writer,
        )
    else:
        with multiprocessing.Pool() as pool:
            molecules = pool.imap_unordered(
                map_conf_gen,
                prep_input(search_config.mols_to_align, num_conformers),
                chunksize=10,
            )
            align_n_score_mols(molecules, href_mol, output_writer)

    if output_writer:
        logger.info(
            f"3D overlays & calculated scores stored in {output_file}", colorize=True
        )
        output_writer.close()


def prep_input(
    mol_dict: Dict[str, Chem.Mol], num_conformers: int
) -> Generator[Tuple[Chem.Mol, str, int], None, None]:
    """
    Prepare input generator for conformer generation.

    Parameters
    ----------
    mol_dict : Dict[str, Chem.Mol]
        Dictionary of molecule IDs and molecules.
    num_conformers : int
        Number of conformers to generate for each molecule.

    Yields
    ------
    Generator[Tuple[Chem.Mol, str, int], None, None]
        Generator of molecule, molecule ID, and number of conformers tuples.
    """
    for mol_id, mol in tqdm(mol_dict.items(), desc="Conformer generation"):
        yield mol, mol_id, num_conformers
