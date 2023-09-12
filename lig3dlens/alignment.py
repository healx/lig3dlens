import multiprocessing
from pathlib import Path
from typing import Optional

from loguru import logger
from rdkit import Chem
from tqdm import tqdm

from lig3dlens import parsers
from lig3dlens.align3D_score import score_alignment
from lig3dlens.gen_conformers import generate_conformers, map_conf_gen


def run_alignment(
    search_config: parsers.SearchConfig,
    num_conformers: int,
    output_file: Optional[str] = None,
) -> None:
    """
    Run alignment virtual screening pipeline

    Parameters
    ----------
    search_config: A SearchConfig with the reference molecule and the molecules to be
        aligned
    num_conformers: Number of conformers
    output_file: Optional file path

    Returns
    -------
    None, writes output to file
    """

    # Generate reference mol
    href_mol = generate_conformers(
        search_config.ref_mol, num_conformers=1, optimize=True
    )

    # Store the 3D conformer of the reference compound
    # But first get the output file's path
    output_path = Path(output_file).resolve().parent
    Chem.SDWriter(f"{output_path}/ref_mol-3D.sdf").write(href_mol)

    logger.info(
        "3D conformation of the reference compound has been stored in current working"
        " dir."
    )

    # output file to store the best overlays
    if output_file is not None:
        output_file = Path(output_file)
        w = Chem.SDWriter(str(output_file))

    num_processes = multiprocessing.cpu_count() - 2
    with multiprocessing.Pool(num_processes) as pool:
        try:
            for mol_id, mol in pool.imap_unordered(
                map_conf_gen,
                prep_input(search_config.mols_to_align, num_conformers),
                chunksize=10,
            ):
                if mol is not None:
                    # Retrieve scores and conformer with the best 3D shape score
                    alignment_score = score_alignment(mol, href_mol)

                    best_mol = alignment_score.best_mol

                    # Scoring
                    best_mol.SetProp("CPD_ID", mol_id)
                    best_mol.SetProp(
                        "SC_RDKit_score", f"{alignment_score.rdkit_score:.4f}"
                    )
                    best_mol.SetProp("ESPSim", f"{alignment_score.esp_score:.4f}")
                    best_mol.SetProp("ShapeSim", f"{alignment_score.shape_score:.4f}")

                    # Write the conformer with the highest 3D shape score to the SD file
                    if output_file is not None:
                        w.write(best_mol)

        finally:
            if output_file is not None:
                logger.info(
                    f"3D overlays & calculated scores stored in {output_file}",
                    colorize=True,
                )
                w.close()


def prep_input(mol_dict: dict[str, Chem.Mol], num_conformers: int):
    for mol_id, mol in tqdm(mol_dict.items(), desc="Conformer generation"):
        yield mol, mol_id, num_conformers
