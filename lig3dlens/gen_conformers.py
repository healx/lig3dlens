# ----------------------------------------------------------------------------
# Authors: Bill Tatsis, Matt Seddon, Dan Mason, and Dan O'Donovan
# Company: Healx
# Date: August 17, 2023
#
# Description: This script generates 3D conformers for a given Mol.
# ----------------------------------------------------------------------------

from typing import Tuple

from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem, Mol
from rdkit.Chem.rdDistGeom import ETKDGv3


def generate_conformers(
    molecule: Mol,
    num_conformers: int,
    prune_rms_threshold: float = 0.5,
    add_hydrogens: bool = True,
    optimize: bool = True,  # this option needs to be ON!
) -> Mol:
    """
    Generate RDKit conformers for a Mol object.

    Parameters
    ----------
    molecule: RDKit Mol object
        The molecule for which conformers need to be generated.
    num_conformers: int
        Number of conformers to generate.
    prune_rms_threshold: float, optional
        RMS threshold for pruning, default 0.5.
    add_hydrogens: bool, optional
        Whether to add Hydrogen atoms to the Mol, default True.
    optimize: bool, optional
        Whether to optimize conformer generation, default False.

    Returns
    -------
    Mol
        Mutated Mol object with conformers generated.
    """

    params = ETKDGv3()
    # Configure embedding parameters explicitly rather than mixing signatures
    params.pruneRmsThresh = prune_rms_threshold
    params.randomSeed = 42
    params.useRandomCoords = True
    params.enforceChirality = True

    mol = Mol(molecule)
    mol.RemoveAllConformers()

    # Add "explicit" hydrogens
    if add_hydrogens:
        mol = Chem.AddHs(mol, addCoords=True)

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

    if len(cids) == 0:
        raise ValueError(
            "Failed to generate any 3D conformers for "
            f"{Chem.MolToSmiles(Chem.RemoveHs(mol))}"
        )

    if optimize:
        try:
            # Energy-minimize the conformers
            AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s", numThreads=0)
        except Exception as e:
            logger.warning(f"Could not optimize conformers due to an error: {e}")

    return mol


def conf_gen(molecule: Mol, mol_id: str, num_conformers: int) -> Tuple[str, Mol]:
    """
    Generate conformers for an RDKit Mol.

    Parameters
    ----------
    molecule: RDKit Mol object
        The molecule for which conformers need to be generated.
    mol_id: str
        A string representing the molecule ID.
    num_conformers: int
        Number of conformers to be generated.

    Returns
    -------
    Tuple[str, Mol]
        A tuple containing the molecule ID and the Mol object mutated to have conformers
        computed.
    """
    if isinstance(molecule, Mol):
        h_mol = generate_conformers(
            molecule, num_conformers, prune_rms_threshold=0.5, add_hydrogens=True
        )

        return mol_id, h_mol


# Helper function
# We can also use functools.partial to create a partial function of conf_gen
# map_conf_gen = partial(conf_gen)


def map_conf_gen(args):
    return conf_gen(*args)
