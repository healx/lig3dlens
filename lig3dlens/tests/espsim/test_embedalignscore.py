# Copied and modified from https://github.com/hesther/espsim/blob/master/scripts/

from espsim import EmbedAlignConstrainedScore
from rdkit import Chem
from rdkit.Chem import AllChem


def test_embed_and_align_score(train_and_test_smiles, train_and_test_mols):
    # This isn't a good test as it just outputs a score. How do we know it is what we expect?

    # Load molecules (9 reference molecules, 3 probe molecules)
    train_smiles, test_smiles = train_and_test_smiles
    train_mols, test_mols = train_and_test_mols

    # Embed first training molecule and optimize structure to extract the 3D pattern to get a 3D core molecule:
    AllChem.EmbedMolecule(train_mols[0], AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(train_mols[0])
    patt = Chem.MolFromSmiles("[H]OC([H])(C)C(=O)O[H]", sanitize=False)
    core = AllChem.DeleteSubstructs(
        AllChem.ReplaceSidechains(train_mols[0], patt), Chem.MolFromSmiles("*")
    )
    core.UpdatePropertyCache()
    train_mols[0].RemoveAllConformers()

    # For each probe molecule, call EmbedAlignConstrainedScore() and print all scores.
    for j, test_mol in enumerate(test_mols):
        results_shape, results_electrostatics = EmbedAlignConstrainedScore(
            test_mol, train_mols, core
        )
        print("Results for probe molecule", test_smiles[j])
        print("%40s %8s %8s" % ("Reference", "Shape", "ESP"))
        for i in range(len(results_shape)):
            print(
                "%40s %8.2f %8.2f"
                % (train_smiles[i], results_shape[i], results_electrostatics[i])
            )
