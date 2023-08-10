import pytest
from rdkit import Chem


@pytest.fixture(scope="function")
def train_and_test_smiles():
    train_smiles = [
        "C1=CC=C(C=C1)C(C(=O)O)O",
        "CCC(C(=O)O)O",
        "OC(C(O)=O)c1ccc(Cl)cc1",
        "C1=CC(=CC=C1C(C(=O)O)O)O",
        "COc1ccc(cc1)C(O)C(O)=O",
        "OC(C(O)=O)c1ccc(cc1)[N+]([O-])=O",
        "CCCC(C(=O)O)O",
        "CCC(C)C(C(=O)O)O",
        "CC(C(=O)O)O",
    ]
    test_smiles = ["c1c(C)cc(cc1)C(O)C(O)=O", "Cc1ccc(cc1)C(O)C(O)=O", "C(C(C(=O)O)O)O"]
    return train_smiles, test_smiles


@pytest.fixture(scope="function")
def train_and_test_mols(train_and_test_smiles):
    train_smiles, test_smiles = train_and_test_smiles
    train_mols = [Chem.AddHs(Chem.MolFromSmiles(x)) for x in train_smiles]
    test_mols = [Chem.AddHs(Chem.MolFromSmiles(x)) for x in test_smiles]
    return train_mols, test_mols
