import pytest
from rdkit import Chem

from lig3dlens.gen_conformers import conf_gen, generate_conformers


def _count_explicit_hydrogens(mol: Chem.Mol) -> int:
    return sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")


def test_generate_conformers_with_no_pruning_produces_requested_count():
    mol = Chem.MolFromSmiles("CCO")

    # No pruning ensures we get exactly the requested number of conformers
    out = generate_conformers(
        mol,
        num_conformers=3,
        prune_rms_threshold=0.0,
        add_hydrogens=True,
        optimize=False,
    )

    assert out.GetNumConformers() == 3
    # Adding hydrogens should make explicit H atoms present
    assert _count_explicit_hydrogens(out) > 0


def test_generate_conformers_zero_requested_raises():
    mol = Chem.MolFromSmiles("CC")

    with pytest.raises(ValueError, match="Failed to generate any 3D conformers"):
        generate_conformers(
            mol,
            num_conformers=0,
            prune_rms_threshold=0.0,
            add_hydrogens=True,
            optimize=False,
        )


def test_generate_conformers_respects_add_hydrogens_flag():
    mol = Chem.MolFromSmiles("CCC")

    no_h = generate_conformers(
        mol,
        num_conformers=1,
        prune_rms_threshold=0.0,
        add_hydrogens=False,
        optimize=False,
    )
    with_h = generate_conformers(
        mol,
        num_conformers=1,
        prune_rms_threshold=0.0,
        add_hydrogens=True,
        optimize=False,
    )

    assert _count_explicit_hydrogens(no_h) == 0
    assert _count_explicit_hydrogens(with_h) > 0


def test_conf_gen_returns_tuple_and_conformers():
    mol = Chem.MolFromSmiles("CCO")
    mol_id, out = conf_gen(mol, "mol-001", 2)

    assert mol_id == "mol-001"
    assert isinstance(out, Chem.Mol)
    assert out.GetNumConformers() >= 1
