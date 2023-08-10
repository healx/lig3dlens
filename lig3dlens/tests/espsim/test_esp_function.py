# Copied and modified from https://github.com/hesther/espsim/blob/master/scripts/

import os

from espsim import GetEspSim, readMol2File, readMolFile, readSdfFile


def test_mol_esp_sim():
    # The following block of code reads in prealigned molecules in mol format and calculates their ESP similarity:
    mol1 = readMolFile(
        os.path.abspath(
            os.path.join(os.path.dirname(__file__), "test_data/prbmol1.mol")
        )
    )
    mol2 = readMolFile(
        os.path.abspath(
            os.path.join(os.path.dirname(__file__), "test_data/refmol1.mol")
        )
    )
    sim_esp = GetEspSim(mol1, mol2)
    print("%15s %5.2f" % ("ESP similarity (mols read from mol file):", sim_esp))


def test_mol2_esp_sim():
    # The following block of code reads in prealigned molecules in mol2 format with custom charges and calculates their ESP similarity:
    mol1, charge1 = readMol2File(
        os.path.abspath(
            os.path.join(os.path.dirname(__file__), "test_data/prbmol1.mol2")
        )
    )
    mol2, charge2 = readMol2File(
        os.path.abspath(
            os.path.join(os.path.dirname(__file__), "test_data/refmol1.mol2")
        )
    )
    sim_esp = GetEspSim(mol1, mol2, prbCharge=charge1, refCharge=charge2)
    print("%15s %5.2f" % ("ESP similarity (mols read from mol2 file):", sim_esp))


def test_sdf_esp_sim():
    # The following block of code reads in prealigned molecules in sdf format with custom charges as a comma-separated list
    # in the SDF file and calculates their ESP similarity:
    mol1, charge1 = readSdfFile(
        os.path.abspath(
            os.path.join(os.path.dirname(__file__), "test_data/prbmol1.sdf")
        )
    )
    mol2, charge2 = readSdfFile(
        os.path.abspath(
            os.path.join(os.path.dirname(__file__), "test_data/prbmol1.sdf")
        )
    )
    sim_esp = GetEspSim(mol1, mol2, prbCharge=charge1, refCharge=charge2)
    print("%15s %5.2f" % ("ESP similarity (mols read from sdf file):", sim_esp))
