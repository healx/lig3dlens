import pytest

from lig3dlens.parsers import MolFileReader


@pytest.fixture
def smiles_path(tmp_path):
    fh = tmp_path / "smiles.smi"
    fh.write_text("CN1C=NC2=C1C(=O)N(C(=O)N2C)C, caffeine")
    return fh


@pytest.fixture
def sdf_path(tmp_path):
    fh = tmp_path / "structure.sdf"
    fh.write_text(
        """\

     RDKit          3D

 24 25  0  0  0  0  0  0  0  0999 V2000
    3.1885   -0.9504   -0.2719 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1389    0.0635   -0.2674 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2940    1.4000   -0.4343 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1226    2.0856   -0.3746 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2155    1.1111   -0.1594 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8168   -0.0984   -0.0927 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0677   -1.2461    0.1310 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6254   -2.3770    0.1708 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2912   -1.1124    0.3062 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8872    0.1350    0.2110 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1389    0.2442    0.3474 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1309    1.2617   -0.0325 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7496    2.5929   -0.0888 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1235   -2.3072    0.5258 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9981   -1.6831   -1.0838 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.1809   -0.4811   -0.4418 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.2061   -1.4756    0.7060 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.2513    1.8777   -0.5989 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7163    2.5509   -0.6348 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9275    2.9640    0.9424 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1039    3.3192   -0.6270 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4969   -2.6845   -0.4494 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5517   -3.1135    1.0328 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9884   -2.0765    1.1838 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  7  9  1  0
  9 10  1  0
 10 11  2  0
 10 12  1  0
 12 13  1  0
  9 14  1  0
  6  2  1  0
 12  5  1  0
  1 15  1  0
  1 16  1  0
  1 17  1  0
  3 18  1  0
 13 19  1  0
 13 20  1  0
 13 21  1  0
 14 22  1  0
 14 23  1  0
 14 24  1  0
M  END
$$$$
    """
    )
    return fh


def test_smiles_parser_reads_smiles(smiles_path):
    parser = MolFileReader("", "")
    smiles = parser.smiles_parser(smiles_path)

    assert (
        "".join(map(lambda atom: atom.GetSymbol(), smiles.GetAtoms()))
        == "CNCNCCCONCONCC"
    )


def test_sdf_parser_reads_sdf(sdf_path):
    parser = MolFileReader("", "")
    sdf = parser.sdf_parser(sdf_path)

    assert (
        "".join(map(lambda atom: atom.GetSymbol(), sdf.GetAtoms()))
        == "C[n]1c(C(N(C)C(N2C)=O)=O)c2nc1"
    )
