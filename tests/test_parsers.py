from textwrap import dedent

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
CT1001987571


 24 25  0  3  0               999 V2000
   -0.0171    1.4073    0.0098 C   0 00  0  0  0  0  0  0  0  0  0  0
    0.0021   -0.0041    0.0020 C   0 00  0  0  0  0  0  0  0  0  0  0
    1.1868    2.1007    0.0020 C   0 00  0  0  0  0  0  0  0  0  0  0
   -1.0133    2.3630    0.0190 N   0 00  0  0  0  0  0  0  0  0  0  0
    2.3717    1.3829   -0.0136 N   0 00  0  0  0  0  0  0  0  0  0  0
    0.8932    3.4034    0.0118 N   0 00  0  0  0  0  0  0  0  0  0  0
    1.1884   -0.6467   -0.0128 N   0 00  0  0  0  0  0  0  0  0  0  0
   -1.0401   -0.6344    0.0090 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3458    0.0368   -0.0214 C   0 00  0  0  0  0  0  0  0  0  0  0
    3.6549    2.0897   -0.0220 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2155   -2.1115   -0.0209 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3959   -0.5761   -0.0355 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4053    3.5654    0.0231 C   0 00  0  0  0  0  0  0  0  0  0  0
   -2.4574    2.1166    0.0226 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9831    2.2592    1.0035 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.3975    1.4884   -0.5465 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5388    3.0475   -0.5293 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124   -2.4692   -1.0505 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1169   -2.4610    0.4825 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3373   -2.4940    0.4993 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9129    4.5186    0.0303 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8119    2.0494    1.0512 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9671    2.9358   -0.4846 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6677    1.1812   -0.4960 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2 01  0  1  0  0
  1  3 02  0  2  0  0
  1  4 01  0  1  0  0
  2  7 01  0  1  0  0
  2  8  2  0  0  0  0
  3  5 01  0  1  0  0
  3  6 01  0  1  0  0
  4 13 01  0  1  0  0
  4 14  1  0  0  0  0
  5  9 01  0  1  0  0
  5 10  1  0  0  0  0
  6 13 02  0  1  0  0
  7  9 01  0  1  0  0
  7 11  1  0  0  0  0
  9 12  2  0  0  0  0
 10 15  1  0  0  0  0
 10 16  1  0  0  0  0
 10 17  1  0  0  0  0
 11 18  1  0  0  0  0
 11 19  1  0  0  0  0
 11 20  1  0  0  0  0
 13 21  1  0  0  0  0
 14 22  1  0  0  0  0
 14 23  1  0  0  0  0
 14 24  1  0  0  0  0
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
        "".join(map(lambda atom: atom.GetSymbol(), sdf.GetAtoms())) == "CCCNNNNOCCCOCC"
    )
