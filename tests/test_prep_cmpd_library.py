import pandas as pd
import pytest

from lig3dlens.prep_cmpd_library import _calculate_descriptors


class TestCalculateDescriptors:
    def test_invalid_smiles_raises_value_error(self):
        sr = pd.Series({"test-col": "O=C([N-]c1cn+no1)C1CCCCC1"})

        with pytest.raises(ValueError, match="Invalid SMILES"):
            _calculate_descriptors(0, sr, smiles_column="test-col")

    def test_mesoionic_smiles(self):
        sr = pd.Series({"test-col": "O=[S](=O)([O-])[NH3+]"})

        desc = _calculate_descriptors(0, sr, smiles_column="test-col")

        # assert that we get a series with expected columns
        assert isinstance(desc, pd.Series)
        assert desc.shape == (12,)
        assert desc.index.to_list() == [
            "test-col",
            "mw",
            "hba",
            "hbd",
            "rot_bonds",
            "arom_rings",
            "fsp3",
            "tpsa",
            "clogp",
            "qed",
            "heavy_atms",
            "stereo_cntrs",
        ]
