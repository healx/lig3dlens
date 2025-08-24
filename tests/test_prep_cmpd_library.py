import pandas as pd

from lig3dlens.prep_cmpd_library import _calculate_descriptors


class TestCalculateDescriptors:
    def test_invalid_smiles_returns_none(self):
        sr = pd.Series({"test-col": "O=C([N-]c1cn+no1)C1CCCCC1"})

        # Should return None instead of raising exception
        result = _calculate_descriptors(0, sr, smiles_column="test-col")
        assert result is None

    def test_complex_heterocyclic_smiles_handling(self):
        """Test handling of complex heterocyclic SMILES that may cause RDKit parsing issues

        This tests a specific SMILES structure (thiazolo-pyrimidine nucleoside) that was
        reported to cause the prep_cmpd_library script to stall during processing.
        The function should gracefully handle any parsing failures by returning None
        rather than raising exceptions that halt the entire workflow.
        """
        # SMILES for a thiazolo-pyrimidine nucleoside that caused processing issues
        problematic_smiles = "O=c1cc(=O)n2ccsc2n1C1OC(CO)C(O)C1O"
        sr = pd.Series({"smiles_sdt": problematic_smiles})

        # This should not crash and should return None for invalid SMILES
        result = _calculate_descriptors(0, sr, smiles_column="smiles_sdt")
        # If the SMILES is actually valid, result should be a Series; if invalid, should be None
        assert result is None or isinstance(result, pd.Series)

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
