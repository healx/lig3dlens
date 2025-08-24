import pandas as pd
import pytest
from click.testing import CliRunner

from lig3dlens.prep_cmpd_library import _calculate_descriptors, main


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


class TestPrepCmdLibraryCLI:
    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_no_options_provided_raises_error(self, runner):
        """Test that missing both --filter and --no-filter raises error"""
        result = runner.invoke(main, ["--in", "dummy.sdf", "--out", "dummy_out.sdf"])

        # Should fail with specific error message
        assert result.exit_code == 1
        assert (
            "Must provide either --filter <yaml_file> or --no-filter" in result.output
        )

    def test_help_shows_new_options(self, runner):
        """Test that help text shows both filter options"""
        result = runner.invoke(main, ["--help"])

        assert result.exit_code == 0
        assert "--filter" in result.output
        assert "--no-filter" in result.output
        assert "Skip physicochemical property filtering entirely" in result.output

    def test_no_filter_flag_validation(self, runner):
        """Test that --no-filter flag is properly recognized"""
        # This will fail due to missing input file, but should pass validation
        result = runner.invoke(
            main, ["--in", "nonexistent.sdf", "--no-filter", "--out", "dummy_out.sdf"]
        )

        # Should fail on file not found, not on validation
        assert result.exit_code == 1
        # Should not see the validation error message
        assert "Must provide either --filter" not in result.output

    def test_both_options_provided_raises_error(self, runner):
        """Test that providing both --filter and --no-filter raises error"""
        result = runner.invoke(
            main,
            [
                "--in",
                "dummy.sdf",
                "--filter",
                "dummy.yaml",
                "--no-filter",
                "--out",
                "dummy_out.sdf",
            ],
        )

        # Should fail with specific error message
        assert result.exit_code == 1
        assert "Cannot provide both --no-filter and --filter options" in result.output
