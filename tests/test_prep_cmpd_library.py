import pandas as pd
import pytest
from click.testing import CliRunner

from lig3dlens.prep_cmpd_library import _calculate_descriptors, main


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


class TestPrepCmdLibraryCLI:
    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_no_options_provided_raises_error(self):
        """Test that missing both --filter and --no-filter raises error"""
        runner = CliRunner()
        result = runner.invoke(main, ["--in", "dummy.sdf", "--out", "dummy_out.sdf"])

        # Should fail with specific error message
        assert result.exit_code == 1
        assert (
            "Must provide either --filter <yaml_file> or --no-filter" in result.output
        )

    def test_help_shows_new_options(self):
        """Test that help text shows both filter options"""
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])

        assert result.exit_code == 0
        assert "--filter" in result.output
        assert "--no-filter" in result.output
        assert "Skip physicochemical property filtering entirely" in result.output

    def test_no_filter_flag_validation(self):
        """Test that --no-filter flag is properly recognized"""
        runner = CliRunner()
        # This will fail due to missing input file, but should pass validation
        result = runner.invoke(
            main, ["--in", "nonexistent.sdf", "--no-filter", "--out", "dummy_out.sdf"]
        )

        # Should fail on file not found, not on validation
        assert result.exit_code == 1
        # Should not see the validation error message
        assert "Must provide either --filter" not in result.output
