import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest
from rdkit.Chem import PandasTools


@pytest.fixture(scope="function")
def setup_teardown_test():
    # Set up
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    # Teardown
    shutil.rmtree(temp_dir)


@pytest.mark.skip("not a unit tests")
def test_lig3dlens_out(setup_teardown_test):
    temp_dir = setup_teardown_test

    lig3dlens_py = Path("./main.py")
    ref_mol = Path("./tests/lig3dlens/test_data/ref_mol.smi")
    input_file = Path(
        "./tests/lig3dlens/test_data/Enamine_hts_collection_202303_first500.smi"
    )
    expected_output = Path(
        "./tests/lig3dlens/test_data/Enamine_hts_collection_202303_first500_VS_results.sdf"
    )

    # Run the main script
    subprocess.run(
        [
            "python",
            str(lig3dlens_py),
            "--ref",
            str(ref_mol),
            "--lib",
            str(input_file),
            "--out",
            str(temp_dir / "test_out.sdf"),
        ],
        stdout=subprocess.PIPE,
        text=True,
    )

    # Check if the scores in the temp. output file are the same as in the expected
    # output file
    temp_df = PandasTools.LoadSDF(temp_dir / "test_out.sdf")
    expected_output_df = PandasTools.LoadSDF(expected_output)

    assert temp_df["SC_RDKit_score"].equals(expected_output_df["SC_RDKit_score"])
    assert temp_df["ESPSim"].equals(expected_output_df["ESPSim"])
    assert temp_df["ShapeSim"].equals(expected_output_df["ShapeSim"])
