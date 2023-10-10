# ----------------------------------------------------------------------------
# Authors: Bill Tatsis, Matt Seddon, Dan Mason, and Dan O'Donovan
# Company: Healx
# Date: August 17, 2023
#
# Description: This is the main script used to run the ligand-based 3D VS
# ----------------------------------------------------------------------------

import click
from loguru import logger

from lig3dlens.alignment import run_alignment
from lig3dlens.parsers import MolFileReader


@click.command(name="alignment")
@click.option(
    "--ref",
    "input_refmol",
    type=str,
    required=True,
    help="Input reference cmpd file name",
)
@click.option(
    "--lib",
    "input_mols_to_align",
    type=str,
    required=True,
    help="Input cmpd library file name",
)
@click.option(
    "--out",
    "output_file",
    type=str,
    required=True,
    help="Output .sd file with Mol blocks, ID columns, and 3D scores",
)
@click.option(
    "--conf",
    "num_conformers",
    type=int,
    default=10,
    help="Number of conformers for each of the cmpds in the library",
)
def main(input_refmol, input_mols_to_align, output_file, num_conformers):
    mol_file_reader = MolFileReader(input_refmol, input_mols_to_align)
    search_config = mol_file_reader.process()

    logger.info(
        "3D VS initialising...",
        colorize=True,
        format="<green>{time}</green> <level>{message}</level>",
    )

    run_alignment(search_config, num_conformers, output_file)

    logger.success("3D VS successfully terminated!", colorize=True)


if __name__ == "__main__":
    main()
