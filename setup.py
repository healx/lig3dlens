from setuptools import find_packages, setup

setup(
    name="lig3dlens",
    version="0.0.1",
    url="https://github.com/healx/lig3dlens.git",
    description="Open source ligand-based 3D VS toolbox",
    packages=find_packages(include=["lig3dlens", "lig3dlens.*"]),
    entry_points={"console_scripts": ["lig3dlens = lig3dlens.main:main"]}
    # install_requires=['numpy >= 1.11.1', 'matplotlib >= 1.5.1'],
)
