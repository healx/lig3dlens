from setuptools import find_packages, setup

setup(
    name="lig3dlens",
    version="0.0.3",
    url="https://github.com/healx/lig3dlens.git",
    description="Open source ligand-based 3D VS toolbox",
    packages=find_packages(include=["lig3dlens", "lig3dlens.*"]),
    install_requires=[
        # pin RDKit to reduce breakages
        "rdkit~=2024.9",
        # undocumented RDKit pins
        "pandas~=1.5",
        "numpy~=1.21",
    ],
    entry_points={
        "console_scripts": [
            "lig3dlens-prepare = lig3dlens.prep_cmpd_library:main",
            "lig3dlens-align = lig3dlens.main:main",
            "lig3dlens-cluster = lig3dlens.kmeans_clustering:main",
        ]
    },
)
