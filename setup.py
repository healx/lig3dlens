from setuptools import setup, find_packages

setup(
    name='lig3dlens',
    version='0.0.1',
    url='https://github.com/healx/lig3dlens.git',
    # description='Description of my package',
    packages=find_packages(include=["lig3dlens", "lig3dlens.*"]),
    # install_requires=['numpy >= 1.11.1', 'matplotlib >= 1.5.1'],
)
