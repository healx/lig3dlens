import collections
from loguru import logger
from rdkit import Chem

SearchConfig = collections.namedtuple("SearchConfig", ["ref_mol", "mols_to_align"])

class MolFileReader:
    def __init__(self, ref_mol_file_path: str, mols_to_align_file_path: str):
        self.ref_mol_file_path = ref_mol_file_path
        self.mols_to_align_file_path = mols_to_align_file_path

    def smiles_parser(self, file_path):
        with open(file_path, "r") as f:
            # Smiles files contains SMILES and title/ID, separated by commas
            ref_smi = f.readline().split(',')[0]
            ref_mol = Chem.MolFromSmiles(ref_smi)
        return ref_mol

    def sdf_parser(self, file_path):
        ref_mol = next(Chem.SDMolSupplier(file_path))
        # remove all conformers and return the Mol object
        ref_mol.RemoveAllConformers()
        return ref_mol

    def smiles_supplier_parser(self, file_path):
        suppl = Chem.SmilesMolSupplier(
            file_path, smilesColumn=0, nameColumn=1, delimiter=",", titleLine=True
        )
        return {m.GetProp("_Name"): m for m in suppl if m is not None}

    def sdf_supplier_parser(self, file_path):
        return {
            m.GetProp("Catalog ID"): m
            for m in Chem.SDMolSupplier(file_path)
            if m is not None
        }

    def csv_supplier_parser(self, file_path):
        suppl = Chem.SmilesMolSupplier(
            file_path, smilesColumn=0, nameColumn=1, delimiter=",", titleLine=True
        )
        return {m.GetProp("_Name"): m for m in suppl if m is not None}

    def determine_file_type(self, file_path):
        if file_path.endswith(".smi") or file_path.endswith(".smiles"):
            return "smiles"
        elif file_path.endswith(".sdf") or file_path.endswith(".sd"):
            return "sdf"
        elif file_path.endswith(".csv"):
            return "csv"
        else:
            raise ValueError(f"Unsupported file type: {file_path}")

    def parse_files(self):
        ref_file_type = self.determine_file_type(self.ref_mol_file_path)
        align_file_type = self.determine_file_type(self.mols_to_align_file_path)

        if ref_file_type == "smiles":
            ref_mol = self.smiles_parser(self.ref_mol_file_path)
        elif ref_file_type == "sdf":
            ref_mol = self.sdf_parser(self.ref_mol_file_path)
        else:
            raise ValueError("Only `.smi` and `.sdf` files are supported for reference mol.")

        if align_file_type == "smiles":
            mols_to_align = self.smiles_supplier_parser(self.mols_to_align_file_path)
        elif align_file_type == "sdf":
            mols_to_align = self.sdf_supplier_parser(self.mols_to_align_file_path)
        elif align_file_type == "csv":
            mols_to_align = self.csv_supplier_parser(self.mols_to_align_file_path)
        else:
            raise ValueError("Only `csv`, `.smi` and `.sdf` files are supported.")

        logger.debug(f"#{len(mols_to_align)} compounds were parsed")
        return ref_mol, mols_to_align
        

    def process(self) -> SearchConfig:
        ref_mol, mols_to_align = self.parse_files()
        return SearchConfig(ref_mol, mols_to_align)
