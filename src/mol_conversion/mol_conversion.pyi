"""
Type stubs for mol_conversion package
"""

import numpy as np
from rdkit.Chem import Mol

class FileConversion:
    """
    File conversion utilities for molecular formats.
    All methods are static and use Open Babel command-line tool.
    """
    
    @staticmethod
    def xyz_to_sdf(xyz_file: str, sdf_file: str) -> None: ...
    
    @staticmethod
    def xyz_to_inchi(xyz_file: str, fixed_h: bool = ...) -> str: ...
    
    @staticmethod
    def xyz_to_inchikey(xyz_file: str) -> str: ...
    
    @staticmethod
    def xyz_to_smiles(xyz_file: str) -> str: ...
    
    @staticmethod
    def xyz_to_pdb(xyz_file: str, pdb_file: str) -> None: ...
    
    @staticmethod
    def xyz_to_mol2(xyz_file: str, mol2_file: str) -> None: ...
    
    @staticmethod
    def xyz_to_cif(xyz_file: str, cif_file: str) -> None: ...
    
    @staticmethod
    def xyz_to_pdb_string(xyz_file: str) -> str: ...
    
    @staticmethod
    def xyz_to_mol2_string(xyz_file: str) -> str: ...
    
    @staticmethod
    def xyz_to_cif_string(xyz_file: str) -> str: ...
    
    @staticmethod
    def merge_sdf_files(sdf_list: list[str], output_path: str, need_remove: bool = ...) -> str: ...


class MemoryConversion:
    """
    Memory conversion utilities for molecular formats.
    All methods are static and use Open Babel Python API.
    """
    
    @staticmethod
    def xyz_to_mol_string(xyz_str: str) -> str: ...
    
    @staticmethod
    def xyz_to_rdkit_mol(xyz_str: str) -> Mol: ...
    
    @staticmethod
    def xyz_to_bond_order_matrix(xyz_str: str) -> np.ndarray: ...
    
    @staticmethod
    def xyz_to_inchi_string(xyz_str: str, fixed_h: bool = ...) -> str: ...
    
    @staticmethod
    def xyz_to_inchikey_string(xyz_str: str) -> str: ...
    
    @staticmethod
    def xyz_to_smiles_string(xyz_str: str) -> str: ...
    
    @staticmethod
    def xyz_to_sdf_string(xyz_str: str) -> str: ...
    
    @staticmethod
    def xyz_to_pdb_string(xyz_str: str) -> str: ...
    
    @staticmethod
    def xyz_to_mol2_string(xyz_str: str) -> str: ...
    
    @staticmethod
    def xyz_to_cif_string(xyz_str: str) -> str: ...


class GCNEncoding:
    """
    Graph Convolutional Network encoding for molecules.
    Generates three levels of atom encodings: GCN0, GCN1, GCN2.
    """
    
    def __init__(self, xyz_str: str) -> None: ...
    
    @property
    def encodings(self) -> dict[str, list[str]]: ...
    
    @property
    def atoms(self) -> list[str]: ...
    
    @property
    def adjacency(self) -> np.ndarray: ...
    
    @property
    def gcn0(self) -> list[str]: ...
    
    @property
    def gcn1(self) -> list[str]: ...
    
    @property
    def gcn2(self) -> list[str]: ...


def xyz_string_to_bond_order_matrix(xyz_str: str) -> np.ndarray: ...

def split_multiframe_xyz(xyz_str: str) -> list[str]: ...

def split_multiframe_xyz_with_comments(xyz_str: str) -> tuple[list[str], list[str]]: ...


__all__ = [
    "FileConversion",
    "MemoryConversion",
    "GCNEncoding",
    "xyz_string_to_bond_order_matrix",
    "split_multiframe_xyz",
    "split_multiframe_xyz_with_comments",
]
