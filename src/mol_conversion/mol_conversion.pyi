"""
Type stubs for mol_conversion package
"""

from typing import Any, Dict, List, Optional, Union, Tuple

class FileConversion:
    """File conversion utilities for molecular formats"""
    
    @staticmethod
    def xyz_to_inchi(xyz_file: str) -> str: ...
    
    @staticmethod
    def xyz_to_inchikey(xyz_file: str) -> str: ...
    
    @staticmethod
    def xyz_to_smiles(xyz_file: str) -> str: ...
    
    @staticmethod
    def xyz_to_sdf(xyz_file: str, output_file: str) -> None: ...

class MemoryConversion:
    """Memory conversion utilities for molecular formats"""
    
    def __init__(self, xyz_string: str) -> None: ...
    
    def xyz_to_inchi(self) -> str: ...
    
    def xyz_to_inchikey(self) -> str: ...
    
    def xyz_to_smiles(self) -> str: ...
    
    def xyz_to_mol_string(self) -> str: ...

class GCNEncoding:
    """Graph Convolutional Network encoding for molecules"""
    
    def __init__(self, xyz_string: str) -> None: ...
    
    @property
    def encodings(self) -> Dict[str, List[str]]: ...

class XYZStandardizer:
    """XYZ coordinate standardizer"""
    
    def __init__(self, xyz_string: str) -> None: ...
    
    def to_standard_xyz(self) -> str: ...

def split_multiframe_xyz(xyz_string: str) -> List[str]: ...

def quick_xyz_to_inchi(xyz_file: str) -> str: ...

def quick_xyz_to_inchi_from_string(xyz_str: str) -> str: ...

def quick_split_multiframe(xyz_str: str) -> List[str]: ...

__all__ = [
    "FileConversion",
    "MemoryConversion", 
    "GCNEncoding",
    "XYZStandardizer",
    "split_multiframe_xyz",
    "quick_xyz_to_inchi",
    "quick_xyz_to_inchi_from_string", 
    "quick_split_multiframe"
]