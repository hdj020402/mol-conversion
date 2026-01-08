"""
Molecular Conversion Tools

A comprehensive toolkit for molecular structure format conversion,
focusing on XYZ format processing with support for multiple output formats.

Features:
- File-to-file conversion using Open Babel subprocess
- Memory-to-memory conversion using Open Babel Python API
- GCN-based molecular encoding
- Molecular coordinate standardization
- Multi-frame XYZ string parsing
"""

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"
__license__ = "MIT"

from typing import TYPE_CHECKING

# Core classes that don't require Open Babel
from .file_conversion import FileConversion

# Open Babel-dependent imports - these will be imported lazily
_memory_conversion_imported = False
_utils_imported = False

# Type checking imports
if TYPE_CHECKING:
    from .memory_conversion import MemoryConversion
    from .utils import GCNEncoding

def _import_memory_conversion():
    """Lazy import of MemoryConversion class"""
    global _memory_conversion_imported, MemoryConversion
    if not _memory_conversion_imported:
        try:
            from .memory_conversion import MemoryConversion
            _memory_conversion_imported = True
        except ImportError as e:
            raise ImportError(
                f"MemoryConversion requires Open Babel. Please install with: "
                f"conda install openbabel -c conda-forge. Original error: {e}"
            )
    return MemoryConversion

def _import_utils():
    """Lazy import of utility functions"""
    global _utils_imported, GCNEncoding, split_multiframe_xyz
    if not _utils_imported:
        try:
            from .utils import GCNEncoding, split_multiframe_xyz
            _utils_imported = True
        except ImportError as e:
            raise ImportError(
                f"Utils module requires Open Babel. Please install with: "
                f"conda install openbabel -c conda-forge. Original error: {e}"
            )
    return GCNEncoding, split_multiframe_xyz

def __getattr__(name):
    """Lazy attribute access for Open Babel-dependent classes and functions"""
    if name == 'MemoryConversion':
        return _import_memory_conversion()
    elif name == 'GCNEncoding':
        return _import_utils()[0]
    elif name == 'split_multiframe_xyz':
        return _import_utils()[1]
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

__all__ = [
    "FileConversion",
    "MemoryConversion", 
    "GCNEncoding",
    "split_multiframe_xyz"
]

# Convenience imports
def quick_xyz_to_inchi(xyz_file: str) -> str:
    """Quick conversion from XYZ file to InChI"""
    return FileConversion.xyz_to_inchi(xyz_file)

def quick_xyz_to_inchi_from_string(xyz_str: str) -> str:
    """Quick conversion from XYZ string to InChI"""
    return _import_memory_conversion().xyz_to_inchi_string(xyz_str)

def quick_split_multiframe(xyz_str: str):
    """Quick multi-frame XYZ splitting"""
    return _import_utils()[1](xyz_str)

# Expose quick functions in __all__
__all__.extend([
    "quick_xyz_to_inchi",
    "quick_xyz_to_inchi_from_string", 
    "quick_split_multiframe"
])