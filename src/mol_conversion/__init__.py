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

__version__ = "0.2.0"
__author__ = "Dejun Hu"
__email__ = "hudejun2002@gmail.com"
__license__ = "MIT"

from typing import TYPE_CHECKING, Literal

# Core classes that don't require Open Babel
from .file_conversion import FileConverter
from .memory_conversion import MemoryConverter
from .gcn_encoding import GCNEncoder
from .utils import split_multiframe_xyz, split_multiframe_xyz_with_comments

# Log level type
LogLevel = Literal["none", "error", "warning", "info", "debug"]

# Current log level (module-level state)
_current_log_level: LogLevel = "error"


def set_log_level(level: LogLevel) -> None:
    """
    Set the log level for Open Babel and RDKit libraries.
    
    Args:
        level: Log level, one of 'none', 'error', 'warning', 'info', 'debug'.
            - 'none': Disable all logging output
            - 'error': Show only errors
            - 'warning': Show warnings and errors
            - 'info': Show info, warnings, and errors
            - 'debug': Show all log messages
    
    Example:
        >>> from mol_conversion import set_log_level
        >>> set_log_level('none')  # Silence all logs
        >>> set_log_level('error')  # Show only errors (default)
    """
    global _current_log_level
    
    level = level.lower()
    valid_levels = ("none", "error", "warning", "info", "debug")
    if level not in valid_levels:
        raise ValueError(f"Invalid log level '{level}'. Must be one of: {valid_levels}")
    
    _current_log_level = level
    
    # Configure RDKit log level
    from rdkit import RDLogger
    logger = RDLogger.logger()
    
    if level == "none":
        RDLogger.DisableLog('rdApp.*')
    elif level == "error":
        logger.setLevel(RDLogger.ERROR)
    elif level == "warning":
        logger.setLevel(RDLogger.WARNING)
    elif level == "info":
        logger.setLevel(RDLogger.INFO)
    elif level == "debug":
        logger.setLevel(RDLogger.DEBUG)
    
    # Configure Open Babel log level
    from openbabel import openbabel
    
    # Open Babel log levels: 0=Error, 1=Warning, 2=Info, 3=Debug
    level_map = {
        "error": 0,
        "warning": 1,
        "info": 2,
        "debug": 3,
    }
    
    if level == "none":
        openbabel.obErrorLog.StopLogging()
    else:
        openbabel.obErrorLog.SetOutputLevel(level_map[level])


def get_log_level() -> LogLevel:
    """
    Get the current log level.
    
    Returns:
        Current log level string.
    """
    return _current_log_level


# Set default log level on import
set_log_level("error")

# Open Babel-dependent imports - these will be imported lazily
_memory_conversion_imported = False
_utils_imported = False
_gcn_encoding_imported = False


def _import_memory_conversion():
    """Lazy import of MemoryConverter class"""
    global _memory_conversion_imported, MemoryConverter
    if not _memory_conversion_imported:
        try:
            from .memory_conversion import MemoryConverter
            _memory_conversion_imported = True
        except ImportError as e:
            raise ImportError(
                f"MemoryConverter requires openbabel-wheel. Please install with: "
                f"pip install openbabel-wheel. Original error: {e}"
            )
    return MemoryConverter

def _import_utils():
    """Lazy import of utility functions"""
    global _utils_imported, split_multiframe_xyz, split_multiframe_xyz_with_comments
    if not _utils_imported:
        try:
            from .utils import split_multiframe_xyz, split_multiframe_xyz_with_comments
            _utils_imported = True
        except ImportError as e:
            raise ImportError(
                f"Utils module requires openbabel-wheel. Please install with: "
                f"pip install openbabel-wheel. Original error: {e}"
            )
    return split_multiframe_xyz, split_multiframe_xyz_with_comments

def _import_gcn_encoding():
    """Lazy import of GCNEncoder class"""
    global _gcn_encoding_imported, GCNEncoder
    if not _gcn_encoding_imported:
        try:
            from .gcn_encoding import GCNEncoder
            _gcn_encoding_imported = True
        except ImportError as e:
            raise ImportError(
                f"GCNEncoder requires openbabel-wheel. Please install with: "
                f"pip install openbabel-wheel. Original error: {e}"
            )
    return GCNEncoder

def __getattr__(name):
    """Lazy attribute access for Open Babel-dependent classes and functions"""
    if name == 'MemoryConverter':
        return _import_memory_conversion()
    elif name == 'GCNEncoder':
        return _import_gcn_encoding()
    elif name == 'split_multiframe_xyz':
        return _import_utils()[0]
    elif name == 'split_multiframe_xyz_with_comments':
        return _import_utils()[1]
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

__all__ = [
    "FileConverter",
    "MemoryConverter",
    "GCNEncoder",
    "split_multiframe_xyz",
    "split_multiframe_xyz_with_comments",
    "set_log_level",
    "get_log_level",
]
