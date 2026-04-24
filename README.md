# mol-conversion

A comprehensive Python package for molecular structure format conversion, focusing on XYZ format processing with support for multiple output formats.

## Features

- **File-to-file conversion** using Open Babel subprocess
- **Memory-to-memory conversion** using Open Babel Python API
- **GCN-based molecular encoding** with three levels of complexity
- **Molecular coordinate standardization**
- **Multi-frame XYZ string parsing**
- **Support for 8+ formats**: XYZ, SDF, InChI, InChIKey, SMILES, PDB, MOL2, and CIF
- **Lazy loading** for Open Babel-dependent features
- **Comprehensive test suite** with real molecular data

## Installation

```bash
pip install mol-conversion
```

Or for development:

```bash
git clone https://github.com/hdj020402/mol-conversion.git
cd mol-conversion
pip install -e .
```

## Architecture

The package uses a modern `src` layout with lazy loading:

```
src/mol_conversion/
├── __init__.py          # Package entry point with lazy loading
├── file_conversion.py   # File-to-file conversions
├── memory_conversion.py # Memory-to-memory conversions
├── gcn_encoding.py      # GCN molecular encoding
└── utils.py             # Utility functions (multi-frame XYZ parsing, etc.)
```

## Quick Start

### File Conversion

```python
from mol_conversion import FileConverter

# Convert XYZ to various formats
FileConverter.xyz_to_sdf("molecule.xyz", "molecule.sdf")
FileConverter.xyz_to_pdb("molecule.xyz", "molecule.pdb")

# Get string representations
inchi = FileConverter.xyz_to_inchi("molecule.xyz")
inchikey = FileConverter.xyz_to_inchikey("molecule.xyz")
smiles = FileConverter.xyz_to_smiles("molecule.xyz")
```

### Memory Conversion

```python
from mol_conversion import MemoryConverter

# Convert XYZ string to various formats
xyz_string = """70
Energy:-1025.075684 a.u. s1:  1 s2:  1 rot_angle:260 escan:0.0044
C      1.381157     0.209710     1.978440
C      1.326816     0.000000     0.469020
... (remaining atoms)
"""

mol_string = MemoryConverter.xyz_to_mol_string(xyz_string)
inchi = MemoryConverter.xyz_to_inchi_string(xyz_string)
bond_matrix = MemoryConverter.xyz_to_bond_order_matrix(xyz_string)
```

### GCN Encoding

```python
from mol_conversion import GCNEncoder

# Generate hierarchical GCN encodings from various formats
encoder = GCNEncoder.from_xyz(xyz_string)
encoder = GCNEncoder.from_smiles("c1ccccc1")
encoder = GCNEncoder.from_inchi("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")

# Access different levels of encoding
encodings = encoder.encodings
gcn0 = encodings["gcn0"]  # Basic atom + neighbor count encoding
gcn1 = encodings["gcn1"]  # First-order neighbor encoding
gcn2 = encodings["gcn2"]  # Second-order neighbor encoding

# Hydrogen atoms return empty strings
```

### Multi-frame XYZ Processing

```python
from mol_conversion import split_multiframe_xyz, split_multiframe_xyz_with_comments

# Split multi-frame XYZ into individual frames
frames = split_multiframe_xyz(multiframe_xyz_string)

# Split with comments for additional metadata
frames, comments = split_multiframe_xyz_with_comments(multiframe_xyz_string)
```

### Log Level Configuration

```python
from mol_conversion import set_log_level, get_log_level

# Get current log level (default: 'error')
print(get_log_level())  # 'error'

# Set log level
set_log_level('none')     # Disable all logging
set_log_level('error')    # Show only errors (default)
set_log_level('warning')  # Show warnings and errors
set_log_level('info')     # Show info, warnings, and errors
set_log_level('debug')    # Show all log messages

# Case-insensitive
set_log_level('ERROR')    # Also works
```

## API Reference

### FileConverter Class

Static methods for file-to-file conversions:

- `xyz_to_sdf(xyz_file, sdf_file)` - Convert XYZ to SDF
- `xyz_to_inchi(xyz_file)` - Convert XYZ to InChI string
- `xyz_to_inchikey(xyz_file)` - Convert XYZ to InChIKey string
- `xyz_to_smiles(xyz_file)` - Convert XYZ to SMILES string
- `xyz_to_pdb(xyz_file, pdb_file)` - Convert XYZ to PDB
- `xyz_to_mol2(xyz_file, mol2_file)` - Convert XYZ to MOL2
- `xyz_to_cif(xyz_file, cif_file)` - Convert XYZ to CIF
- `merge_sdf_files(sdf_list, output_path, need_remove=False)` - Merge multiple SDF files

### MemoryConverter Class

Static methods for memory-to-memory conversions:

- `xyz_to_mol_string(xyz_string)` - Convert XYZ to MOL string
- `xyz_to_rdkit_mol(xyz_string)` - Convert XYZ to RDKit Mol object
- `xyz_to_bond_order_matrix(xyz_string)` - Generate bond order matrix
- `xyz_to_inchi_string(xyz_string)` - Convert XYZ to InChI string
- `xyz_to_inchikey_string(xyz_string)` - Convert XYZ to InChIKey string
- `xyz_to_smiles_string(xyz_string)` - Convert XYZ to SMILES string
- `xyz_to_sdf_string(xyz_string)` - Convert XYZ to SDF string
- `xyz_to_pdb_string(xyz_string)` - Convert XYZ to PDB string
- `xyz_to_mol2_string(xyz_string)` - Convert XYZ to MOL2 string
- `xyz_to_cif_string(xyz_string)` - Convert XYZ to CIF string

### GCNEncoder Class

Generate atom-level GCN-style encodings from multiple molecular formats:

- `from_xyz(xyz_string)` - Initialize from XYZ format
- `from_smiles(smiles_string)` - Initialize from SMILES format
- `from_inchi(inchi_string)` - Initialize from InChI format
- `from_sdf(sdf_string)` - Initialize from SDF format
- `from_mol(mol_string)` - Initialize from MOL format
- `from_mol2(mol2_string)` - Initialize from MOL2 format
- `from_pdb(pdb_string)` - Initialize from PDB format
- `from_cif(cif_string)` - Initialize from CIF format
- `encodings` property - Returns dictionary with 'gcn0', 'gcn1', 'gcn2' encodings

### Utility Functions

- `split_multiframe_xyz(xyz_string)` - Split multi-frame XYZ into list of frames
- `split_multiframe_xyz_with_comments(xyz_string)` - Split with comments
- `xyz_string_to_bond_order_matrix(xyz_string)` - Generate bond order matrix

### Log Level Functions

- `set_log_level(level)` - Set log level for Open Babel and RDKit
  - `level`: One of `'none'`, `'error'`, `'warning'`, `'info'`, `'debug'`
- `get_log_level()` - Get current log level

## Testing

The package includes comprehensive tests using real molecular data:

```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=mol_conversion

# Run specific module tests
pytest tests/test_file_conversion.py
pytest tests/test_memory_conversion.py
pytest tests/test_utils.py

# Run integration tests
pytest tests/test_integration.py
```

### Test Features

- **Real molecular data** - Tests use actual 70-atom molecules
- **Dependency detection** - Automatically skips tests if Open Babel unavailable
- **Error handling** - Comprehensive error condition testing
- **Performance testing** - Large molecule handling verification
- **Integration testing** - Cross-module consistency validation

## Dependencies

- **Python 3.8+**
- **openbabel-wheel** - Core molecular conversion (pip installable)
- **NumPy** - Numerical computations
- **RDKit** - Chemical informatics
- **tqdm** - Progress bars

## Project Structure

```
mol_conversion/
├── src/mol_conversion/          # Source code
├── tests/                       # Test suite
├── tests/data/                  # Test data files
├── docs/                        # Documentation
├── pyproject.toml              # Modern packaging config
└── README.md                   # This file
```

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Changelog

### v0.1.0
- Initial release with core conversion functionality
- GCN encoding support
- Multi-frame XYZ processing
- Comprehensive test suite
- Modern src layout packaging