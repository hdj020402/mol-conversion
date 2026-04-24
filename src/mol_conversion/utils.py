import numpy as np
from openbabel import openbabel

def xyz_string_to_bond_order_matrix(xyz_str):
    """
    Build bond order matrix from XYZ string using Open Babel low-level API (no file I/O)
    
    Args:
        xyz_str (str): XYZ format molecular structure string
        
    Returns:
        np.ndarray: Bond order matrix with shape (n_atoms, n_atoms)
    """
    # Initialize converter and molecule object
    ob_conv = openbabel.OBConversion()
    ob_mol = openbabel.OBMol()
    
    # Set input format to "xyz"
    if not ob_conv.SetInFormat("xyz"):
        raise ValueError("XYZ format not supported")
    
    # Read from string (key! no file I/O)
    success = ob_conv.ReadString(ob_mol, xyz_str)
    if not success:
        raise ValueError("Failed to parse XYZ string")

    n_atoms = ob_mol.NumAtoms()
    bond_matrix = np.zeros((n_atoms, n_atoms), dtype=int)

    # Iterate through all bonds
    for ob_bond in openbabel.OBMolBondIter(ob_mol):
        i = ob_bond.GetBeginAtomIdx() - 1  # OB index starts from 1
        j = ob_bond.GetEndAtomIdx() - 1
        bo = ob_bond.GetBondOrder()
        bond_matrix[i, j] = bo
        bond_matrix[j, i] = bo

    return bond_matrix

def _parse_multiframe_xyz_core(xyz_str: str):
    """
    Internal core parser: single pass through, returns both frame string list and comment list.
    Not exposed externally, only called by public functions.
    
    Args:
        xyz_str (str): Multi-frame XYZ string
        
    Returns:
        tuple: (frames list, comments list)
    """
    lines = xyz_str.strip().splitlines()
    frames = []
    comments = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n_atoms = int(lines[i].strip())
        except ValueError:
            raise ValueError(f"Expected atom count at line {i+1}, got: {lines[i]}")

        if i + n_atoms + 2 > len(lines):
            raise ValueError(f"Incomplete frame starting at line {i+1}")

        # Extract comment (second line)
        comment = lines[i + 1] if i + 1 < len(lines) else ""

        # Construct frame string
        frame_lines = lines[i : i + n_atoms + 2]
        frame_str = "\n".join(frame_lines) + "\n"

        frames.append(frame_str)
        comments.append(comment)

        i += n_atoms + 2

    return frames, comments


# Public interface 1: split frames only (maintain original function behavior)
def split_multiframe_xyz(xyz_str: str):
    """
    Split multi-frame XYZ string into individual frame strings
    
    Args:
        xyz_str (str): Multi-frame XYZ string
        
    Returns:
        list: List of individual frame strings
    """
    frames, _ = _parse_multiframe_xyz_core(xyz_str)
    return frames


# Public interface 2: get frames and comments
def split_multiframe_xyz_with_comments(xyz_str: str):
    """
    Split multi-frame XYZ string and return both frames and comments
    
    Args:
        xyz_str (str): Multi-frame XYZ string
        
    Returns:
        tuple: (frames list, comments list)
    """
    return _parse_multiframe_xyz_core(xyz_str)
