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

class GCNEncoding:
    """
    Generate atom-level GCN-style encodings (GCN0, GCN1, GCN2) from an XYZ string.
    Hydrogen atoms are encoded as empty strings ('').
    """

    def __init__(self, xyz_str: str):
        self.atoms, self.adjacency = self._parse_xyz(xyz_str)
        self._compute_encodings()

    def _parse_xyz(self, xyz_str: str):
        """Parse XYZ string and extract atom symbols and adjacency matrix"""
        ob_conv = openbabel.OBConversion()
        if not ob_conv.SetInFormat("xyz"):
            raise ValueError("Unsupported format: xyz")

        mol = openbabel.OBMol()
        if not ob_conv.ReadString(mol, xyz_str):
            raise ValueError("Failed to parse XYZ string")

        n = mol.NumAtoms()
        atoms = []
        for i in range(1, n + 1):
            atom = mol.GetAtom(i)
            atoms.append(openbabel.GetSymbol(atom.GetAtomicNum()))

        adjacency = np.zeros((n, n), dtype=bool)
        for bond in openbabel.OBMolBondIter(mol):
            i = bond.GetBeginAtomIdx() - 1
            j = bond.GetEndAtomIdx() - 1
            adjacency[i, j] = True
            adjacency[j, i] = True

        return atoms, adjacency

    def _compute_encodings(self):
        """Compute all three levels of GCN encodings"""
        self._compute_gcn0()
        self._compute_gcn1()
        self._compute_gcn2()

    @property
    def encodings(self):
        """Return dictionary containing all GCN encodings"""
        return {"gcn0": self.gcn0, "gcn1": self.gcn1, "gcn2": self.gcn2}

    def _compute_gcn0(self):
        """Compute GCN0 encodings based on atom type and neighbor count"""
        self.gcn0 = []
        for i, symbol in enumerate(self.atoms):
            if symbol == "H":
                self.gcn0.append("")
            else:
                neighbor_symbols = [
                    self.atoms[j] for j, connected in enumerate(self.adjacency[i]) if connected
                ]
                num_h = neighbor_symbols.count("H")
                num_heavy = len(neighbor_symbols) - num_h
                self.gcn0.append(f"{symbol}{num_heavy}{num_h}")

    def _compute_gcn1(self):
        """Compute GCN1 encodings including first-neighbor information"""
        self.gcn1 = []
        for i, symbol in enumerate(self.atoms):
            if symbol == "H":
                self.gcn1.append("")
            else:
                neighbor_indices = np.where(self.adjacency[i])[0]
                neighbor_encodings = [self.gcn0[j] for j in neighbor_indices if self.gcn0[j] != ""]
                if neighbor_encodings:
                    neighbor_encodings.sort(reverse=True)
                    encoding = f"{self.gcn0[i]}({','.join(neighbor_encodings)})"
                else:
                    encoding = self.gcn0[i]
                self.gcn1.append(encoding)

    def _compute_gcn2(self):
        """Compute GCN2 encodings including second-neighbor information"""
        self.gcn2 = []
        for i, symbol in enumerate(self.atoms):
            if symbol == "H":
                self.gcn2.append("")
            else:
                neighbor_indices = np.where(self.adjacency[i])[0]
                neighbor_encodings = [self.gcn1[j] for j in neighbor_indices if self.gcn1[j] != ""]
                if neighbor_encodings:
                    neighbor_encodings.sort(reverse=True)
                    encoding = f"{self.gcn0[i]}[{','.join(neighbor_encodings)}]"
                else:
                    encoding = self.gcn0[i]
                self.gcn2.append(encoding)

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
