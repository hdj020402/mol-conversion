from __future__ import annotations

import numpy as np
from openbabel import openbabel, pybel


class GCNEncoder:
    """
    Generate atom-level GCN-style encodings (GCN0, GCN1, GCN2) from molecular structures.
    Hydrogen atoms are encoded as empty strings ('').

    Use from_* class methods to create instances from various molecular formats:
        - from_xyz(xyz_str)
        - from_smiles(smiles_str)
        - from_inchi(inchi_str)
        - from_sdf(sdf_str)
        - from_mol(mol_str)
        - from_mol2(mol2_str)
        - from_pdb(pdb_str)
        - from_cif(cif_str)
    """

    def __init__(self, atoms: list, adjacency: np.ndarray):
        """
        Internal constructor. Use from_* class methods to create instances.

        Args:
            atoms: List of atom symbols.
            adjacency: Boolean adjacency matrix with shape (n_atoms, n_atoms).
        """
        self.atoms = atoms
        self.adjacency = adjacency
        self._compute_encodings()

    @classmethod
    def from_xyz(cls, xyz_str: str) -> GCNEncoder:
        """Create GCNEncoder from XYZ format string."""
        return cls._from_format(xyz_str, "xyz")

    @classmethod
    def from_smiles(cls, smiles_str: str) -> GCNEncoder:
        """Create GCNEncoder from SMILES format string."""
        return cls._from_format(smiles_str, "smi")

    @classmethod
    def from_inchi(cls, inchi_str: str) -> GCNEncoder:
        """Create GCNEncoder from InChI format string."""
        return cls._from_format(inchi_str, "inchi")

    @classmethod
    def from_sdf(cls, sdf_str: str) -> GCNEncoder:
        """Create GCNEncoder from SDF format string."""
        return cls._from_format(sdf_str, "sdf")

    @classmethod
    def from_mol(cls, mol_str: str) -> GCNEncoder:
        """Create GCNEncoder from MOL format string."""
        return cls._from_format(mol_str, "mol")

    @classmethod
    def from_mol2(cls, mol2_str: str) -> GCNEncoder:
        """Create GCNEncoder from MOL2 format string."""
        return cls._from_format(mol2_str, "mol2")

    @classmethod
    def from_pdb(cls, pdb_str: str) -> GCNEncoder:
        """Create GCNEncoder from PDB format string."""
        return cls._from_format(pdb_str, "pdb")

    @classmethod
    def from_cif(cls, cif_str: str) -> GCNEncoder:
        """Create GCNEncoder from CIF format string."""
        return cls._from_format(cif_str, "cif")

    @classmethod
    def _from_format(cls, mol_str: str, fmt: str) -> GCNEncoder:
        """
        Internal method: parse molecular string and create GCNEncoder instance.

        Args:
            mol_str: Molecular structure string.
            fmt: Open Babel format identifier (e.g. "xyz", "smi", "inchi").
        """
        mol = pybel.readstring(fmt, mol_str)
        atoms, adjacency = cls._extract_atoms_and_adjacency(mol.OBMol)
        return cls(atoms, adjacency)

    @staticmethod
    def _extract_atoms_and_adjacency(ob_mol: openbabel.OBMol):
        """
        Extract atom symbols and boolean adjacency matrix from an OBMol object.

        Args:
            ob_mol: Open Babel OBMol object.

        Returns:
            tuple: (atoms list, adjacency ndarray)
        """
        n = ob_mol.NumAtoms()
        atoms = []
        for i in range(1, n + 1):
            atom = ob_mol.GetAtom(i)
            atoms.append(openbabel.GetSymbol(atom.GetAtomicNum()))

        adjacency = np.zeros((n, n), dtype=bool)
        for bond in openbabel.OBMolBondIter(ob_mol):
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
