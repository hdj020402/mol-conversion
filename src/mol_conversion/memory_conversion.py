import numpy as np
from openbabel import pybel, openbabel
from rdkit import Chem

__all__ = ["MemoryConverter"]


class MemoryConverter:
    """
    Memory conversion methods collection class.

    Provides memory-to-memory molecular format conversions across:
        - XYZ (no native charge info; charged molecules need ``formal_charges``)
        - MOL block, SMILES (charge info is native to the format)
        - InChI / InChIKey
        - RDKit Mol / Open Babel pybel.Molecule
    """

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _xyz_to_pybel_mol(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> pybel.Molecule:
        """
        Parse an XYZ string into a ``pybel.Molecule``, optionally setting
        formal charges on individual atoms (and adjusting total charge).

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional dict mapping 0-based atom index to formal
                charge. Example: ``{2: 1}`` means atom at index 2 has charge +1.

        Returns:
            pybel.Molecule with charges applied to its OBMol if provided.
        """
        mol = pybel.readstring("xyz", xyz_str)
        if formal_charges:
            obmol = mol.OBMol
            for atom_idx, charge in formal_charges.items():
                obmol.GetAtom(atom_idx + 1).SetFormalCharge(charge)  # OB uses 1-based index
            obmol.SetTotalCharge(sum(formal_charges.values()))
        return mol

    # ------------------------------------------------------------------
    # XYZ -> other formats
    # XYZ does not carry charge information; pass ``formal_charges`` for
    # charged molecules.
    # ------------------------------------------------------------------
    @staticmethod
    def xyz_to_mol_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> str:
        """
        Convert XYZ format string to MOL block string (memory-to-memory).

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                When provided, M CHG lines reflecting the charges are
                emitted in the MOL block.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges)
        return mol.write("mol")

    @staticmethod
    def xyz_to_rdkit_mol(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> Chem.Mol:
        """
        Convert XYZ string to RDKit Mol object (memory-only, no file I/O).

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                When provided, charges are baked into the intermediate MOL
                block as ``M CHG`` lines, and RDKit reads them natively
                during sanitization.
        """
        molblock = MemoryConverter.xyz_to_mol_string(xyz_str, formal_charges)
        return MemoryConverter.molblock_to_rdkit_mol(molblock)

    @staticmethod
    def xyz_to_bond_order_matrix(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> np.ndarray:
        """
        Build bond order matrix from XYZ string using Open Babel.

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                Affects Open Babel's perceived bond orders for charged
                atoms (e.g. distinguishing ``[NH3+]`` from radical ``NH3``).

        Returns:
            np.ndarray of shape ``(n_atoms, n_atoms)`` with integer bond orders.
        """
        ob_mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges).OBMol
        n_atoms = ob_mol.NumAtoms()
        bond_matrix = np.zeros((n_atoms, n_atoms), dtype=int)
        for ob_bond in openbabel.OBMolBondIter(ob_mol):
            i = ob_bond.GetBeginAtomIdx() - 1  # OB index starts from 1
            j = ob_bond.GetEndAtomIdx() - 1
            bo = ob_bond.GetBondOrder()
            bond_matrix[i, j] = bo
            bond_matrix[j, i] = bo
        return bond_matrix

    @staticmethod
    def xyz_to_inchi_string(
        xyz_str: str,
        fixed_h: bool = False,
        formal_charges: dict[int, int] | None = None,
    ) -> str:
        """
        Convert XYZ string to InChI format (memory conversion).

        Args:
            xyz_str: XYZ format molecular structure string.
            fixed_h: Whether to generate InChI with fixed hydrogen layer.
            formal_charges: Optional 0-based atom-index -> charge dict.
                Required to obtain a correct InChI for charged molecules
                since XYZ alone does not carry charge information.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges)
        ob_conv = openbabel.OBConversion()
        ob_conv.SetOutFormat("inchi")
        ob_conv.AddOption("readconformer", openbabel.OBConversion.INOPTIONS)
        if fixed_h:
            ob_conv.AddOption("F", openbabel.OBConversion.OUTOPTIONS)
        return ob_conv.WriteString(mol.OBMol).strip()

    @staticmethod
    def xyz_to_inchikey_string(
        xyz_str: str,
        fixed_h: bool = False,
        formal_charges: dict[int, int] | None = None,
    ) -> str:
        """
        Convert XYZ string to InChIKey format.

        Mirrors ``xyz_to_inchi_string`` so that the resulting key is
        consistent with the InChI that would be produced under the same
        flags.

        Args:
            xyz_str: XYZ format molecular structure string.
            fixed_h: Whether to derive the key from a fixed-H InChI.
            formal_charges: Optional 0-based atom-index -> charge dict.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges)
        ob_conv = openbabel.OBConversion()
        ob_conv.SetOutFormat("inchikey")
        ob_conv.AddOption("readconformer", openbabel.OBConversion.INOPTIONS)
        if fixed_h:
            ob_conv.AddOption("F", openbabel.OBConversion.OUTOPTIONS)
        return ob_conv.WriteString(mol.OBMol).strip()

    @staticmethod
    def xyz_to_smiles_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> str:
        """
        Convert XYZ string to canonical SMILES.

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                Required to produce charge-aware SMILES (e.g. ``[NH3+]``)
                for charged molecules.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges)
        # mol.write("can") returns "<canonical_smiles>\t<title>\n";
        # split()[0] robustly extracts the SMILES regardless of the title.
        return mol.write("can").split()[0]

    @staticmethod
    def xyz_to_sdf_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> str:
        """Convert XYZ string to SDF format string."""
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges)
        return mol.write("sdf")

    @staticmethod
    def xyz_to_pdb_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> str:
        """Convert XYZ string to PDB format string."""
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges)
        return mol.write("pdb")

    @staticmethod
    def xyz_to_mol2_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
    ) -> str:
        """Convert XYZ string to MOL2 format string."""
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges)
        return mol.write("mol2")

    @staticmethod
    def xyz_to_cif_string(xyz_str: str) -> str:
        """Convert XYZ string to CIF format string."""
        mol = pybel.readstring("xyz", xyz_str)
        return mol.write("cif")

    # ------------------------------------------------------------------
    # MOL block -> other formats
    # MOL blocks carry charge info natively (M CHG lines), so these methods
    # do not accept ``formal_charges``.
    # ------------------------------------------------------------------
    @staticmethod
    def molblock_to_inchi_string(molblock_str: str, fixed_h: bool = False) -> str:
        """
        Convert MOL block string to InChI via Open Babel (memory conversion).

        Unlike ``xyz_to_inchi_string``, MOL block format carries formal charge
        information (M CHG lines), so charged molecules are handled correctly
        without explicit charge parameters.

        Args:
            molblock_str: MOL block format string.
            fixed_h: Whether to generate InChI with fixed hydrogen layer.
        """
        mol = pybel.readstring("mol", molblock_str)
        ob_conv = openbabel.OBConversion()
        ob_conv.SetOutFormat("inchi")
        if fixed_h:
            ob_conv.AddOption("F", openbabel.OBConversion.OUTOPTIONS)
        return ob_conv.WriteString(mol.OBMol).strip()

    @staticmethod
    def molblock_to_inchikey_string(molblock_str: str, fixed_h: bool = False) -> str:
        """
        Convert MOL block string to InChIKey (consistent with the matching
        ``molblock_to_inchi_string`` flag set).
        """
        mol = pybel.readstring("mol", molblock_str)
        ob_conv = openbabel.OBConversion()
        ob_conv.SetOutFormat("inchikey")
        if fixed_h:
            ob_conv.AddOption("F", openbabel.OBConversion.OUTOPTIONS)
        return ob_conv.WriteString(mol.OBMol).strip()

    @staticmethod
    def molblock_to_smiles_string(molblock_str: str) -> str:
        """Convert MOL block string to canonical SMILES."""
        mol = pybel.readstring("mol", molblock_str)
        return mol.write("can").split()[0]

    @staticmethod
    def molblock_to_rdkit_mol(
        molblock_str: str,
        removeHs: bool = False,
    ) -> Chem.Mol:
        """
        Convert MOL block string to RDKit Mol.

        Charges are read natively from M CHG lines; no extra parameter needed.

        Args:
            molblock_str: MOL block format string.
            removeHs: Whether RDKit should strip explicit hydrogens.
        """
        rdkit_mol = Chem.MolFromMolBlock(molblock_str, removeHs=removeHs, sanitize=True)
        if rdkit_mol is None:
            raise ValueError("RDKit failed to parse the MOL block.")
        return rdkit_mol

    @staticmethod
    def molblock_to_xyz_string(molblock_str: str) -> str:
        """
        Convert MOL block string to XYZ format string.

        Note: XYZ does not carry bond/charge information; only atom positions
        and elements are preserved.
        """
        mol = pybel.readstring("mol", molblock_str)
        return mol.write("xyz")
