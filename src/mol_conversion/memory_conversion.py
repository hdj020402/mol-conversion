import numpy as np
from openbabel import pybel, openbabel
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

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
    def _infer_formal_charges(xyz_str: str, total_charge: int) -> dict[int, int]:
        """
        Infer per-atom formal charges from total charge using RDKit's xyz2mol
        algorithm (``rdDetermineBonds.DetermineBonds``).

        ``embedChiral=False`` is passed to skip the downstream ``sanitizeMol``
        call: this helper only reads formal charges, which are assigned during
        bond-order determination (before sanitize). Skipping sanitize avoids
        unnecessary work and a class of failures (aromaticity / Kekule /
        valence checks) that would otherwise abort the call after the charges
        we want are already in place.

        Args:
            xyz_str: XYZ format molecular structure string. Multi-frame XYZ
                trajectories are accepted; only the first frame is used (this
                matches Open Babel's default behavior elsewhere in this class).
            total_charge: Total formal charge of the molecule.

        Returns:
            dict mapping 0-based atom index to formal charge, excluding atoms
            whose formal charge is zero. Atom indices match the order in the
            XYZ block (the same order Open Babel uses when parsing it).

        Raises:
            ValueError: if RDKit cannot find a bond ordering consistent with
                the requested total charge (raised by ``DetermineBonds`` and
                propagated unchanged).
        """
        lines = xyz_str.splitlines()
        n_atoms = int(lines[0].strip())
        first_frame = "\n".join(lines[: n_atoms + 2]) + "\n"
        mol = Chem.MolFromXYZBlock(first_frame)
        rdDetermineBonds.DetermineBonds(mol, charge=total_charge, embedChiral=False)
        return {
            a.GetIdx(): a.GetFormalCharge()
            for a in mol.GetAtoms()
            if a.GetFormalCharge() != 0
        }

    @staticmethod
    def _xyz_to_pybel_mol(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
        total_charge: int | None = None,
    ) -> pybel.Molecule:
        """
        Parse an XYZ string into a ``pybel.Molecule``, optionally setting
        formal charges on individual atoms (and adjusting total charge).

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional dict mapping 0-based atom index to formal
                charge. Example: ``{2: 1}`` means atom at index 2 has charge +1.
            total_charge: Optional molecule-level total charge. When given
                alone, the per-atom assignment is inferred via RDKit's
                xyz2mol algorithm. When given together with ``formal_charges``,
                the two must agree (``sum(formal_charges.values()) ==
                total_charge``); otherwise a ``ValueError`` is raised.

        Returns:
            pybel.Molecule with charges applied to its OBMol if provided.
        """
        if formal_charges is not None and total_charge is not None:
            observed = sum(formal_charges.values())
            if observed != total_charge:
                raise ValueError(
                    f"total_charge={total_charge} contradicts "
                    f"sum(formal_charges)={observed}"
                )
        elif formal_charges is None and total_charge is not None:
            formal_charges = MemoryConverter._infer_formal_charges(xyz_str, total_charge)

        mol = pybel.readstring("xyz", xyz_str)
        if formal_charges:
            obmol = mol.OBMol
            for atom_idx, charge in formal_charges.items():
                obmol.GetAtom(atom_idx + 1).SetFormalCharge(charge)  # OB uses 1-based index
            obmol.SetTotalCharge(sum(formal_charges.values()))
        return mol

    # ------------------------------------------------------------------
    # XYZ -> other formats
    # XYZ does not carry charge information; pass either ``total_charge``
    # (to let RDKit's xyz2mol infer the per-atom assignment) or
    # ``formal_charges`` (to specify it explicitly) for charged molecules.
    # ------------------------------------------------------------------
    @staticmethod
    def xyz_to_mol_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
        total_charge: int | None = None,
    ) -> str:
        """
        Convert XYZ format string to MOL block string (memory-to-memory).

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                When provided, M CHG lines reflecting the charges are
                emitted in the MOL block.
            total_charge: Optional molecule total charge; when given alone,
                per-atom charges are inferred via RDKit's xyz2mol. See
                ``_xyz_to_pybel_mol`` for the full semantics.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge)
        return mol.write("mol")

    @staticmethod
    def xyz_to_rdkit_mol(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
        total_charge: int | None = None,
    ) -> Chem.Mol:
        """
        Convert XYZ string to RDKit Mol object (memory-only, no file I/O).

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                When provided, charges are baked into the intermediate MOL
                block as ``M CHG`` lines, and RDKit reads them natively
                during sanitization.
            total_charge: Optional molecule total charge; see
                ``_xyz_to_pybel_mol`` for the full semantics.
        """
        molblock = MemoryConverter.xyz_to_mol_string(xyz_str, formal_charges, total_charge)
        return MemoryConverter.molblock_to_rdkit_mol(molblock)

    @staticmethod
    def xyz_to_bond_order_matrix(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
        total_charge: int | None = None,
    ) -> np.ndarray:
        """
        Build bond order matrix from XYZ string using Open Babel.

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                Affects Open Babel's perceived bond orders for charged
                atoms (e.g. distinguishing ``[NH3+]`` from radical ``NH3``).
            total_charge: Optional molecule total charge; see
                ``_xyz_to_pybel_mol`` for the full semantics.

        Returns:
            np.ndarray of shape ``(n_atoms, n_atoms)`` with integer bond orders.
        """
        ob_mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge).OBMol
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
        total_charge: int | None = None,
    ) -> str:
        """
        Convert XYZ string to InChI format (memory conversion).

        Args:
            xyz_str: XYZ format molecular structure string.
            fixed_h: Whether to generate InChI with fixed hydrogen layer.
            formal_charges: Optional 0-based atom-index -> charge dict.
                Required to obtain a correct InChI for charged molecules
                since XYZ alone does not carry charge information.
            total_charge: Optional molecule total charge; see
                ``_xyz_to_pybel_mol`` for the full semantics.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge)
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
        total_charge: int | None = None,
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
            total_charge: Optional molecule total charge; see
                ``_xyz_to_pybel_mol`` for the full semantics.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge)
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
        total_charge: int | None = None,
    ) -> str:
        """
        Convert XYZ string to canonical SMILES.

        Args:
            xyz_str: XYZ format molecular structure string.
            formal_charges: Optional 0-based atom-index -> charge dict.
                Required to produce charge-aware SMILES (e.g. ``[NH3+]``)
                for charged molecules.
            total_charge: Optional molecule total charge; see
                ``_xyz_to_pybel_mol`` for the full semantics.
        """
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge)
        # mol.write("can") returns "<canonical_smiles>\t<title>\n";
        # split()[0] robustly extracts the SMILES regardless of the title.
        return mol.write("can").split()[0]

    @staticmethod
    def xyz_to_sdf_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
        total_charge: int | None = None,
    ) -> str:
        """Convert XYZ string to SDF format string."""
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge)
        return mol.write("sdf")

    @staticmethod
    def xyz_to_pdb_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
        total_charge: int | None = None,
    ) -> str:
        """Convert XYZ string to PDB format string."""
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge)
        return mol.write("pdb")

    @staticmethod
    def xyz_to_mol2_string(
        xyz_str: str,
        formal_charges: dict[int, int] | None = None,
        total_charge: int | None = None,
    ) -> str:
        """Convert XYZ string to MOL2 format string."""
        mol = MemoryConverter._xyz_to_pybel_mol(xyz_str, formal_charges, total_charge)
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
