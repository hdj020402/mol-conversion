"""Tests for memory conversion module"""

import pytest
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from mol_conversion import MemoryConverter


class TestMemoryConverter:
    """Test class for MemoryConverter methods"""
    
    @pytest.fixture
    def test_xyz_string(self):
        """Real XYZ string from test data file"""
        test_file = Path(__file__).parent / "data" / "test.xyz"
        if test_file.exists():
            return test_file.read_text()
        else:
            pytest.skip(f"Test file {test_file} not found")
    
    def test_xyz_to_mol_string(self, test_xyz_string):
        """Test XYZ to MOL string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        mol_string = MemoryConverter.xyz_to_mol_string(test_xyz_string)
        assert isinstance(mol_string, str)
        assert len(mol_string) > 0
        # Basic MOL file validation
        assert "V2000" in mol_string or "M  END" in mol_string
    
    def test_xyz_to_rdkit_mol(self, test_xyz_string):
        """Test XYZ to RDKit Mol object conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        # Skip test if RDKit is not available
        if not self._is_rdkit_available():
            pytest.skip("RDKit not available")
        
        mol = MemoryConverter.xyz_to_rdkit_mol(test_xyz_string)
        assert mol is not None
        # Real file has 70 atoms
        assert mol.GetNumAtoms() == 70
        assert mol.GetNumBonds() > 0
    
    def test_xyz_to_rdkit_mol_invalid_xyz(self):
        """Test XYZ to RDKit Mol object conversion with invalid XYZ"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        # Skip test if RDKit is not available
        if not self._is_rdkit_available():
            pytest.skip("RDKit not available")
        
        invalid_xyz = """2
Invalid molecule
X    0.000000    0.000000    0.000000
Y    1.000000    0.000000    0.000000
"""
        
        # Should not raise an exception, but might return None
        mol = MemoryConverter.xyz_to_rdkit_mol(invalid_xyz)
        # The behavior depends on Open Babel and RDKit versions
        # Either it returns None or creates a molecule with unknown atoms
        assert mol is None or mol.GetNumAtoms() == 2
    
    def test_xyz_to_bond_order_matrix(self, test_xyz_string):
        """Test XYZ to bond order matrix conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        bond_matrix = MemoryConverter.xyz_to_bond_order_matrix(test_xyz_string)
        assert isinstance(bond_matrix, np.ndarray)
        # Real file has 70 atoms
        assert bond_matrix.shape == (70, 70)
        
        # Check diagonal is zero
        assert np.all(np.diag(bond_matrix) == 0)
        
        # Check matrix is symmetric
        assert np.allclose(bond_matrix, bond_matrix.T)
        
        # Should have some bonds (non-zero entries)
        assert np.sum(bond_matrix > 0) > 0
    
    def test_xyz_to_inchi_string(self, test_xyz_string):
        """Test XYZ to InChI string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        inchi = MemoryConverter.xyz_to_inchi_string(test_xyz_string)
        inchi_fixedh = MemoryConverter.xyz_to_inchi_string(test_xyz_string, fixed_h=True)
        assert isinstance(inchi, str)
        assert len(inchi) > 0
        assert inchi.startswith("InChI=1S/")
        assert isinstance(inchi_fixedh, str)
        assert len(inchi_fixedh) > 0
        assert inchi_fixedh.startswith("InChI=1/")
    
    def test_xyz_to_inchikey_string(self, test_xyz_string):
        """Test XYZ to InChIKey string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        inchikey = MemoryConverter.xyz_to_inchikey_string(test_xyz_string)
        assert isinstance(inchikey, str)
        assert len(inchikey) == 27  # Standard InChIKey length
    
    def test_xyz_to_smiles_string(self, test_xyz_string):
        """Test XYZ to SMILES string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        smiles = MemoryConverter.xyz_to_smiles_string(test_xyz_string)
        assert isinstance(smiles, str)
        assert len(smiles) > 0
    
    def test_xyz_to_sdf_string(self, test_xyz_string):
        """Test XYZ to SDF string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        sdf_string = MemoryConverter.xyz_to_sdf_string(test_xyz_string)
        assert isinstance(sdf_string, str)
        assert len(sdf_string) > 0
        # Basic SDF validation
        assert "V2000" in sdf_string or "M  END" in sdf_string
    
    def test_xyz_to_pdb_string(self, test_xyz_string):
        """Test XYZ to PDB string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        pdb_string = MemoryConverter.xyz_to_pdb_string(test_xyz_string)
        assert isinstance(pdb_string, str)
        assert len(pdb_string) > 0
        # Basic PDB validation
        assert "ATOM" in pdb_string or "HETATM" in pdb_string
    
    def test_xyz_to_mol2_string(self, test_xyz_string):
        """Test XYZ to MOL2 string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        mol2_string = MemoryConverter.xyz_to_mol2_string(test_xyz_string)
        assert isinstance(mol2_string, str)
        assert len(mol2_string) > 0
        # Basic MOL2 validation
        assert "@<TRIPOS>" in mol2_string
    
    def test_xyz_to_cif_string(self, test_xyz_string):
        """Test XYZ to CIF string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        cif_string = MemoryConverter.xyz_to_cif_string(test_xyz_string)
        assert isinstance(cif_string, str)
        assert len(cif_string) > 0
        # Basic CIF validation
        assert "data_" in cif_string
    
    def test_all_conversion_methods(self, test_xyz_string):
        """Test all conversion methods"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        # Test all conversion methods
        mol_string = MemoryConverter.xyz_to_mol_string(test_xyz_string)
        assert isinstance(mol_string, str) and len(mol_string) > 0
        
        bond_matrix = MemoryConverter.xyz_to_bond_order_matrix(test_xyz_string)
        assert isinstance(bond_matrix, np.ndarray) and bond_matrix.size > 0
        
        inchi = MemoryConverter.xyz_to_inchi_string(test_xyz_string)
        assert isinstance(inchi, str) and len(inchi) > 0
        
        inchikey = MemoryConverter.xyz_to_inchikey_string(test_xyz_string)
        assert isinstance(inchikey, str) and len(inchikey) == 27
        
        smiles = MemoryConverter.xyz_to_smiles_string(test_xyz_string)
        assert isinstance(smiles, str) and len(smiles) > 0
        
        sdf_string = MemoryConverter.xyz_to_sdf_string(test_xyz_string)
        assert isinstance(sdf_string, str) and len(sdf_string) > 0
        
        pdb_string = MemoryConverter.xyz_to_pdb_string(test_xyz_string)
        assert isinstance(pdb_string, str) and len(pdb_string) > 0
        
        mol2_string = MemoryConverter.xyz_to_mol2_string(test_xyz_string)
        assert isinstance(mol2_string, str) and len(mol2_string) > 0
        
        cif_string = MemoryConverter.xyz_to_cif_string(test_xyz_string)
        assert isinstance(cif_string, str) and len(cif_string) > 0
    
    def test_invalid_xyz_string(self):
        """Test conversion methods with invalid XYZ string"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        invalid_xyz = """1
Invalid atom
X    0.000000    0.000000    0.000000
"""

        # Should not raise exceptions for invalid input
        # The behavior depends on Open Babel's handling of unknown atoms
        try:
            mol_string = MemoryConverter.xyz_to_mol_string(invalid_xyz)
            assert isinstance(mol_string, str)
        except Exception:
            # It's acceptable if Open Babel raises an exception for invalid input
            pass

    # ------------------------------------------------------------------
    # Charged-molecule fixtures (cation + anion)
    # XYZ generated via RDKit EmbedMolecule with fixed seed for reproducibility.
    # ------------------------------------------------------------------
    @pytest.fixture
    def cation_data(self):
        """Ethylammonium [CC-NH3]+ as (xyz_str, formal_charges_dict)."""
        if not self._is_openbabel_available() or not self._is_rdkit_available():
            pytest.skip("Open Babel or RDKit not available")
        mol = Chem.AddHs(Chem.MolFromSmiles("CC[NH3+]"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        xyz = Chem.MolToXYZBlock(mol)
        n_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N")
        return xyz, {n_idx: 1}

    @pytest.fixture
    def anion_data(self):
        """Ethoxide CH3-CH2-O- as (xyz_str, formal_charges_dict).

        Chosen over carboxylates (e.g. acetate) because the C-O bond is
        unambiguous: there is no resonance for OB's geometry-based bond
        perception to misassign.
        """
        if not self._is_openbabel_available() or not self._is_rdkit_available():
            pytest.skip("Open Babel or RDKit not available")
        mol = Chem.AddHs(Chem.MolFromSmiles("CC[O-]"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        xyz = Chem.MolToXYZBlock(mol)
        o_idx = next(a.GetIdx() for a in mol.GetAtoms()
                     if a.GetSymbol() == "O" and a.GetFormalCharge() == -1)
        return xyz, {o_idx: -1}

    # ------------------------------------------------------------------
    # Group A: XYZ + formal_charges
    # ------------------------------------------------------------------
    def test_xyz_to_inchi_string_cation(self, cation_data):
        xyz, charges = cation_data
        inchi = MemoryConverter.xyz_to_inchi_string(xyz, formal_charges=charges)
        assert inchi == "InChI=1S/C2H7N/c1-2-3/h2-3H2,1H3/p+1"

    def test_xyz_to_inchi_string_cation_fixed_h(self, cation_data):
        xyz, charges = cation_data
        inchi = MemoryConverter.xyz_to_inchi_string(xyz, fixed_h=True, formal_charges=charges)
        assert inchi == "InChI=1/C2H7N/c1-2-3/h2-3H2,1H3/p+1/fC2H8N/h3H/q+1"

    def test_xyz_to_inchikey_string_cation(self, cation_data):
        xyz, charges = cation_data
        key = MemoryConverter.xyz_to_inchikey_string(xyz, formal_charges=charges)
        assert key == "QUSNBJAOOMFDIB-UHFFFAOYSA-O"

    def test_xyz_to_smiles_string_cation(self, cation_data):
        """Cation SMILES is canonical and exactly [NH3+]."""
        xyz, charges = cation_data
        smiles = MemoryConverter.xyz_to_smiles_string(xyz, formal_charges=charges)
        assert smiles == "CC[NH3+]"

    def test_xyz_to_rdkit_mol_cation(self, cation_data):
        xyz, charges = cation_data
        mol = MemoryConverter.xyz_to_rdkit_mol(xyz, formal_charges=charges)
        assert mol is not None
        assert Chem.GetFormalCharge(mol) == 1

    def test_xyz_to_mol_string_cation(self, cation_data):
        """MOL block must carry the charge as an M CHG line."""
        xyz, charges = cation_data
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        assert "M  CHG" in molblock

    def test_xyz_to_inchi_string_anion(self, anion_data):
        xyz, charges = anion_data
        inchi = MemoryConverter.xyz_to_inchi_string(xyz, formal_charges=charges)
        assert inchi == "InChI=1S/C2H5O/c1-2-3/h2H2,1H3/q-1"

    def test_xyz_to_inchi_string_anion_fixed_h(self, anion_data):
        xyz, charges = anion_data
        inchi = MemoryConverter.xyz_to_inchi_string(xyz, fixed_h=True, formal_charges=charges)
        assert inchi == "InChI=1/C2H5O/c1-2-3/h2H2,1H3/q-1"

    def test_xyz_to_inchikey_string_anion(self, anion_data):
        xyz, charges = anion_data
        key = MemoryConverter.xyz_to_inchikey_string(xyz, formal_charges=charges)
        assert key == "HHFAWKCIHAUFRX-UHFFFAOYSA-N"

    def test_xyz_to_smiles_string_anion(self, anion_data):
        xyz, charges = anion_data
        smiles = MemoryConverter.xyz_to_smiles_string(xyz, formal_charges=charges)
        assert smiles == "[O-]CC"

    def test_xyz_to_rdkit_mol_anion(self, anion_data):
        xyz, charges = anion_data
        mol = MemoryConverter.xyz_to_rdkit_mol(xyz, formal_charges=charges)
        assert mol is not None
        assert Chem.GetFormalCharge(mol) == -1

    # ------------------------------------------------------------------
    # Group B: molblock_to_* (new methods, cation only)
    # ------------------------------------------------------------------
    def test_molblock_to_inchi_string(self, cation_data):
        xyz, charges = cation_data
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        inchi = MemoryConverter.molblock_to_inchi_string(molblock)
        assert inchi == "InChI=1S/C2H7N/c1-2-3/h2-3H2,1H3/p+1"

    def test_molblock_to_inchikey_string(self, cation_data):
        xyz, charges = cation_data
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        assert MemoryConverter.molblock_to_inchikey_string(molblock) == "QUSNBJAOOMFDIB-UHFFFAOYSA-O"

    def test_molblock_to_smiles_string(self, cation_data):
        xyz, charges = cation_data
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        smiles = MemoryConverter.molblock_to_smiles_string(molblock)
        rdmol = Chem.MolFromSmiles(smiles)
        assert rdmol is not None and Chem.GetFormalCharge(rdmol) == 1

    def test_molblock_to_rdkit_mol(self, cation_data):
        xyz, charges = cation_data
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        mol = MemoryConverter.molblock_to_rdkit_mol(molblock)
        assert mol is not None
        assert Chem.GetFormalCharge(mol) == 1

    def test_molblock_to_xyz_string(self, cation_data):
        xyz, charges = cation_data
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        xyz_back = MemoryConverter.molblock_to_xyz_string(molblock)
        assert int(xyz_back.strip().split('\n')[0]) == 11

    # ------------------------------------------------------------------
    # Group C: cross-path consistency (the original bug's regression test)
    # ------------------------------------------------------------------
    def test_xyz_path_matches_molblock_path_cation(self, cation_data):
        xyz, charges = cation_data
        via_xyz = MemoryConverter.xyz_to_inchi_string(xyz, fixed_h=True, formal_charges=charges)
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        via_mol = MemoryConverter.molblock_to_inchi_string(molblock, fixed_h=True)
        assert via_xyz == via_mol

    def test_xyz_path_matches_molblock_path_anion(self, anion_data):
        xyz, charges = anion_data
        via_xyz = MemoryConverter.xyz_to_inchi_string(xyz, fixed_h=True, formal_charges=charges)
        molblock = MemoryConverter.xyz_to_mol_string(xyz, formal_charges=charges)
        via_mol = MemoryConverter.molblock_to_inchi_string(molblock, fixed_h=True)
        assert via_xyz == via_mol

    # ------------------------------------------------------------------
    # Group E: XYZ + total_charge (RDKit xyz2mol inference)
    # ------------------------------------------------------------------
    @pytest.fixture
    def zwitterion_data(self):
        """Glycine zwitterion [NH3+]CC(=O)[O-] as (xyz_str, expected_per_atom_dict).

        Total charge is 0 but two atoms carry non-zero formal charges; this is
        the case where total_charge mode shines, since specifying the per-atom
        dict requires the user to know the N and O indices ahead of time.
        """
        if not self._is_openbabel_available() or not self._is_rdkit_available():
            pytest.skip("Open Babel or RDKit not available")
        mol = Chem.AddHs(Chem.MolFromSmiles("[NH3+]CC(=O)[O-]"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        xyz = Chem.MolToXYZBlock(mol)
        expected = {
            a.GetIdx(): a.GetFormalCharge()
            for a in mol.GetAtoms() if a.GetFormalCharge() != 0
        }
        return xyz, expected

    def test_total_charge_cation_matches_formal_charges(self, cation_data):
        """xyz_to_inchi(total_charge=+1) should match the formal_charges path."""
        xyz, charges = cation_data
        via_total = MemoryConverter.xyz_to_inchi_string(xyz, total_charge=1)
        via_dict = MemoryConverter.xyz_to_inchi_string(xyz, formal_charges=charges)
        assert via_total == via_dict == "InChI=1S/C2H7N/c1-2-3/h2-3H2,1H3/p+1"

    def test_total_charge_anion_matches_formal_charges(self, anion_data):
        xyz, charges = anion_data
        via_total = MemoryConverter.xyz_to_inchi_string(xyz, total_charge=-1)
        via_dict = MemoryConverter.xyz_to_inchi_string(xyz, formal_charges=charges)
        assert via_total == via_dict == "InChI=1S/C2H5O/c1-2-3/h2H2,1H3/q-1"

    def test_total_charge_zwitterion(self, zwitterion_data):
        """Glycine zwitterion: total=0 but per-atom non-zero. Inference must
        find both the +1 N and -1 O without the caller specifying indices."""
        xyz, expected = zwitterion_data
        inferred = MemoryConverter._infer_formal_charges(xyz, total_charge=0)
        assert inferred == expected
        # Round-trip through OB: SMILES should carry both [N+] and [O-]
        smiles = MemoryConverter.xyz_to_smiles_string(xyz, total_charge=0)
        assert "[NH3+]" in smiles and "[O-]" in smiles

    def test_total_charge_smiles_cation(self, cation_data):
        xyz, _ = cation_data
        smiles = MemoryConverter.xyz_to_smiles_string(xyz, total_charge=1)
        assert smiles == "CC[NH3+]"

    def test_total_charge_rdkit_mol_cation(self, cation_data):
        xyz, _ = cation_data
        mol = MemoryConverter.xyz_to_rdkit_mol(xyz, total_charge=1)
        assert Chem.GetFormalCharge(mol) == 1

    def test_total_charge_and_formal_charges_consistent(self, cation_data):
        """Passing both with matching sum should succeed (no spurious error)."""
        xyz, charges = cation_data
        inchi = MemoryConverter.xyz_to_inchi_string(
            xyz, formal_charges=charges, total_charge=1
        )
        assert inchi == "InChI=1S/C2H7N/c1-2-3/h2-3H2,1H3/p+1"

    def test_total_charge_and_formal_charges_contradiction_raises(self, cation_data):
        """sum(formal_charges) != total_charge is a genuine input error."""
        xyz, charges = cation_data  # sum == +1
        with pytest.raises(ValueError, match="contradicts"):
            MemoryConverter.xyz_to_inchi_string(
                xyz, formal_charges=charges, total_charge=2
            )

    def test_wrong_total_charge_raises(self, cation_data):
        """RDKit's DetermineBonds raises when total_charge is inconsistent
        with what the geometry can support; we propagate that unchanged."""
        xyz, _ = cation_data
        with pytest.raises(ValueError):
            MemoryConverter.xyz_to_inchi_string(xyz, total_charge=0)

    def test_total_charge_zero_neutral_equiv_no_charge(self, test_xyz_string):
        """For a neutral molecule, total_charge=0 should yield the same
        InChI as not specifying any charge info."""
        if not self._is_openbabel_available() or not self._is_rdkit_available():
            pytest.skip("Open Babel or RDKit not available")
        inchi_none = MemoryConverter.xyz_to_inchi_string(test_xyz_string)
        inchi_zero = MemoryConverter.xyz_to_inchi_string(test_xyz_string, total_charge=0)
        assert inchi_none == inchi_zero

    # ------------------------------------------------------------------
    # Group D: edge cases / backward compat
    # ------------------------------------------------------------------
    def test_empty_formal_charges_equiv_none(self, test_xyz_string):
        """formal_charges={} should produce the same output as formal_charges=None."""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        inchi_none = MemoryConverter.xyz_to_inchi_string(test_xyz_string)
        inchi_empty = MemoryConverter.xyz_to_inchi_string(test_xyz_string, formal_charges={})
        assert inchi_none == inchi_empty

    def test_neutral_methods_backward_compat(self, test_xyz_string):
        """All XYZ-source methods callable with positional XYZ arg only."""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        MemoryConverter.xyz_to_mol_string(test_xyz_string)
        MemoryConverter.xyz_to_rdkit_mol(test_xyz_string)
        MemoryConverter.xyz_to_inchi_string(test_xyz_string)
        MemoryConverter.xyz_to_inchikey_string(test_xyz_string)
        MemoryConverter.xyz_to_smiles_string(test_xyz_string)
        MemoryConverter.xyz_to_bond_order_matrix(test_xyz_string)
        MemoryConverter.xyz_to_sdf_string(test_xyz_string)
        MemoryConverter.xyz_to_pdb_string(test_xyz_string)
        MemoryConverter.xyz_to_mol2_string(test_xyz_string)

    def _is_openbabel_available(self):
        """Check if Open Babel Python bindings are available"""
        try:
            from openbabel import pybel
            return True
        except ImportError:
            return False
    
    def _is_rdkit_available(self):
        """Check if RDKit is available"""
        try:
            from rdkit import Chem
            return True
        except ImportError:
            return False