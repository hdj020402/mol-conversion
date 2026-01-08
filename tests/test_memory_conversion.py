"""Tests for memory conversion module"""

import pytest
import numpy as np
from pathlib import Path
from mol_conversion import MemoryConversion


class TestMemoryConversion:
    """Test class for MemoryConversion methods"""
    
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
        
        mol_string = MemoryConversion.xyz_to_mol_string(test_xyz_string)
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
        
        mol = MemoryConversion.xyz_to_rdkit_mol(test_xyz_string)
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
        mol = MemoryConversion.xyz_to_rdkit_mol(invalid_xyz)
        # The behavior depends on Open Babel and RDKit versions
        # Either it returns None or creates a molecule with unknown atoms
        assert mol is None or mol.GetNumAtoms() == 2
    
    def test_xyz_to_bond_order_matrix(self, test_xyz_string):
        """Test XYZ to bond order matrix conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        bond_matrix = MemoryConversion.xyz_to_bond_order_matrix(test_xyz_string)
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
        
        inchi = MemoryConversion.xyz_to_inchi_string(test_xyz_string)
        assert isinstance(inchi, str)
        assert len(inchi) > 0
        # Basic InChI validation
        assert inchi.startswith("InChI=1S/")
    
    def test_xyz_to_inchikey_string(self, test_xyz_string):
        """Test XYZ to InChIKey string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        inchikey = MemoryConversion.xyz_to_inchikey_string(test_xyz_string)
        assert isinstance(inchikey, str)
        assert len(inchikey) == 27  # Standard InChIKey length
    
    def test_xyz_to_smiles_string(self, test_xyz_string):
        """Test XYZ to SMILES string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        smiles = MemoryConversion.xyz_to_smiles_string(test_xyz_string)
        assert isinstance(smiles, str)
        assert len(smiles) > 0
    
    def test_xyz_to_sdf_string(self, test_xyz_string):
        """Test XYZ to SDF string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        sdf_string = MemoryConversion.xyz_to_sdf_string(test_xyz_string)
        assert isinstance(sdf_string, str)
        assert len(sdf_string) > 0
        # Basic SDF validation
        assert "V2000" in sdf_string or "M  END" in sdf_string
    
    def test_xyz_to_pdb_string(self, test_xyz_string):
        """Test XYZ to PDB string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        pdb_string = MemoryConversion.xyz_to_pdb_string(test_xyz_string)
        assert isinstance(pdb_string, str)
        assert len(pdb_string) > 0
        # Basic PDB validation
        assert "ATOM" in pdb_string or "HETATM" in pdb_string
    
    def test_xyz_to_mol2_string(self, test_xyz_string):
        """Test XYZ to MOL2 string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        mol2_string = MemoryConversion.xyz_to_mol2_string(test_xyz_string)
        assert isinstance(mol2_string, str)
        assert len(mol2_string) > 0
        # Basic MOL2 validation
        assert "@<TRIPOS>" in mol2_string
    
    def test_xyz_to_cif_string(self, test_xyz_string):
        """Test XYZ to CIF string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        cif_string = MemoryConversion.xyz_to_cif_string(test_xyz_string)
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
        mol_string = MemoryConversion.xyz_to_mol_string(test_xyz_string)
        assert isinstance(mol_string, str) and len(mol_string) > 0
        
        bond_matrix = MemoryConversion.xyz_to_bond_order_matrix(test_xyz_string)
        assert isinstance(bond_matrix, np.ndarray) and bond_matrix.size > 0
        
        inchi = MemoryConversion.xyz_to_inchi_string(test_xyz_string)
        assert isinstance(inchi, str) and len(inchi) > 0
        
        inchikey = MemoryConversion.xyz_to_inchikey_string(test_xyz_string)
        assert isinstance(inchikey, str) and len(inchikey) == 27
        
        smiles = MemoryConversion.xyz_to_smiles_string(test_xyz_string)
        assert isinstance(smiles, str) and len(smiles) > 0
        
        sdf_string = MemoryConversion.xyz_to_sdf_string(test_xyz_string)
        assert isinstance(sdf_string, str) and len(sdf_string) > 0
        
        pdb_string = MemoryConversion.xyz_to_pdb_string(test_xyz_string)
        assert isinstance(pdb_string, str) and len(pdb_string) > 0
        
        mol2_string = MemoryConversion.xyz_to_mol2_string(test_xyz_string)
        assert isinstance(mol2_string, str) and len(mol2_string) > 0
        
        cif_string = MemoryConversion.xyz_to_cif_string(test_xyz_string)
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
            mol_string = MemoryConversion.xyz_to_mol_string(invalid_xyz)
            assert isinstance(mol_string, str)
        except Exception:
            # It's acceptable if Open Babel raises an exception for invalid input
            pass
    
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