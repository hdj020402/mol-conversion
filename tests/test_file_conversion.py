"""Tests for file conversion module"""

import pytest
import os
import tempfile
import subprocess
from pathlib import Path
from mol_conversion import FileConverter


class TestFileConverter:
    """Test class for FileConverter methods"""
    
    @pytest.fixture
    def test_xyz_file(self):
        """Path to the real test XYZ file"""
        test_dir = Path(__file__).parent / "data"
        return str(test_dir / "test.xyz")
    
    def test_xyz_to_inchi(self, test_xyz_file):
        """Test XYZ to InChI conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        inchi = FileConverter.xyz_to_inchi(test_xyz_file)
        assert isinstance(inchi, str)
        assert len(inchi) > 0
        # Basic InChI validation
        assert inchi.startswith("InChI=1S/")
    
    def test_xyz_to_inchikey(self, test_xyz_file):
        """Test XYZ to InChIKey conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        inchikey = FileConverter.xyz_to_inchikey(test_xyz_file)
        assert isinstance(inchikey, str)
        assert len(inchikey) == 27  # Standard InChIKey length
    
    def test_xyz_to_smiles(self, test_xyz_file):
        """Test XYZ to SMILES conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        smiles = FileConverter.xyz_to_smiles(test_xyz_file)
        assert isinstance(smiles, str)
        assert len(smiles) > 0
    
    def test_xyz_to_sdf(self, test_xyz_file):
        """Test XYZ to SDF conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False) as f:
            sdf_file = f.name
        
        try:
            FileConverter.xyz_to_sdf(test_xyz_file, sdf_file)
            assert os.path.exists(sdf_file)
            assert os.path.getsize(sdf_file) > 0
            # Basic SDF validation
            with open(sdf_file, 'r') as f:
                content = f.read()
                assert "V2000" in content or "END" in content
        finally:
            if os.path.exists(sdf_file):
                os.unlink(sdf_file)
    
    def test_xyz_to_pdb(self, test_xyz_file):
        """Test XYZ to PDB conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as f:
            pdb_file = f.name
        
        try:
            FileConverter.xyz_to_pdb(test_xyz_file, pdb_file)
            assert os.path.exists(pdb_file)
            assert os.path.getsize(pdb_file) > 0
        finally:
            if os.path.exists(pdb_file):
                os.unlink(pdb_file)
    
    def test_xyz_to_pdb_string(self, test_xyz_file):
        """Test XYZ to PDB string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        pdb_string = FileConverter.xyz_to_pdb_string(test_xyz_file)
        assert isinstance(pdb_string, str)
        assert len(pdb_string) > 0
        # Basic PDB validation
        assert "ATOM" in pdb_string or "HETATM" in pdb_string
    
    def test_xyz_to_mol2_string(self, test_xyz_file):
        """Test XYZ to MOL2 string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        mol2_string = FileConverter.xyz_to_mol2_string(test_xyz_file)
        assert isinstance(mol2_string, str)
        assert len(mol2_string) > 0
        # Basic MOL2 validation
        assert "@<TRIPOS>" in mol2_string
    
    def test_xyz_to_cif_string(self, test_xyz_file):
        """Test XYZ to CIF string conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        cif_string = FileConverter.xyz_to_cif_string(test_xyz_file)
        assert isinstance(cif_string, str)
        assert len(cif_string) > 0
        # Basic CIF validation
        assert "data_" in cif_string
    
    def test_xyz_to_inchi_file_not_found(self):
        """Test XYZ to InChI conversion with non-existent file"""
        with pytest.raises(FileNotFoundError):
            FileConverter.xyz_to_inchi("nonexistent.xyz")
    
    def _is_openbabel_available(self):
        """Check if Open Babel is available on the system"""
        try:
            result = subprocess.run(['obabel', '-V'], 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=5)
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False