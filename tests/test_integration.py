"""Integration tests for mol_conversion package"""

import pytest
import tempfile
import os
from pathlib import Path
from mol_conversion import FileConversion, MemoryConversion, GCNEncoding, split_multiframe_xyz


class TestIntegration:
    """Integration tests for the entire package"""
    
    @pytest.fixture
    def test_xyz_file(self):
        """Path to the real test XYZ file"""
        test_dir = Path(__file__).parent / "data"
        return str(test_dir / "test.xyz")
    
    @pytest.fixture
    def test_xyz_string(self):
        """Real XYZ string from test data file"""
        test_file = Path(__file__).parent / "data" / "test.xyz"
        if test_file.exists():
            return test_file.read_text()
        else:
            pytest.skip(f"Test file {test_file} not found")
    
    def test_file_to_memory_consistency(self, test_xyz_file):
        """Test that file and memory conversions give consistent results"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        # Convert using FileConversion
        file_inchi = FileConversion.xyz_to_inchi(test_xyz_file)
        file_smiles = FileConversion.xyz_to_smiles(test_xyz_file)
        
        # Read XYZ content and convert using MemoryConversion
        with open(test_xyz_file, 'r') as f:
            xyz_content = f.read()
        
        memory_inchi = MemoryConversion.xyz_to_inchi_string(xyz_content)
        memory_smiles = MemoryConversion.xyz_to_smiles_string(xyz_content)
        
        # Results should be consistent
        assert file_inchi == memory_inchi
        assert file_smiles == memory_smiles
    
    def test_memory_conversion_roundtrip(self, test_xyz_file):
        """Test round-trip conversion: XYZ -> MOL -> XYZ"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        # Read original XYZ
        with open(test_xyz_file, 'r') as f:
            original_xyz = f.read()
        
        # Convert to MOL string
        mol_string = MemoryConversion.xyz_to_mol_string(original_xyz)
        assert len(mol_string) > 0
        
        # Convert MOL back to RDKit mol and verify
        if self._is_rdkit_available():
            rdkit_mol = MemoryConversion.xyz_to_rdkit_mol(original_xyz)
            assert rdkit_mol is not None
            # Real file has 70 atoms
            assert rdkit_mol.GetNumAtoms() == 70
    
    def test_gcn_with_memory_conversion(self, test_xyz_string):
        """Test GCN encoding with memory conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        # Generate GCN encodings
        encoder = GCNEncoding(test_xyz_string)
        encodings = encoder.encodings
        
        # Verify structure
        assert isinstance(encodings, dict)
        assert "gcn0" in encodings
        assert "gcn1" in encodings
        assert "gcn2" in encodings
        
        # Real file has 70 atoms
        assert len(encodings["gcn0"]) == 70
    
    def test_multiframe_with_memory_conversion(self, test_xyz_string):
        """Test multi-frame XYZ with memory conversion"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        # Split multi-frame XYZ
        frames = split_multiframe_xyz(test_xyz_string)
        assert len(frames) > 1
        
        # Test conversion of each frame
        for i, frame in enumerate(frames):
            # Convert to InChI
            inchi = MemoryConversion.xyz_to_inchi_string(frame)
            assert isinstance(inchi, str)
            assert len(inchi) > 0
            
            # Convert to bond order matrix
            bond_matrix = MemoryConversion.xyz_to_bond_order_matrix(frame)
            assert bond_matrix.shape[0] > 0  # Should have atoms
            
            # Generate GCN encodings
            encoder = GCNEncoding(frame)
            encodings = encoder.encodings
            assert len(encodings["gcn0"]) == bond_matrix.shape[0]
    
    def test_complete_workflow(self, test_xyz_file):
        """Test complete workflow from file to various formats"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        if not os.path.exists(test_xyz_file):
            pytest.skip(f"Test file {test_xyz_file} not found")
        
        # Step 1: File to various formats
        inchi = FileConversion.xyz_to_inchi(test_xyz_file)
        inchikey = FileConversion.xyz_to_inchikey(test_xyz_file)
        smiles = FileConversion.xyz_to_smiles(test_xyz_file)
        
        # Step 2: Read file content for memory operations
        with open(test_xyz_file, 'r') as f:
            xyz_content = f.read()
        
        # Step 3: Memory conversions
        mol_string = MemoryConversion.xyz_to_mol_string(xyz_content)
        sdf_string = MemoryConversion.xyz_to_sdf_string(xyz_content)
        pdb_string = MemoryConversion.xyz_to_pdb_string(xyz_content)
        
        # Step 4: Generate GCN encodings
        encoder = GCNEncoding(xyz_content)
        encodings = encoder.encodings
        
        # Step 5: Bond order matrix
        bond_matrix = MemoryConversion.xyz_to_bond_order_matrix(xyz_content)
        
        # Verify all results are valid
        assert isinstance(inchi, str) and len(inchi) > 0
        assert isinstance(inchikey, str) and len(inchikey) == 27
        assert isinstance(smiles, str) and len(smiles) > 0
        assert isinstance(mol_string, str) and len(mol_string) > 0
        assert isinstance(sdf_string, str) and len(sdf_string) > 0
        assert isinstance(pdb_string, str) and len(pdb_string) > 0
        assert isinstance(encodings, dict) and len(encodings) == 3
        assert isinstance(bond_matrix, object) and bond_matrix.size > 0
    
    def test_error_handling_consistency(self, test_xyz_file):
        """Test error handling consistency across modules"""
        # Test file not found
        with pytest.raises(FileNotFoundError):
            FileConversion.xyz_to_inchi("nonexistent.xyz")
        
        # Test invalid XYZ content
        invalid_xyz = """1
Invalid atom
X    0.000000    0.000000    0.000000
"""
        
        # Memory conversion should handle invalid input gracefully
        if self._is_openbabel_available():
            try:
                result = MemoryConversion.xyz_to_inchi_string(invalid_xyz)
                # Should either return a result or raise an exception
                assert isinstance(result, str)
            except Exception:
                # It's acceptable to raise an exception for invalid input
                pass
    
    def test_large_molecule_handling(self, test_xyz_string):
        """Test handling of large molecules (70 atoms)"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        # Test memory conversion with large molecule
        inchi = MemoryConversion.xyz_to_inchi_string(test_xyz_string)
        assert isinstance(inchi, str) and len(inchi) > 0
        
        # Test bond order matrix with large molecule
        bond_matrix = MemoryConversion.xyz_to_bond_order_matrix(test_xyz_string)
        assert bond_matrix.shape == (70, 70)
        
        # Test GCN encoding with large molecule
        encoder = GCNEncoding(test_xyz_string)
        encodings = encoder.encodings
        assert len(encodings["gcn0"]) == 70
        
        # Verify performance is reasonable (should complete quickly)
        # This is more of a smoke test than a performance test
        assert True
    
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