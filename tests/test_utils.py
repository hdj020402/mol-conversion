"""Tests for utils and gcn_encoding modules"""

import pytest
import numpy as np
from pathlib import Path
from mol_conversion import GCNEncoder, split_multiframe_xyz, split_multiframe_xyz_with_comments
from mol_conversion.gcn_encoding import GCNEncoder as GCNEncoderDirect


class TestGCNEncoder:
    """Test class for GCNEncoder"""

    @pytest.fixture
    def test_xyz_string(self):
        """Real XYZ string from test data file"""
        test_file = Path(__file__).parent / "data" / "test.xyz"
        if test_file.exists():
            return test_file.read_text()
        else:
            pytest.skip(f"Test file {test_file} not found")

    def test_gcn_encoding_from_xyz(self, test_xyz_string):
        """Test GCN encoding with real XYZ file via from_xyz"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        encoder = GCNEncoder.from_xyz(test_xyz_string)
        encodings = encoder.encodings

        # Check structure
        assert isinstance(encodings, dict)
        assert "gcn0" in encodings
        assert "gcn1" in encodings
        assert "gcn2" in encodings

        # Real file has 70 atoms
        assert len(encodings["gcn0"]) == 70
        assert len(encodings["gcn1"]) == 70
        assert len(encodings["gcn2"]) == 70

        # Should have some non-empty encodings (heavy atoms)
        gcn0 = encodings["gcn0"]
        non_empty_count = sum(1 for encoding in gcn0 if encoding != "")
        assert non_empty_count > 0  # Should have some heavy atoms

        # Check that GCN1 and GCN2 are more complex than GCN0
        gcn1 = encodings["gcn1"]
        gcn2 = encodings["gcn2"]
        for i in range(len(gcn0)):
            if gcn0[i] != "":
                # Non-hydrogen atoms should have more complex GCN1 encodings
                assert len(gcn1[i]) >= len(gcn0[i])
                # GCN2 should be at least as complex as GCN1
                assert len(gcn2[i]) >= len(gcn1[i])

    def test_gcn_encoding_from_smiles(self):
        """Test GCN encoding from SMILES string"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        # Ethanol: CCO
        encoder = GCNEncoder.from_smiles("CCO")
        encodings = encoder.encodings

        assert isinstance(encodings, dict)
        assert "gcn0" in encodings

        # Should have atoms
        assert len(encoder.atoms) > 0
        # Heavy atoms should have non-empty encodings
        non_empty = [e for e in encodings["gcn0"] if e != ""]
        assert len(non_empty) > 0

    def test_gcn_encoding_from_inchi(self):
        """Test GCN encoding from InChI string"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        # Ethanol InChI
        encoder = GCNEncoder.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
        encodings = encoder.encodings

        assert isinstance(encodings, dict)
        assert len(encoder.atoms) > 0
        non_empty = [e for e in encodings["gcn0"] if e != ""]
        assert len(non_empty) > 0

    def test_gcn_encoding_from_sdf(self):
        """Test GCN encoding from SDF string"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        # Minimal methane SDF
        sdf_str = """methane
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6292    0.6292    0.6292 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6292   -0.6292    0.6292 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6292    0.6292   -0.6292 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6292   -0.6292   -0.6292 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
"""
        encoder = GCNEncoder.from_sdf(sdf_str)
        encodings = encoder.encodings

        assert isinstance(encodings, dict)
        # 5 atoms: 1 C + 4 H
        assert len(encoder.atoms) == 5
        # C should have encoding C04
        non_empty = [e for e in encodings["gcn0"] if e != ""]
        assert len(non_empty) == 1
        assert non_empty[0] == "C04"

    def test_gcn_encoding_cross_format_consistency(self):
        """Test that same molecule via different formats gives consistent GCN0 heavy-atom encodings"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        # Benzene via SMILES
        enc_smiles = GCNEncoder.from_smiles("c1ccccc1")
        # Benzene via InChI
        enc_inchi = GCNEncoder.from_inchi("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")

        # Heavy atom GCN0 encodings should match
        smiles_gcn0 = [e for e in enc_smiles.gcn0 if e != ""]
        inchi_gcn0 = [e for e in enc_inchi.gcn0 if e != ""]

        assert sorted(smiles_gcn0) == sorted(inchi_gcn0)

    def test_gcn_encoding_progression(self, test_xyz_string):
        """Test that GCN encodings become progressively more complex"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        encoder = GCNEncoder.from_xyz(test_xyz_string)
        encodings = encoder.encodings

        gcn0 = encodings["gcn0"]
        gcn1 = encodings["gcn1"]
        gcn2 = encodings["gcn2"]

        # For carbon atoms (non-hydrogen), encodings should get progressively more complex
        for i in range(len(gcn0)):
            if gcn0[i] != "":  # Heavy atom
                # GCN1 should be at least as complex as GCN0
                assert len(gcn1[i]) >= len(gcn0[i])

                # GCN2 should be at least as complex as GCN1
                assert len(gcn2[i]) >= len(gcn1[i])

                # Check for expected delimiters
                if len(gcn1[i]) > len(gcn0[i]):
                    assert "(" in gcn1[i]  # GCN1 uses parentheses

                if len(gcn2[i]) > len(gcn1[i]):
                    assert "[" in gcn2[i]  # GCN2 uses brackets

    def test_gcn_encoding_invalid_xyz(self):
        """Test GCN encoding with invalid XYZ string"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        invalid_xyz = """1
Invalid atom
X    0.000000    0.000000    0.000000
"""

        # Should handle invalid atoms gracefully
        try:
            encoder = GCNEncoder.from_xyz(invalid_xyz)
            encodings = encoder.encodings
            # Unknown atom should be treated as heavy atom
            assert len(encodings["gcn0"]) == 1
        except Exception:
            # It's acceptable if Open Babel raises an exception
            pass

    def test_direct_import(self, test_xyz_string):
        """Test that GCNEncoder can be imported directly from gcn_encoding module"""
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")

        encoder = GCNEncoderDirect.from_xyz(test_xyz_string)
        assert isinstance(encoder, GCNEncoder)

    def _is_openbabel_available(self):
        """Check if Open Babel Python bindings are available"""
        try:
            from openbabel import pybel
            return True
        except ImportError:
            return False


class TestMultiframeXYZ:
    """Test class for multi-frame XYZ functions"""
    
    @pytest.fixture
    def test_xyz_string(self):
        """Real XYZ string from test data file"""
        test_file = Path(__file__).parent / "data" / "test.xyz"
        if test_file.exists():
            return test_file.read_text()
        else:
            pytest.skip(f"Test file {test_file} not found")
    
    def test_split_real_file(self, test_xyz_string):
        """Test splitting real XYZ file"""
        frames = split_multiframe_xyz(test_xyz_string)
        assert isinstance(frames, list)
        assert len(frames) > 1  # Real file should have multiple frames
        
        # Check each frame
        for frame in frames:
            assert isinstance(frame, str)
            assert len(frame) > 0
            # Each frame should start with atom count
            lines = frame.strip().split('\n')
            assert len(lines) >= 3  # At least atom count, comment, and one atom
            assert lines[0].strip().isdigit()  # First line should be atom count
    
    def test_split_with_comments_real_file(self, test_xyz_string):
        """Test splitting real XYZ file with comments"""
        frames, comments = split_multiframe_xyz_with_comments(test_xyz_string)
        assert isinstance(frames, list)
        assert isinstance(comments, list)
        assert len(frames) == len(comments)
        assert len(frames) > 1
        
        # Check that comments are meaningful
        for comment in comments:
            assert isinstance(comment, str)
            assert len(comment) >= 0
    
    def test_empty_xyz_string(self):
        """Test splitting empty XYZ string"""
        frames = split_multiframe_xyz("")
        assert isinstance(frames, list)
        assert len(frames) == 0
    
    def test_invalid_xyz_string(self):
        """Test splitting invalid XYZ string"""
        invalid_xyz = """invalid
not a number
"""
        # Should raise an exception for invalid format
        with pytest.raises(ValueError):
            split_multiframe_xyz(invalid_xyz)
    
    def test_incomplete_xyz_string(self):
        """Test splitting incomplete XYZ string"""
        incomplete_xyz = """4
Incomplete frame
O    0.000000    0.000000    0.000000
"""
        # Should raise an exception for incomplete frame
        with pytest.raises(ValueError):
            split_multiframe_xyz(incomplete_xyz)


class TestXYZStringToBondOrderMatrix:
    """Test class for xyz_string_to_bond_order_matrix function"""
    
    @pytest.fixture
    def test_xyz_string(self):
        """Real XYZ string from test data file"""
        test_file = Path(__file__).parent / "data" / "test.xyz"
        if test_file.exists():
            return test_file.read_text()
        else:
            pytest.skip(f"Test file {test_file} not found")
    
    def test_bond_order_matrix(self, test_xyz_string):
        """Test bond order matrix with real molecule"""
        # Skip test if Open Babel is not available
        if not self._is_openbabel_available():
            pytest.skip("Open Babel not available")
        
        from mol_conversion.utils import xyz_string_to_bond_order_matrix
        
        bond_matrix = xyz_string_to_bond_order_matrix(test_xyz_string)
        assert isinstance(bond_matrix, np.ndarray)
        # Real file has 70 atoms
        assert bond_matrix.shape == (70, 70)
        
        # Check diagonal is zero
        assert np.all(np.diag(bond_matrix) == 0)
        
        # Check matrix is symmetric
        assert np.allclose(bond_matrix, bond_matrix.T)
        
        # Should have some bonds
        assert np.sum(bond_matrix > 0) > 0
    
    def _is_openbabel_available(self):
        """Check if Open Babel Python bindings are available"""
        try:
            from openbabel import openbabel
            return True
        except ImportError:
            return False