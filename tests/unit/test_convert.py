"""
Unit tests for the convert module.
"""

from gfftk.convert import gff2gtf, gff2proteins, gff2transcripts


class TestConvert:
    """Tests for conversion functions."""

    def test_convert_functions_exist(self):
        """Test that the convert functions exist and can be imported."""
        # Just check that the functions exist
        assert callable(gff2gtf)
        assert callable(gff2proteins)
        assert callable(gff2transcripts)
