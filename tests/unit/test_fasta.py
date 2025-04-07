"""
Unit tests for the fasta module.
"""

import os
import tempfile

from gfftk.fasta import RevComp, fasta2dict, fasta2headers, getSeqRegions, softwrap, translate


class TestFasta:
    """Tests for FASTA functions."""

    def test_fasta2dict(self):
        """Test converting FASTA to dictionary."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp:
            temp.write(">seq1\n")
            temp.write("ATGCATGCATGC\n")
            temp.write(">seq2\n")
            temp.write("GTATCGATCGAT\n")
            temp_name = temp.name

        try:
            # Convert FASTA to dictionary
            result = fasta2dict(temp_name)

            # Check the result
            assert len(result) == 2
            assert "seq1" in result
            assert "seq2" in result
            assert result["seq1"] == "ATGCATGCATGC"
            assert result["seq2"] == "GTATCGATCGAT"
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_softwrap(self):
        """Test the softwrap function."""
        # Test with a short string
        string = "ATGCATGCATGC"
        result = softwrap(string, every=4)
        assert result == "ATGC\nATGC\nATGC"

        # Test with a string that's not a multiple of the wrap length
        string = "ATGCATGCATG"
        result = softwrap(string, every=4)
        assert result == "ATGC\nATGC\nATG"

        # Test with an empty string
        string = ""
        result = softwrap(string)
        assert result == ""

    def test_fasta2headers(self):
        """Test extracting headers from FASTA."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp:
            temp.write(">seq1 description 1\n")
            temp.write("ATGCATGCATGC\n")
            temp.write(">seq2 description 2\n")
            temp.write("GTATCGATCGAT\n")
            temp_name = temp.name

        try:
            # Extract headers from FASTA
            result = fasta2headers(temp_name)

            # Check the result
            # fasta2headers returns a set of sequence IDs
            assert len(result) == 2
            assert "seq1" in result
            assert "seq2" in result
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_RevComp(self):
        """Test reverse complementing DNA sequences."""
        # Test with a simple sequence
        seq = "ATGCATGCATGC"
        result = RevComp(seq)
        assert result == "GCATGCATGCAT"

        # Test with a sequence containing N
        seq = "ATGCNNGCATGC"
        result = RevComp(seq)
        assert result == "GCATGCNNGCAT"

        # Test with an empty sequence
        seq = ""
        result = RevComp(seq)
        assert result == ""

    def test_translate(self):
        """Test translating DNA to protein."""
        # Test with a simple sequence
        seq = "ATGCATGCATGC"
        result = translate(seq, "+", 0)
        assert result == "MHAC"

        # Test with a sequence containing N
        seq = "ATGCNNGCATGC"
        result = translate(seq, "+", 0)
        assert "X" in result  # N translates to X

        # Test with a sequence that's not a multiple of 3
        seq = "ATGCATGCATG"  # 11 bp
        result = translate(seq, "+", 0)
        assert len(result) == 3  # 11 // 3 = 3 complete codons

        # Test with an empty sequence
        seq = ""
        result = translate(seq, "+", 0)
        assert result == ""

    def test_getSeqRegions(self):
        """Test extracting regions from sequences."""
        # Create a dictionary of sequences
        fasta_dict = {
            "seq1": "ATGCATGCATGC",  # 12 bp
            "seq2": "GTATCGATCGAT",  # 12 bp
        }

        # Define regions to extract
        regions = [
            (1, 6),  # ATGCAT
            (7, 12),  # GCATGC
        ]

        # Extract regions
        result = getSeqRegions(fasta_dict, "seq1", regions)

        # Check the result
        assert isinstance(result, str)
        assert result == "ATGCATGCATGC"  # ATGCAT + GCATGC

        # Test with regions that are out of bounds
        regions = [
            (1, 20),  # Out of bounds
        ]

        # Extract regions
        result = getSeqRegions(fasta_dict, "seq1", regions)

        # Check the result
        assert isinstance(result, str)
        # Out of bounds returns the sequence up to the end
        assert result == "ATGCATGCATGC"

        # Test with non-existent sequence
        try:
            result = getSeqRegions(fasta_dict, "seq3", [(1, 6)])
            assert False, "Should have raised KeyError"
        except KeyError:
            pass
