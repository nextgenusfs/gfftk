"""
Fixed unit tests for the consensus module.
"""

import os
import tempfile

from gfftk.consensus import (
    contained,
    fasta_length,
    get_overlap,
    getAED,
    map_coords,
    reasonable_model,
    score_evidence,
)


class TestConsensusFixed:
    """Fixed tests for consensus functions."""

    def test_getAED_identical(self):
        """Test the getAED function with identical transcripts."""
        # Test with identical transcripts
        query = [
            [1, 100],
            [200, 300],
        ]

        reference = [
            [1, 100],
            [200, 300],
        ]

        aed = getAED(query, reference)
        assert aed == 0.0  # Perfect match should have AED of 0

    def test_getAED_different(self):
        """Test the getAED function with different transcripts."""
        # Test with different transcripts
        query = [
            [1, 100],
            [200, 300],
        ]

        reference = [
            [1, 100],
            [250, 350],  # Different second exon
        ]

        aed = getAED(query, reference)
        assert aed > 0.0  # Different transcripts should have AED > 0
        assert aed <= 1.0  # AED should be <= 1

    def test_getAED_no_overlap(self):
        """Test the getAED function with no overlap."""
        # Test with completely different transcripts
        query = [
            [1, 100],
            [200, 300],
        ]

        reference = [
            [500, 600],
            [700, 800],  # Completely different exons
        ]

        aed = getAED(query, reference)
        assert aed == 1.0  # Completely different transcripts should have AED of 1

    def test_score_evidence_identical(self):
        """Test the score_evidence function with identical coordinates."""
        # Test with perfect overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[1, 100], [200, 300]]

        score = score_evidence(g_coords, e_coords)
        assert score == 20  # Perfect overlap with default weight of 2 (10 * 2)

    def test_score_evidence_partial(self):
        """Test the score_evidence function with partial overlap."""
        # Test with partial overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[50, 150], [250, 350]]

        score = score_evidence(g_coords, e_coords)
        assert score == 6  # Partial overlap with coverage adjustment

    def test_score_evidence_no_overlap(self):
        """Test the score_evidence function with no overlap."""
        # Test with no overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[400, 500], [600, 700]]

        score = score_evidence(g_coords, e_coords)
        assert score == 0  # No overlap should have score of 0

    def test_score_evidence_custom_weight(self):
        """Test the score_evidence function with custom weight."""
        # Test with custom weight
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[1, 100], [200, 300]]

        score = score_evidence(g_coords, e_coords, weight=5)
        assert score == 50  # Perfect overlap with weight of 5 (10 * 5)

    def test_map_coords_identical(self):
        """Test the map_coords function with identical coordinates."""
        # Test with simple coordinates
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[1, 100], [200, 300]]

        result = map_coords(g_coords, e_coords)

        # Check the result
        assert isinstance(result, list)
        assert len(result) == 2

        # The result should be a list of lists with differences
        # For identical coordinates, the differences should be [0, 0]
        assert result[0] == [0, 0]
        assert result[1] == [0, 0]

    def test_map_coords_partial(self):
        """Test the map_coords function with partial overlap."""
        # Test with partial overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[10, 90], [210, 290]]

        result = map_coords(g_coords, e_coords)

        # Check the result
        assert isinstance(result, list)
        assert len(result) == 2

        # The result should be a list of lists with differences
        # For partial overlap, the differences should be [9, -10] and [10, -10]
        assert result[0] == [9, -10]
        assert result[1] == [10, -10]

    def test_map_coords_no_overlap(self):
        """Test the map_coords function with no overlap."""
        # Test with no overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[400, 500], [600, 700]]

        result = map_coords(g_coords, e_coords)

        # Check the result
        assert isinstance(result, list)
        assert len(result) == 2

        # The result should be a list of empty lists for no overlap
        assert result[0] == []
        assert result[1] == []

    def test_reasonable_model(self):
        """Test the reasonable_model function."""
        # Test with a reasonable model
        coords = [[1, 100], [200, 300]]

        result = reasonable_model(coords)

        # Check the result
        assert result is True  # The model should be reasonable

        # Test with an unreasonable model (exon too short)
        coords = [[1, 2], [200, 300]]

        result = reasonable_model(coords)

        # Check the result
        # The function returns a string with the reason for failure
        assert isinstance(result, str)
        assert "exon < 3" in result  # The model should fail because exon is too short

    def test_fasta_length(self):
        """Test the fasta_length function."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp:
            temp.write(">contig1\n")
            temp.write("ATGCATGCATGCATGCATGC\n")  # 20 bases
            temp.write(">contig2\n")
            temp.write("ATGCATGCATGC\n")  # 12 bases
            temp_name = temp.name

        try:
            # Get the length of the FASTA file
            length = fasta_length(temp_name)

            # Check that the length is correct
            assert length == 32  # 20 + 12 = 32
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_get_overlap(self):
        """Test the get_overlap function."""
        # Test with overlapping ranges
        a = [10, 20]
        b = [15, 25]
        overlap = get_overlap(a, b)
        assert overlap == 5  # Overlap is 15-20 = 5

        # Test with non-overlapping ranges
        a = [10, 20]
        b = [30, 40]
        overlap = get_overlap(a, b)
        assert overlap == 0

        # Test with one range contained in the other
        a = [10, 30]
        b = [15, 25]
        overlap = get_overlap(a, b)
        assert overlap == 10  # Overlap is 15-25 = 10

        # Test with ranges touching but not overlapping
        a = [10, 20]
        b = [20, 30]
        overlap = get_overlap(a, b)
        assert overlap == 0

    def test_contained(self):
        """Test the contained function."""
        # Test with one range contained in the other
        a = [15, 25]
        b = [10, 30]
        result = contained(a, b)
        assert result is True

        # Test with ranges that overlap but neither contains the other
        a = [10, 20]
        b = [15, 25]
        result = contained(a, b)
        assert result is False

        # Test with non-overlapping ranges
        a = [10, 20]
        b = [30, 40]
        result = contained(a, b)
        assert result is False

        # Test with identical ranges
        a = [10, 20]
        b = [10, 20]
        result = contained(a, b)
        assert result is False  # The function returns False for identical ranges
