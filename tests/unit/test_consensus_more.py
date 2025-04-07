"""
Additional unit tests for the consensus module.
"""

from gfftk.consensus import getAED, map_coords, reasonable_model, score_evidence


class TestConsensusMore:
    """Additional tests for consensus functions."""

    def test_getAED(self):
        """Test the getAED function."""
        # Test with identical transcripts
        query = [[1, 100], [200, 300]]
        reference = [[1, 100], [200, 300]]

        aed = getAED(query, reference)
        assert aed == 0.0  # Perfect match should have AED of 0

        # Test with different transcripts
        query = [[1, 100], [200, 300]]
        reference = [[1, 100], [250, 350]]  # Different second exon

        aed = getAED(query, reference)
        assert aed > 0.0  # Different transcripts should have AED > 0
        assert aed <= 1.0  # AED should be <= 1

        # Test with completely different transcripts
        query = [[1, 100], [200, 300]]
        reference = [[500, 600], [700, 800]]  # Completely different exons

        aed = getAED(query, reference)
        assert aed == 1.0  # Completely different transcripts should have AED of 1

    def test_score_evidence(self):
        """Test the score_evidence function."""
        # Test with perfect overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[1, 100], [200, 300]]

        score = score_evidence(g_coords, e_coords)
        assert score == 20  # Perfect overlap with default weight of 2 (10 * 2)

        # Test with partial overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[50, 150], [250, 350]]

        score = score_evidence(g_coords, e_coords)
        assert score > 0  # Partial overlap should have score > 0
        assert score < 20  # Partial overlap should have score < 20

        # Test with no overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[400, 500], [600, 700]]

        score = score_evidence(g_coords, e_coords)
        assert score == 0  # No overlap should have score of 0

        # Test with custom weight
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[1, 100], [200, 300]]

        score = score_evidence(g_coords, e_coords, weight=5)
        assert score == 50  # Perfect overlap with weight of 5 (10 * 5)

    def test_map_coords(self):
        """Test the map_coords function."""
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

        # Test with partial overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[10, 90], [210, 290]]

        result = map_coords(g_coords, e_coords)

        # Check the result
        assert isinstance(result, list)
        assert len(result) == 2

        # The result should be a list of lists with differences
        # For partial overlap, the differences should be [9, -10] and [10, -10]
        assert result[0][0] > 0  # Start difference is positive
        assert result[0][1] < 0  # End difference is negative
        assert result[1][0] > 0  # Start difference is positive
        assert result[1][1] < 0  # End difference is negative

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
