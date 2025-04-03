"""
Unit tests for the gff module functions.
"""

import os
import tempfile
from gfftk.gff import (
    _detect_format,
    simplifyGO,
)


class TestGFFFunctions:
    """Tests for GFF functions."""

    def test_detect_format(self):
        """Test the _detect_format function."""
        # Create a GFF3-like file for testing
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("##gff-version 3\n")
            temp.write(
                "contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene\n"
            )
            temp_name = temp.name

        try:
            # Test with a GFF3 file
            format_type = _detect_format(temp_name)
            # The function returns a tuple with a parser function and a format name
            assert isinstance(format_type, tuple)
            assert format_type[1] in ["gff3", "default"]
        finally:
            # Clean up
            os.unlink(temp_name)

        # Create a GTF-like file for testing
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write(
                'contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tgene_id "gene1"; transcript_id "mRNA1";\n'
            )
            temp_name = temp.name

        try:
            # Test with a GTF file
            format_type = _detect_format(temp_name)
            # The function returns a tuple with a parser function and a format name
            assert isinstance(format_type, tuple)
            assert format_type[1] in ["gtf", "default"]
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_simplify_go(self):
        """Test the simplifyGO function."""
        # Test with a simple GO term
        go_terms = ["GO:0005634"]
        result = simplifyGO(go_terms)
        assert result == ["GO:0005634"]

        # Test with multiple GO terms
        go_terms = ["GO:0005634", "GO:0003677", "GO:0006355"]
        result = simplifyGO(go_terms)
        assert set(result) == set(["GO:0005634", "GO:0003677", "GO:0006355"])

        # Test with an empty list
        go_terms = []
        result = simplifyGO(go_terms)
        assert result == []
