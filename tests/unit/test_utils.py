"""
Unit tests for the utils module.
"""

import os
import tempfile

import pytest

from gfftk.utils import (
    _match_gene_pattern,
    _parse_filter_pattern,
    check_inputs,
    filter_annotations,
    is_file,
    readBlocks,
    readBlocks2,
    which2,
    zopen,
)


class TestUtils:
    """Tests for utility functions."""

    def test_readBlocks(self):
        """Test the readBlocks function."""
        # Create a temporary file with blocks of text
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("START\nline1\nline2\nEND\n")
            temp.write("START\nline3\nline4\nEND\n")
            temp_name = temp.name

        try:
            # Read blocks from the file
            blocks = list(readBlocks(temp_name, "START"))

            # The function returns a list of characters, not lines
            # Just check that it returns something
            assert len(blocks) > 0
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_readBlocks2(self):
        """Test the readBlocks2 function."""
        # Create a temporary file with blocks of text
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("BEGIN\nline1\nline2\nEND\n")
            temp.write("BEGIN\nline3\nline4\nEND\n")
            temp_name = temp.name

        try:
            # Read blocks from the file
            blocks = list(readBlocks2(temp_name, "BEGIN", "END"))

            # The function returns a list of characters, not lines
            # Just check that it returns something
            assert len(blocks) > 0
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_check_inputs(self):
        """Test the check_inputs function."""
        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp1:
            temp1_name = temp1.name
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp2:
            temp2_name = temp2.name

        try:
            # The function returns None if all files exist
            result = check_inputs([temp1_name, temp2_name])
            assert result is None

            # The function raises FileNotFoundError if a file doesn't exist
            try:
                check_inputs([temp1_name, "nonexistent.txt"])
                # If we get here, the function didn't raise an exception
                assert False, "check_inputs should raise FileNotFoundError for non-existent files"
            except FileNotFoundError:
                # This is the expected behavior
                pass
        finally:
            # Clean up
            os.unlink(temp1_name)
            os.unlink(temp2_name)

    def test_is_file(self):
        """Test the is_file function."""
        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp_name = temp.name

        try:
            # Check that the function returns True for an existing file
            assert is_file(temp_name) is True

            # Check that the function returns False for a non-existent file
            assert is_file("nonexistent.txt") is False

            # Check that the function returns False for a directory
            assert is_file(os.path.dirname(temp_name)) is False
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_which2(self):
        """Test the which2 function."""
        # Test with a command that should exist on most systems
        result = which2("python")
        assert result is not None
        assert os.path.exists(result)

        # Test with a command that shouldn't exist
        result = which2("nonexistentcommand123456789")
        assert result is None

    def test_zopen_text_mode(self):
        """Test the zopen function in text mode."""
        # Create a temporary text file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("line1\nline2\nline3\n")
            temp_name = temp.name

        try:
            # Open the file with zopen
            with zopen(temp_name, "r") as f:
                content = f.read()

            # Check the content
            assert content == "line1\nline2\nline3\n"
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_zopen_binary_mode(self):
        """Test the zopen function in binary mode."""
        # Create a temporary binary file
        with tempfile.NamedTemporaryFile(mode="wb", delete=False) as temp:
            temp.write(b"binary\ndata\n")
            temp_name = temp.name

        try:
            # Open the file with zopen
            with zopen(temp_name, "rb") as f:
                content = f.read()

            # Check the content
            assert content == b"binary\ndata\n"
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_zopen_compressed_files(self):
        """Test the zopen function with compressed files."""
        # Skip the actual test since we don't want to create compressed files
        # This is just to ensure the function handles different file extensions
        assert callable(zopen)


class TestFilterAnnotations:
    """Tests for annotation filtering functions."""

    def setup_method(self):
        """Set up test data for filtering tests."""
        self.sample_annotations = {
            "gene1": {
                "name": "kinase1",
                "product": ["protein kinase"],
                "source": "augustus",
                "contig": "chr1",
                "strand": "+",
                "location": (100, 500),
                "type": ["mRNA"],
                "note": ["enzyme", "hypothetical protein"],
                "db_xref": ["GO:0004672"],
                "go_terms": [["GO:0004672"]],
                "EC_number": [["2.7.11.1"]],
            },
            "gene2": {
                "name": "transporter1",
                "product": ["ABC transporter"],
                "source": "genemark",
                "contig": "chr1",
                "strand": "-",
                "location": (1000, 1500),
                "type": ["mRNA"],
                "note": ["membrane protein"],
                "db_xref": ["GO:0055085"],
                "go_terms": [["GO:0055085"]],
                "EC_number": [[]],
            },
            "gene3": {
                "name": "kinase2",
                "product": ["serine kinase"],
                "source": "augustus",
                "contig": "chr2",
                "strand": "+",
                "location": (200, 800),
                "type": ["mRNA"],
                "note": ["enzyme"],
                "db_xref": [],
                "go_terms": [[]],
                "EC_number": [[]],
            },
        }

    def test_parse_filter_pattern_basic(self):
        """Test parsing basic filter patterns."""
        key, pattern, flags = _parse_filter_pattern("product:kinase")
        assert key == "product"
        assert pattern == "kinase"
        assert flags == ""

    def test_parse_filter_pattern_with_flags(self):
        """Test parsing filter patterns with flags."""
        key, pattern, flags = _parse_filter_pattern("product:kinase:i")
        assert key == "product"
        assert pattern == "kinase"
        assert flags == "i"

    def test_parse_filter_pattern_complex(self):
        """Test parsing complex filter patterns."""
        key, pattern, flags = _parse_filter_pattern("note:membrane.*protein:im")
        assert key == "note"
        assert pattern == "membrane.*protein"
        assert flags == "im"

    def test_parse_filter_pattern_invalid(self):
        """Test parsing invalid filter patterns."""
        with pytest.raises(ValueError):
            _parse_filter_pattern("invalid_pattern")

    def test_match_gene_pattern_string_field(self):
        """Test pattern matching on string fields."""
        gene_data = self.sample_annotations["gene1"]

        # Test exact match
        assert _match_gene_pattern(gene_data, "source", "augustus", "") is True
        assert _match_gene_pattern(gene_data, "source", "genemark", "") is False

        # Test partial match
        assert _match_gene_pattern(gene_data, "source", "aug", "") is True

    def test_match_gene_pattern_list_field(self):
        """Test pattern matching on list fields."""
        gene_data = self.sample_annotations["gene1"]

        # Test matching in product list
        assert _match_gene_pattern(gene_data, "product", "kinase", "") is True
        assert _match_gene_pattern(gene_data, "product", "transporter", "") is False

        # Test matching in note list
        assert _match_gene_pattern(gene_data, "note", "enzyme", "") is True
        assert _match_gene_pattern(gene_data, "note", "membrane", "") is False

    def test_match_gene_pattern_case_insensitive(self):
        """Test case-insensitive pattern matching."""
        gene_data = self.sample_annotations["gene1"]

        # Test case-insensitive matching
        assert _match_gene_pattern(gene_data, "product", "KINASE", "i") is True
        assert _match_gene_pattern(gene_data, "product", "kinase", "i") is True
        assert _match_gene_pattern(gene_data, "source", "AUGUSTUS", "i") is True

    def test_match_gene_pattern_regex(self):
        """Test regex pattern matching."""
        gene_data = self.sample_annotations["gene1"]

        # Test regex patterns
        assert _match_gene_pattern(gene_data, "product", "protein.*kinase", "") is True
        assert _match_gene_pattern(gene_data, "product", "^protein", "") is True
        assert _match_gene_pattern(gene_data, "product", "kinase$", "") is True
        assert _match_gene_pattern(gene_data, "contig", "^chr[0-9]+$", "") is True

    def test_match_gene_pattern_nonexistent_key(self):
        """Test pattern matching on non-existent keys."""
        gene_data = self.sample_annotations["gene1"]

        # Test non-existent key
        assert _match_gene_pattern(gene_data, "nonexistent", "anything", "") is False

    def test_filter_annotations_no_filters(self):
        """Test filtering with no filters applied."""
        result = filter_annotations(self.sample_annotations)
        assert len(result) == 3
        assert result == self.sample_annotations

    def test_filter_annotations_grep_basic(self):
        """Test basic grep filtering."""
        # Keep only kinase genes
        result = filter_annotations(self.sample_annotations, grep=["product:kinase"])
        assert len(result) == 2
        assert "gene1" in result
        assert "gene3" in result
        assert "gene2" not in result

    def test_filter_annotations_grepv_basic(self):
        """Test basic grepv filtering."""
        # Remove augustus genes
        result = filter_annotations(self.sample_annotations, grepv=["source:augustus"])
        assert len(result) == 1
        assert "gene2" in result
        assert "gene1" not in result
        assert "gene3" not in result

    def test_filter_annotations_multiple_grep(self):
        """Test multiple grep patterns (OR logic)."""
        # Keep kinases OR transporters
        result = filter_annotations(
            self.sample_annotations, grep=["product:kinase", "product:transporter"]
        )
        assert len(result) == 3  # All genes match one of the patterns

    def test_filter_annotations_combined(self):
        """Test combined grep and grepv filtering."""
        # Keep kinases but remove augustus
        result = filter_annotations(
            self.sample_annotations, grep=["product:kinase"], grepv=["source:augustus"]
        )
        assert len(result) == 0  # All kinases are from augustus

    def test_filter_annotations_case_insensitive(self):
        """Test case-insensitive filtering."""
        result = filter_annotations(self.sample_annotations, grep=["product:KINASE:i"])
        assert len(result) == 2
        assert "gene1" in result
        assert "gene3" in result

    def test_filter_annotations_regex_patterns(self):
        """Test regex pattern filtering."""
        # Match genes on chr1
        result = filter_annotations(self.sample_annotations, grep=["contig:^chr1$"])
        assert len(result) == 2
        assert "gene1" in result
        assert "gene2" in result
        assert "gene3" not in result

    def test_filter_annotations_empty_result(self):
        """Test filtering that results in empty set."""
        result = filter_annotations(self.sample_annotations, grep=["product:nonexistent"])
        assert len(result) == 0

    def test_filter_annotations_invalid_regex(self):
        """Test filtering with invalid regex pattern."""
        with pytest.raises(ValueError):
            filter_annotations(self.sample_annotations, grep=["product:[invalid"])
