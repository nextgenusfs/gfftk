"""
Unit tests for the convert module.
"""

import os
import tempfile
from unittest.mock import patch

from gfftk.convert import gff2gtf, gff2proteins, gff2transcripts


class TestConvert:
    """Tests for conversion functions."""

    def test_convert_functions_exist(self):
        """Test that the convert functions exist and can be imported."""
        # Just check that the functions exist
        assert callable(gff2gtf)
        assert callable(gff2proteins)
        assert callable(gff2transcripts)


class TestConvertFiltering:
    """Tests for convert functions with filtering."""

    def setup_method(self):
        """Set up test data for convert filtering tests."""
        # Mock annotation data that would be returned by gff2dict
        self.mock_annotations = {
            "gene1": {
                "name": "kinase1",
                "product": ["protein kinase"],
                "source": "augustus",
                "contig": "chr1",
                "strand": "+",
                "location": (100, 500),
                "type": ["mRNA"],
                "note": ["enzyme"],
                "ids": ["gene1-T1"],
                "CDS": [[(100, 200), (300, 500)]],
                "mRNA": [[(100, 500)]],
                "protein": ["MKTEST*"],
                "transcript": ["ATGAAATESTTAGTAG"],
                "cds_transcript": ["ATGAAATEST"],
                "partialStart": [False],
                "partialStop": [False],
                "pseudo": False,
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
                "ids": ["gene2-T1"],
                "CDS": [[(1000, 1200), (1300, 1500)]],
                "mRNA": [[(1000, 1500)]],
                "protein": ["MTRANS*"],
                "transcript": ["ATGTRANSTAGCTAG"],
                "cds_transcript": ["ATGTRANS"],
                "partialStart": [False],
                "partialStop": [False],
                "pseudo": False,
            },
        }

    @patch("gfftk.convert.gff2dict")
    @patch("gfftk.convert.dict2gtf")
    def test_gff2gtf_with_grep_filtering(self, mock_dict2gtf, mock_gff2dict):
        """Test gff2gtf with grep filtering."""
        # Setup mocks
        mock_gff2dict.return_value = self.mock_annotations

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as gff_file:
            gff_file.write("##gff-version 3\n")
            gff_file_path = gff_file.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_file:
            fasta_file.write(">chr1\nATGCGCGCGCGC\n")
            fasta_file_path = fasta_file.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False) as output_file:
            output_file_path = output_file.name

        try:
            # Test with grep filter
            gff2gtf(
                gff_file_path,
                fasta_file_path,
                output=output_file_path,
                grep=["product:kinase"],
            )

            # Verify gff2dict was called
            mock_gff2dict.assert_called_once()

            # Verify dict2gtf was called with filtered data
            mock_dict2gtf.assert_called_once()
            filtered_data = mock_dict2gtf.call_args[0][0]

            # Should only contain kinase gene
            assert len(filtered_data) == 1
            assert "gene1" in filtered_data
            assert "gene2" not in filtered_data

        finally:
            # Clean up
            for path in [gff_file_path, fasta_file_path, output_file_path]:
                if os.path.exists(path):
                    os.unlink(path)

    @patch("gfftk.convert.gff2dict")
    @patch("gfftk.convert.dict2gtf")
    def test_gff2gtf_with_grepv_filtering(self, mock_dict2gtf, mock_gff2dict):
        """Test gff2gtf with grepv filtering."""
        # Setup mocks
        mock_gff2dict.return_value = self.mock_annotations

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as gff_file:
            gff_file.write("##gff-version 3\n")
            gff_file_path = gff_file.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_file:
            fasta_file.write(">chr1\nATGCGCGCGCGC\n")
            fasta_file_path = fasta_file.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False) as output_file:
            output_file_path = output_file.name

        try:
            # Test with grepv filter
            gff2gtf(
                gff_file_path,
                fasta_file_path,
                output=output_file_path,
                grepv=["source:augustus"],
            )

            # Verify dict2gtf was called with filtered data
            mock_dict2gtf.assert_called_once()
            filtered_data = mock_dict2gtf.call_args[0][0]

            # Should only contain non-augustus gene
            assert len(filtered_data) == 1
            assert "gene2" in filtered_data
            assert "gene1" not in filtered_data

        finally:
            # Clean up
            for path in [gff_file_path, fasta_file_path, output_file_path]:
                if os.path.exists(path):
                    os.unlink(path)

    @patch("gfftk.convert.gff2dict")
    @patch("gfftk.convert._dict2transcripts")
    def test_gff2transcripts_with_filtering(self, mock_dict2transcripts, mock_gff2dict):
        """Test gff2transcripts with filtering."""
        # Setup mocks
        mock_gff2dict.return_value = self.mock_annotations

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as gff_file:
            gff_file.write("##gff-version 3\n")
            gff_file_path = gff_file.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_file:
            fasta_file.write(">chr1\nATGCGCGCGCGC\n")
            fasta_file_path = fasta_file.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as output_file:
            output_file_path = output_file.name

        try:
            # Test with combined filtering
            gff2transcripts(
                gff_file_path,
                fasta_file_path,
                output=output_file_path,
                grep=["product:transporter"],
                grepv=["source:augustus"],
            )

            # Verify _dict2transcripts was called with filtered data
            mock_dict2transcripts.assert_called_once()
            filtered_data = mock_dict2transcripts.call_args[0][0]

            # Should only contain transporter gene (not augustus)
            assert len(filtered_data) == 1
            assert "gene2" in filtered_data

        finally:
            # Clean up
            for path in [gff_file_path, fasta_file_path, output_file_path]:
                if os.path.exists(path):
                    os.unlink(path)

    def test_convert_functions_accept_filter_params(self):
        """Test that convert functions accept grep and grepv parameters."""
        # Test that functions accept the parameters without error
        # (We can't easily test the full functionality without complex setup)

        # Check function signatures include grep and grepv
        import inspect

        # Test gff2gtf signature
        sig = inspect.signature(gff2gtf)
        assert "grep" in sig.parameters
        assert "grepv" in sig.parameters
        assert sig.parameters["grep"].default == []
        assert sig.parameters["grepv"].default == []

        # Test gff2transcripts signature
        sig = inspect.signature(gff2transcripts)
        assert "grep" in sig.parameters
        assert "grepv" in sig.parameters
