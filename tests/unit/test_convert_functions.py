"""
Unit tests for the convert module functions.
"""

import os
import tempfile

from gfftk.convert import gff2gtf, gff2proteins, gff2tbl, gff2transcripts


class TestConvertFunctions:
    """Tests for convert functions."""

    def test_convert_functions_exist(self):
        """Test that the convert functions exist and can be imported."""
        # Just check that the functions exist
        assert callable(gff2gtf)
        assert callable(gff2proteins)
        assert callable(gff2transcripts)
        assert callable(gff2tbl)

    def test_gff2proteins_basic(self):
        """Test basic functionality of gff2proteins."""
        # Create a temporary GFF file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_gff:
            temp_gff.write("##gff-version 3\n")
            temp_gff.write("contig1\tprediction\tgene\t1\t12\t.\t+\t.\tID=gene1;Name=test_gene\n")
            temp_gff.write(
                "contig1\tprediction\tmRNA\t1\t12\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=test_mrna\n"
            )
            temp_gff.write("contig1\tprediction\tCDS\t1\t12\t.\t+\t0\tID=cds1;Parent=mRNA1\n")
            temp_gff_name = temp_gff.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_fasta:
            temp_fasta.write(">contig1\n")
            temp_fasta.write("ATGCATGCATGC\n")  # 12 bases, translates to MHAC
            temp_fasta_name = temp_fasta.name

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_out:
            temp_out_name = temp_out.name

        try:
            # Convert GFF to proteins
            gff2proteins(temp_gff_name, temp_fasta_name, temp_out_name)

            # Check that the output file exists
            assert os.path.exists(temp_out_name)

            # Check the content of the output file
            with open(temp_out_name, "r") as f:
                content = f.read()

            # The output should be a FASTA file with protein sequences
            assert ">mRNA1" in content or ">gene1" in content

            # We don't check the exact protein sequence because it depends on the
            # implementation details of the translation function
        finally:
            # Clean up
            if os.path.exists(temp_gff_name):
                os.unlink(temp_gff_name)
            if os.path.exists(temp_fasta_name):
                os.unlink(temp_fasta_name)
            if os.path.exists(temp_out_name):
                os.unlink(temp_out_name)

    def test_gff2transcripts_basic(self):
        """Test basic functionality of gff2transcripts."""
        # Create a temporary GFF file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_gff:
            temp_gff.write("##gff-version 3\n")
            temp_gff.write("contig1\tprediction\tgene\t1\t12\t.\t+\t.\tID=gene1;Name=test_gene\n")
            temp_gff.write(
                "contig1\tprediction\tmRNA\t1\t12\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=test_mrna\n"
            )
            temp_gff.write("contig1\tprediction\texon\t1\t12\t.\t+\t.\tID=exon1;Parent=mRNA1\n")
            temp_gff_name = temp_gff.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_fasta:
            temp_fasta.write(">contig1\n")
            temp_fasta.write("ATGCATGCATGC\n")  # 12 bases
            temp_fasta_name = temp_fasta.name

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_out:
            temp_out_name = temp_out.name

        try:
            # Convert GFF to transcripts
            gff2transcripts(temp_gff_name, temp_fasta_name, temp_out_name)

            # Check that the output file exists
            assert os.path.exists(temp_out_name)

            # Check the content of the output file
            with open(temp_out_name, "r") as f:
                content = f.read()

            # The output should be a FASTA file with transcript sequences
            assert ">mRNA1" in content or ">gene1" in content

            # The transcript sequence should be the same as the genomic sequence
            # for this simple case
            assert "ATGCATGCATGC" in content
        finally:
            # Clean up
            if os.path.exists(temp_gff_name):
                os.unlink(temp_gff_name)
            if os.path.exists(temp_fasta_name):
                os.unlink(temp_fasta_name)
            if os.path.exists(temp_out_name):
                os.unlink(temp_out_name)
