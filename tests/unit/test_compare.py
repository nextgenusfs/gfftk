"""
Unit tests for the compare module.
"""

import os
import tempfile
from gfftk.compare import (
    compareAnnotations,
    pairwiseAED,
)


class TestCompare:
    """Tests for compare functions."""

    def test_pairwiseAED(self):
        """Test calculating pairwise AED."""
        # Create sample data
        query = [
            {
                "exon": [[1, 100], [200, 300]],
                "CDS": [[1, 100], [200, 300]],
            }
        ]

        reference = [
            {
                "exon": [[1, 100], [200, 300]],
                "CDS": [[1, 100], [200, 300]],
            }
        ]

        # Calculate pairwise AED
        result = pairwiseAED(query, reference)

        # Check the result
        assert isinstance(result, str)
        # Perfect match should have AED of 0
        assert float(result) == 0.0

    def test_compareAnnotations(self):
        """Test comparing annotations."""
        # Create a reference GFF file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as ref_gff:
            ref_gff.write("##gff-version 3\n")
            ref_gff.write(
                "contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene1\n"
            )
            ref_gff.write(
                "contig1\tprediction\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=test_mrna1\n"
            )
            ref_gff.write(
                "contig1\tprediction\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=mRNA1\n"
            )
            ref_gff.write(
                "contig1\tprediction\texon\t200\t300\t.\t+\t.\tID=exon2;Parent=mRNA1\n"
            )
            ref_gff.write(
                "contig1\tprediction\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1\n"
            )
            ref_gff.write(
                "contig1\tprediction\tCDS\t200\t300\t.\t+\t0\tID=cds2;Parent=mRNA1\n"
            )
            ref_gff_name = ref_gff.name

        # Create a prediction GFF file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as pred_gff:
            pred_gff.write("##gff-version 3\n")
            pred_gff.write(
                "contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene2;Name=test_gene2\n"
            )
            pred_gff.write(
                "contig1\tprediction\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA2;Parent=gene2;Name=test_mrna2\n"
            )
            pred_gff.write(
                "contig1\tprediction\texon\t1\t100\t.\t+\t.\tID=exon3;Parent=mRNA2\n"
            )
            pred_gff.write(
                "contig1\tprediction\texon\t200\t300\t.\t+\t.\tID=exon4;Parent=mRNA2\n"
            )
            pred_gff.write(
                "contig1\tprediction\tCDS\t1\t100\t.\t+\t0\tID=cds3;Parent=mRNA2\n"
            )
            pred_gff.write(
                "contig1\tprediction\tCDS\t200\t300\t.\t+\t0\tID=cds4;Parent=mRNA2\n"
            )
            pred_gff_name = pred_gff.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_fasta:
            temp_fasta.write(">contig1\n")
            temp_fasta.write("A" * 1000 + "\n")
            temp_fasta_name = temp_fasta.name

        try:
            # Compare the annotations
            result = compareAnnotations(ref_gff_name, pred_gff_name, temp_fasta_name)

            # Check the result
            assert isinstance(result, list)
            # The first element should be the number of shared genes
            assert isinstance(result[0], int)
        finally:
            # Clean up
            for filename in [ref_gff_name, pred_gff_name, temp_fasta_name]:
                if os.path.exists(filename):
                    os.unlink(filename)
