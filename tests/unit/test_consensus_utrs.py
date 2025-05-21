"""
Unit tests for the UTR extension functionality in the consensus module.
"""

from unittest.mock import MagicMock, patch

from gfftk.consensus import check_intron_compatibility, extend_utrs, select_best_utrs


class TestUTRExtension:
    """Tests for the UTR extension functionality in the consensus module."""

    def setup_method(self):
        """Set up test data for each test method."""
        # Create a simple consensus model
        self.consensus_models = {
            "gene1": {
                "contig": "contig1",
                "strand": "+",
                "location": (1000, 2000),
                "coords": [(1000, 1200), (1400, 1600), (1800, 2000)],
                "score": 10,
                "source": "test_source",
            },
            "gene2": {
                "contig": "contig1",
                "strand": "-",
                "location": (3000, 4000),
                "coords": [(3000, 3200), (3400, 3600), (3800, 4000)],
                "score": 10,
                "source": "test_source",
            },
        }

        # Create transcript alignments
        self.transcript_alignments = {
            "transcript1": {
                "contig": "contig1",
                "strand": "+",
                "coords": [(900, 1200), (1400, 1600), (1800, 2100)],
                "source": "transcript",
            },
            "transcript2": {
                "contig": "contig1",
                "strand": "+",
                "coords": [(950, 1200), (1400, 1600), (1800, 2050)],
                "source": "transcript",
            },
            "transcript3": {
                "contig": "contig1",
                "strand": "-",
                "coords": [(2900, 3200), (3400, 3600), (3800, 4100)],
                "source": "transcript",
            },
        }

        # Mock genome FASTA
        self.genome_fasta = "/Users/jon/software/gfftk/tests/unit/test_data/test_genome.fa"

    def test_check_intron_compatibility(self):
        """Test the check_intron_compatibility function."""
        # Test case 1: Compatible introns
        model_coords = [(1000, 1200), (1400, 1600), (1800, 2000)]
        transcript_coords = [(900, 1200), (1400, 1600), (1800, 2100)]
        assert check_intron_compatibility(model_coords, transcript_coords, "+") is True

        # Test case 2: Incompatible introns (different boundaries)
        model_coords = [(1000, 1200), (1400, 1600), (1800, 2000)]
        transcript_coords = [(900, 1200), (1450, 1600), (1800, 2100)]
        assert check_intron_compatibility(model_coords, transcript_coords, "+") is False

        # Test case 3: No introns in model
        model_coords = [(1000, 2000)]
        transcript_coords = [(900, 2100)]
        assert check_intron_compatibility(model_coords, transcript_coords, "+") is True

        # Test case 4: No introns in transcript but introns in model
        model_coords = [(1000, 1200), (1400, 1600), (1800, 2000)]
        transcript_coords = [(900, 2100)]
        assert check_intron_compatibility(model_coords, transcript_coords, "+") is False

    def test_select_best_utrs(self):
        """Test the select_best_utrs function."""
        # Test case 1: Single UTR set
        utr_exons = [[(900, 950), (960, 990)]]
        best_utrs, method = select_best_utrs(utr_exons, "+")
        assert best_utrs == [(900, 950), (960, 990)]
        assert method == "single_transcript"

        # Test case 2: Multiple UTR sets with same exon count
        utr_exons = [[(900, 950), (960, 990)], [(880, 940), (950, 980)]]
        best_utrs, method = select_best_utrs(utr_exons, "+")
        assert method in ["median_length", "common_exon_count", "consensus"]

        # Test case 3: Empty UTR list
        best_utrs, method = select_best_utrs([], "+")
        assert best_utrs is None
        assert method == "none"

        # Test case 4: UTRs that don't meet length constraints
        utr_exons = [[(900, 905)]]  # 6 bp, less than default min_length of 10
        best_utrs, method = select_best_utrs(utr_exons, "+")
        assert best_utrs is None
        assert method == "none"

    @patch("gfftk.consensus.FASTA")
    def test_extend_utrs_plus_strand(self, mock_fasta):
        """Test UTR extension for plus strand genes."""
        # Mock the FASTA class
        mock_fasta_instance = MagicMock()
        mock_fasta.return_value = mock_fasta_instance

        # Run the function
        result = extend_utrs(
            self.consensus_models,
            self.transcript_alignments,
            "/Users/jon/software/gfftk/tests/unit/test_data/test_genome.fa",
            min_utr_length=10,
            max_utr_length=2000,
            log=lambda _: None,
        )

        # Check that UTRs were added to gene1
        assert "gene1" in result
        assert "five_prime_utr" in result["gene1"]
        assert "three_prime_utr" in result["gene1"]

        # Check that the gene boundaries were extended or UTRs were added
        assert result["gene1"]["location"][0] <= self.consensus_models["gene1"]["location"][0]
        assert result["gene1"]["location"][1] >= self.consensus_models["gene1"]["location"][1]

        # Check that CDS coordinates exist
        assert "cds" in result["gene1"]
        # The CDS coordinates might be adjusted based on transcript evidence
        # Just check that they exist and have the same number of exons
        assert len(result["gene1"]["cds"]) == len(self.consensus_models["gene1"]["coords"])

        # Check that the mRNA coordinates include UTRs
        assert len(result["gene1"]["coords"]) >= len(self.consensus_models["gene1"]["coords"])

    @patch("gfftk.consensus.FASTA")
    def test_extend_utrs_minus_strand(self, mock_fasta):
        """Test UTR extension for minus strand genes."""
        # Mock the FASTA class
        mock_fasta_instance = MagicMock()
        mock_fasta.return_value = mock_fasta_instance

        # Run the function
        result = extend_utrs(
            self.consensus_models,
            self.transcript_alignments,
            "/Users/jon/software/gfftk/tests/unit/test_data/test_genome.fa",
            min_utr_length=10,
            max_utr_length=2000,
            log=lambda _: None,
        )

        # Check that UTRs were added to gene2
        assert "gene2" in result
        assert "five_prime_utr" in result["gene2"] or "three_prime_utr" in result["gene2"]

        # Check that the gene boundaries were extended or UTRs were added
        assert (
            result["gene2"]["location"][0] <= self.consensus_models["gene2"]["location"][0]
            or result["gene2"]["location"][1] >= self.consensus_models["gene2"]["location"][1]
        )

        # Check that CDS coordinates exist
        assert "cds" in result["gene2"]
        # The CDS coordinates might be adjusted based on transcript evidence
        # Just check that they exist and have the same number of exons
        assert len(result["gene2"]["cds"]) == len(self.consensus_models["gene2"]["coords"])

    @patch("gfftk.consensus.FASTA")
    def test_extend_utrs_no_transcripts(self, mock_fasta):
        """Test UTR extension with no matching transcripts."""
        # Mock the FASTA class
        mock_fasta_instance = MagicMock()
        mock_fasta.return_value = mock_fasta_instance

        # Create a gene with no matching transcripts
        models = {
            "gene3": {
                "contig": "contig2",  # Different contig
                "strand": "+",
                "location": (1000, 2000),
                "coords": [(1000, 1200), (1400, 1600), (1800, 2000)],
                "score": 10,
                "source": "test_source",
            }
        }

        # Run the function
        result = extend_utrs(
            models,
            self.transcript_alignments,
            "/Users/jon/software/gfftk/tests/unit/test_data/test_genome.fa",
            min_utr_length=10,
            max_utr_length=2000,
            log=lambda _: None,
        )

        # Check that no UTRs were added
        assert "gene3" in result
        assert "five_prime_utr" not in result["gene3"]
        assert "three_prime_utr" not in result["gene3"]

        # Check that CDS coordinates were copied
        assert "cds" in result["gene3"]
        assert result["gene3"]["cds"] == models["gene3"]["coords"]

    @patch("gfftk.consensus.FASTA")
    def test_extend_utrs_incompatible_transcripts(self, mock_fasta):
        """Test UTR extension with incompatible transcripts."""
        # Mock the FASTA class
        mock_fasta_instance = MagicMock()
        mock_fasta.return_value = mock_fasta_instance

        # Add an incompatible transcript
        self.transcript_alignments["incompatible"] = {
            "contig": "contig1",
            "strand": "+",
            "coords": [
                (900, 1200),
                (1450, 1600),
                (1800, 2100),
            ],  # Different intron boundary
            "source": "transcript",
        }

        # Run the function
        result = extend_utrs(
            self.consensus_models,
            self.transcript_alignments,
            "/Users/jon/software/gfftk/tests/unit/test_data/test_genome.fa",
            min_utr_length=10,
            max_utr_length=2000,
            log=lambda _: None,
        )

        # Check that UTRs were still added to gene1 from compatible transcripts
        assert "gene1" in result
        assert "five_prime_utr" in result["gene1"]
        assert "three_prime_utr" in result["gene1"]
