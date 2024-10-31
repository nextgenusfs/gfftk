from .gff import gff2dict, dict2gff3


def sanitize(args):
    """
    Sanitize GFF3 and FASTA files by converting them to a standardized format and back.

    This function parses GFF3 and FASTA files using `gff2dict` to convert them into a standardized dictionary format. It then converts the dictionary back to GFF3 format using `dict2gff3`, effectively sanitizing the input files.

    Parameters:
        args (argparse.Namespace): Command-line arguments containing paths to GFF3 and FASTA files, output file path, and an optional debug flag.

    Returns:
        None
    """
    Genes = gff2dict(args.gff3, args.fasta, debug=args.debug)
    dict2gff3(Genes, output=args.out, debug=args.debug)
