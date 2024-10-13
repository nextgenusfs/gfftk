from .gff import gff2dict, dict2gff3


def sanitize(args):
    Genes = gff2dict(args.gff3, args.fasta, debug=args.debug)
    dict2gff3(Genes, output=args.out, debug=args.debug)
