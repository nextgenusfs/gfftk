#!/usr/bin/env python3

import sys
import os
import argparse
from .__version__ import __version__
from .help_formatter import MyParser, MyHelpFormatter
from .consensus import consensus
from .convert import convert
from .sort import sort
from .sanitize import sanitize
from .rename import rename
from .stats import stats
from .compare import compare


def main():
    args = parse_args(sys.argv[1:])
    if args.subparser_name == 'consensus':
        consensus(args)
    elif args.subparser_name == 'convert':
        convert(args)
    elif args.subparser_name == 'sort':
        sort(args)
    elif args.subparser_name == 'sanitize':
        sanitize(args)
    elif args.subparser_name == 'rename':
        rename(args)
    elif args.subparser_name == 'stats':
        stats(args)
    elif args.subparser_name == 'compare':
        compare(args)


def parse_args(args):
    description = 'GFFtk: tool kit for GFF3 genome annotation manipulation'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)
    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    # add subtools here
    consensus_subparser(subparsers)
    convert_subparser(subparsers)
    sort_subparser(subparsers)
    sanitize_subparser(subparsers)
    rename_subparser(subparsers)
    stats_subparser(subparsers)
    compare_subparser(subparsers)

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def consensus_subparser(subparsers):
    group = subparsers.add_parser('consensus',
                                  description='''EvidenceModeler-like tool to generate consensus gene predictions. All gene models are loaded and sorted into loci based on genomic location, in each locus the gene models are scored based on:|n\
    1) overlap with protein/transcript evidence|n\
    2) AED score in relation to all models at that locus|n\
    3) user supplied source weights|n\
    4) de novo source weights estimated from evidence overlaps|n\
    5) for each locus the gene model with highest score is retained.''',
                                  help='EvidenceModeler-like tool to generate consensus gene predictions.',
                                  formatter_class=MyHelpFormatter,
                                  add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-f', '--fasta', required=True, help='genome in FASTA format', metavar='')
    required_args.add_argument('-g', '--genes', required=True, nargs='+',
                               help='gene model predictions in GFF3 format [accepts multiple files: space separated]',
                               metavar='')
    required_args.add_argument('-o', '--out', required=True, help='output in GFF3 format', metavar='')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('-p', '--proteins', nargs='+', default=[],
                        help='protein alignments in GFF3 format [accepts multiple files: space separated]', metavar='')
    optional_args.add_argument('-t', '--transcripts', nargs='+', default=[],
                        help='transcripts alignments in GFF3 format [accepts multiple files: space separated]', metavar='')
    optional_args.add_argument('-r', '--repeats', nargs='+', default=[],
                        help='repeat alignments in BED or GFF3 format [accepts multiple files: space separated]', metavar='')
    optional_args.add_argument('-w', '--weights', nargs='+', default=[],
                        help='user supplied source weights [accepts multiple: space separated source:weight]', metavar='')
    optional_args.add_argument('-m', '--minscore', type=int,
                        help='minimum score to retain gene model (default: auto)', metavar='')
    optional_args.add_argument('--repeat-overlap', dest='repeat_overlap', default=90, type=int, choices=range(0, 101),
                        help='percent gene model overlap with repeats to remove', metavar='')
    optional_args.add_argument('-l', '--logfile',
                        help='write logs to file', metavar='')
    optional_args.add_argument('--silent', action='store_true',
                        help='do not write anything to terminal/stderr')
    optional_args.add_argument('--debug', action='store_true',
                        help='write/keep intermediate files')

    other_args = group.add_argument_group('Other arguments')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")


def convert_subparser(subparsers):
    group = subparsers.add_parser('convert',
                                  description='convert GFF3/tbl format into another format [output gff3, gtf, tbl, protein fasta]. ',
                                  help='convert GFF3/tbl format into another format [output gff3, gtf, tbl, protein fasta].',
                                  formatter_class=MyHelpFormatter,
                                  add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-f', '--fasta', required=True, help='genome in FASTA format', metavar='')
    required_args.add_argument('-i', '--input', required=True, help='annotation in GFF3 or TBL format', metavar='')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--out', help='write converted output to file (default: stdout)', metavar='')
    optional_args.add_argument('--input-format', dest='input_format', help='format of input file [gff3, tbl]. (default: auto)',
                               choices=['gff3', 'tbl'], metavar='')
    optional_args.add_argument('--output-format', dest='output_format', help='format of output file [gff3, gtf, tbl, proteins, nucleotides]. (default: auto)',
                               choices=['gff3', 'gtf', 'tbl', 'proteins', 'nucleotides'], metavar='')
    optional_args.add_argument('-n', '--no-stop', dest='strip_stop', action='store_true',
                        help='for proteins output, do not write stop codons (*)')
    optional_args.add_argument('--debug', action='store_true',
                        help='write parsing errors to stderr')

    other_args = group.add_argument_group('Other arguments')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")


def sort_subparser(subparsers):
    group = subparsers.add_parser('sort',
                                  description='sort GFF3 file properly [maintain feature order: gene, mrna, exon, cds].',
                                  help='sort GFF3 file properly [maintain feature order: gene, mrna, exon, cds].',
                                  formatter_class=MyHelpFormatter,
                                  add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-g', '--gff3', required=True, help='GFF3 file to sort', metavar='')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--out', help='write sorted output to file (default: stdout)', metavar='')

    other_args = group.add_argument_group('Other arguments')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")

def sanitize_subparser(subparsers):
    group = subparsers.add_parser('sanitize',
                                  description='sanitize GFF3 file, load GFF3 and output cleaned up GFF3 output.',
                                  help='sanitize GFF3 file, load GFF3 and output cleaned up GFF3 output.',
                                  formatter_class=MyHelpFormatter,
                                  add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-f', '--fasta', required=True, help='genome in FASTA format', metavar='')
    required_args.add_argument('-g', '--gff3', required=True, help='annotation in GFF3 format', metavar='')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--out', help='write santized GFF3 output to file (default: stdout)', metavar='')
    optional_args.add_argument('--debug', action='store_true',
                        help='write/keep intermediate files')

    other_args = group.add_argument_group('Other arguments')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")

def rename_subparser(subparsers):
    group = subparsers.add_parser('rename',
                                  description='rename gene models in GFF3 annotation file. Script will sort genes by contig and location and then rename using --locus-tag, ie PREFIX_000001.',
                                  help='rename gene models in GFF3 annotation file.',
                                  formatter_class=MyHelpFormatter,
                                  add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-f', '--fasta', required=True, help='genome in FASTA format', metavar='')
    required_args.add_argument('-g', '--gff3', required=True, help='annotation in GFF3 format', metavar='')
    required_args.add_argument('-l', '--locus-tag', dest='locus_tag', required=True,
                               help='Locus tag for gene names', metavar='')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--out', help='write santized GFF3 output to file (default: stdout)', metavar='')
    optional_args.add_argument('-n', '--numbering', default=1, type=int, help='start locus tag numbering', metavar='')

    other_args = group.add_argument_group('Other arguments')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")


def stats_subparser(subparsers):
    group = subparsers.add_parser('stats',
                                  description='parse annotation GFF3/tbl and output summary statistics.',
                                  help='parse annotation GFF3/tbl and output summary statistics.',
                                  formatter_class=MyHelpFormatter,
                                  add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-f', '--fasta', required=True, help='genome in FASTA format', metavar='')
    required_args.add_argument('-i', '--input', required=True, help='annotation in GFF3 or TBL format', metavar='')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--out', help='write converted output to file (default: stdout)', metavar='')
    optional_args.add_argument('--input-format', dest='input_format', help='format of input file [gff3, tbl]. (default: auto)',
                               choices=['gff3', 'tbl'], metavar='')
    optional_args.add_argument('--debug', action='store_true',
                        help='write parsing errors to stderr')

    other_args = group.add_argument_group('Other arguments')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")


def compare_subparser(subparsers):
    group = subparsers.add_parser('compare',
                                  description='compare two GFF3 annotations of a genome.',
                                  help='compare two GFF3 annotations of a genome.',
                                  formatter_class=MyHelpFormatter,
                                  add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-f', '--fasta', required=True, help='genome in FASTA format', metavar='')
    required_args.add_argument('-q', '--query', required=True, help='query annotation in GFF3 format', metavar='')
    required_args.add_argument('-r', '--reference', required=True, help='query annotation in GFF3 format', metavar='')

    optional_args = group.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--out', help='write converted output to file (default: stdout)', metavar='')
    optional_args.add_argument('--debug', action='store_true',
                        help='write parsing errors to stderr')

    other_args = group.add_argument_group('Other arguments')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")


if __name__ == '__main__':
    main()