from collections import OrderedDict
from .utils import zopen
import io


codon_table = {
    1: {
        "start": ["ATG"],
        "table": {
            "TTT": "F",
            "TTC": "F",
            "TTA": "L",
            "TTG": "L",
            "TCT": "S",
            "TCC": "S",
            "TCA": "S",
            "TCG": "S",
            "TAT": "Y",
            "TAC": "Y",
            "TGT": "C",
            "TGC": "C",
            "TGG": "W",
            "CTT": "L",
            "CTC": "L",
            "CTA": "L",
            "CTG": "L",
            "CCT": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "CAT": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGT": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "ATT": "I",
            "ATC": "I",
            "ATA": "I",
            "ATG": "M",
            "ACT": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "AAT": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "AGT": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GTT": "V",
            "GTC": "V",
            "GTA": "V",
            "GTG": "V",
            "GCT": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",
            "GAT": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",
            "GGT": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
            "TAA": "*",
            "TAG": "*",
            "TGA": "*",
        },
    },
    11: {
        "start": ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
        "table": {
            "TTT": "F",
            "TTC": "F",
            "TTA": "L",
            "TTG": "L",
            "TCT": "S",
            "TCC": "S",
            "TCA": "S",
            "TCG": "S",
            "TAT": "Y",
            "TAC": "Y",
            "TGT": "C",
            "TGC": "C",
            "TGG": "W",
            "CTT": "L",
            "CTC": "L",
            "CTA": "L",
            "CTG": "L",
            "CCT": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "CAT": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGT": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "ATT": "I",
            "ATC": "I",
            "ATA": "I",
            "ATG": "M",
            "ACT": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "AAT": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "AGT": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GTT": "V",
            "GTC": "V",
            "GTA": "V",
            "GTG": "V",
            "GCT": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",
            "GAT": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",
            "GGT": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
            "TAA": "*",
            "TAG": "*",
            "TGA": "*",
        },
    },
}


def RevComp(s):
    rev_comp_lib = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "U": "A",
        "M": "K",
        "R": "Y",
        "W": "W",
        "S": "S",
        "Y": "R",
        "K": "M",
        "V": "B",
        "H": "D",
        "D": "H",
        "B": "V",
        "X": "X",
        "N": "N",
    }
    cseq = ""
    n = len(s)
    s = s.upper()
    for i in range(0, n):
        c = s[n - i - 1]
        cseq += rev_comp_lib[c]
    return cseq


def translate(dna, strand, phase, table=1):
    """Translates DNA sequence into proteins.

    Takes DNA (or rather cDNA sequence) and translates to proteins/amino acids.
    It requires the DNA sequence, the strand, translation phase, and translation
    table.

    Parameters
    ----------
    dna : str
        DNA (cDNA) sequence as nucleotides
    strand : str, (+/-)
        strand to translate (+ or -)
    phase : int
        phase to start translation [0,1,2]
    table : int, default=1
        translation table [1]

    Returns
    -------
    protSeq : str
        string of translated amino acid sequence

    """

    def _split(str, num):
        return [str[start : start + num] for start in range(0, len(str), num)]

    if strand == "-" or strand == -1:
        seq = RevComp(dna)
    else:
        seq = dna
    seq = seq[phase:]
    # map seq to proteins
    protSeq = []
    for c, i in enumerate(_split(seq, 3)):
        if len(i) == 3:
            iSeq = i.upper()
            if c == 0:  # first codon
                if iSeq in codon_table[table]["start"]:
                    aa = "M"
                    protSeq.append(aa)
                else:
                    if iSeq in codon_table[table]["table"]:
                        aa = codon_table[table]["table"][iSeq]
                        protSeq.append(aa)
                    else:
                        protSeq.append("X")
            else:
                if iSeq in codon_table[table]["table"]:
                    aa = codon_table[table]["table"][iSeq]
                    protSeq.append(aa)
                else:
                    protSeq.append("X")
    return "".join(protSeq)


def fastaparser(handle):
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if isinstance(line, bytes):
            line = line.decode()
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return
    lines = []
    for line in handle:
        if isinstance(line, bytes):
            line = line.decode()
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())
    yield title, "".join(lines).replace(" ", "").replace("\r", "")


def fasta2dict(fasta, full_header=False):
    """Read FASTA file to dictionary.

    This is same as biopython SeqIO.to_dict(), return dictionary
    keyed by contig name and value is the sequence string.

    Parameters
    ----------
    fasta : filename
        FASTA input file (can be gzipped)
    full_header : bool, default=False
        return full header for contig names, default is split at first space

    Returns
    -------
    seqs : dict
        returns OrderedDict() of header: seq

    """
    seqs = OrderedDict()
    if isinstance(fasta, io.BytesIO):
        fasta.seek(0)
        infile = fasta
    else:
        infile = zopen(fasta)
    for title, seq in fastaparser(infile):
        if full_header:
            title = title
        else:
            title = title.split()[0]
        seqs[title] = seq
    if not isinstance(fasta, io.BytesIO):
        infile.close()
    return seqs


def fasta2headers(fasta, full_header=False):
    """Read FASTA file set of headers.

    Simple function to read FASTA file and return set of contig names

    Parameters
    ----------
    fasta : filename
        FASTA input file (can be gzipped)
    full_header : bool, default=False
        return full header for contig names, default is split at first space

    Returns
    -------
    headers : set
        returns set() of header names

    """
    # generate a set of the contig/scaffold names
    headers = set()
    if isinstance(fasta, io.BytesIO):
        fasta.seek(0)
        infile = fasta
    else:
        infile = zopen(fasta)
    for title, seq in fastaparser(infile):
        if full_header:
            title = title
        else:
            title = title.split()[0]
        headers.add(title)
    if not isinstance(fasta, io.BytesIO):
        infile.close()
    return headers


def fasta2lengths(fasta, full_header=False):
    """Read FASTA file to dictionary of sequence lengths.

    Reads FASTA file (optionally gzipped) and returns dictionary of
    contig header names as keys with length of sequences as values

    Parameters
    ----------
    fasta : filename
        FASTA input file (can be gzipped)
    full_header : bool, default=False
        return full header for contig names, default is split at first space

    Returns
    -------
    seqs : dict
        returns dictionary of header: len(seq)

    """
    seqs = {}
    if isinstance(fasta, io.BytesIO):
        fasta.seek(0)
        infile = fasta
    else:
        infile = zopen(fasta)
    for title, seq in fastaparser(infile):
        if full_header:
            title = title
        else:
            title = title.split()[0]
        seqs[title] = len(seq)
    if not isinstance(fasta, io.BytesIO):
        infile.close()
    return seqs


def getSeqRegions(seqs, header, coordinates, coords=False):
    """From sequence dictionary return spliced coordinates.

    Takes a sequence dictionary (ie from fasta2dict), the contig name (header)
    and the coordinates to fetch (list of tuples)

    Parameters
    ----------
    seqs : dict
        dictionary of sequences keyed by contig name/ header
    header : str
        contig name (header) for sequence in seqs dictionary
    coordinates : list of tuples
        list of tuples of sequence coordinates to return [(1,10), (20,30)]

    Returns
    -------
    result : str
        returns spliced DNA sequence

    """
    # takes SeqRecord dictionary or Index, returns sequence string
    # coordinates is a list of tuples [(1,10), (20,30)]
    result = []
    sorted_coordinates = sorted(coordinates, key=lambda tup: tup[0])
    for x in sorted_coordinates:
        result.append(seqs[header][x[0] - 1 : x[1]])
    if coords:
        return "".join(result), result
    else:
        return "".join(result)


def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i : i + every])
    return "\n".join(lines)
