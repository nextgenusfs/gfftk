import sys
from natsort import natsorted
from .utils import zopen


def sort(args):
    sortGFF3(args.gff3, output=args.out)


def sortGFF3(input, output=False):
    # function to sort GFF3 file but maintain gene, mrna, exon, cds order
    data = []
    features = set()
    comments = []
    with zopen(input) as infile:
        for line in infile:
            if line.startswith("\n"):
                continue
            if line.startswith("#"):
                comments.append(line)
                continue
            line = line.rstrip()
            cols = line.split("\t")
            data.append(cols)
            features.add(cols[2])
    # build sort order dictionary for features
    order_map = {
        "gene": 0,
        "mRNA": 1,
        "transcript": 2,
        "tRNA": 3,
        "ncRNA": 4,
        "rRNA": 5,
        "pseudogene": 6,
        "five_prime_utr": 7,
        "five_prime_UTR": 8,
        "exon": 9,
        "CDS": 10,
        "three_prime_utr": 11,
        "three_prime_UTR": 12,
    }
    idx = len(order_map)
    for x in features:
        if x not in order_map:
            order_map[x] = idx
            idx += 1
    # we can now sort
    sort_data = natsorted(data, key=lambda x: (x[0], int(x[3]), order_map[x[2]]))
    # now we can write back out to file
    if output:
        outfile = zopen(output, mode="w")
    else:
        outfile = sys.stdout
    for y in comments:
        outfile.write(y)
    for x in sort_data:
        outfile.write("{}\n".format("\t".join(x)))
    if output:
        outfile.close()
