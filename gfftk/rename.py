from .gff import gff2dict, dict2gff3
from natsort import natsorted
from collections import OrderedDict


def rename(args):
    # load into dictionary
    Genes = gff2dict(args.gff3, args.fasta)

    # now create ordered dictionary and sort by contig and position
    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    sGenes = natsorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    counter = args.numbering
    args.locus_tag = args.locus_tag.rstrip("_")
    transcripts = 0
    for k, v in list(sortedGenes.items()):
        locusTag = args.locus_tag + "_" + str(counter).zfill(6)
        renamedGenes[locusTag] = v
        renamedGenes[locusTag]["gene_synonym"].append(k)
        newIds = []
        for i in range(0, len(v["ids"])):
            newIds.append("{}-T{}".format(locusTag, i + 1))
            transcripts += 1
        renamedGenes[locusTag]["ids"] = newIds
        counter += 1

    # write to gff3
    dict2gff3(renamedGenes, output=args.out)
