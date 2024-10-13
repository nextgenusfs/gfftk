import sys
from .fasta import RevComp
from .utils import zopen


'''
# + strand example
OPO1_006208-T1	2883	0	2883	+	scaffold_60	41779	19205	22545	2883	2883	60	NM:i:0	ms:i:2883	AS:i:2623	nn:i:0	ts:A:+	tp:A:P	cm:i:940	s1:i:2816	s2:i:0	de:f:0	rl:i:0	cs:Z::60~gt63ag:365~gt49ag:520~gt53ag:228~gt51ag:560~gt57ag:97~gt67ag:266~gt64ag:660~gt53ag:127
# cs string cs:Z::60~gt63ag:365~gt49ag:520~gt53ag:228~gt51ag:560~gt57ag:97~gt67ag:266~gt64ag:660~gt53ag:127
scaffold_60	funannotate	exon	19206	19265	.	+	.	ID=OPO1_006208-T1.exon1;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	19329	19693	.	+	.	ID=OPO1_006208-T1.exon2;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	19743	20262	.	+	.	ID=OPO1_006208-T1.exon3;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	20316	20543	.	+	.	ID=OPO1_006208-T1.exon4;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	20595	21154	.	+	.	ID=OPO1_006208-T1.exon5;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	21212	21308	.	+	.	ID=OPO1_006208-T1.exon6;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	21376	21641	.	+	.	ID=OPO1_006208-T1.exon7;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	21706	22365	.	+	.	ID=OPO1_006208-T1.exon8;Parent=OPO1_006208-T1;
scaffold_60	funannotate	exon	22419	22545	.	+	.	ID=OPO1_006208-T1.exon9;Parent=OPO1_006208-T1;

# - strand example
OPO1_006778-T1	3150	0	3150	-	scaffold_82	14472	8938	12444	3150	3152	60	NM:i:2	ms:i:3146	AS:i:2910	nn:i:0	ts:A:+	tp:A:P	cm:i:1015	s1:i:3101	s2:i:0	de:f:0.0003	rl:i:0	cs:Z::377~ct45ac:171~ct46ac:72-ct~ct42ac:546~ct60ac:381~ct64ac:695~ct54ac:830~ct43ac:78
# cs string: cs:Z::377~ct45ac:171~ct46ac:72-ct~ct42ac:546~ct60ac:381~ct64ac:695~ct54ac:830~ct43ac:78
scaffold_82	funannotate	exon	12367	12444	.	-	.	ID=OPO1_006778-T1.exon1;Parent=OPO1_006778-T1;
scaffold_82	funannotate	exon	11494	12323	.	-	.	ID=OPO1_006778-T1.exon2;Parent=OPO1_006778-T1;
scaffold_82	funannotate	exon	10745	11439	.	-	.	ID=OPO1_006778-T1.exon3;Parent=OPO1_006778-T1;
scaffold_82	funannotate	exon	10300	10680	.	-	.	ID=OPO1_006778-T1.exon4;Parent=OPO1_006778-T1;
scaffold_82	funannotate	exon	9694	10239	.	-	.	ID=OPO1_006778-T1.exon5;Parent=OPO1_006778-T1;
scaffold_82	funannotate	exon	9578	9649	.	-	.	ID=OPO1_006778-T1.exon6;Parent=OPO1_006778-T1;
scaffold_82	funannotate	exon	9361	9531	.	-	.	ID=OPO1_006778-T1.exon7;Parent=OPO1_006778-T1;
scaffold_82	funannotate	exon	8939	9315	.	-	.	ID=OPO1_006778-T1.exon8;Parent=OPO1_006778-T1;
'''

def paf2dict(paf, fasta, annotation={}, min_mapq=2):
    # minimap2 paf output needs to be run with --cs flag
    with zopen(paf) as infile:
        for line in infile:
            line = line.rstrip()
            cols = line.split('\t')
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, matches, alen, mapq = cols[:12]
            if int(mapq) <= min_mapq:
                continue
            tags = cols[12:]
            cs = None
            for tag in tags:
                if tag.startswith('cs:Z:'):
                    cs = tag
            if not cs:
                sys.stderr.write('Error: {} file does not contain cs:Z: tags, run minimap2 with the -x splice --cs flags\n'.format(paf))
                return {}
            ref_coords, query_coords, properSplice = cs2coords(int(tstart), qstart, int(qlen), strand, cs)
            annotation[qname] = (ref_coords, query_coords, properSplice, cs)
    return annotation


def cs2coords(start, qstart, length, strand, cs, offset=1,
              splice_donor=['gt', 'at'], splice_acceptor=['ag', 'ac']):
    """
    # From minimap2 manual this is the cs flag definitions
    Op	Regex	Description
    =	[ACGTN]+	Identical sequence (long form)
    :	[0-9]+	Identical sequence length
    *	[acgtn][acgtn]	Substitution: ref to query
    +	[acgtn]+	Insertion to the reference
    -	[acgtn]+	Deletion from the reference
    ~	[acgtn]{2}[0-9]+[acgtn]{2}	Intron length and splice signal
    """
    cs = cs.replace('cs:Z:', '')
    ProperSplice = True
    exons = [int(start)]
    position = int(start)
    query = [int(qstart)]
    querypos = 0
    num_exons = 1
    gaps = 0
    indels = []
    if strand == '+':
        sp_donor = splice_donor
        sp_acceptor = splice_acceptor
        sort_orientation = False
    elif strand == '-':
        # rev comp and swap donor/acceptor
        sp_donor = [RevComp(x).lower() for x in splice_acceptor]
        sp_acceptor = [RevComp(x).lower() for x in splice_donor]
        sort_orientation = True
    for s, value in cs2tuples(cs):
        if s == ':':
            position += int(value)
            querypos += int(value)
            indels.append(0)
        elif s == '-':
            gaps += 1
            position += len(value)
            querypos += len(value)
            indels.append(-len(value))
        elif s == '+':
            gaps += 1
            position += len(value)
            querypos += len(value)
            indels.append(len(value))
        elif s == '~':
            if value.startswith(tuple(sp_donor)) and value.endswith(tuple(sp_acceptor)):
                ProperSplice = True
            else:
                ProperSplice = False
            num_exons += 1
            exons.append(position+indels[-1])
            query.append(querypos)
            query.append(querypos+1)
            intronLen = int(value[2:-2])
            position += intronLen
            exons.append(position)
            indels.append(0)
    # add last Position
    exons.append(position)
    query.append(int(length))
    # convert exon list into list of exon tuples
    exontmp = list(zip(exons[0::2], exons[1::2]))
    queryList = list(zip(query[0::2], query[1::2]))
    exonList = []
    for x in sorted(exontmp, key=lambda tup: tup[0], reverse=sort_orientation):
        exonList.append((x[0]+offset, x[1]))
    return exonList, queryList, ProperSplice


def cs2tuples(cs, separators=[':', '*', '+', '-', '~']):
    # separators is an array of strings that are being used to split the the string.
    tmpList = []
    i = 0
    while i < len(cs):
        theSeparator = ""
        for current in separators:
            if current == cs[i:i+len(current)]:
                theSeparator = current
        if theSeparator != "":
            tmpList += [theSeparator]
            i = i + len(theSeparator)
        else:
            if tmpList == []:
                tmpList = [""]
            if(tmpList[-1] in separators):
                tmpList += [""]
            tmpList[-1] += cs[i]
            i += 1
    tupList = list(zip(tmpList[::2], tmpList[1::2]))
    return tupList
