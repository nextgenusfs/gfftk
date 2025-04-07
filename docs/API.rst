API
===

.. toctree::
   :maxdepth: 2

   API/convert
   API/gff
   API/genbank
   API/fasta
   API/consensus
   API/stats
   API/go
   API/utils



GFFtk works by parsing annotation files and storing in a python dictionary. After initial parsing
the records are sorted by contig and start position, translated into protein space to test complete
gene models or not, and then output in an a python OrderedDict(). The structure looks like this:

.. code-block:: none

   locustag: {
      'contig': contigName,  #string
      'type': [],  # list of str one for each transcript mRNA/rRNA/tRNA/ncRNA
      'location': (start, end), #integer tuple
      'strand': +/-, #string
      'ids': [transcript/protein IDs], #list
      'mRNA':[[(ex1,ex1),(ex2,ex2)]], #list of lists of tuples (start, end)
      'CDS':[[(cds1,cds1),(cds2,cds2)]], #list of lists of tuples (start, end)
      'transcript': [seq1, seq2], #list of mRNA trnascripts
      'cds_transcript': [seq1, seq2], #list of mRNA trnascripts (no UTRs)
      'protein': [protseq1,protseq2], #list of CDS translations
      'codon_start': [1,1], #codon start for translations
      'note': [[first note, second note], [first, second, etc]], #list of lists
      'name': genename, # str common gene name
      'product': [hypothetical protein, velvet complex], #list of product definitions
      'gene_synonym': [], # list of gene name Aliases
      'EC_number': [[ec number]], # list of lists
      'go_terms': [[GO:0000001,GO:0000002]],  #list of lists
      'db_xref': [[InterPro:IPR0001,PFAM:004384]],  #list of lists
      'partialStart': [bool], # list of True/False for each transcript
      'partialStop': [bootl], $ list of True/False for each transcript
      'source': source, # string annotation source
      'phase': [[0,2,1]], list of lists
      '5UTR': [[(),()]], #list of lists of tuples (start, end)
      '3UTR': [[(),()]] #list of lists of tuples (start, end)
      }
   }
