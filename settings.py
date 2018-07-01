from collections import namedtuple


# https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0
BAM_CSOFT_CLIP = 4              # used to identify soft clip

HEADER = [
    'seqname', 'strand', 'clv',
    'evidence_type', 'contig_id', 'contig_length', 'contig_mapq',
    'num_suffix_reads', 'suffix_contig_tail_length',    # suffix
    'num_bridge_reads', 'max_bridge_read_tail_length',  # bridge
    'num_link_reads',                                   # link
    'num_blank_contigs',                                # blank

]


ClvRecord = namedtuple('ClvRecord', HEADER)
