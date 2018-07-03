from collections import namedtuple


# http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0                  # M
BAM_CINS = 1                    # I
BAM_CDEL = 2                    # D
BAM_CREF_SKIP = 3               # N
BAM_CSOFT_CLIP = 4              # S
BAM_CHARD_CLIP = 5              # H
BAM_CPAD = 6                    # P
BAM_CEQUAL = 7                  # =
BAM_CDIFF = 8                   # X
BAM_CBACK = 9                   # B


HEADER = [
    'seqname', 'strand', 'clv',
    'evidence_type', 'contig_id', 'contig_len', 'contig_mapq',
    'num_suffix_reads', 'suffix_contig_tail_len',    # suffix
    'num_bridge_reads', 'max_bridge_read_tail_len',  # bridge
    'num_link_reads',                                # link
    'num_blank_contigs',                             # blank
]


ClvRecord = namedtuple('ClvRecord', HEADER)
