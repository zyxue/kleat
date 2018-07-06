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

    'ctg_hex', 'ctg_hex_id', 'ctg_hex_pos',
    'ref_hex', 'ref_hex_id', 'ref_hex_pos'
]


ClvRecord = namedtuple('ClvRecord', HEADER)


CANDIDATE_HEXAMERS = [
    ('AATAAA', 1),
    ('ATTAAA', 2),
    ('AGTAAA', 3),
    ('TATAAA', 4),
    ('CATAAA', 5),
    ('GATAAA', 6),
    ('AATATA', 7),
    ('AATACA', 8),
    ('AATAGA', 9),
    ('AAAAAG', 10),
    ('ACTAAA', 11),
    ('AAGAAA', 12),
    ('AATGAA', 13),
    ('TTTAAA', 14),
    ('AAAACA', 15),
    ('GGGGCT', 16),

    # 6 more from QAPA paper: but not sure of their strength (TODO)
    # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1414-4
    # ('AACAAA', 17),
    # ('AACAAG', 18),
    # ('AATAAG', 19),
    # ('AATAAT', 20),
    # ('ATTACA', 21),
    # ('ATTATA', 22)
]


UCSC_SEQNAMES = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22',
    'chrX', 'chrY', 'chrM'
]

ENSEMBL_SEQNAMES = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
    '21', '22',
    'X', 'Y', 'MT'
]

UCSC_TO_ENSEMBL_SEQNAME = dict(zip(UCSC_SEQNAMES, ENSEMBL_SEQNAMES))
ENSEMBL_TO_UCSC_SEQNAME = dict(zip(ENSEMBL_SEQNAMES, UCSC_SEQNAMES))


# # kept for reference
# COMPLEMENT_DICT = str.maketrans("ACTG", "TGAC")


# For later use
# # used when merging polyA evidence
# COLUMS_TO_SUM = [
#     'num_suffix_reads',
#     'num_bridge_reads',
#     'num_link_reads',
#     'num_blank_contigs'
# ]

# COLUMS_TO_MAX = [
#     'contig_len',
#     'contig_mapq',
#     'suffix_contig_tail_len',
#     'max_bridge_read_tail_len'
# ]
