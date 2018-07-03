"""
Utilities for analyzing suffix-specific evidence. Contigs here are always
assumed to be suffix contigs
"""

from kleat.misc import apautils
from kleat.misc.settings import ClvRecord


def calc_num_suffix_reads(r2c_bam, suffix_contig, ref_clv):
    """
    calculate the number of reads aligned to the cleavage site

    https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
    """
    num_tail_reads = r2c_bam.count(
        # region is half-open
        # https://pysam.readthedocs.io/en/latest/glossary.html#term-region
        suffix_contig.query_name, ref_clv, ref_clv + 1
    )
    return num_tail_reads


def gen_clv_record(suffix_contig, r2c_bam, tail_side):
    ref_clv = apautils.calc_ref_clv(suffix_contig, tail_side)
    tail_len = apautils.calc_tail_length(suffix_contig, tail_side)

    num_suffix_reads = calc_num_suffix_reads(r2c_bam, suffix_contig, ref_clv)

    return ClvRecord(
        suffix_contig.reference_name,
        apautils.calc_strand(tail_side),
        ref_clv,

        'suffix',
        suffix_contig.query_name,
        suffix_contig.query_length,
        suffix_contig.mapq,

        num_suffix_reads,
        tail_len,

        # other types of evidence are left empty
        num_bridge_reads=0,
        max_bridge_read_tail_length=0,
        num_link_reads=0,
        num_blank_contigs=0
    )
