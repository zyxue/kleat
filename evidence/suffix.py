"""
Utilities for analyzing suffix-specific evidence. Contigs here are always
assumed to be suffix contigs
"""


from utils.tail import calc_strand, calc_ref_clv, calc_tail_length
from settings import ClvRecord


def calc_num_suffix_reads(r2c_sam, suffix_contig, ref_clv=None):
    """
    calculate the number of reads aligned to the cleavage site

    https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
    """
    if ref_clv is None:
        ref_clv = calc_ref_clv(suffix_contig)

    num_tail_reads = r2c_sam.count(
        # region is half-open
        # https://pysam.readthedocs.io/en/latest/glossary.html#term-region
        suffix_contig.query_name, ref_clv, ref_clv + 1
    )
    return num_tail_reads


def gen_clv_record(suffix_contig, r2c_bam):
    ref_clv = calc_ref_clv(suffix_contig)
    num_tail_reads = calc_num_suffix_reads(r2c_bam, suffix_contig, ref_clv)
    tail_length = calc_tail_length(suffix_contig)

    return ClvRecord(
        suffix_contig.reference_name,
        calc_strand(suffix_contig),
        ref_clv,

        'tail_contig',
        suffix_contig.query_name,
        suffix_contig.query_length,
        suffix_contig.mapq,

        # bridge_contig evidence is left empty
        num_tail_reads, tail_length,

        num_bridge_reads=0,
        max_bridge_read_tail_length=0,
        num_link_reads=0,
        num_blank_contigs=0
    )
