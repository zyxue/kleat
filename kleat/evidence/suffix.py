"""
Utilities for analyzing suffix-specific evidence. Contigs here are always
assumed to be suffix contigs
"""

from kleat.misc import apautils
from kleat.misc.settings import ClvRecord
from kleat.misc.search_hexamer import (
    gen_contig_hexamer_tuple, gen_reference_hexamer_tuple
)




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


def gen_clv_record(contig, r2c_bam, tail_side, ref_fa=None):
    """
    :param contig: suffix contig
    :param r2c_bam: pysam instance of read2genome alignment BAM
    :param tail_side (TODO, rename to tail_direction): 'left' or 'right'
    :param ref_fa: pysam instance of reference genome fasta. if provided,
                   will also search PAS hexamer on reference genome.
    """
    strand = apautils.calc_strand(tail_side)
    ref_clv = apautils.calc_ref_clv(contig, tail_side)
    tail_len = apautils.calc_tail_length(contig, tail_side)

    num_suffix_reads = calc_num_suffix_reads(r2c_bam, contig, ref_clv)

    ctg_hex, ctg_hex_id, ctg_hex_pos = gen_contig_hexamer_tuple(
        contig, strand, ref_clv)

    ref_hex, ref_hex_id, ref_hex_pos = gen_reference_hexamer_tuple(
        ref_fa, contig.reference_name, strand, ref_clv)

    return ClvRecord(
        contig.reference_name,
        strand,
        ref_clv,

        'suffix',
        contig.query_name,
        contig.query_length,
        contig.mapq,

        num_suffix_reads,
        tail_len,

        # other types of evidence are left empty
        0,    # num_bridge_reads
        0,    # max_bridge_read_tail_len

        0,    # num_link_reads

        0,    # num_blank_contigs

        ctg_hex, ctg_hex_id, ctg_hex_pos,
        ref_hex, ref_hex_id, ref_hex_pos
    )
