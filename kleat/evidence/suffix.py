"""
Utilities for analyzing suffix-specific evidence. Contigs here are always
assumed to be suffix contigs
"""

from kleat.misc import apautils
from kleat.misc.settings import ClvRecord
from kleat.hexamer.hexamer import (
    gen_contig_hexamer_tuple,
    gen_reference_hexamer_tuple
)


def calc_ref_clv(suffix_segment, tail_side=None):
    """
    Calculate cleavage site position wst the reference

    :param tail_sideed: pass to avoid redundant checking of tail
    """
    if tail_side is None:
        tail_side = apautils.has_tail(suffix_segment)

    # the coordinates (+1 or not) are verified against visualization on IGV
    if tail_side == 'left':
        return suffix_segment.reference_start
    elif tail_side == 'right':
        return suffix_segment.reference_end - 1
    else:
        raise ValueError('{0} is not a suffix segment'.format(suffix_segment))


def calc_strand(tail_side):
    if tail_side == 'left':
        return '-'
    elif tail_side == 'right':
        return '+'
    else:
        raise ValueError('tail_side must be "left" or "right", '
                         'but {0} passed'.format(tail_side))


def is_a_suffix_read(read, contig, ctg_tail_len=None):
    """
    :param read: read
    :param contig: a suffix contig
    :param tail_len: the tail length of the suffix_contig
    """
    if ctg_tail_len is None:
        ctg_tail_len = apautils.calc_tail_length(contig)
    if apautils.left_tail(read):
        return read.reference_start == ctg_tail_len
    elif apautils.right_tail(read):
        return read.reference_end == contig.infer_query_length(always=True) - 1 - ctg_tail_len
    else:
        return False


def analyze_suffix_reads(r2c_bam, suffix_contig, ctg_clv, ctg_tail_len):
    """
    calculate the number of reads aligned to the cleavage site

    https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
    """
    num_suffix_reads, max_tail_len = 0, 0
    for read in r2c_bam.fetch(suffix_contig.query_name, ctg_clv, ctg_clv + 1):
        if read.is_unmapped:
            continue
        if is_a_suffix_read(read, suffix_contig, ctg_tail_len):
            max_tail_len = max(max_tail_len, apautils.calc_tail_length(read))
            num_suffix_reads += 1
    return num_suffix_reads, max_tail_len


def gen_clv_record(contig, r2c_bam, tail_side, ref_fa):
    """
    :param contig: suffix contig
    :param r2c_bam: pysam instance of read2genome alignment BAM
    :param tail_side (TODO, rename to tail_direction): 'left' or 'right'
    :param ref_fa: pysam instance of reference genome fasta. if provided,
                   will also search PAS hexamer on reference genome.
    """
    strand = calc_strand(tail_side)
    ref_clv = calc_ref_clv(contig, tail_side)
    ctg_tail_len = apautils.calc_tail_length(contig, tail_side)

    if strand == '-':
        ctg_clv = ctg_tail_len
    else:
        ctg_seq_len = contig.infer_query_length(always=True)
        ctg_clv = ctg_seq_len - ctg_tail_len - 1

    num_suffix_reads, max_suffix_read_tail_len = analyze_suffix_reads(
        r2c_bam, contig, ctg_clv, ctg_tail_len)

    ctg_hex, ctg_hex_id, ctg_hex_pos = gen_contig_hexamer_tuple(
        contig, strand, ref_clv, ref_fa, ctg_clv)

    ref_hex, ref_hex_id, ref_hex_pos = gen_reference_hexamer_tuple(
        ref_fa, contig.reference_name, strand, ref_clv)

    return ClvRecord(
        contig.reference_name,
        strand,
        ref_clv,

        'suffix',
        '{0}@{1}'.format(contig.query_name, ctg_clv),
        contig.query_length,
        contig.mapq,
        apautils.is_hardclipped(contig),

        num_suffix_reads,
        max_suffix_read_tail_len,
        ctg_tail_len,
        1,                      # num_suffix_contigs

        # other types of evidence are left empty
        0,                      # num_bridge_reads
        0,                      # max_bridge_read_tail_len
        0,                      # num_bridge_contigs

        0,                      # num_link_reads
        0,                      # num_link_contigs

        0,                      # num_blank_contigs

        ctg_hex, ctg_hex_id, ctg_hex_pos,
        ref_hex, ref_hex_id, ref_hex_pos
    )
