"""
Utilities for analyzing suffix-specific evidence. Contigs here are always
assumed to be suffix contigs
"""

from kleat.misc import apautils
import kleat.misc.search_hexamer as srch_hex
import kleat.misc.settings as S


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


def extract_seq(contig):
    """remove clipped ends before searching for hexamer, the clipped ends would
    affect calculation of genomics coordinates of the found PAS hexamer"""
    beg_offset = 0
    end_offset = -1
    for key, val in contig.cigartuples:
        if key == S.BAM_CSOFT_CLIP:
            beg_offset = val
        if key == S.BAM_CHARD_CLIP:
            end_offset = val
    return contig.query_sequence[beg_offset: end_offset]


def gen_clv_record(contig, r2c_bam, tail_side, ref_fa=None):
    """
    :param contig: suffix contig
    :param r2c_bam: pysam read2genome alignment BAM
    :param tail_side (TODO, rename to tail_direction): 'left' or 'right'
    :param ref_fa: if provided, search PAS hexamer on reference genome, too
    """
    strand = apautils.calc_strand(tail_side)
    ref_clv = apautils.calc_ref_clv(contig, tail_side)
    tail_len = apautils.calc_tail_length(contig, tail_side)

    num_suffix_reads = calc_num_suffix_reads(r2c_bam, contig, ref_clv)

    ctg_hex, ctg_hex_id, ctg_hex_pos = 'NA', -1, -1
    ref_hex, ref_hex_id, ref_hex_pos = 'NA', -1, -1

    # no need to reverse_complement the seq as the hexamer search function is
    # designed to search reference genome sequence, rev_comp is taken care of
    # within the function
    ctg_seq = extract_seq(contig)
    ctg_hex_tuple = srch_hex.search(strand, ref_clv, ctg_seq)

    if ctg_hex_tuple is not None:
        ctg_hex, ctg_hex_id, ctg_hex_pos = ctg_hex_tuple

    if ref_fa is not None:
        ref_hex_tuple = srch_hex.search_reference_genome(
            ref_fa, contig.reference_name, ref_clv, strand)
        if ref_hex_tuple is not None:
            ref_hex, ref_hex_id, ref_hex_pos = ref_hex_tuple

    return S.ClvRecord(
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
