"""
Utilities for analyzing suffix evidence
"""


# https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0
BAM_CSOFT_CLIP = 4              # used to identify soft clip 


__all__ = ['is_suffix_contig', 'analyze_suffix_config']


def calc_strand(contig):
    """
    calculate the strand of clv (hence the corresponding gene) this contig may
    support
    """
    if is_right_tail_segment(contig):
        return '+'
    elif is_left_tail_segment(contig):
        return '-'
    else:
        raise ValueError(f'cannot infer strand for contig: {contig.query_name}')


def gen_suffix_contig_info(contig):
    return [
        contig.reference_name,
        calc_strand(contig),
        contig.query_name,
        contig.query_length,
        contig.mapq
    ]


def analyze_suffix_contig(contig, r2c_bam):
    contig_info = gen_suffix_contig_info(contig)
    ref_clv = calc_ref_clv(contig)
    num_tail_reads, _ = calc_num_tail_reads(contig, r2c_bam)
    tail_length = calc_tail_length(contig)

    clv_record = (
        *contig_info, ref_clv, 'tail_contig',
        # bridge_contig evidence is left empty
        num_tail_reads, tail_length, 0, 0, 0, 0
    )
    return clv_record


def is_suffix_contig(segment):
    return is_right_tail_segment(segment) or is_left_tail_segment(segment)


def is_right_tail_segment(segment, tail_base='A'):
    """
    tseg: tail segment

    default tail_base applies main to alignment to genome, where the polyA tail
    strand is known
    """
    seq = segment.query_sequence
    last_cigar = segment.cigartuples[-1]
    # potential test case == "A0.R100820":
    return (
        seq.endswith(tail_base)
        # right clipped
        and last_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all As
        and set(seq[-last_cigar[1]:]) == {tail_base}
    )


def is_left_tail_segment(segment, tail_base='T'):
    seq = segment.query_sequence
    first_cigar = segment.cigartuples[0]
    # potential test case A0.S36384
    return (
        seq.startswith(tail_base)
        # left clipped
        and first_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all Ts
        and set(seq[:first_cigar[1]]) == {tail_base}
    )


def calc_ref_clv(segment):
    """calculate cleavage site position wst the reference"""
    if is_right_tail_segment(segment):      # + strand
        return segment.reference_end + 1
    elif is_left_tail_segment(segment):  # - strand
        return segment.reference_start
    else:
        raise


def calc_tail_length(segment):
    """
    calculate A/T length of a tail contig, this information is extracted from
    softclip in the CIGAR
    """
    if is_right_tail_segment(segment):
        return segment.cigartuples[-1][1]
    elif is_left_tail_segment(segment):
        return segment.cigartuples[0][1]
    else:
        return '?'


def calc_num_tail_reads(segment, r2c_sam):
    """
    calculate the number of reads aligned to the cleavage site

    https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
    """
    if is_right_tail_segment(segment):
        # clv position wst. the segment
        seg_clv = segment.query_alignment_end
    elif is_left_tail_segment(segment):
        seg_clv = segment.query_alignment_start - 1
    else:
        raise

    num_tail_reads = r2c_sam.count(
        # region is half-open
        # https://pysam.readthedocs.io/en/latest/glossary.html#term-region
        segment.query_name, seg_clv, seg_clv + 1
    )
    return num_tail_reads, seg_clv
