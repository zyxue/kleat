import utils as U
from settings import ClvRecord, BAM_CMATCH, BAM_CREF_SKIP, BAM_CDEL


def calc_offset(contig, match_len_cutoff):
    """
    Calculate the offset caused by skipped region (e.g. intron) until stop_pos
    relative to contig
    """
    match_len = 0
    offset = 0
    for key, val in contig.cigartuples:
        if key in [BAM_CMATCH, BAM_CDEL]:
            match_len += val
            if match_len >= match_len_cutoff:
                delta = val - (match_len - match_len_cutoff)
                offset += delta
                break
            offset += val
        if key == BAM_CREF_SKIP:
            offset += val
    return offset


def analyze_bridge_read(contig, read):
    # beginning and end wst to genome
    gnm_beg = U.infer_contig_abs_ref_start(contig)

    seqname = contig.reference_name
    if not contig.is_reverse:
        if U.left_tail(read, 'T'):
            strand = '-'
            match_len_cutoff = read.reference_start
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset + 1
            tail_len = read.cigartuples[0][1]
        elif U.right_tail(read, 'A'):
            strand = '+'
            match_len_cutoff = read.reference_end - 1
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset
            tail_len = read.cigartuples[-1][1]
        else:
            raise ValueError(f'no tail found for read {read}')
    else:
        if U.left_tail(read, 'T'):
            strand = '+'
            match_len_cutoff = contig.query_length - read.reference_start
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset + 1
            tail_len = read.cigartuples[0][1]
        elif U.right_tail(read, 'A'):
            strand = '-'
            match_len_cutoff = contig.query_length - read.reference_end
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset + 1
            tail_len = read.cigartuples[-1][1]
        else:
            raise ValueError(f'no tail found for read {read}')
    return seqname, strand, gnm_clv, tail_len


def gen_clv_record(
        bridge_contig, clv_key_tuple, num_bridge_reads, max_bridge_tail_len):
    seqname, strand, gnm_clv = clv_key_tuple
    return ClvRecord(
        seqname, strand, gnm_clv,

        'bridge',
        bridge_contig.query_name,
        bridge_contig.query_length,
        bridge_contig.mapq,

        0,                      # num_tail_reads
        0,                      # tail_length

        num_bridge_reads,
        max_bridge_tail_len,

        num_link_reads=0,
        num_blank_contigs=0
    )
