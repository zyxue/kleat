import utils as U
from settings import ClvRecord


def analyze_bridge_read(contig, read):
    ctg_beg = U.infer_contig_abs_ref_start(contig)
    ctg_end = U.infer_contig_abs_ref_end(contig)

    seqname = contig.reference_name
    if not contig.is_reverse:
        if U.left_tail(read, 'T'):
            strand = '-'
            contig_clv = read.reference_start  # in contig coordinate
            ref_clv = ctg_beg + contig_clv     # convert to genome coordinate
            tail_len = read.cigartuples[0][1]
        elif U.right_tail(read, 'A'):
            strand = '+'
            contig_clv = read.reference_end
            ref_clv = ctg_beg + contig_clv
            tail_len = read.cigartuples[-1][1]
        else:
            raise ValueError(f'no tail found for read {read}')
    else:
        if U.left_tail(read, 'T'):
            strand = '+'
            contig_clv = read.reference_start
            ref_clv = ctg_end - contig_clv
            tail_len = read.cigartuples[0][1]
        elif U.right_tail(read, 'A'):
            strand = '-'
            contig_clv = read.reference_end
            ref_clv = ctg_end - contig_clv
            tail_len = read.cigartuples[-1][1]
        else:
            raise ValueError(f'no tail found for read {read}')
    return seqname, strand, ref_clv, tail_len


def gen_clv_record(
        bridge_contig, clv_key_tuple, num_bridge_reads, max_bridge_tail_len):
    seqname, strand, ref_clv = clv_key_tuple
    return ClvRecord(
        seqname, strand, ref_clv,

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
