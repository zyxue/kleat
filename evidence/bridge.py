import utils as U


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
