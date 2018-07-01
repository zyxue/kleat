import utils as U


def allN(seq, N):
    return set(seq) == {N}


def calc_next_reference_end(read):
    return read.next_reference_start + read.query_length


def analyze_link(contig, poly_A_or_T_read):
    """poly_read refers to the read with all A or T rather than its mate"""
    ctg_beg = U.infer_contig_abs_ref_start(contig)
    ctg_end = U.infer_contig_abs_ref_end(contig)

    read = poly_A_or_T_read
    seq = read.query_sequence

    seqname = contig.reference_name

    if not contig.is_reverse:
        if allN(seq, 'T'):
            strand = '-'
            contig_clv = read.next_reference_start
            ref_clv = ctg_beg + contig_clv     # convert to genome coordinate
        elif allN(seq, 'A'):
            strand = '+'
            # oddly, next_reference_end API doesn't exist, just compute it
            contig_clv = calc_next_reference_end(read)
            ref_clv = ctg_beg + contig_clv
        else:
            raise ValueError(f'NOT a polyA/T read: {read}')
    else:
        if allN(seq, 'T'):
            strand = '+'
            contig_clv = read.next_reference_start
            ref_clv = ctg_end - contig_clv
        elif allN(seq, 'A'):
            strand = '-'
            contig_clv = calc_next_reference_end(read)
            ref_clv = ctg_end - contig_clv
        else:
            raise ValueError(f'NOT a polyA/T read: {read}')
    return seqname, strand, ref_clv
