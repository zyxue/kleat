def analyze_bridge_candidate(contig, read):
    pass
    # if is_left_tail_segment(read):
    #     ref_clv = calc_ref_clv_from_r2c_alignment(contig, read.reference_end)
    #     strand = '-'
    #     tail_len = calc_tail_length(read)
    #     num_bdg_reads_dd[ref_clv] = num_bdg_reads_dd.get(ref_clv, 0) + 1
    #     max_bdg_tail_len_dd[ref_clv] = max(max_bdg_tail_len_dd.get(ref_clv, 0), tail_len)
    # if is_right_tail_segment(read):
    #     ref_clv = calc_ref_clv_from_r2c_alignment(contig, read.reference_start)
    #     strand = '+'
    #     tail_len = calc_tail_length(read)
    #     num_bdg_reads_dd[ref_clv] = num_bdg_reads_dd.get(ref_clv, 0) + 1
    #     max_bdg_tail_len_dd[ref_clv] = max(max_bdg_tail_len_dd.get(ref_clv, 0), tail_len)

        # if read.is_reverse:
        #     # read2contig alignment cannot be reverse to be a read
        #     # supporting polyA tail
        #     continue
