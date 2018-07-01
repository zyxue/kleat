def analyze_link_read_candidate(contig, poly_read):
    """poly_read refers to the read with all A or T rather than its mate"""
    pass
    # if read.query_sequence.startswith('T' * 50) and contig.query_length > 150:
    #     import pdb; pdb.set_trace()


    #         # assume the start/end of the paired read is the
    #         # cleavage site
    #         ref_clv, strand = calc_ref_clv_from_r2c_alignment(
    #             contig, read.next_reference_start)
