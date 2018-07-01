from settings import ClvRecord


def gen_two_clv_records(blank_contig):
    """
    Assume there is still a clv at the 3' end of the contig even without any
    polya evidence, in thus case, there is no direction, so either end of the
    contig could be a clv. Hence add two to the function name explicitly
    """

    strands = ['+', '-']
    ref_clv_candidates = [
        blank_contig.reference_end + 1,
        blank_contig.reference_start
    ]

    for strand, ref_clv in zip(strands, ref_clv_candidates):
        yield ClvRecord(
            blank_contig.reference_name, strand, ref_clv,

            'blank',
            blank_contig.query_name,
            blank_contig.query_length,
            blank_contig.mapq,

            0,                   # num_tail_reads
            0,                   # tail_length
            0,                   # num_bridge_reads
            0,                   # max_bridge_tail_len

            0,                   # num_link_reads
            num_blank_contigs=1
        )
