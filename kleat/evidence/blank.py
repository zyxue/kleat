from kleat.misc import apautils
from kleat.misc.settings import ClvRecord

from kleat.hexamer.hexamer import (
    gen_contig_hexamer_tuple,
    gen_reference_hexamer_tuple
)


def gen_two_clv_records(contig, ref_fa, already_supported_clv_keys):
    """
    Assume there is still a clv at the 3' end of the contig even without any
    polya evidence, in thus case, there is no direction, so either end of the
    contig could be a clv. Hence add two to the function name explicitly
    """
    seqname = contig.reference_name
    strands = ['+', '-']
    ref_clv_candidates = [
        contig.reference_end - 1,
        contig.reference_start
    ]
    is_hardclipped = apautils.is_hardclipped(contig)

    for strand, ref_clv in zip(strands, ref_clv_candidates):
        clv_key = apautils.gen_clv_key_tuple(seqname, strand, ref_clv)
        if clv_key in already_supported_clv_keys:
            continue

        if strand == '-':
            ctg_clv = 0
        else:
            ctg_seq_len = contig.infer_query_length(always=True)
            ctg_clv = ctg_seq_len - 1

        ctg_hex, ctg_hex_id, ctg_hex_pos = gen_contig_hexamer_tuple(
            contig, strand, ref_clv, ref_fa, ctg_clv)

        ref_hex, ref_hex_id, ref_hex_pos = gen_reference_hexamer_tuple(
            ref_fa, contig.reference_name, strand, ref_clv)

        yield ClvRecord(
            contig.reference_name,
            strand,
            ref_clv,

            'blank',
            '{0}@{1}'.format(contig.query_name, ctg_clv),
            contig.query_length,
            contig.mapq,
            is_hardclipped,

            0,                   # num_tail_reads
            0,                   # tail_length
            0,                   # num_suffix_contigs

            0,                   # num_bridge_reads
            0,                   # max_bridge_tail_len
            0,                   # num_bridge_contigs

            0,                   # num_link_reads
            0,                   # num_link_contigs

            1,                   # num_blank_contigs

            ctg_hex, ctg_hex_id, ctg_hex_pos,
            ref_hex, ref_hex_id, ref_hex_pos

        )
