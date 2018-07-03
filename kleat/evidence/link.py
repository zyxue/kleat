from collections import defaultdict

from kleat.misc.settings import ClvRecord
from kleat.misc import apautils


def write_evidence(dd_link, contig, csvwriter):
    for clv_key in dd_link['num_reads']:
        clv_record = gen_clv_record(
            contig, clv_key, dd_link['num_reads'][clv_key])
        apautils.write_row(clv_record, csvwriter)


def update_evidence(evid_tuple, evid_holder):
    seqname, strand, ref_clv = evid_tuple
    clv_key = apautils.gen_clv_key_tuple(seqname, strand, ref_clv)
    evid_holder['num_reads'][clv_key] += 1


def is_a_link_read(read):
    # Here we focused on the unmapped all A/T read, but in
    # principle, we could also check from the perspecitve and
    # the mate of a link read, but it would be harder to verify
    # the sequence composition of the link read
    return (not read.mate_is_unmapped
            # assume reference_id comparison is faster than
            # reference_name
            and read.reference_id == read.next_reference_id
            and set(read.query_sequence) in [{'A'}, {'T'}])


def init_evidence_holder():
    return {
        'num_reads': defaultdict(int)
    }


def allN(seq, N):
    return set(seq) == {N}


def calc_next_reference_end(read):
    return read.next_reference_start + read.query_length


def analyze_link(contig, poly_A_or_T_read):
    """poly_read refers to the read with all A or T rather than its mate"""
    ctg_beg = apautils.infer_abs_ref_start(contig)
    ctg_end = apautils.infer_abs_ref_end(contig)

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
            raise ValueError('NOT a polyA/T read: {0}'.format(read))
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
            raise ValueError('NOT a polyA/T read: {0}'.format(read))
    return seqname, strand, ref_clv


def gen_clv_record(
        link_contig, clv_key_tuple, num_link_reads):
    seqname, strand, ref_clv = clv_key_tuple
    return ClvRecord(
        seqname, strand, ref_clv,

        'link',
        link_contig.query_name,
        link_contig.query_length,
        link_contig.mapq,

        0,                      # num_tail_reads
        0,                      # tail_length

        0,                      # num_bridge_reads
        0,                      # max_bridge_tail_len

        num_link_reads,
        num_blank_contigs=0
    )
