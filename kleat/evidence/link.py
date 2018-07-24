from collections import defaultdict

from kleat.misc.settings import ClvRecord
from kleat.misc import apautils
from kleat.hexamer.hexamer import (
    gen_contig_hexamer_tuple,
    gen_reference_hexamer_tuple
)


def write_evidence(dd_link, contig, ref_fa, csvwriter):
    res = []
    for clv_key in dd_link['num_reads']:
        clv_record = gen_clv_record(
            contig, clv_key, dd_link['num_reads'][clv_key], ref_fa)
        apautils.write_row(clv_record, csvwriter)
        res.append(clv_record)
    return res


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


def analyze_forward_link(contig, read):
    """
    :param read: a pysam.libcalignedsegment.AlignedSegment instance for a
                 polyA/T read
    :rtype: A tuple of (strand, reference_clv)
    """
    seq = read.query_sequence
    if allN(seq, 'T'):
        return '-', contig.reference_start
    elif allN(seq, 'A'):
        return '+', contig.reference_end - 1
    else:
        raise ValueError('NOT a polyA/T read: {0}'.format(read))


def analyze_reverse_link(contig, read):
    """
    :param read: a pysam.libcalignedsegment.AlignedSegment instance for a
                 polyA/T read
    :rtype: A tuple of (strand, reference_clv)
    """
    seq = read.query_sequence
    if allN(seq, 'T'):
        return '+', contig.reference_end - 1
    elif allN(seq, 'A'):
        return '-', contig.reference_start
    else:
        raise ValueError('NOT a polyA/T read: {0}'.format(read))


def analyze_link(contig, polyA_or_T_read):
    """poly_read refers to the read with all A or T rather than its mate"""
    seqname = contig.reference_name

    if contig.is_reverse:
        strand, ref_clv = analyze_reverse_link(contig, polyA_or_T_read)
    else:
        strand, ref_clv = analyze_forward_link(contig, polyA_or_T_read)

    return seqname, strand, ref_clv


def calc_ctg_clv(strand, ctg_seq_len):
    if strand == '-':
        ctg_clv = 0
    else:
        ctg_clv = ctg_seq_len - 1
    return ctg_clv


def gen_clv_record(contig, clv_key_tuple, num_link_reads, ref_fa):
    """
    :param contig: link contig
    :clv_key_tuple: a tuple of (seqname, strand, cleavage_site_position)
    :param ref_fa: if provided, search PAS hexamer on reference genome, too
    """
    seqname, strand, ref_clv = clv_key_tuple
    ctg_seq_len = contig.infer_query_length(always=True)
    ctg_clv = calc_ctg_clv(strand, ctg_seq_len)

    ctg_hex, ctg_hex_id, ctg_hex_pos = gen_contig_hexamer_tuple(
        contig, strand, ref_clv, ref_fa, ctg_clv)

    ref_hex, ref_hex_id, ref_hex_pos = gen_reference_hexamer_tuple(
        ref_fa, contig.reference_name, strand, ref_clv)

    return ClvRecord(
        seqname,
        strand,
        ref_clv,

        ctg_hex,
        ctg_hex_id,
        ctg_hex_pos,

        ref_hex,
        ref_hex_id,
        ref_hex_pos,

        evidence_type='link',
        contig_id_at_pos='{0}@{1}'.format(contig.query_name, ctg_clv),
        contig_len=contig.query_length,
        contig_mapq=contig.mapq,
        contig_is_hardclipped=apautils.is_hardclipped(contig),

        num_suffix_reads=0,
        max_suffix_read_tail_len=0,
        suffix_contig_tail_len=0,
        num_suffix_contigs=0,

        num_bridge_reads=0,
        max_bridge_read_tail_len=0,
        num_bridge_contigs=0,

        num_link_reads=num_link_reads,
        num_link_contigs=1,

        num_blank_contigs=0,
    )
