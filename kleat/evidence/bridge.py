from collections import defaultdict

from kleat.misc import apautils
import kleat.misc.settings as S
from kleat.hexamer.search import (
    search,  # do PAS hexamer search in a more customized way instead of using
             # gen_contig_hexamer_tuple
    gen_reference_hexamer_tuple
)
from kleat.hexamer.extract_seq import extract_seq


def write_evidence(dd_bridge, contig, ref_fa, csvwriter):
    res = []
    for clv_key in dd_bridge['num_reads']:
        clv_record = gen_clv_record(
            contig, clv_key,
            dd_bridge['num_reads'][clv_key],
            dd_bridge['max_tail_len'][clv_key],
            dd_bridge['hexamer_tuple'][clv_key],
            ref_fa
        )
        apautils.write_row(clv_record, csvwriter)
        res.append(clv_record)
    return res


def update_evidence(evid_tuple, evid_holder):
    """
    update the bridge evidence holder with extracted bridge evidence

    :param evid_tuple: evidence tuple
    :param evid_holder: A dict holder bridge evidence for a given contig
    """
    seqname, strand, ref_clv, tail_len, hex_tuple = evid_tuple
    clv_key = apautils.gen_clv_key_tuple(seqname, strand, ref_clv)
    evid_holder['num_reads'][clv_key] += 1
    evid_holder['max_tail_len'][clv_key] = max(
        evid_holder['max_tail_len'][clv_key], tail_len)
    if hex_tuple is not None:
        evid_holder['hexamer_tuple'][clv_key] = hex_tuple


def is_a_bridge_read(read):
    return not read.is_unmapped and apautils.has_tail(read)


def init_evidence_holder():
    """
    initialize holders for bridge and link evidence of a given contig
    """
    return {
        'num_reads': defaultdict(int),
        'max_tail_len': defaultdict(int),
        'hexamer_tuple': defaultdict(lambda: None),  # TODO: maybe defaultdict(tuple)
    }


def do_fwd_ctg_lt_bdg(read):
    """fwd: forwad, ctg: contig, lt: left-tailed, bdg: bridge"""
    return '-', read.reference_start, read.cigartuples[0][1]


def do_fwd_ctg_rt_bdg(read):
    """rt: right-tailed"""
    return '+', read.reference_end - 1, read.cigartuples[-1][1]


def do_rev_ctg_lt_bdg(read, ctg_len):
    # rev (reverse), opposite of fwd (forward)
    return '+', ctg_len - read.reference_start - 1, read.cigartuples[0][1]


def do_rev_ctg_rt_bdg(read, ctg_len):
    return '-', ctg_len - read.reference_end, read.cigartuples[-1][1]


def do_forward_contig(contig, read):
    if apautils.left_tail(read, 'T'):
        return do_fwd_ctg_lt_bdg(read) + ('left', )
    elif apautils.right_tail(read, 'A'):
        return do_fwd_ctg_rt_bdg(read) + ('right', )
    else:
        raise ValueError('no tail found for read {0}'.format(read))


def do_reverse_contig(contig, read):
    # set always=True to include hard-clipped bases
    # https://pysam.readthedocs.io/en/latest/api.html?highlight=parse_region#pysam.AlignedSegment.infer_query_length
    ctg_len = contig.infer_query_length(always=True)
    # TODO: avoid multiple calling of left_tail/right_tail
    if apautils.left_tail(read, 'T'):
        # it's right because it should be reversed again to match the forward
        # direction
        return do_rev_ctg_lt_bdg(read, ctg_len) + ('right',)
    elif apautils.right_tail(read, 'A'):
        return do_rev_ctg_rt_bdg(read, ctg_len) + ('left',)
    else:
        raise ValueError('no tail found for read {0}'.format(read))


def do_bridge(contig, read):
    if contig.is_reverse:
        return do_reverse_contig(contig, read)
    else:
        return do_forward_contig(contig, read)


def gen_hex_tuple(contig, strand, ref_clv, ref_fa, ctg_clv, dd_bridge):
    # TODO: the returning of None is pretty ugly, to refactor
    seqname = contig.reference_name
    clv_key = apautils.gen_clv_key_tuple(seqname, strand, ref_clv)
    if dd_bridge['hexamer_tuple'][clv_key] is None:  # do search
        hex_src_seq = extract_seq(
            contig, strand, ref_clv, ref_fa, ctg_clv)

        ctg_hex_tuple = search(strand, ref_clv, hex_src_seq)
        if ctg_hex_tuple is None:
            ctg_hex_tuple = ('NA', -1, -1)
    else:
        ctg_hex_tuple = None
    return ctg_hex_tuple


def analyze_bridge(contig, read, ref_fa, dd_bridge):
    """
    :param dd_bridge: holds bridge_evidence for a given contig, here it's just
    used to check if hexamer_search has already been done for a given ref_clv
    """
    seqname = contig.reference_name

    strand, ctg_offset, tail_len, tail_direction = do_bridge(contig, read)

    offset = apautils.calc_genome_offset(
        contig.cigartuples, ctg_offset, tail_direction)

    ref_clv = contig.reference_start + offset

    # ctg_clv is ctg_offset, i.e. the pos of clv in contig coordinate
    ctg_hex_tuple = gen_hex_tuple(
        contig, strand, ref_clv, ref_fa, ctg_offset, dd_bridge)

    return seqname, strand, ref_clv, tail_len, ctg_hex_tuple


def gen_clv_record(contig, clv_key_tuple,
                   num_bridge_reads, max_bridge_tail_len,
                   ctg_hex_tuple, ref_fa):
    """
    :param contig: bridge contig
    :clv_key_tuple: a tuple of (seqname, strand, cleavage_site_position)
    :param ref_fa: if provided, search PAS hexamer on reference genome, too
    """
    seqname, strand, ref_clv = clv_key_tuple

    ctg_hex, ctg_hex_id, ctg_hex_pos = ctg_hex_tuple
    ref_hex, ref_hex_id, ref_hex_pos = gen_reference_hexamer_tuple(
        ref_fa, contig.reference_name, strand, ref_clv)

    return S.ClvRecord(
        seqname,
        strand,
        ref_clv,

        'bridge',
        contig.query_name,
        contig.infer_query_length(True),
        contig.mapq,
        apautils.is_hardclipped(contig),

        0,                      # num_tail_reads
        0,                      # tail_length
        0,                      # num_suffix_contigs

        num_bridge_reads,
        max_bridge_tail_len,
        1,                      # num_bridge_contigs

        0,                      # num_link_reads,
        0,                      # num_link_contigs

        0,                      # num_blank_contigs

        ctg_hex, ctg_hex_id, ctg_hex_pos,
        ref_hex, ref_hex_id, ref_hex_pos
    )
