from collections import defaultdict

from kleat.evidence.do_bridge import do_bridge
from kleat.misc import apautils
from kleat.misc.calc_genome_offset import calc_genome_offset
import kleat.misc.settings as S

# for bridge PAS hexamer search needs done in a more customized way than just
# using gen_contig_hexamer_tuple
from kleat.hexamer.search import search
from kleat.hexamer.hexamer import extract_seq, gen_reference_hexamer_tuple


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
    seqname, strand, ref_clv, ctg_clv, tail_len, hex_tuple = evid_tuple
    clv_key = apautils.gen_clv_key_tuple_with_ctg_clv(seqname, strand, ref_clv, ctg_clv)
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


def gen_hex_tuple(contig, strand, ref_clv, ref_fa, ctg_clv, dd_bridge):
    # TODO: the returning of None is pretty ugly, to refactor
    seqname = contig.reference_name
    clv_key = apautils.gen_clv_key_tuple_with_ctg_clv(seqname, strand, ref_clv, ctg_clv)
    if dd_bridge['hexamer_tuple'][clv_key] is None:  # do search
        hex_src_seq = extract_seq(contig, strand, ref_clv, ref_fa, ctg_clv)

        ctg_hex_tuple = search(strand, ref_clv, hex_src_seq)
        if ctg_hex_tuple is None:
            ctg_hex_tuple = ('NA', -1, -1)
    else:
        ctg_hex_tuple = None
    return ctg_hex_tuple


def analyze_bridge(contig, read, ref_fa, dd_bridge, bridge_skip_check_size):
    """
    :param dd_bridge: holds bridge_evidence for a given contig, here it's just
    used to check if hexamer_search has already been done for a given ref_clv
    """
    seqname = contig.reference_name

    bdg_support = do_bridge(contig, read)
    if bdg_support is None:     # likely a chimeric contig
        return

    strand, ctg_clv, tail_len, tail_direction = bdg_support

    offset = calc_genome_offset(contig.cigartuples, ctg_clv, tail_direction, skip_check_size=3)

    if offset < 0:             # meaning the clv is on soft/hard clipped region
        return

    ref_clv = contig.reference_start + offset

    ctg_hex_tuple = gen_hex_tuple(
        contig, strand, ref_clv, ref_fa, ctg_clv, dd_bridge)

    return seqname, strand, ref_clv, ctg_clv, tail_len, ctg_hex_tuple


def gen_clv_record(contig, clv_key_tuple,
                   num_bridge_reads, max_bridge_read_tail_len,
                   ctg_hex_tuple, ref_fa):
    """
    :param contig: bridge contig
    :clv_key_tuple: a tuple of (seqname, strand, cleavage_site_position)
    :param ref_fa: if provided, search PAS hexamer on reference genome, too
    """
    seqname, strand, ref_clv, ctg_clv = clv_key_tuple

    ctg_hex, ctg_hex_id, ctg_hex_pos = ctg_hex_tuple
    ref_hex, ref_hex_id, ref_hex_pos = gen_reference_hexamer_tuple(
        ref_fa, contig.reference_name, strand, ref_clv)

    return S.ClvRecord(
        seqname,
        strand,
        ref_clv,

        ctg_hex,
        ctg_hex_id,
        ctg_hex_pos,

        ref_hex,
        ref_hex_id,
        ref_hex_pos,

        evidence_type='bridge',
        contig_id_at_pos='{0}@{1}'.format(contig.query_name, ctg_clv),
        contig_len=contig.infer_query_length(True),
        contig_mapq=contig.mapq,
        contig_is_hardclipped=apautils.is_hardclipped(contig),

        num_suffix_reads=0,
        max_suffix_read_tail_len=0,
        suffix_contig_tail_len=0,
        num_suffix_contigs=0,

        num_bridge_reads=num_bridge_reads,
        max_bridge_read_tail_len=max_bridge_read_tail_len,
        num_bridge_contigs=1,

        num_link_reads=0,
        num_link_contigs=0,

        num_blank_contigs=0,
    )
