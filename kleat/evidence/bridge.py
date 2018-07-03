from collections import defaultdict

from kleat.misc import apautils
from kleat.misc.settings import ClvRecord, BAM_CMATCH, BAM_CREF_SKIP, BAM_CDEL


def write_evidence(dd_bridge, contig, csvwriter):
    for clv_key in dd_bridge['num_reads']:
        clv_record = gen_clv_record(
            contig, clv_key,
            dd_bridge['num_reads'][clv_key],
            dd_bridge['max_tail_len'][clv_key]
        )
        apautils.write_row(clv_record, csvwriter)


def update_evidence(evid_tuple, evid_holder):
    """
    update the bridge evidence holder with extracted bridge evidence

    :param evid_tuple: evidence tuple
    :param evid_holder: A dict holder bridge evidence for a given contig
    """
    seqname, strand, ref_clv, tail_len = evid_tuple
    clv_key = apautils.gen_clv_key_tuple(seqname, strand, ref_clv)
    evid_holder['num_reads'][clv_key] += 1
    evid_holder['max_tail_len'][clv_key] = max(
        evid_holder['max_tail_len'][clv_key], tail_len)


def is_a_bridge_read(read):
    return not read.is_unmapped and apautils.has_tail(read)


def init_evidence_holder():
    """
    initialize holders for bridge and link evidence of a given contig
    """
    return {
        'num_reads': defaultdict(int),
        'max_tail_len': defaultdict(int),
    }


def calc_offset(contig, match_len_cutoff):
    """
    Calculate the offset caused by skipped region (e.g. intron) until stop_pos
    relative to contig
    """
    match_len = 0
    offset = 0
    for key, val in contig.cigartuples:
        if key in [BAM_CMATCH, BAM_CDEL]:
            match_len += val
            if match_len >= match_len_cutoff:
                delta = val - (match_len - match_len_cutoff)
                offset += delta
                break
            offset += val
        if key == BAM_CREF_SKIP:
            offset += val
    return offset


def analyze_bridge_read(contig, read):
    # beginning and end wst to genome
    gnm_beg = apautils.infer_contig_abs_ref_start(contig)

    seqname = contig.reference_name
    if not contig.is_reverse:
        if apautils.left_tail(read, 'T'):
            strand = '-'
            match_len_cutoff = read.reference_start
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset + 1
            tail_len = read.cigartuples[0][1]
        elif apautils.right_tail(read, 'A'):
            strand = '+'
            match_len_cutoff = read.reference_end - 1
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset
            tail_len = read.cigartuples[-1][1]
        else:
            raise ValueError(f'no tail found for read {read}')
    else:
        # set always=True to include hard-clipped bases
        # https://pysam.readthedocs.io/en/latest/api.html?highlight=parse_region#pysam.AlignedSegment.infer_query_length
        ctg_len = contig.infer_query_length(always=True)
        if apautils.left_tail(read, 'T'):
            strand = '+'
            match_len_cutoff = ctg_len - read.reference_start
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset + 1
            tail_len = read.cigartuples[0][1]
        elif apautils.right_tail(read, 'A'):
            strand = '-'
            match_len_cutoff = ctg_len - read.reference_end
            offset = calc_offset(contig, match_len_cutoff)
            gnm_clv = gnm_beg + offset + 1
            tail_len = read.cigartuples[-1][1]
        else:
            raise ValueError(f'no tail found for read {read}')
    return seqname, strand, gnm_clv, tail_len


def gen_clv_record(
        bridge_contig, clv_key_tuple, num_bridge_reads, max_bridge_tail_len):
    seqname, strand, gnm_clv = clv_key_tuple
    return ClvRecord(
        seqname, strand, gnm_clv,

        'bridge',
        bridge_contig.query_name,
        bridge_contig.infer_query_length(True),
        bridge_contig.mapq,

        0,                      # num_tail_reads
        0,                      # tail_length

        num_bridge_reads,
        max_bridge_tail_len,

        num_link_reads=0,
        num_blank_contigs=0
    )
