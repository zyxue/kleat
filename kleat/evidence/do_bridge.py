from kleat.misc import apautils
import kleat.misc.settings as S


def calc_hardclips(cigartuples):
    """
    calculate the number of bases hard clipped at both ends, given cigartuples
    of a contig
    """
    first_idx, last_idx = 0, len(cigartuples) - 1
    left_hc, right_hc = 0, 0    # hc: hardclip
    for k, (key, val) in enumerate(cigartuples):
        if key == S.BAM_CHARD_CLIP:
            if k == first_idx:
                left_hc += val
            elif k == last_idx:
                right_hc += val
    return left_hc, right_hc


def do_fwd_ctg_lt_bdg(read, contig):
    """
    fwd: forwad, ctg: contig, lt: left-tailed, bdg: bridge
    """
    ctg_len_with_hc = contig.infer_query_length(always=True)
    left_hc, right_hc = calc_hardclips(contig.cigartuples)

    pre_ctg_offset = read.reference_start
    if pre_ctg_offset < left_hc:
        # meaning clv is within left hardclip
        return
    elif pre_ctg_offset >= ctg_len_with_hc - right_hc:
        # meaning clv is within right hardclip
        return
    else:
        ctg_offset = pre_ctg_offset - left_hc
        tail_len = read.cigartuples[0][1]
        return '-', ctg_offset, tail_len


def do_fwd_ctg_rt_bdg(read, contig):
    """rt: right-tailed"""
    ctg_len_with_hc = contig.infer_query_length(always=True)
    left_hc, right_hc = calc_hardclips(contig.cigartuples)

    pre_ctg_offset = read.reference_end - 1
    if pre_ctg_offset < left_hc:
        # meaning clv is within left hardclip
        return
    elif pre_ctg_offset >= ctg_len_with_hc - right_hc:
        # meaning clv is within right hardclip
        return
    else:
        ctg_offset = pre_ctg_offset - left_hc
        return '+', ctg_offset, read.cigartuples[-1][1]


def do_rev_ctg_lt_bdg(read, contig):
    ctg_len_with_hc = contig.infer_query_length(always=True)
    left_hc, right_hc = calc_hardclips(list(reversed(contig.cigartuples)))

    pre_ctg_offset = read.reference_start
    if pre_ctg_offset < left_hc:
        # meaning clv is within left hardclip
        return
    elif pre_ctg_offset >= ctg_len_with_hc - right_hc:
        # meaning clv is within right hardclip
        return
    else:
        ctg_offset = ctg_len_with_hc - pre_ctg_offset - 1 - right_hc
        return '+', ctg_offset, read.cigartuples[0][1]


def do_rev_ctg_rt_bdg(read, contig):
    ctg_len_with_hc = contig.infer_query_length(always=True)
    left_hc, right_hc = calc_hardclips(list(reversed(contig.cigartuples)))

    pre_ctg_offset = read.reference_end - 1
    if pre_ctg_offset < left_hc:
        # meaning clv is within left hardclip
        return
    elif pre_ctg_offset >= ctg_len_with_hc - right_hc:
        # meaning clv is within right hardclip
        return
    else:
        ctg_offset = ctg_len_with_hc - pre_ctg_offset - 1 - right_hc
        return '-', ctg_offset, read.cigartuples[-1][1]


def do_forward_contig(contig, read):
    if apautils.left_tail(read, 'T'):
        res = do_fwd_ctg_lt_bdg(read, contig)
        if res is not None:
            return res + ('left', )
    elif apautils.right_tail(read, 'A'):
        res = do_fwd_ctg_rt_bdg(read, contig)
        if res is not None:
            return res + ('right', )
    else:
        raise ValueError('no tail found for read {0}'.format(read))


def do_reverse_contig(contig, read):
    # set always=True to include hard-clipped bases
    # https://pysam.readthedocs.io/en/latest/api.html?highlight=parse_region#pysam.AlignedSegment.infer_query_length
    # TODO: avoid multiple calling of left_tail/right_tail
    if apautils.left_tail(read, 'T'):
        # it's right because it should be reversed again to match the forward
        # direction
        res = do_rev_ctg_lt_bdg(read, contig)
        if res is not None:
            return res + ('right',)
    elif apautils.right_tail(read, 'A'):
        res = do_rev_ctg_rt_bdg(read, contig)
        if res is not None:
            return res + ('left',)
    else:
        raise ValueError('no tail found for read {0}'.format(read))


def do_bridge(contig, read):
    if contig.is_reverse:
        return do_reverse_contig(contig, read)
    else:
        return do_forward_contig(contig, read)
