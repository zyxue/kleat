import kleat.misc.settings as S
from kleat.misc.apautils import (
    calc_genome_offset,
    fetch_seq                   # deal with circular DNA, too, e.g. chrM
)

"""xseq: extract_seq"""


def do_match(ctg_e, ref_e, cigar_val, ctg_seq, ctg_clv):
    """
    :param ctg_e: contig end index
    :param ref_e: reference end index
    :param ctg_clv: cleavage site in contig coordinate
    """
    ctg_b = ctg_e - cigar_val
    ref_b = ref_e - cigar_val

    if ctg_e >= ctg_clv:
        if ctg_b < ctg_clv:
            seq_to_add = ctg_seq[ctg_b: ctg_clv + 1]
        else:
            seq_to_add = ''
    else:
        seq_to_add = ctg_seq[ctg_b: ctg_e]

    next_ctg_e = ctg_b
    next_ref_e = ref_b
    return next_ctg_e, next_ref_e, seq_to_add


def do_skip(ref_e, cigar_val, ref_fa, seqname, ref_clv):
    """
    :param ref_e: reference end index
    :param ref_fa: a pysam.libcalignmentfile.AlignmentFile instance for reference genome
    :param ref_clv: cleavage site in reference genome coordinate
    """
    ref_b = ref_e - cigar_val

    if ref_e >= ref_clv:
        if ref_b < ref_clv:
            seq_to_add = fetch_seq(ref_fa, seqname, ref_b, ref_clv + 1)
        else:
            seq_to_add = ''
    else:
        seq_to_add = fetch_seq(ref_fa, seqname, ref_b, ref_e)

    next_ref_e = ref_b
    return next_ref_e, seq_to_add


def init_ctg_end(ctg_seq):
    return len(ctg_seq)


def init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq):
    """
    Initialize the end index in genome coordinate by calculating the offset
    from right, using `calc_genome_offset`

    comparing to the minus corresponding function, there is one additional
    argument neeed, i.e. `ctg_seq`

    :param ctg_seq: should include soft/hardclipped region if it's clipped
    """
    cigartuples = list(reversed(cigartuples))
    ctg_clv = len(ctg_seq) - ctg_clv  # from_the_right

    # TODO: left may not matter in such case
    offset = calc_genome_offset(cigartuples, ctg_clv, 'left')

    cgr = cigartuples[0]
    # needs to add back clipped regions here
    if cgr[0] == S.BAM_CSOFT_CLIP or cgr[0] == S.BAM_CHARD_CLIP:
        offset += cgr[1]
    return ref_clv + offset


def extract(cigartuples, ctg_seq, seqname, strand, ctg_clv, ref_clv, ref_fa, window):
    """
    scan from Right to Left,

    compared to minus strand, initialize contig end (ce) and reference end (fe)
    instead of cb and fb

    use fe instead of re because re is a Python module for regex
    """
    ce = init_ctg_end(ctg_seq)
    fe = init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq)
    target_fe = max(0, ref_clv - window)

    res_seq = ''
    for idx, (key, val) in enumerate(reversed(cigartuples)):
        if key == S.BAM_CSOFT_CLIP or key == S.BAM_CHARD_CLIP:
            if ce <= ctg_clv:
                seq = ctg_seq[ce - val:ce]
                res_seq = seq + res_seq
            ce -= val
            fe -= val

        elif key == S.BAM_CMATCH:
            ce, fe, seq = do_match(ce, fe, val, ctg_seq, ctg_clv)
            res_seq = seq + res_seq

        elif key == S.BAM_CREF_SKIP:
            fe, seq = do_skip(fe, val, ref_fa, seqname, ref_clv)
            res_seq = seq + res_seq

        elif key == S.BAM_CDEL:
            res_seq = ' ' * val + res_seq  # placeholders for deletion
            fe -= val

        elif key == S.BAM_CINS:
            cb = ce - val
            res_seq = ctg_seq[cb: ce] + res_seq
            ce = cb

        else:
            err = ("cigar '{0}' hasn't been delta properly "
                   "for '{1}' strand, please report".format(key, strand))
            raise NotImplementedError(err)

        if fe <= target_fe:
            res_seq = res_seq[target_fe - fe + 1:]
            break

    # remove placeholders for deletion
    res_seq = res_seq.replace(' ', '')
    return res_seq
