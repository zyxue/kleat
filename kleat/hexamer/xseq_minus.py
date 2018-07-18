import kleat.misc.settings as S
from kleat.misc.apautils import (
    calc_genome_offset,
    fetch_seq                   # deal with circular DNA, too, e.g. chrM
)

"""xseq: extract_seq"""


def do_match(ctg_b, ref_b, cigar_val, ctg_seq, ctg_clv):
    ctg_e = ctg_b + cigar_val
    ref_e = ref_b + cigar_val

    if ctg_b <= ctg_clv:
        if ctg_e > ctg_clv:
            seq_to_add = ctg_seq[ctg_clv: ctg_e]
        else:
            seq_to_add = ''
    else:
        seq_to_add = ctg_seq[ctg_b: ctg_e]

    next_ctg_b = ctg_e
    next_ref_b = ref_e
    return next_ctg_b, next_ref_b, seq_to_add


def do_skip(ref_b, cigar_val, ref_fa, seqname, ref_clv):
    """
    Handle BAM CREF_SKIP CIGAR
    """
    ref_e = ref_b + cigar_val

    if ref_b <= ref_clv:
        if ref_e > ref_clv:
            seq_to_add = fetch_seq(ref_fa, seqname, ref_clv, ref_e)
        else:           # still before clv
            seq_to_add = ''
    else:
        seq_to_add = fetch_seq(ref_fa, seqname, ref_b, ref_e)

    next_ref_b = ref_e
    return next_ref_b, seq_to_add


def init_ctg_beg(ctg_seq):
    return 0


def init_ref_beg(ref_clv, cigartuples, ctg_clv):
    """
    Initialize the beginning index in genome coordinate by calculating the
    offset from left, using `calc_genome_offset`
    """
    # TODO: left may not matter in such case
    offset = calc_genome_offset(cigartuples, ctg_clv, 'left')

    # ref_clv = contig.reference_start + offset, so
    # ref_clv - offset = contig.reference_start,
    # so still needs to subtract clipped region, which isn't in reference_start

    cgr = cigartuples[0]
    if cgr[0] == S.BAM_CSOFT_CLIP or cgr[0] == S.BAM_CHARD_CLIP:
        offset += cgr[1]
    return ref_clv - offset


def extract(cigartuples, ctg_seq, seqname, strand, ctg_clv, ref_clv, ref_fa, window):
    """
    scan from Left => Right

    compared to plus strand, initialize contig beginning (cb) and reference
    beginning (fb) instead of ce and fe

    use fe instead of re because re is a Python module for regex
    """
    # aka. ctg_b: contig beginning
    cb = init_ctg_beg(ctg_seq)
    # aka. ref_b: reference genome ending
    fb = init_ref_beg(ref_clv, cigartuples, ctg_clv)

    res_seq = ''
    for idx, (key, val) in enumerate(cigartuples):
        if key == S.BAM_CSOFT_CLIP or key == S.BAM_CHARD_CLIP:
            ce = cb + val
            if cb <= ctg_clv:
                if ce > ctg_clv:  # meaning clip spans the clv
                    seq = ctg_seq[ctg_clv: ce]
                    res_seq += seq
                cb += val
                fb += val
            else:
                seq = ctg_seq[cb: cb + val]
                res_seq += seq
                cb += val
                fb += val

        elif key == S.BAM_CMATCH:
            cb, fb, seq = do_match(cb, fb, val, ctg_seq, ctg_clv)
            res_seq += seq

        elif key == S.BAM_CREF_SKIP:
            fb, seq = do_skip(fb, val, ref_fa, seqname, ref_clv)
            res_seq += seq

        elif key == S.BAM_CDEL:
            fb += val

        elif key == S.BAM_CINS:
            res_seq += ctg_seq[cb: cb + val]
            cb += val

        else:
            err = ("cigar '{0}' hasn't been delta properly "
                   "for '{1}' strand, please report".format(key, strand))
            raise NotImplementedError(err)

        if len(res_seq) >= window:
            res_seq = res_seq[:window]
            break

    return res_seq
