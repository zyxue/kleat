import kleat.misc.settings as S
from kleat.misc.apautils import calc_genome_offset

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
            seq_to_add = ref_fa.fetch(seqname, ref_b, ref_clv + 1)
        else:
            seq_to_add = ''
    else:
        seq_to_add = ref_fa.fetch(seqname, ref_b, ref_e)

    next_ref_e = ref_b
    return next_ref_e, seq_to_add


def init_ref_end(ref_clv, cigartuples, ctg_seq, ctg_clv):
    """
    Taking advantage of `calc_genome_offset`, initialize the reference end
    index by calculating the genome offset from right
    """
    rev_cigartuples = reversed(cigartuples)
    ctg_clv_from_the_right = len(ctg_seq) - ctg_clv
    tail_side = 'left'           # reverse complement the + strand

    ref_end = ref_clv + calc_genome_offset(
        rev_cigartuples,
        ctg_clv_from_the_right,
        tail_side
    )
    return ref_end


def xseq(cigartuples, ctg_seq, seqname, strand, ctg_clv, ref_clv, ref_fa, window):
    """
    Scan from Right to Left,

    compared to minus strand, initialize contig end (ce) and reference end (fe)
    instead of cb and fb

    use fe instead of re because re is a Python module for regex
    """

    ce = len(ctg_seq)
    fe = init_ref_end(ref_clv, cigartuples, ctg_seq, ctg_clv)

    res_seq = ''
    for idx, (key, val) in enumerate(reversed(cigartuples)):
        if key in [S.BAM_CSOFT_CLIP, S.BAM_CHARD_CLIP]:
            if idx == 0:        # meaning its the upstream clip
                ce -= val
            # else, downstream clip should just be ignored

        elif key == S.BAM_CMATCH:
            ce, fe, seq = do_match(ce, fe, val, ctg_seq, ctg_clv)
            res_seq = seq + res_seq

        elif key == S.BAM_CREF_SKIP:
            fe, seq = do_skip(fe, val, ref_fa, seqname, ref_clv)
            res_seq = seq + res_seq

        elif key == S.BAM_CDEL:
            fe -= val

        elif key == S.BAM_CINS:
            cb = ce - val
            res_seq = ctg_seq[cb: ce] + res_seq
            ce = cb

        else:
            err = ("cigar '{0}' hasn't been delta properly "
                   "for '{1}' strand, please report".format(key, strand))
            raise NotImplementedError(err)

        if len(res_seq) >= window:
            res_seq = res_seq[-window:]
            break
    return res_seq
