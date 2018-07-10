import kleat.misc.settings as S
from kleat.misc.apautils import calc_genome_offset


"""Extract relevatn sequence, in which PAS hexamer is searched"""


def extract_seq_for_plus_strand(cigartuples, ctg_seq, seqname, strand,
                                ctg_clv, ref_clv, ref_fa, window):
    """scan from Right to Left, init ctg_end and ref_end first"""
    ctg_e = len(ctg_seq)
    # take advantage of calc_genome_offset to calculate the offset from right
    ref_e = ref_clv + calc_genome_offset(
        reversed(cigartuples), len(ctg_seq) - ctg_clv, 'left')
    print(ref_e)

    res_seq = ''
    for idx, (key, val) in enumerate(reversed(cigartuples)):
        print(res_seq)
        if key in [S.BAM_CSOFT_CLIP, S.BAM_CHARD_CLIP]:
            if idx == 0:
                # meaning it's the upstream clip, downstream clip should
                # just be ignored
                ctg_e -= val
        elif key in [S.BAM_CMATCH]:
            ctg_b = ctg_e - val
            if ctg_e >= ctg_clv:
                if ctg_b < ctg_clv:
                    _seq = ctg_seq[ctg_b: ctg_clv + 1]
                else:
                    _seq = ''
            else:
                _seq = ctg_seq[ctg_b: ctg_e]

            res_seq = _seq + res_seq

            ctg_e = ctg_b
            ref_e -= val
        elif key in [S.BAM_CREF_SKIP]:
            ref_b = ref_e - val
            if ref_e >= ref_clv:
                if ref_b < ref_clv:
                    _seq = ref_fa.fetch(seqname, ref_b, ref_clv + 1)
                else:
                    _seq = ''
            else:
                _seq = ref_fa.fetch(seqname, ref_b, ref_e)
            res_seq = _seq + res_seq

            ref_e -= val
        elif key in [S.BAM_CDEL]:
            ref_e -= val
        elif key in [S.BAM_CINS]:
            ctg_b = ctg_e - val
            _seq = ctg_seq[ctg_b: ctg_e]
            res_seq = _seq + res_seq
            ctg_e -= val
        else:
            err = ("cigar '{0}' hasn't been delta properly "
                   "for '{1}' strand, please report".format(key, strand))
            raise NotImplementedError(err)
        if len(res_seq) >= window:
            res_seq = res_seq[-window:]
            break
    return res_seq


def calc_ref_beg(ref_clv, strand, cigartuples, ctg_clv):
    """calculate the beginning index in ref cooridnate"""
    tail_dir = 'left' if strand == '-' else 'right'
    ref_offset = calc_genome_offset(cigartuples, ctg_clv, tail_dir)
    ref_beg = ref_clv - ref_offset
    return ref_beg


def extract_seq_for_minus_strand(cigartuples, ctg_seq, seqname, strand,
                                 ctg_clv, ref_clv, ref_fa, window):
    """scan from Left => Right"""
    ctg_b = 0
    ref_b = calc_ref_beg(ref_clv, strand, cigartuples, ctg_clv)

    res_seq = ''
    for idx, (key, val) in enumerate(cigartuples):
        if key in [S.BAM_CSOFT_CLIP, S.BAM_CHARD_CLIP]:
            if idx == 0:
                ctg_b += val
        elif key in [S.BAM_CMATCH]:
            ctg_e = ctg_b + val
            if ctg_b <= ctg_clv:
                if ctg_e > ctg_clv:
                    res_seq += ctg_seq[ctg_clv: ctg_e]
                else:
                    pass
            else:
                res_seq += ctg_seq[ctg_b: ctg_e]
            ctg_b = ctg_e
            ref_b += val
        elif key in [S.BAM_CREF_SKIP]:
            ref_e = ref_b + val
            if ref_b <= ref_clv:
                if ref_e > ref_clv:
                    res_seq += ref_fa.fetch(seqname, ref_clv, ref_e)
                else:           # still before clv
                    pass
            else:
                res_seq += ref_fa.fetch(seqname, ref_b, ref_e)
            ref_b = ref_e
        elif key in [S.BAM_CDEL]:
            ref_b += val
        elif key in [S.BAM_CINS]:
            res_seq += ctg_seq[ctg_b: ctg_b + val]
            ctg_b += val
        else:
            err = ("cigar '{0}' hasn't been delta properly "
                   "for '{1}' strand, please report".format(key, strand))
            raise NotImplementedError(err)
        if len(res_seq) >= window:
            res_seq = res_seq[:window]
            break
    return res_seq


def extract_seq(contig, strand, ref_clv, ref_fa, ctg_clv=0, window=50):
    """remove clipped ends before searching for hexamer, the clipped ends would
    affect calculation of genomics coordinates of the found PAS hexamer

    :param ref_clv: would be necessary needs to combine sequence from both
                    contig and reference genome
    :param ctg_clv: the position of clv in contig coordinate.
    """
    seqname = contig.reference_name
    ctg_seq = contig.query_sequence
    cigartuples = contig.cigartuples

    if strand == '+':
        return extract_seq_for_plus_strand(
            cigartuples, ctg_seq, seqname, strand,
            ctg_clv, ref_clv, ref_fa, window
        )
    elif strand == '-':
        return extract_seq_for_minus_strand(
            cigartuples, ctg_seq, seqname, strand,
            ctg_clv, ref_clv, ref_fa, window
        )
    else:
        raise ValueError('unknown strand: "{0}"'.format(strand))

