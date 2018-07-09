import kleat.misc.settings as S


def extract_seq_for_plus_strand(cigartuples, ctg_seq, seqname, strand,
                                ctg_idx, ref_idx, ref_fa, window):
    res_seq = ''
    for idx, (key, val) in enumerate(reversed(cigartuples)):
        if key in [S.BAM_CSOFT_CLIP, S.BAM_CHARD_CLIP]:
            if idx == 0:
                # meaning it's the upstream clip, downstream clip should
                # just be ignored
                ctg_idx -= val
                ref_idx -= val
        elif key in [S.BAM_CMATCH]:
            res_seq = ctg_seq[ctg_idx - val + 1: ctg_idx + 1] + res_seq
            ctg_idx -= val
            ref_idx -= val
        elif key in [S.BAM_CREF_SKIP]:
            ref_seq = ref_fa.fetch(seqname, ref_idx - val + 1, ref_idx + 1)
            res_seq = ref_seq + res_seq
            ref_idx -= val
        elif key in [S.BAM_CDEL]:
            ref_idx -= val
        else:
            err = ("cigar '{0}' hasn't been delta properly "
                   "for '{1}' strand, please report".format(key, strand))
            raise NotImplementedError(err)
        if len(res_seq) >= window:
            res_seq = res_seq[-window:]
            break
    return res_seq


def extract_seq_for_minus_strand(cigartuples, ctg_seq, seqname, strand,
                                 ctg_idx, ref_idx, ref_fa, window):
    res_seq = ''
    for idx, (key, val) in enumerate(cigartuples):
        if key in [S.BAM_CSOFT_CLIP, S.BAM_CHARD_CLIP]:
            if idx == 0:
                ctg_idx += val
                ref_idx += val
        elif key in [S.BAM_CMATCH]:
            res_seq += ctg_seq[ctg_idx: ctg_idx + val]
            ctg_idx += val
            ref_idx += val
        elif key in [S.BAM_CREF_SKIP]:
            print('ref_idx', ref_idx)
            ref_seq = ref_fa.fetch(seqname, ref_idx, ref_idx + val)
            res_seq += ref_seq
        else:
            err = ("cigar '{0}' hasn't been delta properly "
                   "for '{1}' strand, please report".format(key, strand))
            raise NotImplementedError(err)
        if len(res_seq) >= window:
            res_seq = res_seq[:window]
            break
    return res_seq


def extract_seq(contig, strand, ref_clv, ref_fa, window=50, ctg_clv=0):
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
        init_ctg_idx = len(ctg_seq) - 1
        init_ref_idx = ref_clv + (len(ctg_seq) - ctg_clv - 1)
        func = extract_seq_for_plus_strand
    elif strand == '-':
        init_ctg_idx = 0
        init_ref_idx = ref_clv - (len(ctg_seq) - ctg_clv)
        func = extract_seq_for_minus_strand
    else:
        raise ValueError('unknown strand: "{0}"'.format(strand))

    return func(cigartuples, ctg_seq, seqname, strand,
                init_ctg_idx, init_ref_idx, ref_fa, window)
