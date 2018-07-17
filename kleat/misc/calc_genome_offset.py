import kleat.misc.settings as S


CIGARS_FOR_CONTIG_SEQ_LENGTH = {
    S.BAM_CSOFT_CLIP,
    S.BAM_CHARD_CLIP,

    S.BAM_CMATCH,
    S.BAM_CINS,

    S.BAM_CDIFF,
    S.BAM_CEQUAL,
}


CIGARS_FOR_GENOME_SEQ_LENGTH = {
    S.BAM_CMATCH,
    S.BAM_CREF_SKIP,
    S.BAM_CDEL,

    S.BAM_CDIFF,
    S.BAM_CEQUAL,
}


def calc_contig_seq_len(cigartuples):
    return sum(_[1] for _ in cigartuples
               if _[0] in CIGARS_FOR_CONTIG_SEQ_LENGTH)


def calc_genome_seq_len(cigartuples):
    return sum(_[1] for _ in cigartuples
               if _[0] in CIGARS_FOR_GENOME_SEQ_LENGTH)


def calc_offset(cigartuples, ctg_clv, tail_side, skip_check_size):
    # current position in contig coordinate and genome offset coordinates
    ctg_pos, ref_pos = 0, 0

    last_cigar_is_skip = None
    last_skip_size = 0

    for key, val in cigartuples:
        if key == S.BAM_CSOFT_CLIP or key == S.BAM_CHARD_CLIP:
            ctg_clv -= val

        elif key == S.BAM_CMATCH:
            ctg_pos += val

            if ctg_pos > ctg_clv:
                delta = ctg_pos - ctg_clv
                next_ref_pos = ref_pos + val - delta

                if last_cigar_is_skip:
                    if skip_check_size is not None and next_ref_pos - ref_pos < skip_check_size:
                        ref_pos = ref_pos - last_skip_size - 1
                    else:
                        ref_pos = next_ref_pos
                else:
                    ref_pos = next_ref_pos
                break
            else:
                ref_pos = ref_pos + val

        elif key == S.BAM_CREF_SKIP or key == S.BAM_CDEL:
            ref_pos += val
            last_skip_size = val

        elif key == S.BAM_CINS:
            ctg_pos += val
            if ctg_pos > ctg_clv:
                if tail_side == 'left':
                    ref_pos -= 1
                    break
                else:
                    break
        else:
            err_msg = (' S.BAM_CEQUA, S.BAM_CDIFF & S.BAM_CPAD & BAM_CBACK '
                       'cigar value are note implemented yet. Please report. '
                       'Your cigar: ({0}, {1})\n'
                       '{2}'.format(key, val, S.CIGAR_TABLE))
            raise NotImplementedError(err_msg)

        if key == S.BAM_CREF_SKIP:
            last_cigar_is_skip = True
        else:
            last_cigar_is_skip = False
    return ref_pos


def calc_genome_offset(cigartuples, ctg_clv, tail_side='left', skip_check_size=None):
    """Calculate the genome offset based on clv in the contig coordinate and the
    contig cigars when it's aligned to the reference genome.

    Note the ctg_clv and tail_base should be in contig coordinate in the same
    direction of the reference genome. If the contig is reversed when aligned
    to the genome, the contig should be flipped to derive its corresponding
    ctg_clv and tail_base serving this function. see corresponding test cases
    for detail.

    ctg_clv is independent of the reference genome except that its direction
    needs to be the same as the reference genome. When the contig is aligned to
    it may include softclip, hardclip, skips (e.g introns) and indels, these
    regions are properly dealt by this function to obtain an genome offset that
    is directly addable to the genome coordinate.

    For example, This genome offset is needed for inferring the clv in genomic
    coordinate,

    i.e. gnm_clv = ctg.reference_start + gnm_offset

    :param ctg_clv: clv in contig coordinate.
    :param tail_side: T or A, meaning the tail is polyA or polyT. If not
    specified, default to T. This parameter is only used when the clv happens
    to within a insertion. In such case, the exact coordinate of the clv in the
    genome is not available, but trying to be as accurate as possible, one
    strategy is to use the position of leftmost genome base around the
    insertion when it's polyT; otherwise, use the position of the rightmost
    base.

    """
    if tail_side == 'left':     # flip
        ctg_seq_len = calc_contig_seq_len(cigartuples)
        cigartuples = list(reversed(cigartuples))
        ctg_clv = ctg_seq_len - ctg_clv - 1

    ref_offset = calc_offset(cigartuples, ctg_clv, tail_side, skip_check_size)

    if tail_side == 'left':     # flip back
        ref_seq_len = calc_genome_seq_len(cigartuples)
        ref_offset = ref_seq_len - ref_offset - 1

    return ref_offset
