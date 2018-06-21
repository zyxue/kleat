import pysam


# https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
BAM_CSOFT_CLIP = 4              # used to identify soft clip 


def scan_c2g_alignments(samfile):
    """classify contigs into tail and potential bridge"""
    pass


def scan_r2c_alignments():
    """identify tail reads and bridge reads"""


def is_tail_contig(segment):
    seg = segment
    if seg.is_reverse:
        return is_reverse_tc(seg)
    else:
        return is_forward_tc(seg)


def is_reverse_tc(segment):
    """tc: tail contig"""
    seq = segment.query_sequence
    last_cigar = seg.cigartuples[-1]
    # potential test case A0.S36384
    return (
        seq.endswith('A')
        # right clipped
        and last_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all As
        and set(seq[-last_cigar[1]:]) == {'A'}
    )


def is_forward_tc(segment):
    seq = segment.query_sequence
    first_cigar = seg.cigartuples[0]
    # potential test case == "A0.R100820":
    return (
        seq.startswith('T')
        # left clipped
        and first_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all Ts
        and set(seq[:first_cigar[1]]) == {'T'}
    )


if __name__ == "__main__":

    # samfile = pysam.AlignmentFile('./test-dataset/assembly/c2g.bam')
    samfile = pysam.AlignmentFile('./cba.sorted.bam')

    for k, seg in enumerate(samfile):
        if seg.is_unmapped:
            continue

        tail_contigs = []
        if is_tail_contig(seg):
            print(seg.query_sequence, seg.reference_name, seg.reference_start, seg.reference_end)
