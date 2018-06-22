import pysam


# https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0
BAM_CSOFT_CLIP = 4              # used to identify soft clip 


def scan_c2g_alignments(samfile):
    """classify contigs into tail and potential bridge"""
    pass


def scan_r2c_alignments():
    """identify tail reads and bridge reads"""


def is_tail_segment(segment):
    seg = segment
    if seg.is_reverse:
        return is_reverse_tseg(seg)
    else:
        return is_forward_tseg(seg)


def is_reverse_tseg(segment):
    """tseg: tail segment"""
    seq = segment.query_sequence
    last_cigar = segment.cigartuples[-1]
    # potential test case == "A0.R100820":
    return (
        seq.endswith('A')
        # right clipped
        and last_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all As
        and set(seq[-last_cigar[1]:]) == {'A'}
    )


def is_forward_tseg(segment):
    seq = segment.query_sequence
    first_cigar = segment.cigartuples[0]
    # potential test case A0.S36384
    return (
        seq.startswith('T')
        # left clipped
        and first_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all Ts
        and set(seq[:first_cigar[1]]) == {'T'}
    )


def calc_ref_clv(segment):
    """calculate cleavage site position wst the reference"""
    if segment.is_reverse:
        return segment.reference_end + 1
    else:
        return segment.reference_start


def calc_tail_length(segment):
    """
    calculate A/T length of a tail contig, this information is extracted from
    softclip in the CIGAR
    """
    if segment.is_reverse:
        return segment.cigartuples[-1][1]
    else:
        return segment.cigartuples[0][1]


def calc_num_tail_reads(segment, r2c_sam):
    """
    calculate the number of reads aligned to the cleavage site

    https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
    """
    if segment.is_reverse:
        # clv position wst. the segment
        seg_clv = segment.query_alignment_end
    else:
        seg_clv = segment.query_alignment_start - 1
    num_tail_reads = r2c_sam.count(
        # region is half-open
        # https://pysam.readthedocs.io/en/latest/glossary.html#term-region
        segment.query_name, seg_clv, seg_clv + 1
    )
    return num_tail_reads, seg_clv


def infer_contig_abs_ref_start(contig):
    """
    infer the absolute reference starting position taking into consideration
    the non-M bases (esp. softclipped bases)
    """
    pos = contig.reference_start
    for key, val in contig.cigartuples:
        if key != BAM_CMATCH:
            pos -= val
        break
    return pos


def infer_contig_abs_ref_end(contig):
    """
    infer the absolute reference starting position taking into consideration
    the non-M bases (esp. softclipped bases)
    """
    pos = contig.reference_end
    for key, val in reversed(contig.cigartuples):
        if key != BAM_CMATCH:
            pos += val
        break
    return pos


def calc_ref_clv_from_r2c_alignment(contig, read):
    """calculate
    cleavage site position wst the reference based on bridge read, and
    read2contig and contig2genome alignments
    """
    if contig.is_reverse:
        abs_ref_end = infer_contig_abs_ref_end(contig)
        ref_clv = abs_ref_end - read.reference_start
    else:
        abs_ref_beg = infer_contig_abs_ref_start(contig)
        ref_clv = abs_ref_beg + read.reference_start
    return ref_clv



if __name__ == "__main__":
    c2g_sam = pysam.AlignmentFile('../kleat3-test-data/tasrkleat-results/align_contigs2genome/cba.sorted.bam')
    r2c_sam = pysam.AlignmentFile('../kleat3-test-data/tasrkleat-results/align_reads2contigs/cba.sorted.bam')

    tmp_dd = {}

    print('identifying tail contigs...')
    for k, contig in enumerate(c2g_sam):
        if contig.is_unmapped:
            continue

        tmp_dd[contig.query_name] = contig

        if is_tail_segment(contig):
            ref_clv = calc_ref_clv(contig)
            tail_length = calc_tail_length(contig)
            num_tail_reads, contig_clv = calc_num_tail_reads(contig, r2c_sam)

            # below are for debugging purpose
            if contig.query_name == "A0.R100820":
                test_A0_R100820(
                    contig, ref_clv, tail_length, num_tail_reads, contig_clv)

            if contig.query_name == "A1.R26141":
                test_A1_R26141(
                    contig, ref_clv, tail_length, num_tail_reads, contig_clv)
        else:
            # potential test case for bridge read
            # PTEN	ENST00000371953	+	yes	A0.R100710	chr10	89725287, ref_clv: 89725516, 89725287
            # KRAS	ENST00000256078	-	yes	A1.S26245	chr12	25362769, ref_clv: 25362769

            # loop through all reads that are aligned to this contig looking
            # for bridge reads

            if not contig.query_name == "A0.R100710":
            # if not contig.query_name == "A1.S26245":
                continue
            # else:
            #     pass

            for read in r2c_sam.fetch(
                    contig.query_name, 0, contig.query_length):
                if read.is_unmapped or read.is_reverse:
                    # still possible a read is unmapped even though fetching
                    # used a specific contig, e.g. SN7001282:314:h15b0adxx:1:2206:17178:42842

                    # bridge read should always not be reversed, and its tail
                    # letter is A or T depending on the contig orientation
                    continue

                seq = read.query_sequence
                first_cigar = read.cigartuples[0]
                if (
                        seq.startswith('T')  # cDNA always tails in T
                        and first_cigar[0] == BAM_CSOFT_CLIP
                        and set(seq[:first_cigar[1]]) == {'T'}
                ):
                    ref_clv = calc_ref_clv_from_r2c_alignment(contig, read)

                    if contig.query_name == "A0.R100710":
                        if read.query_name in [
                                'SN7001282:314:h15b0adxx:1:2203:2771:88598',
                                'SN7001282:314:h15b0adxx:2:1105:18463:70356'
                        ]:
                            assert ref_clv == 89725516
                        else:
                            assert ref_clv == 89725287
                    if contig.query_name == "A1.S26245":
                        assert ref_clv == 25362769

                    print(ref_clv)
                    print(f'{read.query_sequence:75s}\t{read.cigarstring:20s}\t{read.reference_start}\t{read.reference_name}\t{read.is_reverse}\t{read.query_name}')


    # print('identifying bridge reads...')
    # for k, read in enumerate(r2c_sam):
    #     if read.is_unmapped:
    #         continue

    #     if is_tail_segment(read) and read.is_reverse:
    #         break
