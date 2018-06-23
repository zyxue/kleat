import csv
import pysam


# https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0
BAM_CSOFT_CLIP = 4              # used to identify soft clip 


def calc_strand(contig):
    return '+' if contig.is_reverse else '-'


def gen_contig_info(contig):
    return [
        contig.reference_name,
        calc_strand(contig),
        contig.query_name,
        contig.query_length,
        contig.mapq
    ]


def analyze_tail_contig(contig):
    contig_info = gen_contig_info(contig)
    ref_clv = calc_ref_clv(contig)
    num_tail_reads, _ = calc_num_tail_reads(contig, r2c_sam)
    tail_length = calc_tail_length(contig)

    clv_record = (
        *contig_info, ref_clv, 'tail_contig',
        # bridge_contig evidence is left empty
        num_tail_reads, tail_length, 0, 0, 0, 0
    )
    return clv_record


def analyze_bridge_contig(contig):
    pass


def is_tail_segment(segment):
    seg = segment
    if seg.is_reverse:
        return is_reverse_tail_segment(seg)
    else:
        return is_forward_tail_segment(seg)


def is_reverse_tail_segment(segment, tail_base='A'):
    """
    tseg: tail segment

    default tail_base applies main to alignment to genome, where the polyA tail
    strand is known
    """
    seq = segment.query_sequence
    last_cigar = segment.cigartuples[-1]
    # potential test case == "A0.R100820":
    return (
        seq.endswith(tail_base)
        # right clipped
        and last_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all As
        and set(seq[-last_cigar[1]:]) == {tail_base}
    )


def is_forward_tail_segment(segment, tail_base='T'):
    seq = segment.query_sequence
    first_cigar = segment.cigartuples[0]
    # potential test case A0.S36384
    return (
        seq.startswith(tail_base)
        # left clipped
        and first_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all Ts
        and set(seq[:first_cigar[1]]) == {tail_base}
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


def calc_ref_clv_from_r2c_alignment(contig, read_reference_start_wst_contig):
    """calculate
    cleavage site position wst the reference based on bridge read, and
    read2contig and contig2genome alignments
    """
    read_start = read_reference_start_wst_contig
    if contig.is_reverse:
        abs_ref_end = infer_contig_abs_ref_end(contig)
        # 0            +---AAA          l              -> contig coordinates
        # <-----------------------------+              -> the contig
        # a                ref_start    b:abs_ref_end  -> genomic coordinates
        ref_clv = abs_ref_end - read_start
    else:
        abs_ref_beg = infer_contig_abs_ref_start(contig)
        # 0              TTT---+        l              -> contig coordinates
        # +----------------------------->              -> the contig
        # a:abs_ref_beg     ref_start   b              -> genomic coordinates
        ref_clv = abs_ref_beg + read_start
    return ref_clv


if __name__ == "__main__":
    c2g_sam = pysam.AlignmentFile('../kleat3-test-data/tasrkleat-results/align_contigs2genome/cba.sorted.bam')
    r2c_sam = pysam.AlignmentFile('../kleat3-test-data/tasrkleat-results/align_reads2contigs/cba.sorted.bam')

    output = './lele.csv'

    # useful for debugging, remove later
    tmp_dd = {}

    with open(output, 'wt') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerow([
            'seqname', 'strand', 'contig_id', 'contig_len', 'contig_mapq',
            'clv', 'polya_evidence_type',
            'num_contig_tail_reads', 'contig_tail_length',
            'num_bridge_reads', 'max_bridge_tail_length',
            'num_link_reads',
            'num_blank_contigs'   # those contigs that provide no polya evidence
        ])

        for k, contig in enumerate(c2g_sam):
            if (k + 1) % 1000 == 0:
                print(f'processed {k + 1} contigs')

            tmp_dd[contig.query_name] = contig

            if (
                    contig.is_unmapped or
                    contig.mapping_quality == 0
            ):
                continue

            if is_tail_segment(contig):
                clv_record = analyze_tail_contig(contig)
                csvwriter.writerow(clv_record)
            else:
                contig_info = gen_contig_info(contig)
                num_bdg_reads_dd = {}
                max_bdg_tail_len_dd = {}

                num_link_reads_dd = {}
                for read in r2c_sam.fetch(
                        contig.query_name, 0, contig.query_length):
                    # if read.query_name == "SN7001282:314:h15b0adxx:1:2102:16317:20542" and read.is_unmapped:
                    #     sys.exit(1)

                    # if read.query_sequence.startswith('T' * 50) and contig.query_length > 150:
                    #     import pdb; pdb.set_trace()

                    if read.is_reverse:
                        # read2contig alignment cannot be reverse to be a read
                        # supporting polyA tail
                        continue

                    if read.is_unmapped:  # a candidate for link read

                        # in principle, could also check from the perspecitve
                        # and the mate of a link read, but it would be harder
                        # to verify the sequence composition of the link read
                        if read.mate_is_unmapped:
                            continue

                        if not set(read.query_sequence) == {'T'}:
                            # needs to be completely T
                            continue

                        if not read.reference_id == read.next_reference_id:
                            # assume reference_id comparison is faster than
                            # reference_name
                            continue

                        # assume the start/end of the paired read is the
                        # cleavage site
                        ref_clv = calc_ref_clv_from_r2c_alignment(contig, read.next_reference_start)

                        # for debug purpose
                        # num_link_reads_dd[ref_clv] = num_link_reads_dd.get(ref_clv, []) + [f'{read.query_sequence}']
                        num_link_reads_dd[ref_clv] = num_link_reads_dd.get(ref_clv, 0) + 1

                    else:       # a candidate for bridge read
                        if is_forward_tail_segment(read):
                            ref_clv = calc_ref_clv_from_r2c_alignment(contig, read.reference_start)
                            tail_len = calc_tail_length(read)
                            num_bdg_reads_dd[ref_clv] = num_bdg_reads_dd.get(ref_clv, 0) + 1
                            max_bdg_tail_len_dd[ref_clv] = max(max_bdg_tail_len_dd.get(ref_clv, 0), tail_len)

                if len(num_bdg_reads_dd) > 0:
                    for ref_clv in num_bdg_reads_dd:
                        clv_record = (
                            *contig_info, ref_clv, 'bridge_contig',
                            # tail_contig evidence is left empty
                            0, 0, num_bdg_reads_dd[ref_clv], max_bdg_tail_len_dd[ref_clv], 0, 1
                        )
                        csvwriter.writerow(clv_record)
                if len(num_link_reads_dd) > 0:
                    for ref_clv in num_link_reads_dd:
                        clv_record = (
                            *contig_info, ref_clv, 'link_contig',
                            0, 0, 0, 0, num_link_reads_dd[ref_clv], 1
                        )
                        csvwriter.writerow(clv_record)
                else:
                    # assume there is still a clv at the 3' end of the contig
                    # even without any polya evidence
                    ref_clv = calc_ref_clv(contig)

                    clv_record = (
                        *contig_info, ref_clv, 'None',
                        0, 0, 0, 0, 0, 1
                    )
                    csvwriter.writerow(clv_record)
