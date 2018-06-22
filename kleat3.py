import csv
import pysam


# https://pysam.readthedocs.io/en/latest/api.html?highlight=AlignmentSegment#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0
BAM_CSOFT_CLIP = 4              # used to identify soft clip 


def calc_strand(contig):
    return '+' if contig.is_reverse else '-'


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
        # 0            +---AAA          l              -> contig coordinates
        # +-----------------------------+              -> the contig
        # a                ref_start    b:abs_ref_end  -> genomic coordinates
        ref_clv = abs_ref_end - read.reference_start
    else:
        abs_ref_beg = infer_contig_abs_ref_start(contig)
        # 0              TTT---+        l              -> contig coordinates
        # +-----------------------------+              -> the contig
        # a:abs_ref_beg     ref_start   b              -> genomic coordinates
        ref_clv = abs_ref_beg + read.reference_start
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
            'seqname', 'strand', 'clv', 'contig_id', 'polya_evidence_type',
            'num_contig_tail_reads', 'contig_tail_length',
            'num_bridge_reads', 'max_bridge_tail_length',
        ])

        for k, contig in enumerate(c2g_sam):
            if (k + 1) % 1000 == 0:
                print(f'processed {k + 1} contigs')

            tmp_dd[contig.query_name] = contig
            if contig.is_unmapped:
                continue

            seqname = contig.reference_name
            strand = calc_strand(contig)
            if is_tail_segment(contig):
                ref_clv = calc_ref_clv(contig)
                num_tail_reads, _ = calc_num_tail_reads(contig, r2c_sam)
                tail_length = calc_tail_length(contig)

                clv_record = (
                    seqname, strand, ref_clv, contig.query_name,
                    'tail_contig',
                    num_tail_reads, tail_length,
                    0, 0      # bridge_contig evidence is left empty
                )
                csvwriter.writerow(clv_record)
            else:
                num_bdg_reads_dd = {}
                max_bdg_tail_len_dd = {}
                for read in r2c_sam.fetch(
                        contig.query_name, 0, contig.query_length):
                    # being unmapped is still possible, see test for a reason
                    # r2c (cDNA) alignment, tail is always T, must not reverse
                    if (read.is_unmapped or read.is_reverse):
                        continue

                    if is_forward_tseg(read):
                        ref_clv = calc_ref_clv_from_r2c_alignment(contig, read)
                        tail_len = calc_tail_length(read)
                        num_bdg_reads_dd[ref_clv] = num_bdg_reads_dd.get(ref_clv, 0) + 1
                        max_bdg_tail_len_dd[ref_clv] = max(max_bdg_tail_len_dd.get(ref_clv, 0), tail_len)
                if len(num_bdg_reads_dd.keys()) > 0:
                    for ref_clv in num_bdg_reads_dd:
                        clv_record = (
                            seqname, strand, ref_clv, contig.query_name,
                            'bridge_contig',
                            0, 0,  # tail_contig evidence is left empty
                            num_bdg_reads_dd[ref_clv], max_bdg_tail_len_dd[ref_clv]
                        )
                        csvwriter.writerow(clv_record)
                else:
                    # assume there is still a clv at the 3' end of the contig
                    # even without any polya evidence
                    ref_clv = calc_ref_clv(contig)
                    clv_record = (
                        seqname, strand, ref_clv, contig.query_name,
                        'None', 0, 0, 0, 0
                    )
                    csvwriter.writerow(clv_record)
