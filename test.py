import pysam

from kleat3 import calc_ref_clv, calc_tail_length, calc_num_tail_reads


C2G_BAM_FILE = '../kleat3-test-data/tasrkleat-results/align_contigs2genome/cba.sorted.bam'
R2C_BAM_FILE = '../kleat3-test-data/tasrkleat-results/align_reads2contigs/cba.sorted.bam'

C2G_BAM = pysam.AlignmentFile(C2G_BAM_FILE)
R2C_BAM = pysam.AlignmentFile(R2C_BAM_FILE)


def test_tail_contigs():
    for k, contig in enumerate(C2G_BAM):
        if contig.query_name not in [
                "A0.R100820",
                "A1.R26141"
        ]:
            continue

        ref_clv = calc_ref_clv(contig)
        tail_length = calc_tail_length(contig)
        num_tail_reads, contig_clv = calc_num_tail_reads(contig, R2C_BAM)

        if contig.query_name == "A0.R100820":
            assert contig.is_reverse is True
            assert contig.reference_name == "chr17"
            assert ref_clv == 37884912
            assert tail_length == 20
            assert num_tail_reads == 2
            assert contig_clv == 4668

        elif contig.query_name == "A1.R26141":
            assert contig.is_reverse is False
            assert contig.reference_name == "chr3"
            assert ref_clv == 52435023
            assert tail_length == 9
            assert num_tail_reads == 1
            assert contig_clv == 8


def test_bridge_reads():
    pass
