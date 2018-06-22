import pysam

from kleat3 import (
    calc_strand,
    calc_ref_clv,
    calc_tail_length,
    calc_num_tail_reads,
    is_forward_tseg,
    calc_ref_clv_from_r2c_alignment,
)


C2G_BAM_FILE = '../kleat3-test-data/tasrkleat-results/align_contigs2genome/cba.sorted.bam'
R2C_BAM_FILE = '../kleat3-test-data/tasrkleat-results/align_reads2contigs/cba.sorted.bam'

C2G_BAM = pysam.AlignmentFile(C2G_BAM_FILE)
R2C_BAM = pysam.AlignmentFile(R2C_BAM_FILE)


def test_tail_contigs():
    for k, contig in enumerate(C2G_BAM):
        if contig.query_name not in [
                "A0.R100820",
                "A1.R26141",

                'A0.R101981',
        ]:
            continue

        strand = calc_strand(contig)
        ref_clv = calc_ref_clv(contig)
        contig_tail_len = calc_tail_length(contig)
        num_contig_tail_reads, contig_clv = calc_num_tail_reads(contig, R2C_BAM)

        # ERBB2
        if contig.query_name == "A0.R100820":
            assert contig.reference_name == "chr17"
            assert strand == '+'
            assert ref_clv == 37884912
            assert contig.is_reverse is True
            assert contig_tail_len == 20
            assert num_contig_tail_reads == 2
            assert contig_clv == 4668

        # BAP1
        if contig.query_name == "A1.R26141":
            assert contig.reference_name == "chr3"
            assert strand == '-'
            assert ref_clv == 52435023
            assert contig.is_reverse is False
            assert contig_tail_len == 9
            assert num_contig_tail_reads == 1
            assert contig_clv == 8

        # An example of short contig supporting a UNconfident CS
        if contig.query_name == 'A0.R101981':
            assert contig.reference_name == "chr10"
            assert strand == '-'
            assert ref_clv == 6132816
            assert num_contig_tail_reads == 1
            assert contig_tail_len == 28
            assert contig.query_length == 80


def test_bridge_reads():
    for k, contig in enumerate(C2G_BAM):
        if contig.query_name not in [
                "A0.R100710",
                "A1.S26245",
        ]:
            continue

        strand = calc_strand(contig)
        for read in R2C_BAM.fetch(
                contig.query_name, 0, contig.query_length):
            if read.is_unmapped or read.is_reverse:
                # still possible a read is unmapped even though fetching
                # used a specific contig, e.g.

# could be an example for link reads
# SN7001282:314:h15b0adxx:1:2206:17178:42842      73      A0.S68856       1       60      75M     =       1       0       GGGGGACAGATCTTCAGTTCTCATGACCACAAAAGAGGATACTAAAGCTCAGACAGGAGAAGAGACGTGGCCAGC     CCCFFFFFHHHHHJJJJHIIIJJIJJJIJIJJJIGIIIICGIJJJJGHHGHBHCHIEHEHIEHCHEHABBEEECC     NM:i:0  MD:Z:75 AS:i:75 XS:i:0
# SN7001282:314:h15b0adxx:1:2206:17178:42842      133     A0.S68856       1       0       *       =       1       0       TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT     +1:A:?)8<CFDHA@/4337=:BBBB6630&07670:7317@B@B@5@88889<<<@57<B7@BBBBB;;>>>5<     AS:i:0  XS:i:0
                continue

            if not is_forward_tseg(read):
                continue

            ref_clv = calc_ref_clv_from_r2c_alignment(contig, read)

            # PTEN
            if contig.query_name == "A0.R100710":
                assert contig.reference_name == "chr10"
                assert strand == '+'

                # there are two clvs
                if read.query_name in [
                        'SN7001282:314:h15b0adxx:1:2203:2771:88598',
                        'SN7001282:314:h15b0adxx:2:1105:18463:70356'
                ]:
                    assert ref_clv == 89725516
                else:
                    assert ref_clv == 89725287

            # KRAS
            if contig.query_name == "A1.S26245":
                assert contig.reference_name == "chr12"
                assert strand == '-'
                assert ref_clv == 25362769

            # potential test case for bridge read
            # PTEN	ENST00000371953	+	yes	A0.R100710	chr10	89725287, ref_clv: 89725516, 89725287
            # KRAS	ENST00000256078	-	yes	A1.S26245	chr12	25362769, ref_clv: 25362769
