import csv
import pysam

from evidence import suffix, bridge, link
from settings import HEADER
import utils as U


if __name__ == "__main__":
    # datadir = '../kleat3-test-data/tasrkleat-results'
    datadir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/UHR/C1/tasrkleat-results'
    # datadir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/UHR/C2/tasrkleat-results'
    # datadir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/HBR/C4/tasrkleat-results'
    # datadir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/HBR/C6/tasrkleat-results'

    c2g_bam = pysam.AlignmentFile(f'{datadir}/align_contigs2genome/cba.sorted.bam')
    r2c_bam = pysam.AlignmentFile(f'{datadir}/align_reads2contigs/cba.sorted.bam')

    import sys
    output = sys.argv[1]

    # useful for debugging, remove later
    tmp_dd = {}

    with open(output, 'wt') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerow(HEADER)

        for k, contig in enumerate(c2g_bam):
            if (k + 1) % 1000 == 0:
                print(f'processed {k + 1} contigs')

            tmp_dd[contig.query_name] = contig

            if (contig.is_unmapped or contig.mapping_quality == 0):
                continue

            # suffix evidence
            if U.has_tail(contig):
                clv_record = suffix.gen_clv_record(contig, r2c_bam)
                csvwriter.writerow([getattr(clv_record, _) for _ in HEADER])
                continue

            # bridge or link evidence
            contig_is_blank = True

            dd_num_bdg_reads = {}
            dd_max_bdg_tail_len = {}

            dd_num_link_reads = {}

            query_region = [0, contig.query_length]
            for read in r2c_bam.fetch(contig.query_name, *query_region):
                # if read.query_name == "SN7001282:314:h15b0adxx:1:2102:16317:20542" and read.is_unmapped:
                #     sys.exit(1)

                if not read.is_unmapped:
                    if U.has_tail(read):
                        seqname, strand, ref_clv, tail_len = \
                            bridge.analyze_bridge_read(contig, read)
                        clv_key = (seqname, strand, ref_clv)
                        dd_num_bdg_reads[clv_key] = \
                            dd_num_bdg_reads.get(clv_key, 0) + 1
                        dd_max_bdg_tail_len[clv_key] = \
                            max(dd_max_bdg_tail_len.get(ref_clv, 0), tail_len)
                else:
                    # in principle, could also check from the perspecitve
                    # and the mate of a link read, but it would be harder
                    # to verify the sequence composition of the link read
                    if (not read.mate_is_unmapped
                        # assume reference_id comparison is faster than
                        # reference_name
                        and read.reference_id == read.next_reference_id
                        and set(read.query_sequence) in [{'A'}, {'T'}]):
                        link.analyze_link(read, contig)

            #         # for debug purpose
            #         # num_link_reads_dd[ref_clv] = num_link_reads_dd.get(ref_clv, []) + [f'{read.query_sequence}']
            #         num_link_reads_dd[ref_clv] = num_link_reads_dd.get(ref_clv, 0) + 1

            # contig_info = gen_contig_info(contig)
            # if len(num_bdg_reads_dd) > 0:
            #     contig_is_blank = False
            #     for ref_clv in num_bdg_reads_dd:
            #         clv_record = (
            #             *contig_info, ref_clv, 'bridge_contig',
            #             # tail_contig evidence is left empty
            #             0, 0, num_bdg_reads_dd[ref_clv], max_bdg_tail_len_dd[ref_clv], 0, 1
            #         )
            #         csvwriter.writerow(clv_record)

            # if len(num_link_reads_dd) > 0:
            #     contig_is_blank = False
            #     for ref_clv in num_link_reads_dd:
            #         clv_record = (
            #             *contig_info, ref_clv, 'link_contig',
            #             0, 0, 0, 0, num_link_reads_dd[ref_clv], 1
            #         )
            #         csvwriter.writerow(clv_record)

            # if contig_is_blank:
            #     # assume there is still a clv at the 3' end of the contig
            #     # even without any polya evidence, in thus case, there is
            #     # no direction, so either end of the contig could be a clv
            #     for strand, ref_clv in zip(['+', '-'], [
            #             contig.reference_end + 1,
            #             contig.reference_start
            #     ]):
            #         # replace strand  in contig_info
            #         ref_name, _, qn, ql, mapq = contig_info
            #         clv_record = (
            #             ref_name, strand, qn, ql, mapq,
            #             ref_clv, 'None',
            #             0, 0, 0, 0, 0, 1
            #         )
            #         csvwriter.writerow(clv_record)
