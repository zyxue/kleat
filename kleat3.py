import csv
from collections import defaultdict

import pysam

from evidence import suffix, bridge, link, blank
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
            contig_is_blank = True
            if (k + 1) % 1000 == 0:
                print(f'processed {k + 1} contigs')

            tmp_dd[contig.query_name] = contig

            if (contig.is_unmapped or contig.mapping_quality == 0):
                continue

            # suffix evidence
            if U.has_tail(contig):
                clv_record = suffix.gen_clv_record(contig, r2c_bam)
                contig_is_blank = False
                U.write_row(clv_record, csvwriter)

            dd_bridge = {
                'num_reads': defaultdict(int),
                'max_tail_len': defaultdict(int),
            }

            dd_link = {
                'num_reads': defaultdict(int)
            }

            query_region = [0, contig.query_length]
            for read in r2c_bam.fetch(contig.query_name, *query_region):
                # if read.query_name == "SN7001282:314:h15b0adxx:1:2102:16317:20542" and read.is_unmapped:
                #     sys.exit(1)

                if not read.is_unmapped:
                    if U.has_tail(read):
                        seqname, strand, ref_clv, tail_len = \
                            bridge.analyze_bridge_read(contig, read)

                        clv_key = U.gen_clv_key_tuple(seqname, strand, ref_clv)
                        dd_bridge['num_reads'][clv_key] += 1
                        dd_bridge['max_tail_len'][clv_key] = max(
                            dd_bridge['max_tail_len'][clv_key], tail_len)
                else:
                    # Here we focused on the unmapped all A/T read, but in
                    # principle, we could also check from the perspecitve and
                    # the mate of a link read, but it would be harder to verify
                    # the sequence composition of the link read
                    if (not read.mate_is_unmapped
                        # assume reference_id comparison is faster than
                        # reference_name
                        and read.reference_id == read.next_reference_id
                        and set(read.query_sequence) in [{'A'}, {'T'}]):
                        seqname, strand, ref_clv = \
                            link.analyze_link(contig, read)

                        clv_key = U.gen_clv_key_tuple(seqname, strand, ref_clv)
                        dd_link['num_reads'][clv_key] += 1

            if len(dd_bridge['num_reads']) > 0:
                contig_is_blank = False
                for clv_key in dd_bridge['num_reads']:
                    clv_record = bridge.gen_clv_record(
                        contig, clv_key,
                        dd_bridge['num_reads'][clv_key],
                        dd_bridge['max_tail_len'][clv_key]
                    )
                    U.write_row(clv_record, csvwriter)

            if len(dd_link['num_reads']) > 0:
                contig_is_blank = False
                for ref_clv in dd_link['num_reads']:
                    clv_record = link.gen_clv_record(
                        contig, clv_key, dd_link['num_reads'][clv_key])
                    U.write_row(clv_record, csvwriter)

            if contig_is_blank:
                for clv_rec in blank.gen_two_clv_records(contig):
                    U.write_row(clv_rec, csvwriter)
