import os
import csv
import logging
import tempfile

import pysam
from tqdm import tqdm

from kleat.misc import apautils
from kleat.proc import process_suffix, process_bridge_and_link, process_blank
from kleat.misc import utils as U
from kleat.misc import settings as S

logger = logging.getLogger(__name__)


def gen_tmp_output(output, path=None):
    if path is None:
        path = os.path.dirname(output)
    return os.path.join(path, '__tmp_{0}.tsv'.format(os.path.basename(output)))


def prepare_args_for_collect_polya_evidence(num_cpus, output, c2g_bam_file, *args):
    c2g_bam = pysam.AlignmentFile(c2g_bam_file)

    args_list = []
    tmpdir = tempfile.gettempdir()
    for k, seqname in tqdm(enumerate(c2g_bam.references),
                           desc='prepared', unit=' seqnames'):

        tmp_output_file = gen_tmp_output(output, tmpdir)
        tmp_output_file += '.{0}'.format(seqname)
        U.backup_file(tmp_output_file)

        the_args = (seqname, tmp_output_file, c2g_bam_file) + args
        args_list.append(the_args)
    return args_list


def collect_polya_evidence(seqname, tmp_output_file, c2g_bam_file,
                           r2c_bam_file, ref_fa_file, bridge_skip_check_size):
    """loop through each contig and collect polyA evidence"""
    logging.info('collecting polyA evidence for {0} to {1} ...'.format(seqname, tmp_output_file))

    c2g_bam = pysam.AlignmentFile(c2g_bam_file)
    r2c_bam = pysam.AlignmentFile(r2c_bam_file)
    ref_fa = pysam.FastaFile(ref_fa_file)

    with open(tmp_output_file, 'wt') as opf:
        csvwriter = csv.writer(opf, delimiter='\t')
        csvwriter.writerow(S.HEADER)
        for contig in c2g_bam.fetch(seqname):
            if contig.is_unmapped:
                continue

            do_collection(contig, r2c_bam, ref_fa, csvwriter, bridge_skip_check_size)

    logging.info('collecting polyA evidence for {0} to {1} is done'.format(seqname, tmp_output_file))
    return tmp_output_file


def do_collection(contig, r2c_bam, ref_fa, csvwriter, bridge_skip_check_size):
    gen_key = apautils.gen_clv_key_tuple_from_clv_record

    ascs = []           # already supported clvs
    rec = process_suffix(
        contig, r2c_bam, ref_fa, csvwriter)
    if rec is not None:
        ascs.append(gen_key(rec))

    for rec in process_bridge_and_link(
            contig, r2c_bam, ref_fa, csvwriter, bridge_skip_check_size):
        # TODO: with either bridge or link, they probably won't support
        # clv of the other strand
        ascs.append(gen_key(rec))

    if not apautils.has_tail(contig):
        process_blank(contig, ref_fa, csvwriter, ascs)


def collect_polya_evidence_wrapper(args):
    return collect_polya_evidence(*args)
