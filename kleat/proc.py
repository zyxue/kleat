from kleat.evidence import suffix, bridge, link, blank
from kleat.misc import apautils


def process_suffix(contig, r2c_bam, ref_fa, csvwriter):
    tail_side = apautils.has_tail(contig)
    if tail_side is not None:
        clv_record = suffix.gen_clv_record(contig, r2c_bam, tail_side, ref_fa)
        apautils.write_row(clv_record, csvwriter)
        return clv_record


def process_bridge_and_link(contig, r2c_bam, ref_fa, csvwriter):
    # bridge & link
    aligned_reads = r2c_bam.fetch(contig.query_name)
    dd_bridge, dd_link = extract_bridge_and_link(contig, aligned_reads, ref_fa)
    bdg_clvs, lnk_clvs = [], []
    if len(dd_bridge['num_reads']) > 0:
        bdg_clvs.extend(
            bridge.write_evidence(dd_bridge, contig, ref_fa, csvwriter))
    if len(dd_link['num_reads']) > 0:
        lnk_clvs.extend(
            link.write_evidence(dd_link, contig, ref_fa, csvwriter))
    return bdg_clvs + lnk_clvs


def process_blank(contig, ref_fa, csvwriter, asp_clv_keys):
    """:param asp_clv_keys: list of already supported clv key tuples"""
    asp_clv_keys = set(asp_clv_keys)
    for clv_rec in blank.gen_two_clv_records(contig, ref_fa, asp_clv_keys):
        apautils.write_row(clv_rec, csvwriter)


def extract_bridge_and_link(contig, aligned_reads, ref_fa):
    """
    bridge and link are processed together by looping throught reads aligned to
    the contig

    :param contig: a pysam.libcalignedsegment.AlignedSegment instance
    :param aligned_reads: a pysam.libcalignmentfile.IteratorRowRegion instance
    """
    dd_bridge = bridge.init_evidence_holder()
    dd_link = link.init_evidence_holder()
    for read in aligned_reads:
        if bridge.is_a_bridge_read(read):
            bdg_evid = bridge.analyze_bridge(contig, read, ref_fa, dd_bridge)
            bridge.update_evidence(bdg_evid, dd_bridge)
        elif link.is_a_link_read(read):
            link_evid = link.analyze_link(contig, read)
            link.update_evidence(link_evid, dd_link)
    return dd_bridge, dd_link
