"""
Utilities for postprocessing results from looping through each contig
individually and collected polyA evidence from them

- aggregate polyA evidence per clv (identified by (seqname, strand, clv) tuple)
- calculate the closest annotated clv for each clv
"""
import logging
import multiprocessing

import pandas as pd
import numpy as np
from tqdm import tqdm

from kleat.misc import utils as U
from kleat.misc import settings as S


logger = logging.getLogger(__name__)


def set_sort_join_strs(vals):
    return '|'.join(sorted(set(vals)))


def prepare_grps_for_agg(df_clv):
    clv_id_cols = ['seqname', 'strand', 'clv']
    clv_id_tuples, grps = [], []
    logger.info('preparing arguments for aggregating polya evidence in parallel...')
    iters = tqdm(df_clv.groupby(clv_id_cols),
                 desc='collected', unit=" cleavage sites")
    for clv_id_tuple, grp in iters:
        clv_id_tuples.append(clv_id_tuple)
        grps.append(grp)
    df_clv_ids = pd.DataFrame(clv_id_tuples, columns=clv_id_cols)
    return df_clv_ids, grps


def aggregate_polya_evidence(df_clv, num_cpus):
    df_clv_ids, grps = prepare_grps_for_agg(df_clv)
    with multiprocessing.Pool(num_cpus) as p:
        logger.info('aggregating (map operation) using {0} CPUs...'.format(num_cpus))
        res = U.timeit(p.map)(agg_polya_evidence_per, grps)
    df_res = pd.concat(res, axis=1).T
    ndf_res = pd.concat([df_clv_ids, df_res], axis=1)
    return ndf_res


def agg_polya_evidence_per(grp):
    sum_cols = grp[S.COLS_TO_SUM].sum()
    max_cols = grp[S.COLS_TO_MAX].max()
    any_cols = grp[S.COLS_TO_ANY].any()
    str_cols = grp[S.COLS_TO_JOIN].apply(set_sort_join_strs)
    # pick the strongest PAS hexamer
    hex_cols = grp[S.COLS_CONTIG_HEXAMERS].loc[grp.ctg_hex_id.idxmax()]
    one_cols = grp[S.COLS_PICK_ONE].iloc[0]
    return pd.concat([sum_cols, max_cols, any_cols,
                      str_cols, hex_cols, one_cols])


def calc_abs_dist_to_annot_clv(grp, annot_clvs):
    aclvs = annot_clvs.loc[grp.name]  # grp.name holds the group key
    bcast = np.broadcast_to(grp.clv.values, (aclvs.shape[0], grp.shape[0])).T

    sgn_dists = bcast - aclvs   # sgn: signed
    abs_dists = np.abs(sgn_dists)
    min_idxes = np.argmin(abs_dists, axis=1)

    nrows = abs_dists.shape[0]
    dists = sgn_dists[np.arange(nrows), min_idxes]

    grp['aclv'] = aclvs[min_idxes]
    grp['signed_dist_to_aclv'] = dists
    return grp


def add_abs_dist_to_annot_clv(df_clv, df_mapping):
    """
    add absolute distance to the closest annotated clv as an addition column

    :param df_mapping: clv-stop codon mapping in dataframe
    """
    # do some checking about which version of seqnames are used, use whatever
    # is used by df_clv as the reference
    mapping_seqname_is_ucsc = df_mapping.seqname.values[0] in S.UCSC_SEQNAMES
    cleavge_seqname_is_ucsc = df_clv.seqname.values[0] in S.UCSC_SEQNAMES

    if mapping_seqname_is_ucsc != cleavge_seqname_is_ucsc:
        if mapping_seqname_is_ucsc:
            df_mapping.seqname = df_mapping.seqname.replace(S.UCSC_TO_ENSEMBL_SEQNAME)
        else:
            df_mapping.seqname = df_mapping.seqname.replace(S.ENSEMBL_TO_UCSC_SEQNAME)

    annot_clvs = df_mapping.groupby(['seqname', 'strand']).apply(
        lambda g: g.clv.sort_values().values)

    # remove patch chromosomes
    if cleavge_seqname_is_ucsc:
        ndf_clv = df_clv.query('seqname in {0}'.format(S.UCSC_SEQNAMES))
    else:
        ndf_clv = df_clv.query('seqname in {0}'.format(S.ENSEMBL_SEQNAMES))

    logger.info('calculating absolute distances to annotated cleavage sites')
    timed = U.timeit(
        lambda _df: _df.groupby(['seqname', 'strand'])
        .apply(calc_abs_dist_to_annot_clv, annot_clvs=annot_clvs)
    )
    out = timed(ndf_clv)
    return out
