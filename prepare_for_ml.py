import os
import sys
import logging

import pandas as pd
import kleat.misc.settings as S
from kleat.misc.utils import timeit

logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def gen_outfile(infile):
    base = '.'.join(os.path.basename(infile).split('.')[:-1])
    return os.path.join(os.path.dirname(infile),
                        '{0}.ml_ready.pkl'.format(base))


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = gen_outfile(infile)

    logging.info('reading {0}'.format(infile))
    adf = pd.read_csv(infile, keep_default_na=False, sep='\t')
    adf['abs_dist_to_aclv'] = adf['signed_dist_to_aclv']

    used_hexamers = [_[0] for _ in S.CANDIDATE_HEXAMERS] + ['NA']
    ctg_dum_hxm = pd.get_dummies(adf.ctg_hex, columns=used_hexamers)
    ctg_dum_hxm.columns = ['ctg_' + _ for _ in ctg_dum_hxm.columns.values]

    ref_dum_hxm = pd.get_dummies(adf.ref_hex, columns=used_hexamers)
    ref_dum_hxm.columns = ['ref_' + _ for _ in ref_dum_hxm.columns.values]

    bdf = pd.concat([adf, ctg_dum_hxm, ref_dum_hxm], axis=1)

    logging.info('writing {0}'.format(outfile))
    timeit(bdf.to_pickle)(outfile.replace('.fea', '.pkl'))
