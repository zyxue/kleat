import os
import logging

from settings import BAM_CMATCH

logger = logging.getLogger(__name__)


def gen_clv_key_tuple(seqname, strand, clv):
    return (seqname, strand, clv)


def gen_clv_key_str(seqname, strand, clv):
    return f'{seqname}|{strand}|{clv}'


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


# Below are very generic functions not related to APA

def backup_file(f):
    """
    Back up a file, old_file will be renamed to #old_file.n#, where n is a
    number incremented each time a backup takes place
    """
    if os.path.exists(f):
        dirname = os.path.dirname(f)
        basename = os.path.basename(f)
        count = 1
        rn_to = os.path.join(
            dirname, '#' + basename + '.{0}#'.format(count))
        while os.path.exists(rn_to):
            count += 1
            rn_to = os.path.join(
                dirname, '#' + basename + '.{0}#'.format(count))
        logger.info("Backing up {0} to {1}".format(f, rn_to))
        os.rename(f, rn_to)
        return rn_to
    else:
        logger.warning('{0} doesn\'t exist'.format(f))

