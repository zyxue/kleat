import os
import logging
import time
from functools import update_wrapper

logger = logging.getLogger(__name__)


"""Below are very generic functions not related to APA"""


def decorator(d):
    "Make function d a decorator: d wraps a function fn."
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d


@decorator
def timeit(f):
    """time a function, used as decorator"""
    def new_f(*args, **kwargs):
        bt = time.time()
        r = f(*args, **kwargs)
        et = time.time()
        logger.info("time spent on {0}: {1:.2f}s".format(f.__name__, et - bt))
        return r
    return new_f


def backup_one_file(f):
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


def backup_file(*files):
    """
    Back up every file in files if it exists, old_file will be renamed to
    #old_file.n#, where n is a number incremented each time a backup takes
    place
    """
    for f in files:
        backup_one_file(f)
