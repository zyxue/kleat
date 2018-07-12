import numpy as np

from scipy.cluster.hierarchy import fclusterdata


def gen_cluster_index(x, window):
    return fclusterdata(
        x.reshape(-1, 1),
        window,
        criterion='distance',
        metric='euclidean',
        method='single'
    )


def cluster(x, window):
    x = np.sort(x)
    # used cluster index as id
    cluster_id = gen_cluster_index(x, window)
    return dict(zip(x, cluster_id))


def select_rep(df):
    """select a representative clv per cluster"""
    mode_clv = df.clv.mode()
    if mode_clv.shape[0] > 0:
        # if one or multiple modes are found, return median
        return mode_clv.median()
    else:
        # if no mode is found, then return median
        return df.clv.median()


def cluster_clv_sites(df, window):
    df = df.copy()
    uniq_clvs = df.clv.unique()
    if uniq_clvs.shape[0] > 1:       # other exception would be raised
        cluster_clv2id_dd = cluster(uniq_clvs, window=window)

        df['cluster_id'] = df.clv.replace(cluster_clv2id_dd)
        # select a representative clv for each cluster in a new column
        cluster_id2rep_dd = df.groupby('cluster_id').apply(select_rep).to_dict()
        df['mode_clv'] = df['cluster_id'].replace(cluster_id2rep_dd).astype(np.int64)
    else:
        df['cluster_id'] = 1   # as its own cluster
        df['mode_clv'] = uniq_clvs[0]
    return df


def _test(x, window):
    # in Python3 dict.values() returns a view instead of a list
    # https://docs.python.org/3/library/stdtypes.html#dict-views
    return np.unique(list(cluster(np.array(x), window).values())).tolist()


# some test for single-linkage clustering
assert _test([1, 2, 3], 20) == [1]
assert _test([1, 2, 3, 4], 20) == [1]
assert _test([1, 2, 3, 30, 45], 20) == [1, 2]
assert _test([10, 30, 50], 20) == [1]
assert _test([10, 31, 50], 20) == [1, 2]
assert _test([10, 30, 50, 69], 20) == [1]
assert _test([10, 30, 50, 70], 20) == [1]
assert _test([10, 30, 50, 71], 20) == [1, 2]
assert _test([10, 30, 49, 69], 20) == [1]
assert _test([1, 10, 24, 28], 20) == [1]
assert _test([25, 26, 30, 34, 35], 20) == [1]
assert _test([25, 26, 30, 34, 35, 36, 38, 39, 43, 45, 48, 53, 55, 57, 61, 62, 63,
       66, 67, 68, 73, 74, 80], 20) == [1]
# # make sure sorting and median works
assert _test([10, 20, 30, 40, 50], 15) == [1]
assert _test([50, 10, 30, 20, 40], 15) == [1]
