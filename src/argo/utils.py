import numpy as np

from scipy.sparse import dok_matrix
from sklearn.preprocessing import normalize


def filter_overlap(file, max_hang=1000, int_frac=0.05):
    '''
    Filter out dangling overlaps of a paf file.
    Reference: https://github.com/lh3/miniasm/blob/master/miniasm.h#L86

    Parmeters:
        max_hang: int
            Max overhang length o.
        int_frac: float
            Max overhang to mapping length ratio r.

    Returns:
        overlaps: list
            Filtered overlapping reads with DV.
    '''
    overlaps = []
    with open(file) as f:
        for line in f:
            ls = line.rstrip().split('\t')
            qseqid, sseqid = ls[0], ls[5]
            qlen, slen = int(ls[1]), int(ls[6])

            qstart, qend = int(ls[2]), int(ls[3])
            sstart, send = int(ls[7]), int(ls[8])
            DV = float(ls[15].split('dv:f:')[-1])

            if ls[4] == '-':
                sstart, send = slen - send, slen - sstart

            # filter out dangling overlaps
            overhang = min(qstart, sstart) + min(qlen - qend, slen - send)
            maplen = max(qend - qstart, send - sstart)

            if overhang <= min(max_hang, maplen * int_frac):
                overlaps.append([qseqid, sseqid, DV])

    return overlaps


def set_cover(elements, subsets, scores):
    '''
    Greedy set cover.

    Parameters:
        elements: set
            Reads to be covered
        subsets: defaultdict(set)
            Species-read sets.
        scores: defaultdict(lambda: defaultdict(dict))
            Alignment scores (AS) of each species-read pair.

    Returns:
        covers: list
            Species that cover all reads.
    '''
    covers = []
    keys = list(subsets.keys())
    while elements:
        min_cost, min_index = np.inf, None
        for key in keys:
            if not (subset := subsets[key].intersection(elements)):
                continue

            cost = sum(-score for id, score in scores[key].items() if id in subset)
            if cost < min_cost:
                min_cost, min_index = cost, key

        if min_index is None:
            break

        elements -= subsets[min_index]
        keys.remove(min_index)
        covers.append(min_index)

    return covers


def sparse_allclose(matrix, matrix_hist, rtol=1e-5, atol=1e-8):
    '''
    Sparse version of np.allclose.
    Reference: https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
    '''
    return (np.abs(matrix - matrix_hist) - rtol * np.abs(matrix_hist)).max() <= atol


def sparse_normalize(matrix):
    '''
    Normalize a sparse matrix with respect to l1 norm, assuming no negative values.
    '''
    return normalize(matrix, norm='l1', axis=0)


def iterate(matrix, inflation, expansion, pruning_threshold):
    '''
    Perform expand + inflate + prune combo.
    '''
    matrix = sparse_normalize((matrix ** expansion).power(inflation))

    ## prune entries with small values
    matrix_pruned = dok_matrix(matrix.shape)
    matrix_pruned[matrix >= pruning_threshold] = matrix[matrix >= pruning_threshold]

    row_indices = matrix.argmax(axis=0).reshape((matrix.shape[0], ))
    col_indices = np.arange(matrix.shape[0])
    matrix_pruned[row_indices, col_indices] = matrix[row_indices, col_indices]

    return matrix_pruned.tocsc()


def mcl(matrix, max_iterations=1000, inflation=2, expansion=2):
    '''
    Run MCL to get clusters of reads, based on their alignment identities.
    Code modified from: https://github.com/GuyAllard/markov_clustering/blob/master/markov_clustering/mcl.py#L167

    Parmeters:
        max_iterations: int
            Max. number of iterations.
        inflation: float
            Cluster inflation factor.
        expansion: float
            Cluster expansion factor.

    Returns:
        clusters: list
            List of clusters
    '''
    ## Add self-loops to the matrix
    for i in range(matrix.shape[0]):
        matrix[i, i] = 1

    matrix = sparse_normalize(matrix)
    for i in range(max_iterations):
        matrix_hist = matrix.copy()

        ## expansion + inflation + prune
        matrix = iterate(matrix, inflation, expansion, 1 / matrix.shape[0])

        ## check for convergence
        if sparse_allclose(matrix, matrix_hist):
            break

    attractors = matrix.diagonal().nonzero()[0]
    clusters = set()

    for attractor in attractors:
        cluster = tuple(matrix.getrow(attractor).nonzero()[1].tolist())
        clusters.add(cluster)

    return sorted(list(clusters))
