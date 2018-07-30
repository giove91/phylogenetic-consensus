"""
A phylogenetic tree on leaf set X is represented as a sorted tuple of:
- elements of X, or
- (recursively) tuples representing phylogenetic trees.
For example, the three rooted triples (phylogenetic binary trees on 3 leaves)
on X = {1,2,3} are: (1,(2,3)), (2,(1,3)), (3,(1,2)).
"""

import sys
import itertools
from multiprocessing import Pool


def powerset(X):
    "Generate all subsets of X."
    return itertools.chain.from_iterable(itertools.combinations(X, r) for r in range(len(X)+1))


def set_partitions(X):
    """
    Generate all set partitions of X.
    """
    X = list(X)
    if len(X) == 1:
        yield [X]
        return
    
    x = X[0]
    for Y in powerset(X[1:]):
        Y_set = set(Y)
        Z = [z for z in X[1:] if z not in Y_set]
        
        if len(Z) == 0:
            yield [X]
        
        else:
            for p in set_partitions(Z):
                yield [[x] + list(Y)] + p


def all_trees(X):
    """
    Generate all phylogenetic trees on leaf set X.
    """
    if len(X) == 1:
        yield X[0]
    
    else:
        for p in set_partitions(X):
            if len(p) == 1:
                # trivial partition
                continue
            
            for subtrees in itertools.product(*[all_trees(Y) for Y in p]):
                yield tuple(sorted(subtrees))


def leaf_set(t):
    """
    Find leaf set of a phylogenetic tree (as a list).
    """
    if not isinstance(t, tuple):
        # t is a leaf
        return [t]
    
    res = []
    for s in t:
        res += leaf_set(s)
    return sorted(res)


CLUSTERS = {}
def clusters(t):
    """
    Find the clusters of a phylogenetic tree (as a list of tuples).
    """
    if t not in CLUSTERS:
        if not isinstance(t, tuple):
            # t is a leaf
            CLUSTERS[t] = [(t,)]
        
        else:
            res = []
            for s in t:
                res += clusters(s)
            res.append(tuple(leaf_set(t)))
            CLUSTERS[t] = res
    
    return CLUSTERS[t]


def compare(t, s):
    """
    Check if t <= s, in the sense that all the clusters of t are also clusters of s.
    """
    c = set(clusters(s))
    return all(cluster in c for cluster in clusters(t))


def restriction(t, Y):
    """
    Restriction of the tree t to the leaf set Y.
    Return None if Y is disjoint from the leaf set of t.
    """
    if not isinstance(t, tuple):
        # t is a leaf
        return t if t in Y else None
    
    subtrees = [restriction(s, Y) for s in t]
    subtrees = [s for s in subtrees if s is not None]
    
    if len(subtrees) == 0:
        return None
    elif len(subtrees) == 1:
        return subtrees[0]
    else:
        return tuple(sorted(subtrees))


def is_binary(t):
    """
    Check if the tree is binary.
    """
    if not isinstance(t, tuple):
        # t is a leaf
        return True
    
    return len(t) == 2 and all(is_binary(s) for s in t)


def all_permutations(X):
    """
    Generate all permutations of X, as dictionaries.
    For example, the permutation (1234) is given by {1: 2, 2: 3, 3: 4, 4: 1}.
    """
    X = list(X)
    for permutation in itertools.permutations(X):
        yield {X[i]: permutation[i] for i in xrange(len(X))}


def apply_permutation(t, sigma):
    """
    Apply the permutation sigma on the leaves of t.
    """
    if not isinstance(t, tuple):
        # t is a leaf
        return sigma[t]
    
    return tuple(sorted(apply_permutation(s, sigma) for s in t))


NORMAL_TREES = {}
def normalize_tree(t):
    """
    Find normal form of the given tree w.r.t. the symmetric group action.
    The normal form is the smallest element of the orbit.
    """
    if t not in NORMAL_TREES:
        # compute normal form
        X = leaf_set(t)
        normal_form = t
        
        for sigma in all_permutations(X):
            normal_form = min(normal_form, apply_permutation(t, sigma))
        
        for sigma in all_permutations(X):
            NORMAL_TREES[apply_permutation(t, sigma)] = normal_form
    
    return NORMAL_TREES[t]


NORMAL_TUPLES = {}
def normalize_tuple(tup):
    """
    Find normal form of the given tuple of trees w.r.t. the diagonal symmetric group action.
    For triples, we also allow exchanging the second and the third argument.
    The leaf set is changed to [1, ..., n]
    """
    X = leaf_set(tup[0])
    Y = range(1, len(X)+1)
    if X != Y:
        # normalize range
        sigma = {X[i]: Y[i] for i in xrange(len(X))}
        tup = tuple(apply_permutation(t, sigma) for t in tup)
        X = Y
    
    if tup not in NORMAL_TUPLES:
        # compute normal form
        normal_form = tup
        orbit = [tuple(apply_permutation(t, sigma) for t in tup) for sigma in all_permutations(X)]
        
        if len(tup) == 3:
            orbit += [(a,c,b) for (a,b,c) in orbit]
        
        normal_form = min(orbit)
        
        for tup2 in orbit:
            NORMAL_TUPLES[tup2] = normal_form
    
    return NORMAL_TUPLES[tup]


def find_normal_tuples(t):
    """
    Find normal forms for all pairs and triples that begin with t.
    """
    normal_pairs = set()
    normal_triples = set()
    
    X = leaf_set(t)
    for r in all_trees(X):
        normal_pairs.add(normalize_tuple((t,r)))
        
        for s in all_trees(X):
            normal_triples.add(normalize_tuple((t,r,s)))
    
    sys.stdout.write(".")
    return normal_pairs, normal_triples


def find_normal_forms(X, processes=1):
    """
    Find normal forms for all trees, pairs and triples.
    Return sorted lists of normal forms.
    """
    normal_trees = set()
    normal_pairs = set()
    normal_triples = set()
    
    for t in all_trees(X):
        normal_trees.add(normalize_tree(t))
    
    if processes == 1:
        results = [find_normal_tuples(t) for t in normal_trees]
    else:
        process_pool = Pool(processes=processes)
        results = process_pool.map_async(find_normal_tuples, sorted(normal_trees, key=lambda t: len(t))).get(9999999)
    
    for pairs, triples in results:
        normal_pairs |= pairs
        normal_triples |= triples
    
    return list(sorted(normal_trees)), list(sorted(normal_pairs)), list(sorted(normal_triples))


