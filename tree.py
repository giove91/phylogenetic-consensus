"""
A phylogenetic tree on leaf set X is represented as a sorted tuple of:
- elements of X, or
- (recursively) tuples representing phylogenetic trees.
For example, the three rooted triples (phylogenetic binary trees on 3 leaves)
on X = {1,2,3} are: (1,(2,3)), (2,(1,3)), (3,(1,2)).
"""

import itertools


def powerset(X):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
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
        return
    
    for p in set_partitions(X):
        if len(p) == 1:
            # trivial partition
            continue
        
        for subtrees in itertools.product(*[all_trees(Y) for Y in p]):
            yield tuple(sorted(subtrees))


def restriction(t, Y):
    """
    Restriction of the tree t to the leaf set Y.
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


def apply_permutation(t, sigma):
    """
    Apply the permutation sigma on the leaves of t.
    """
    if not isinstance(t, tuple):
        # t is a leaf
        return sigma[t]
    
    return tuple(sorted(apply_permutation(s, sigma) for s in t))


if __name__ == '__main__':
    X = range(1, 5)
    Y = [1,2,3]
    sigma = {1: 2, 2: 3, 3: 4, 4: 1}
    
    print len(X), sum(1 for t in all_trees(X))
    
    for t in all_trees(X):
        print t, apply_permutation(t, sigma)


