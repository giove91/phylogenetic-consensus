import sys
import itertools
import argparse
from gurobipy import *

from tree import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find a meet-semilattice structure on the set of phylogenetic trees.')
    
    parser.add_argument('n', nargs='?', default=3, type=int, help='number of leaves')
    parser.add_argument('-t', '--threads', default=1, type=int, help='number of threads that can be used')
    args = parser.parse_args()


    n = args.n
    print "n =", n
    X = range(1, n+1)
    print "X =", X
    
    trees = list(all_trees(X))
    
    print "Compute normal forms of tuples of trees"
    normal_trees, normal_pairs, normal_triples = find_normal_forms(X, processes=args.threads)
    print
    
    print "There are %d trees, %d normal trees, %d normal pairs, and %d normal triples" % (len(trees), len(normal_trees), len(normal_pairs), len(normal_triples))
    
    print "Find possible meets"
    possible_meets = []
    for triple in normal_triples:
        t, r, s = triple
        
        if r == s and t != r:
            # (t,r,r) is possible only for t=r
            continue
        
        if t != r and normalize_tree(t) == normalize_tree(r):
            # t < r is not possible if t and r are in the same orbit
            continue
        
        if t != s and normalize_tree(t) == normalize_tree(s):
            # t < s is not possible if t and s are in the same orbit
            continue
        
        possible = True
        for Y in itertools.combinations(X, 3):
            u = restriction(r, Y)
            if len(u) < 3:  # u is not the bottom element
                if restriction(s, Y) == u and restriction(t, Y) != u:
                    possible = False
                    break
        
        if possible:
            possible_meets.append(triple)
    
    """
    # remove the possible meets (t,r,s) for which (t,t,r) and (t,t,s) are not possible meets
    possible_meets = set(possible_meets)
    new_possible_meets = []
    for (t,r,s) in possible_meets:
        if normalize_tuple((t,t,r)) in possible_meets and normalize_tuple((t,t,s)) in possible_meets:
            new_possible_meets.append((t,r,s))
        else:
            print "remove ", (t,r,s), "from possible meets!"
    possible_meets = new_possible_meets
    """
    
    print "There are %d possible meets" % len(possible_meets)
    
    # create optimization model
    model = Model('phylogenetictrees')
    
    print "Create variables"
    # m[t,r,s] == 1 means that t is the meet of r and s
    
    m = model.addVars(possible_meets, name="m", vtype=GRB.BINARY)
    
    # p[t,r] == 1 means that t <= r
    p = {}
    for t, r in normal_pairs:
        triple = normalize_tuple((t,t,r))
        if triple in m:
            p[t,r] = m[triple]
    
    
    print "Add poset constraints"
    print "* Reflexive"
    model.addConstrs((p[t,t] == 1 for t in normal_trees), "refl")
    
    print "* Antisymmetric"
    model.addConstrs((p[t,r] + p[normalize_tuple((r,t))] <= 1 for (t,r) in p if normalize_tuple((r,t)) in p and t != r), "antisym")
    
    print "* Transitive"
    model.addConstrs((
        (p[normalize_tuple((t,s))] if normalize_tuple((t,s)) in p else 0) \
        >= p[t,r] + p[normalize_tuple((r,s))] - 1 \
        for (t,r) in p for s in trees \
        if normalize_tuple((r,s)) in p
    ), "trans")

