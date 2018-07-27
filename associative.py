import sys
import itertools
import argparse
from gurobipy import *

from tree import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find a meet-semilattice structure on the set of phylogenetic trees. Equivalently, find a regular consensus method that is associative and Pareto on rooted triples.')
    
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
    
    
    print "Add meet constraints"
    
    # r^s <= r
    model.addConstrs((
        p[normalize_tuple((t,r))] >= m[t,r,s] \
        for (t,r,s) in possible_meets
    ), "meet1")
    
    # r^s <= s
    model.addConstrs((
        p[normalize_tuple((t,s))] >= m[t,r,s] \
        for (t,r,s) in possible_meets
    ), "meet2")
    
    # t = r^s, u <= r, u <= s imply u <= t
    # the triple (t,r,s) is in the list of possible meets
    model.addConstrs((
        m[t,r,s] + \
        p[normalize_tuple((u,r))] + \
        p[normalize_tuple((u,s))] - \
        (p[normalize_tuple((u,t))] if normalize_tuple((u,t)) in p else 0) \
        <= 2 
        for (t,r,s) in possible_meets for u in trees \
        if normalize_tuple((u,r)) in p and normalize_tuple((u,s)) in p
    ), "meet3")
    
    
    bottom_element = tuple(sorted(X))
    print "Force bottom element", bottom_element
    model.addConstrs((
        m[normalize_tuple((bottom_element, bottom_element, t))] == 1 \
        for t in trees if normalize_tuple((bottom_element, bottom_element, t)) in m
    ), "bottom")
    
    
    print "Force that every (normal) pair has a meet"
    matching = {(r,s): set() for (r,s) in normal_pairs}
    for (t,r,s) in possible_meets:
        matching[normalize_tuple((r,s))].add((t,r,s))
        if r != s:
            matching[normalize_tuple((s,r))].add((t,r,s))
    
    model.addConstrs((
        sum(m[triple] for triple in matching[pair]) == 1 for pair in matching
    ), "meetexists")
    
    
    ### solve ###
    kwargs = {
        "Threads": args.threads,
        "Presolve": 2, # aggressive presolve
    }
    
    for (arg, val) in kwargs.iteritems():
        model.setParam(arg, val)
    
    model.optimize()

    print
    if model.Status == GRB.INFEASIBLE:
        print "There is no meet-semilattice structure for X =", X

    else:
        print "Found meet-semilattice structure for X = %r (normal pairs only):" % X
        for (t,r), v in p.iteritems():
            if v.x > 0.5:
                print "%r <= %r" % (t,r)

