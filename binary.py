import sys
import itertools
import argparse
from gurobipy import *

from tree import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find a regular consensus method that satisfies extension stability on profiles of two binary trees.')
    
    parser.add_argument('n', nargs='?', default=3, type=int, help='number of leaves')
    parser.add_argument('-t', '--threads', default=1, type=int, help='number of threads that can be used')
    args = parser.parse_args()


    n = args.n
    print "n =", n
    X = range(1, n+1)
    print "X =", X
    
    trees = []
    normal_trees = []
    normal_pairs = []
    normal_triples = []
    
    print "Compute normal forms of tuples of trees up to %d leaves" % n
    for i in xrange(3, n+1):
        Y = range(1, i+1)
        
        trees += list(all_trees(Y))
        
        normal_trees_Y, normal_pairs_Y, normal_triples_Y = find_normal_forms(Y, processes=args.threads)
        
        normal_trees += normal_trees_Y
        normal_pairs += normal_pairs_Y
        normal_triples += normal_triples_Y
    
    print
    
    print "There are %d trees, %d normal trees, %d normal pairs, and %d normal triples" % (len(trees), len(normal_trees), len(normal_pairs), len(normal_triples))
    
    
    print "Find possible consensus triples between binary trees"
    possible_triples = []
    for triple in normal_triples:
        t, r, s = triple
        Y = leaf_set(t)
        
        if r == s and t != r:
            # (t,r,r) is possible only for t=r (unanimity)
            continue
        
        if not is_binary(r) or not is_binary(s):
            # we are only interested in pairs of binary trees
            continue
        
        # check Pareto property on rooted triples (for binary trees)
        possible = True
        for Z in itertools.combinations(Y, 3):
            u = restriction(r, Z)
            if len(u) < 3:  # u is not the bottom element
                if restriction(s, Z) == u and restriction(t, Z) != u:
                    possible = False
                    break
        
        if possible:
            possible_triples.append(triple)
    
    
    print "There are %d possible triples" % len(possible_triples)
    
    ### create optimization model ###
    model = Model('phylogenetictrees')
    
    print "Create variables"
    # m[t,r,s] == 1 means that t is the consensus tree of r and s
    
    m = model.addVars(possible_triples, name="m", vtype=GRB.BINARY)
    
    print "Add extension stability constraints on binary trees"
    extension_stability_constraints = []
    for (t, r, s) in possible_triples:
        Y = leaf_set(t)
        for Z in powerset(Y):
            if len(Z) > 3 and len(Z) < len(Y):
                for u in all_trees(Z):
                    if not compare(u, restriction(t, Z)):
                        second_triple = normalize_tuple((u, restriction(r, Z), restriction(s, Z)))
                        if second_triple in m:
                            extension_stability_constraints.append(((t,r,s), second_triple))
    
    model.addConstrs((m[first_triple] + m[second_triple] <= 1 for first_triple, second_triple in extension_stability_constraints), "extstab")
    
    print "Force that every (normal) pair has exactly one consensus tree"
    matching = {(r,s): set() for (r,s) in normal_pairs if is_binary(r) and is_binary(s)}
    for (t,r,s) in possible_triples:
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
        print "There is no valid consensus method for X =", X

    else:
        print "Found consensus method for X = %r:" % X
        print
        for (t,r,s), v in m.iteritems():
            if v.x > 0.5:
                print "%r ^ %r = %r" % (r,s,t)
    


