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
