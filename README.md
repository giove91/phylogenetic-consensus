# Impossibility results on stability of phylogenetic consensus methods

## Requirements

- Python 2.7
- [Gurobi](https://www.gurobi.com) >= 7.5 ([free academic licences](https://www.gurobi.com/academia/for-universities) are available)
- [Gurobipy](https://www.gurobi.com/documentation/8.0/quickstart_linux/the_gurobi_python_interfac.html)

## Usage


### Extension stability on binary trees

`python binary.py [-h] [-t THREADS] [n]`

Find a regular consensus method that satisfies extension stability on profiles
of two binary trees.

Positional arguments:
```
  n                     number of leaves (between 3 and 5, default 3)
```

Optional arguments:
```
  -h, --help            show help message and exit
  -t THREADS, --threads THREADS
                        number of threads that can be used
```

### Associativity

`python associative.py [-h] [-t THREADS] [n]`

Find a regular consensus method that is associative and Pareto on rooted triples.

Positional arguments:
```
  n                     number of leaves (between 3 and 5, default 3)
````

Optional arguments:
```
  -h, --help            show help message and exit
  -t THREADS, --threads THREADS
                        number of threads that can be used
```
