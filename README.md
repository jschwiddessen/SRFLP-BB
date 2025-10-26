# SRFLP-BB
This is a C implementation of a parallel branch-and-bound solver for the single-row facility layout problem using semidefinite programming bounds.

The solver is exact and always finds the optimal solution when terminating. Moreover, the branch-and-bound implementation is parallelized using POSIX Threads.

This code builds upon the _BiqCrunch_ solver available at https://biqcrunch.lipn.univ-paris13.fr/Download. Moreover, the limited-memory quasi Newton code _L-BFGS-B_ https://users.iems.northwestern.edu/~nocedal/lbfgsb.html is used.

## The Single-Row Facility Layout Problem
Given $n$ one-dimensional facilities with positive integer lenghts $\ell_i$ and integer weights $c_{ij}$ for all pairs of facilities, we want to arrange all $n$ facilities on a straight line such that the total weighted sum of center-to-center distances between all pairs of facilities is minimized. Note that any such arrangement can be represented by a permutation $\pi \in \Pi_n$. Thus, the goal is to solve the problem
```math
\displaystyle \min_{\pi \in \Pi_n} \sum c_{ij} d_{ij}^\pi,
```
where $d_{ij}^\pi$ denotes the center-to-center distance of $i$ and $j$ with respect to the permutation $\pi$.

## SDP-based Bounds
The solver uses semidefinite programming (SDP) based bounds within the branch-and-bound approach. Check out my master's thesis https://jan-schwiddessen.com/assets/pdfs/master_thesis.pdf and the references therein. However, note that this repository does not contain the code to produce the results in the master's thesis. It will be published soon.

## Instances
See https://www.philipphungerlaender.com/benchmark-libraries/layout-lib/row-layout-instances/ for the instance format and to download many instances used in the literature.

## Installation

An installation of the Intel® oneAPI Math Kernel Library (oneMKL) is required. After correctly setting (the environment variable) MKLROOT, simply type:

```
mkdir obj
make
```

## Example
To solve an instance, do the following:
```
cd bin/
./srflp_bb example_15.txt
```
Here, "example_15.txt" is a text file that encodes the instance to be solved. The output will look similar to this:
```
Instance name: example_15.txt
Input weights are given as a symmetric matrix.

Number of facilities: 15
Lengths of facilities range: [1, 9]
Weights range: [0, 20]

Parameters:
Max. Threads: 8
Output frequency: 10
With cuts: 0
Gap cuts: -5.000000e-02
Scale factor triangle inequalities: 0.366025
Max. cuts: 1000000
Max. cuts per iteration: 5000
Min. cuts per iteration: 500
Memory corrections: 10
Max. alpha: 1.000000e+00
Min. alpha: 1.000000e-05
Scale alpha: 0.500000
Max. tol: 1.000000e-01
Min. tol: 1.000000e-02
Scale tol: 0.950000
Heuristic level: 3
Heuristic frequency: 3
Max. evals: 10000
Symmetry breaking: 2
Initial subproblems: 0
Branching: 1

Running heuristics ...
Starting with initial primalbound of 22343.5
Starting with initial dualbound of -100781.00

K=36244.0
=============================================================================================================================================
|      time      |     CPU time     |   explored   |  unexplored  |  evaluations  |      dualbound      |   primalbound   |       gap       |
=============================================================================================================================================
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         17550.5 |    44.52215990% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         17359.5 |    42.94934246% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         17346.5 |    42.84229206% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         17103.5 |    40.84127301% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         17092.5 |    40.75069190% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16924.5 |    39.36727132% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16831.5 |    38.60144921% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16577.5 |    36.50984905% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16574.5 |    36.48514511% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16526.5 |    36.08988209% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16518.5 |    36.02400492% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16483.5 |    35.73579230% |
|          0.09s |            0.09s |            0 |            1 |           102 |          12143.8124 |         16474.5 |    35.66168048% |
|          0.16s |            0.16s |            0 |            1 |           203 |          14429.1672 |         16461.5 |    14.08489315% |
|          0.23s |            0.23s |            0 |            1 |           304 |          15473.8333 |         16449.5 |     6.30526838% |
|          0.23s |            0.23s |            0 |            1 |           304 |          15473.8333 |         16448.5 |     6.29880585% |
|          0.29s |            0.29s |            0 |            1 |           405 |          15941.5770 |         16440.5 |     3.12969680% |
|          0.35s |            0.35s |            0 |            1 |           506 |          16154.6926 |         16439.5 |     1.76300113% |
|         10.00s |           42.22s |           15 |           16 |         75543 |          16415.3535 |         16439.5 |     0.14709734% |
|         17.44s |           98.14s |           67 |            0 |        147137 |          16439.5000 |         16439.5 |     0.00000000% |
=============================================================================================================================================

The instance has been solved to global optimality.
Optimal value: 16439.5

Optimal ordering:

2 14 13 12 5 10 1 6 9 11 3 7 4 8 15

Solution file has been written.
```
## Parameters
See the header file "constants_and_macros.h" for all configurable parameters, especially to set the maximum number of threads that should be used.
