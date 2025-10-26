# SRFLP-BB
A parallel branch-and-bound solver for the single-row facility layout problem using semidefinite programming bounds

This code builds upon the _BiqCrunch_ solver available at https://biqcrunch.lipn.univ-paris13.fr/Download. Moreover, the limited-memory quasi Newton code _L-BFGS-B_ is used https://users.iems.northwestern.edu/~nocedal/lbfgsb.html.

## The Single-Row Facility Layout Problem
Given $n$ one-dimensional facilities with positive integer lenghts $\ell_i$ and integer weights $c_{ij}$ for all pairs of facilities, we want to arrange all $n$ facilities on a straight line such that the total weighted sum of center-to-center distances between all pairs of facilities is minimized. Note that any such arrangement can be represented by a permutation $\pi \in \Pi_n$. Thus, the goal is to solve the problem
```math
\displaystyle \min_{\pi \in \Pi_n} \sum c_{ij} d_{ij}^\pi,
```
where $d_{ij}^\pi$ denots the center-to-center distance of $i$ and $j$ with respect to the permutation $\pi$.
