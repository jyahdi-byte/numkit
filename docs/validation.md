# Validation Against Manufactured Solutions

The solvers agree with each other, but solvers sharing a bug would
too. This study grades them against exact answers, and checks that
the error shrinks with grid refinement at the rate theory predicts.

## Method
The trick is to manufacture a problem whose answer is known in
advance. Pick a function T(x,y) whose Laplacian is zero everywhere
(a harmonic function), and set the grid's boundary cells to its
values. Since a settled field is determined entirely by its edges,
the exact interior solution must equal the chosen function, so the
solver's output can be graded pointwise against T. The reported
error is the maximum of |computed - exact| over all interior points.
All runs use SOR with tol = 1e-8 and omega = 1.5.

## Experiment 1: T = x·y

| N   | Max error   |
|-----|-------------|
| 25  | 1.78e-07    |
| 50  | 7.93e-07    |
| 100 | 3.29e-06    |

A hand computation shows the four-neighbor average of exact x·y
values equals x·y exactly, so the discretization error for this
problem is zero. The measured errors are therefore pure iteration
error. That error grows with N because interior points sit more
cells from the boundary, so they are less converged when the
tolerance cutoff fires.

## Experiment 2: T = sin(πx)·sinh(πy)

| N   | Max error   | Ratio vs previous |
|-----|-------------|-------------------|
| 25  | 5.71e-03    | —                 |
| 50  | 1.37e-03    | 4.17              |
| 100 | 3.33e-04    | 4.11              |

This function has nonzero fourth derivatives, so the finite-difference
approximation has genuine discretization error, predicted to shrink as
h². Halving h should divide the error by 4. The measured ratios of
4.17 and 4.11 confirm this, verifying second-order convergence.

## Error budget
At N = 100, iteration error (about 3.3e-6, from Experiment 1) is
roughly 1% of the total measured error, so the convergence measurement
is clean. At larger N the two error sources would cross, since
discretization error falls while iteration error grows. Keeping the
measurement honest would require tightening tol as the grid refines.
