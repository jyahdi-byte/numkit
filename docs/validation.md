# Validation Against Manufactured Solutions

## Method
To measure solver error against exact truth, we pick a function that
satisfies the Laplacian of T(x,y) = 0 for all x & y (i.e., a harmonic function), set the grid's boundary
cells to its values, and solve. Because a settled field is determined
entirely by the edges, the exact interior solution must equal the chosen
function, so solver output can be graded pointwise against the true values of T(x,y).
The reported error is the maximum of |computed - exact| over all interior points. All runs
use SOR with tol = 1e-8 and ω = 1.5.

## Experiment 1: T = x·y

| N   | Max error   |
|-----|-------------|
| 25  | 1.78e-07    |
| 50  | 7.93e-07    |
| 100 | 3.29e-06    |

A hand computation shows the four-neighbor average of exact x·y values
equals x·y exactly, so the discretization error for this problem
is 0. The measured errors are therefore pure iteration error, which
grows with N because as N grows, interior points are more cells from the boundary
so at cutoff they are less converged. 

## Experiment 2: T = sin(πx)·sinh(πy)

| N   | Max error   | Ratio vs previous |
|-----|-------------|-------------------|
| 25  | 5.71e-03    | —                 |
| 50  | 1.37e-03    | 4.17              |
| 100 | 3.33e-04    | 4.11              |

This function has nonzero fourth derivatives, so the finite-difference
approximation has genuine discretization error, predicted to shrink
as h^2. Halving h should divide the error by 4x; the measured
ratios of 4.17 and 4.11 confirm this, verifying second-order
convergence.

## Error budget
At N = 100, iteration error (≈3.3e-6, from Experiment 1) is roughly
1% of the total measured error, so the convergence measurement
is clean. At larger N the two error sources would cross, since
discretization error falls while iteration error increases; keeping the
measurement honest would require a lower tol as the grid refines.