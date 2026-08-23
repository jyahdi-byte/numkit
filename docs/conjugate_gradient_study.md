# Conjugate Gradient: Convergence Study

GS and SOR update one cell at a time and hope the whole grid settles.
CG instead builds up the answer using A-orthogonal search directions,
which for an n-unknown system reaches the exact answer in at most n
steps. This study checks whether that theoretical edge shows up as an
actual reduction in sweeps, and compares it against SOR at its
theoretical-optimal omega, the strongest classical baseline available,
not plain GS.

## Setup
Sweep counts below are produced by `tests/cg_sweep.cpp`. Same N x N grid as red_black_study.md: top edge held at 10, other edges
at 0, tol = 1e-10, max_iter = 20,000. N = 10, 20, 40, 80. SOR uses the
theoretical-optimum omega for each grid size, computed the same way
omega_study.md derived it. Correctness (cg_solve converges to the same
values as gauss_seidel_solve, on both a plain grid and one with a HOLE
region) is verified separately with real asserts in test_solvers.cpp;
this study only measures sweep counts.

## Results
| Grid Size | GS    | SOR | CG  |
|-----------|-------|-----|-----|
| 10        | 179   | 39  | 24  |
| 20        | 756   | 80  | 60  |
| 40        | 2977  | 161 | 124 |
| 80        | 11346 | 323 | 252 |

## Observations
Doubling N roughly quadruples GS's sweep count at every step
(756/179 = 4.22, 2977/756 = 3.94, 11346/2977 = 3.81), converging
toward the theoretical factor of 4 that O(N^2) predicts. SOR and CG
both roughly double instead (SOR: 2.05, 2.01, 2.01; CG: 2.50, 2.07,
2.03), matching O(N). SOR alone already fixes GS's scaling problem;
the comparison worth making is CG against SOR, not CG against GS.

CG wins at every size measured here (24 vs. 39, 60 vs. 80, 124 vs.
161, 252 vs. 323), without ever being told an omega. SOR only hits
those numbers because a closed-form formula for optimal omega exists
on a plain rectangular grid. That formula doesn't exist for an
irregular domain, like the HOLE regions this codebase already
supports. SOR would need omega found by trial and error there; CG
needs no tuning at all.

## Future work
Run this same comparison on a grid with a large HOLE region, where SOR
has no formula to fall back on and has to be tuned by search instead.

## Preconditioned CG

`pcg_solve` adds a diagonal (Jacobi) preconditioner: each cell's
residual is divided by its own diagonal value from `diag()` before
CG picks a search direction. On a plain grid every cell's diagonal
is 4, so this divides everything by the same constant and changes
nothing - preconditioning only has something to work with once
`HOLE` cells make the diagonal vary from cell to cell (`diag()`
subtracts 1 per HOLE neighbor, so it ranges 0-4). Correctness (pcg_solve
converges to the same values as cg_solve, on both a plain grid and one
with HOLE regions) is verified with real asserts in test_solvers.cpp.

## Setup
Sweep counts below are produced by `tests/pcg_sweep.cpp`. Same N x N
grid as the rest of this study: top edge held at 10, other edges at
0, tol = 1e-10, max_iter = 10,000. Two geometries: a single 3x3 HOLE
block placed at roughly the grid's center, and a checkerboard of 3x3
HOLE blocks tiled every 4 cells (so a 1-cell wall separates adjacent
blocks - thin enough that a wall cell touches HOLE on two sides,
dropping its diagonal to 2).

## Results
| Grid Size | CG (single hole) | PCG (single hole) | CG (checkerboard) | PCG (checkerboard) |
|-----------|-------------------|--------------------|--------------------|--------------------|
| 10        | -                 | -                  | 29                 | 26                 |
| 20        | 82                | 80                 | 52                 | 52                 |
| 40        | 164               | 164                | 137                | 135                |
| 80        | 316               | 317                | 266                | 270                |

## PCG observations
Neither geometry shows a real win. Single hole barely moves the
needle (82 vs 80, then a wash, then PCG a sweep *worse* at N=80).
The checkerboard was built specifically to force diagonal variation -
checked directly at N=40, 522 of 693 interior cells sit at diagonal 2
instead of 4 - and still lands in the same place: a small win at
N=10, a wash at N=20 and N=40, PCG four sweeps worse at N=80.

The reason traces back to what a diagonal preconditioner can actually
see. It corrects each cell's own scale, one cell at a time - it has
no information about the grid's overall shape. What makes CG slow on
the checkerboard isn't any single cell's diagonal; it's the maze-like
geometry itself, threading interior cells through a lattice of 1-cell
corridors. That's a global property of the domain, and a diagonal
preconditioner is a purely local tool. It was never going to reach a
problem like this one.

## PCG future work
Give cells their own coefficient (varying material properties instead
of just HOLE/INTERIOR) rather than the fixed 4-per-cell stencil this
solver assumes throughout. That's the kind of variation diagonal
preconditioning is actually built for, and would be a real test of
whether it earns its keep here - a bigger change than anything else in
this study, since apply_A, compute_b, and diag() would all need to
read a per-cell coefficient instead of assuming uniform 4.