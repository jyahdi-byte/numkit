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
Also worth trying: a diagonal (Jacobi) preconditioner on CG, which
should push its sweep count down further without changing correctness.