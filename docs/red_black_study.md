# Red-Black Ordering: Correctness and Convergence Study

Gauss-Seidel and SOR update in place, so naively threading them the way
Jacobi is threaded creates a data race: one thread can read a neighbor
cell while another thread is mid-write to it. Red-black ordering splits
each sweep into two passes over alternating checkerboard colors, so
within a pass no cell being written is ever read by another update in
that same pass. This study checks that the reordering doesn't change
what the solvers converge to, and measures what it costs (or saves) in
sweep count before any threading work is done.

## Setup
All runs solve an N×N grid with the top edge held at 10 degrees and the
other edges at 0, tol = 1e-10, max_iter = 20,000, with a fresh grid
constructed per run. N = 10, 20, 40, 80. SOR and SOR-RB both use the
theoretical-optimum omega for each grid size. Correctness (that
gauss_seidel_rb_solve and sor_rb_solve converge to the same values as
their originals, including that HOLE cells are never updated) is
verified separately with real asserts in test_solvers.cpp; this study
only measures sweep counts.

## Results
| Grid Size | GS    | GS-RB | SOR | SOR-RB |
|-----------|-------|-------|-----|--------|
| 10        | 179   | 181   | 39  | 40     |
| 20        | 756   | 760   | 80  | 82     |
| 40        | 2977  | 2986  | 161 | 162    |
| 80        | 11346 | 11365 | 323 | 321    |

## Observations
Red-black costs almost nothing in sweep count. GS-RB runs 1-2% more
sweeps than plain GS at every size, and SOR-RB is within a sweep or two
of SOR in both directions (it even finishes one sweep faster at N=80).
That's expected: both orderings solve the same fixed point, so the gap
is just visiting order changing the exact convergence path length, not
a difference in the underlying method.

Sweep counts for both GS variants grow close to O(N^2) as N doubles
(756/179 = 4.22, 2977/756 = 3.94, 11346/2977 = 3.81, approaching the
theoretical factor of 4), matching the expected growth as grid spacing
h halves.

Sweep count is essentially a wash. The actual value of red-black is
that each pass is now race-free, so it can be split across threads -
something plain GS/SOR can't do safely.

## Future work
Thread gauss_seidel_rb_solve and sor_rb_solve the way jacobi_mt_solve
threads Jacobi: split each color's pass across a worker pool, syncing
twice per iteration instead of once. Add a CUDA red-black kernel too.
Measuring wall-clock speedup against single-threaded red-black is the
natural next experiment.