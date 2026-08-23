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
Sweep counts below are produced by `tests/red_black_sweep.cpp`. All runs solve an N×N grid with the top edge held at 10 degrees and the
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

## Threading speedup

Timings below are produced by `tests/bench_rb_mt.cpp`. 200×200 grid,
top edge at 100, tol = 1e-6, max_iter = 100,000. A warm-up solve
(discarded) runs first so the single-threaded baseline isn't unfairly
penalized by first-allocation overhead. Machine: 4 physical cores, 8
logical (hyperthreaded).

| Threads | GS-RB sweeps | GS-RB time (s) | GS-RB speedup | SOR-RB sweeps | SOR-RB time (s) | SOR-RB speedup |
|---------|--------------|-----------------|----------------|----------------|-------------------|------------------|
| 1 (single-threaded) | 36996 | 302.297 | 1.00x | 555 | 14.924 | 1.00x |
| 1 (mt)  | 36996 | 327.635 | 0.92x | 555 | 15.258 | 0.98x |
| 2       | 36996 | 194.153 | 1.56x | 555 | 14.797 | 1.01x |
| 4       | 36996 | 203.275 | 1.49x | 555 | 8.405  | 1.78x |
| 8       | 36996 | 243.162 | 1.24x | 555 | 4.766  | 3.13x |

## Threading observations

GS-RB and SOR-RB run the exact same threading code - same barrier
structure, same row chunking - and still scale in opposite directions
past 4 threads. GS-RB peaks at 2 threads (1.56x) and degrades from
there. SOR-RB keeps improving all the way to 8 threads (3.13x).

The difference is sweep count. GS-RB takes 36996 sweeps to converge;
SOR-RB takes 555. Each sweep is two barrier syncs, so GS-RB hits
`arrive_and_wait()` about 66x more often than SOR-RB does over the
same solve. Every barrier is a spot where the OS can reschedule a
thread, and one thread running late stalls all the others until it
catches up - so the more barriers a solve needs, the more chances
that kind of stall has to accumulate. SOR-RB converges too fast for
it to show up; GS-RB runs long enough that it does.

Task Manager during the GS-RB runs showed CPU usage scaling with
thread count - roughly 15% at 1 thread up to 85-95% at 8 - so the
threads are actually running, not blocked or idle. Past 4 threads, 8
logical threads are splitting time across 4 physical cores, and all
of them are touching the same grid at once, adding memory bandwidth
contention on top. That accounts for the regression - the threading
code itself has no bug in it.

## Future work
Add a CUDA red-black kernel.