# Multithreaded Jacobi: Speedup Study

Jacobi's updates within a sweep are independent, so the sweep can be
split across threads. This study measures how much that actually
buys on real hardware, and what limits it.

## Setup
All runs solve a 200×200 grid with the top edge held at 100 degrees,
tol = 1e-6, using jacobi_mt_solve with thread counts 1, 2, 4, and 8
and a fresh grid per run. Hardware: 4 physical cores, 8 logical
processors.

## Results

| Threads | Time (s) | Speedup |
|---------|----------|---------|
| 1       | 152.0    | 1.00x   |
| 2       | 93.0     | 1.63x   |
| 4       | 110.0    | 1.38x   |
| 8       | 107.8    | 1.41x   |

All runs converged in 68,427 sweeps with probe values matching
single-threaded jacobi_solve to every digit. Correctness is verified;
only speed varies.

## Analysis
Speedup peaks at 2 threads and goes backwards at 4. The main limit is
memory bandwidth: Jacobi does little arithmetic per byte read, so the
shared memory bus saturates around 2 workers. The reversal at 4
threads comes from spawn overhead. Threads are re-created every sweep,
so that cost is paid 68,427 times and grows with thread count. Beyond
4 threads, hyperthreading means extra workers share physical cores and
add little. The serial grid copy was measured directly at 1.24e-5 s
per copy (about 0.85 s total, under 1% of runtime) and is negligible.

## Future work
A persistent worker pool using std::barrier would remove the per-sweep
respawn cost. Measuring the before/after difference is the natural
next experiment.
