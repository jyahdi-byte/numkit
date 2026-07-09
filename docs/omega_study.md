# SOR Relaxation Factor Study

SOR speeds up Gauss-Seidel by over-relaxing each update by a factor
omega, and its convergence speed depends sharply on that factor. This
study measures the best omega experimentally and checks it against
theory.

## Setup
All runs solve a 10×10 grid with the top edge held at 100 degrees and
the other edges at 0, using tolerance 1e-6 and max_iter 10,000. Omega
runs from 1.0 to 1.9 in steps of 0.05, with a fresh grid constructed
for each run.

## Results
| omega | sweeps |
|-------|--------|
| 1.00  | 124    |
| 1.05  | 112    |
| 1.10  | 101    |
| 1.15  | 91     |
| 1.20  | 81     |
| 1.25  | 73     |
| 1.30  | 64     |
| 1.35  | 56     |
| 1.40  | 48     |
| 1.45  | 39     |
| **1.50** | **29** |
| 1.55  | 33     |
| 1.60  | 37     |
| 1.65  | 43     |
| 1.70  | 52     |
| 1.75  | 63     |
| 1.80  | 82     |
| 1.85  | 112    |
| 1.90  | 173    |

## Observations
The minimum, 29 sweeps, landed at omega = 1.50, consistent with the
theoretical optimum of about 1.53 for this grid size.

The omega = 1.0 row is a built-in check: SOR at omega = 1 is
mathematically identical to Gauss-Seidel, and it reproduces
Gauss-Seidel's independently measured 124 sweeps, so the experiment
is wired correctly.

The valley is asymmetric. Omega = 1.9 costs 173 sweeps while
omega = 1.0 costs 124, even though 1.9 sits closer to the optimum.
So when tuning, err low. At omega = 2 and beyond the iteration stops
converging entirely.
