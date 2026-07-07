# SOR Relaxation Factor Study

## Setup
All runs solve a 10×10 grid with the top edge held at 100 degrees and the other edges at 0 degrees, using tolerance 1e-6 and max_iter 10,000. Omega was varied from 1 to 1.9 in steps of 0.05, with a fresh grid constructed for each run.

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
- The minimum (29 sweeps) landed at ω = 1.50, consistent with the theoretical optimum of ~1.53 for this grid size.


- The omega = 1.0 is mathematically identical to Gauss-Seidel's method, matching it makes sure the experiment ran correctly.

  
- The valley is asymmetric: ω = 1.9 costs 173 sweeps while ω = 1.0 costs 124, despite 1.9 sitting closer to the optimum, so when tuning, aim low. It can be assumed that at 2 and beyond, the iteration stops converging entirely.