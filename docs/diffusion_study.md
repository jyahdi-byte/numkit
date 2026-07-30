# Numerical Diffusion vs. Courant Number Study

First-order upwind is exact at Courant number 1 — the update
collapses to a pure copy from the upstream neighbor, no averaging at
all. At any Courant number below 1, each step blends the upstream
neighbor's value in, which is unavoidably a weighted average of two
points. Averaging smooths out whatever's different between them, so
away from cn = 1 the scheme blurs a sharp shape instead of translating
it exactly. This study measures how much.

## Setup
50-point periodic grid, dx = 0.02, initial condition a step block at
points 5-8 (value 1, everywhere else 0), 20 timesteps. Courant number
swept from 0.1 to 0.9 in steps of 0.2. Each run is compared against
the exact solution — the same initial block run through the solver at
cn = 1 for the same number of steps, which is exact by construction —
using RMSE (root-mean-square error) between the run's final array and
that exact array.

## Results
| Courant number | RMSE |
|-----------------|------|
| 0.1              | 0.362069 |
| 0.3              | 0.344084 |
| 0.5              | 0.340257 |
| 0.7              | 0.335192 |
| 0.9              | 0.220912 |

[cn = 0.1 space-time diagram — most blurred]
![Courant 0.1](diffusion_cn_0.1.png)

[cn = 0.9 space-time diagram — closest to exact]
![Courant 0.9](diffusion_cn_0.9.png)

## Observations

[TODO — write this part yourself, in your own words. What to cover,
from what you already derived:

1. The overall trend: RMSE falls as Courant number rises toward 1.

2. The per-step mechanism: in the update
   u_i^(n+1) = (1-cn)*u_i^n + cn*u_{i-1}^n, the two coefficients are
   (1-cn) and cn. The gap between them, |2cn - 1|, controls how close
   a single step is to a pure copy versus a 50/50 blend. Compute that
   gap at each of your five Courant numbers and note that it changes
   at a perfectly constant rate (0.4 per 0.2 step in cn) — the
   per-step mechanism alone is linear.

3. Why the measured RMSE isn't linear even though the per-step
   mechanism is: this is 20 repeated applications of that same
   weighted average, not one. A per-step difference that's small
   between adjacent Courant numbers compounds over 20 repetitions,
   which is why the RMSE drop is roughly flat from 0.1 through 0.7
   and then falls sharply between 0.7 and 0.9, close to cn = 1.

4. One sentence on why a curve fit (e.g. a cubic polynomial through
   the five points) isn't used here: it would describe the shape of
   these five numbers without explaining anything about the actual
   compounding mechanism, and wouldn't generalize to a different
   step count or grid size.]