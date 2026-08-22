# Numerical Diffusion vs. Courant Number Study

First-order upwind is exact at Courant number 1 — the update
collapses to a pure copy from the upstream neighbor, no averaging at
all. At any Courant number below 1, each step blends the upstream
neighbor's value in, which is unavoidably a weighted average of two
points. Averaging smooths out whatever's different between them, so
away from cn = 1 the scheme blurs a sharp shape instead of translating
it exactly. This study measures how much.

## Setup
RMSE values below are produced by `tests/diffusion_study.cpp`. 50-point periodic grid, dx = 0.02, initial condition a step block at
points 5-8 (value 1, everywhere else 0), 20 timesteps. Courant number
swept from 0.1 to 0.9 in steps of 0.2. Each run is compared against
the exact solution — the same initial block run through the solver at
cn = 1 for the same number of steps, which is exact by construction —
using RMSE (root-mean-square error) between the run's final array and
that exact array.

## Results
| Courant number | RMSE |
|----------------|------|
| 0.1 | 0.362069 |
| 0.3 | 0.344084 |
| 0.5 | 0.340257 |
| 0.7 | 0.335192 |
| 0.9 | 0.220912 |

![Courant 0.1](diffusion_cn_0.1.png)
*cn = 0.1, most blurred*

![Courant 0.9](diffusion_cn_0.9.png)
*cn = 0.9, closest to exact*

## Observations

RMSE drops as the Courant number rises toward 1. A larger Courant
number means less numerical diffusion per step, so the transported
block stays sharper and closer to the exact solution. Every Courant
number below 1 introduces some blur, but how much depends heavily on
how close cn is to 1.

That falls directly out of the first-order upwind update:

$$u_i^{n+1} = (1-\mathrm{cn})\,u_i^n + \mathrm{cn}\,u_{i-1}^n$$

The coefficients $(1-\mathrm{cn})$ and $\mathrm{cn}$ set how much of
the current cell and the upstream neighbor get blended together each
step. Their gap, $|2\,\mathrm{cn}-1|$, measures how close a single
update is to a pure copy versus a 50/50 average:

| Courant number | $\lvert 2\,\mathrm{cn}-1 \rvert$ |
|----------------|-----------------------:|
| 0.1 | 0.8 |
| 0.3 | 0.4 |
| 0.5 | 0.0 |
| 0.7 | 0.4 |
| 0.9 | 0.8 |

That gap moves by exactly 0.4 for every 0.2 step in cn — the per-step
effect of the Courant number is perfectly linear.

The measured RMSE isn't linear, though. The update doesn't run once;
it runs 20 times in a row. Each step adds a small amount of numerical
diffusion, and those small differences compound over the full run.
That's why RMSE barely moves between cn = 0.1 and cn = 0.7, then drops
sharply between 0.7 and 0.9 — the closer cn gets to 1, the closer each
of those 20 steps gets to a pure copy, and the compounding effect
shrinks fast.

A polynomial curve fit could match these five points well enough, but
it would only be describing the shape of this one dataset. It
wouldn't explain why the trend looks this way — that's the repeated
weighted average, not a coincidence of these particular numbers — and
it wouldn't predict what happens with a different step count, grid
resolution, or initial condition.