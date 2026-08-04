# numkit

[![build](https://github.com/jyahdi-byte/numkit/actions/workflows/build.yml/badge.svg)](https://github.com/jyahdi-byte/numkit/actions/workflows/build.yml)

A numerical PDE library built from scratch in C++ and CUDA, no
external dependencies. All three classical PDE families are
implemented — steady-state (elliptic), spreading (parabolic), and
translating (hyperbolic) — each one derived from the underlying math,
implemented, and checked against an independent exact answer rather
than just run and eyeballed.

## Highlights

* Verified second-order convergence against manufactured analytical
  solutions (measured ratios 4.17 / 4.11 against the theoretical 4.0
  for a 2x grid refinement).
* SOR at its measured optimal relaxation factor is 8x cheaper than
  Jacobi; the closed-form auto-tuned version is 12x cheaper than a
  fixed guess once the grid gets large.
* A CUDA port with a shared-memory tiled kernel that beats the naive
  version, benchmarked across five block sizes on a Tesla T4.
* A second PDE family, hyperbolic advection, built the same way as
  the first: derived from first principles, implemented, validated
  against an exact solution, and wired into CI — not just "it runs."
* A third PDE family, parabolic diffusion, with its stability limit
  derived from a von Neumann (Fourier-mode) analysis rather than
  assumed, and a measured convergence order (ratios of 4.01, 4.00,
  4.00 against a theoretical 4.0) — the same rigor as the elliptic
  study, now on a time-marching scheme.
* Every core numerical claim has a real assertion behind it, not a
  printed number someone has to eyeball: conservation checked to
  `1e-9` on both time-marching solvers, the advection solver checked
  against its exact solution to `1e-9`, the diffusion solver checked
  against its exact solution to `2e-4`, the CUDA kernels checked
  against the CPU solver to `1e-9`.

## Steady-state solvers (elliptic)

**The problem it solves:** predict the steady-state temperature at
every point of a heated plate. Mathematically this means solving
Laplace's equation. numkit discretizes it with finite differences and
solves the resulting sparse linear system (~10⁴–10⁵ unknowns) with
iterative methods, the same pipeline at the core of CFD and finite
element software.

![Steady-state heat solution](docs/heatmap.png)

*Steady-state temperature field on a 300×200 plate with the top edge
held at 100° and the remaining edges at 0°, solved with SOR and
rendered to PPM by the library.*

**What it does**

* **Discretization** — samples the continuous PDE on a uniform grid
  and replaces second derivatives with central finite differences,
  reducing the problem to a sparse linear system with one equation
  per interior point.
* **Iterative solvers** — Jacobi, Gauss-Seidel, and SOR, implemented
  from scratch and validated against each other and against exact
  solutions. SOR also has a closed-form auto-tuned variant that
  computes its relaxation factor from the grid size.
* **Parallelism** — a multithreaded CPU Jacobi solver with a measured,
  diagnosed speedup study, upgraded to a persistent worker pool
  (`std::barrier`) that creates threads once instead of every sweep,
  plus a CUDA port (`cuda/`): naive and shared-memory-tiled kernels,
  both validated and benchmarked.
* **Visualization** — renders solved fields to PPM images with a
  blue-to-red color map.

**Results**

* On a 10×10 test problem, Jacobi / Gauss-Seidel / SOR converge in
  237 / 124 / 29 sweeps. Study: [docs/omega_study.md](docs/omega_study.md).
* Auto-tuned SOR converges in 298 sweeps versus 3,849 for a fixed
  ω = 1.5 on a 100×100 grid — the fixed guess tuned for a small grid
  stops being a good guess as the grid grows.
* Multithreaded Jacobi reaches 1.63x speedup at 2 threads; memory
  bandwidth is the ceiling. Study: [docs/threading_study.md](docs/threading_study.md).
* Persistent worker pool rewrite (`jacobi_mt.hpp`) re-verified
  against the single-threaded solver: identical iteration count, zero
  difference in the converged grid.
* Naive CUDA Jacobi kernel: 93.7 ms mean across 10 trials
  (σ ≈ 1.1 ms) for 10,000 sweeps on a 100×100 grid, Tesla T4.
* Shared-memory tiled kernel, swept across block sizes 8×8 through
  32×32: 8×8 wins at 85.4 ms mean (σ ≈ 0.24 ms), beating the naive
  kernel. Larger blocks were consistently slower — fewer blocks
  resident per SM means less work available to hide memory-fetch
  stalls, which matters more for a memory-bound kernel like Jacobi
  than a compute-bound one.

Full validation methodology: [docs/validation.md](docs/validation.md).

## Time-marching solvers (hyperbolic): 1D advection

The heat solver settles to a steady state. Advection doesn't — a
shape gets carried along by a flow and should translate without
changing, forever. That calls for a different approach: explicit
time-stepping instead of relaxation, an upwind stencil instead of a
centered one, and a CFL stability limit instead of a convergence
tolerance.

![Numerical diffusion at Courant 0.9](docs/diffusion_cn_0.9.png)

*Space-time diagram: each row is a timestep, each column a grid
point. The diagonal streak is the pulse translating correctly; this
run (Courant number 0.9) stays close to the exact shifted solution.*

**What it does**

* **`Grid1D`** — a periodic 1D grid with wraparound indexing built
  into `at()`, so the solver never special-cases the domain's edges.
* **Upwind solver** (`advection.hpp`) — derived from the advection
  equation by hand, not copied from a reference: forward-difference
  in time, backward-difference in space, biased toward the direction
  information actually travels from.
* **CFL as a first-class input** — the Courant number is a direct
  parameter, and `Δt` is derived from it (`compute_dt`/`compute_cn`),
  so a caller can never accidentally pick an unstable step size. A
  violation reports the actual offending value instead of crashing
  silently, and the deliberate CFL > 1 demo lives in its own file
  (`demo_cfl_violation.cpp`) so it can fail on purpose without
  threatening the real test suite.
* **Space-time visualization** — one PPM row per timestep;
  `write_ppm` is templated to take either a 2D `Grid` or a
  `SpaceTimeLog`, so the heat solver's renderer works here with no
  duplicated code.

**Results**

* Exact-solution test: at Courant number 1 the scheme is provably
  exact (a pure copy, no averaging), so the solver's output is
  checked point-by-point against the analytically shifted initial
  condition to `1e-9`, wired into CI on every push.
* Numerical diffusion study: RMSE against the exact solution drops
  from 0.362 at Courant 0.1 to 0.221 at Courant 0.9, and the drop
  isn't linear — it's flat through most of the range and falls
  sharply near Courant 1, because a small per-step blur compounds
  multiplicatively over repeated timesteps. Study:
  [docs/diffusion_study.md](docs/diffusion_study.md).

## Time-marching solvers (parabolic): 1D diffusion

Advection translates a shape without changing it. Diffusion does the
opposite — nothing moves sideways, a sharp pulse just spreads out and
flattens until it's uniform. Same explicit time-marching style as
advection (`Grid1D`, periodic boundaries, no new infrastructure), but
the stability limit comes from a different tool: a von Neumann
(Fourier-mode) analysis rather than a domain-of-dependence argument,
because the failure mode is a point overshooting its own neighbors,
not information arriving from outside the numerical domain.

**What it does**

* **FTCS solver** (`diffusion.hpp`) — forward-difference in time,
  central-difference in space, using both neighbors (`i-1` and
  `i+1`), unlike advection's one-sided upwind stencil.
* **`r` as a first-class input** — same pattern as advection's Courant
  number: `r = αΔt/Δx²` is the parameter, `Δt` is derived from it
  (`compute_r`), and a violation reports the actual offending value.
  The deliberate `r > 1/2` demo lives in its own file
  (`demo_r_violation.cpp`), same reasoning as `demo_cfl_violation.cpp`.
* **Conservation by construction** — summed over the whole periodic
  grid, every correction term cancels exactly, so total heat/mass is
  unchanged by a step. Checked to `1e-9` in `test_diffusion.cpp`.
* **Exact validation** — a Gaussian initial condition stays Gaussian,
  with `σ(t)² = σ₀² + 2αt` (the heat-kernel result) and amplitude
  scaling as `1/σ(t)` to keep the area under the curve fixed. Checked
  pointwise against the closed-form solution.

**Results**

* Exact-solution test: unlike advection at Courant 1, FTCS diffusion
  always carries some discretization error at finite `Δx`, so the
  tolerance (`2e-4`) is set from the scheme's actual measured error at
  this resolution rather than copied from advection's `1e-9`.
* Convergence study: sweeping `Δx` with `r` held fixed (so `Δt` scales
  automatically with `Δx²`) gives error ratios of 4.01, 4.00, and 4.00
  against the theoretical 4.0 for the first three refinements — clean
  second-order convergence. The final ratio drops to 2.55 once
  floating-point roundoff, accumulated over ~128,000 steps at the
  finest resolution, becomes large enough to compete with the
  (very small) discretization error. Study:
  [docs/diffusion_convergence.md](docs/diffusion_convergence.md).

## GPU port

A CUDA port of the Jacobi solver lives in `cuda/`, built and tested
separately from the CPU code since it needs `nvcc` and an NVIDIA
GPU — developed and run on Colab's T4 instances.

* [x] `Grid` moved to device memory and back, verified with real
      assertions (`cuda/grid_transfer_test.cu`).
* [x] Naive Jacobi kernel, validated against `jacobi_solve` to
      `1e-9` (`cuda/jacobi_validate.cu`).
* [x] Benchmarking methodology — `-O3 -arch=sm_75`, warm-up runs
      discarded, 10 timed trials, mean and standard deviation
      reported (`cuda/bench_gpu.cu`, `include/stats.hpp`).
* [x] Shared-memory tiling with halo handling, validated against
      the CPU solver (`include/jacobi_tiled_kernel.cuh`).
* [x] Block-size sweep, 8×8 through 32×32 (`cuda/bench_tiled_sweep.cu`).

## Engineering practices

* **CI on every push** — GitHub Actions builds every CPU target and
  runs the tests that have real assertions behind them
  (`test_solvers.exe`, `test_grid.exe`, `test_jacobi.exe`,
  `test_mt.exe`, `test_omega_auto.exe`, `test_advection.exe`,
  `test_advection_exact.exe`, `test_diffusion.exe`,
  `test_diffusion_exact.exe`, `diffusion_convergence.exe`). The CUDA
  build is compile-checked only, since GitHub's runners have no GPU;
  actual CUDA test execution stays a local Colab run before trusting
  a commit.
* **Assertions over print statements** — a printed number can say
  anything and still exit cleanly; only a failing `assert()` actually
  turns CI red. Every claim above that's marked as tested is backed
  by one, including the original heat-solver tests
  (`test_solvers.cpp`, `test_grid.cpp`, `test_jacobi.cpp`), which
  assert cross-solver agreement, grid indexing, and monotonic
  temperature falloff from a heated edge rather than just printing
  values.

## Building

Requires g++ and make. From the repo root:

    make

builds every CPU target. Run `./heat.exe` to solve the steady-state
demo problem, `./test_advection_exact.exe` to run the advection solver
against its exact solution, or `./test_diffusion_exact.exe` for the
same check on the diffusion solver.

The CUDA targets need `nvcc` and an NVIDIA GPU, which this repo
doesn't assume you have locally — run these on Colab unless you
actually have one:

    make cuda-test              # grid transfer round-trip test
    make jacobi-validate        # naive kernel vs. CPU jacobi_solve
    make bench-gpu              # naive kernel: warm-up + 10 trials, mean/stddev
    make jacobi-tiled-validate  # tiled kernel vs. CPU jacobi_solve
    make bench-tiled-sweep      # block-size sweep: 8x8 through 32x32

## Roadmap

* [ ] Black-Scholes as a heat equation — a change of variables
      (log-price, time-to-maturity, a discounting substitution) turns
      the option-pricing PDE into exactly the constant-coefficient
      heat equation solved above, reusing the FTCS solver unchanged
      with a real (non-periodic) boundary condition for the first
      time, validated against the closed-form Black-Scholes price.

## License

MIT