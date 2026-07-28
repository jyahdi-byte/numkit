# numkit

[![build](https://github.com/jyahdi-byte/numkit/actions/workflows/build.yml/badge.svg)](https://github.com/jyahdi-byte/numkit/actions/workflows/build.yml)

A numerical computing library written from scratch in modern C++, with
no external dependencies.

**The problem it solves:** predict the steady-state temperature at
every point of a heated plate. Mathematically this means solving
Laplace's equation, a partial differential equation. numkit
discretizes it with finite differences and solves the resulting sparse
linear system (~10⁴–10⁵ unknowns) with iterative methods, the same
pipeline at the core of computational fluid dynamics and finite
element software.

![Steady-state heat solution](docs/heatmap.png)

*Steady-state temperature field on a 300×200 plate with the top edge
held at 100° and the remaining edges at 0°, solved with SOR and
rendered to PPM by the library.*

## What it does

* **Discretization** — samples the continuous PDE on a uniform grid
  and replaces second derivatives with central finite differences,
  reducing the problem to a sparse linear system with one equation
  per interior point.
* **Iterative solvers** — Jacobi, Gauss-Seidel, and SOR, implemented
  from scratch and validated against each other and against exact
  solutions. SOR also has a closed-form auto-tuned variant that
  computes its relaxation factor from the grid size.
* **Parallelism** — a multithreaded CPU Jacobi solver with a measured,
  diagnosed speedup study, since upgraded to a persistent worker pool
  (`std::barrier`) that creates threads once instead of every sweep,
  and a CUDA port (`cuda/`) of the same solver: naive and
  shared-memory-tiled kernels, both validated and benchmarked.
* **Visualization** — renders solved fields to PPM images with a
  blue-to-red color map.

## Results

* **Verified second-order convergence** against manufactured
  analytical solutions: refining the grid 2× reduces error 4×
  (measured ratios 4.17 and 4.11 against the theoretical 4.0).
  Study: [docs/validation.md](docs/validation.md).
* On a 10×10 test problem, Jacobi / Gauss-Seidel / SOR converge in
  237 / 124 / 29 sweeps. SOR at its measured optimal relaxation
  factor (ω = 1.50, theory ~1.53) is 8× cheaper than Jacobi.
  Study: [docs/omega_study.md](docs/omega_study.md).
* Auto-tuned SOR (closed-form ω) converges in 298 sweeps versus 3,849
  for a fixed ω = 1.5 on a 100×100 grid, a 12× improvement, showing
  the fixed guess tuned for a small grid stops being a good guess as
  the grid grows.
* Multithreaded Jacobi reaches 1.63× speedup at 2 threads. The study
  identifies memory bandwidth as the ceiling and per-sweep thread
  respawn as the cost of going wider.
  Study: [docs/threading_study.md](docs/threading_study.md).
* Persistent worker pool rewrite (`jacobi_mt.hpp`) closes the
  respawn cost directly — threads are created once and synchronized
  every sweep with `std::barrier` instead of being recreated each
  round. Re-verified against the single-threaded solver: identical
  iteration count, zero difference in the converged grid.
* Naive CUDA Jacobi kernel (`cuda/bench_gpu.cu`, one thread per
  interior point, no shared-memory optimization): 93.7 ms mean
  across 10 trials (σ ≈ 1.1 ms) for 10,000 sweeps on a 100×100 grid,
  Tesla T4.
* Shared-memory tiled kernel with halo handling, swept across block
  sizes 8×8 through 32×32 (`cuda/bench_tiled_sweep.cu`): 8×8 wins at
  85.4 ms mean (σ ≈ 0.24 ms), beating the naive kernel on the same
  workload. Larger blocks were consistently slower (104 ms at
  32×32) — fewer blocks resident per SM means less work available
  to hide memory-fetch stalls, which matters more for a memory-bound
  kernel like Jacobi than a compute-bound one.

## GPU port

A CUDA port of the Jacobi solver lives in `cuda/`. It's built and
tested separately from the CPU code, since it needs `nvcc` and an
NVIDIA GPU — developed and run on Colab's T4 instances. Same stencil
computation as the CPU version, just moved from `std::thread` to GPU
threads.

* [x] Environment check — confirmed `nvcc` compiles and a kernel runs
      end to end.
* [x] `Grid` moved to device memory and back, verified against the
      original with real assertions (`cuda/grid_transfer_test.cu`).
* [x] Naive Jacobi kernel, one thread per interior point, validated
      against the CPU `jacobi_solve` to a `1e-9` tolerance
      (`cuda/jacobi_validate.cu`).
* [x] Benchmarking methodology — optimized build flags
      (`-O3 -arch=sm_75`), warm-up runs discarded before timing,
      10 timed trials, mean and standard deviation reported
      (`cuda/bench_gpu.cu`, stats in `include/stats.hpp`).
* [x] Shared-memory tiling with halo handling, dynamic block size
      via `extern __shared__`, validated against CPU `jacobi_solve`
      (`include/jacobi_tiled_kernel.cuh`, `cuda/jacobi_tiled_validate.cu`).
* [x] Grid/block size sweep study — 8×8 through 32×32, 8×8 wins
      (`cuda/bench_tiled_sweep.cu`).

## Roadmap

**Heat solver (complete)**

* [x] Grid data structure for 2D scalar fields
* [x] Jacobi, Gauss-Seidel, and SOR iterative solvers
* [x] Closed-form auto-tuned SOR
* [x] 2D steady-state heat solver (`apps/heat`) with PPM heatmap output
* [x] Convergence and validation studies in `docs/`
* [x] Multithreaded Jacobi solver
* [x] Persistent worker pool (`std::barrier`), replacing per-sweep
      thread respawn

**Hyperbolic PDEs: 1D advection (in progress)**

The heat solver is elliptic — it settles to a steady state. Hyperbolic
problems move or oscillate instead of settling, which calls for a
different approach: explicit time-stepping instead of relaxation,
upwind stencils instead of centered ones, and a CFL stability limit
in place of a convergence tolerance.

* [x] CFL and upwind derivation from first principles — why the
      numerical domain of dependence has to contain the true one
* [x] `Grid1D`, a periodic 1D grid with wraparound indexing
      (`include/grid1d.hpp`)
* [x] Single-threaded upwind solver (`include/advection.hpp`),
      conservation-checked with a real assertion
      (`tests/test_advection.cpp`)
* [x] Space-time diagram visualization — one PPM row per timestep,
      `write_ppm` templated to take either a `Grid` or a
      `SpaceTimeLog` (`include/space_time_log.hpp`, `include/ppm.hpp`)
* [x] CFL enforcement — Courant number is a direct input, `Δt` is
      derived from it (`compute_dt`/`compute_cn`); a violation
      reports the actual offending value instead of a bare assert,
      and the deliberate CFL > 1 demo lives in its own file
      (`tests/demo_cfl_violation.cpp`) so it can crash on purpose
      without threatening the real test suite
* [ ] Exact-solution test wired into CI
* [ ] Numerical diffusion vs. Courant number study

The wave equation — leapfrog time-stepping, with energy conservation
as the validation method instead of shape-matching — is a deferred
stretch goal once advection is closed out.

**Backlog**

* [ ] Real test assertions for the heat solver — `test_solvers.exe`
      still just prints sweep counts and a sample grid value instead
      of asserting against expected results; CI compiles and runs it
      but doesn't check correctness yet. `test_advection.cpp` already
      proves the pattern works, just not retrofitted here yet.
* [ ] Red-black ordering, a parallelizable reformulation of
      Gauss-Seidel
* [ ] Conjugate gradient solver

## Building

Requires g++ and make. From the repo root:

    make

builds the test programs and the heat app. Run `./heat.exe` to solve
the demo problem and write `heat.ppm`.

The CUDA targets need `nvcc` and an NVIDIA GPU, which this repo
doesn't assume you have locally — run these on Colab, not your own
machine, unless you actually have an NVIDIA GPU:

    make cuda-test              # grid transfer round-trip test
    make jacobi-validate        # naive kernel vs. CPU jacobi_solve
    make bench-gpu              # naive kernel: warm-up + 10 trials, mean/stddev
    make jacobi-tiled-validate  # tiled kernel vs. CPU jacobi_solve
    make bench-tiled-sweep      # block-size sweep: 8x8, 16x16, 24x24, 32x32

Every push is checked by GitHub Actions: the CPU build is compiled
and run for real, the CUDA build is compile-checked only (all five
`cuda/*.cu` files), since GitHub's runners have no GPU — actual test
execution stays a local Colab run before trusting a commit.

## License

MIT