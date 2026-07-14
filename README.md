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
* **Parallelism** — a multithreaded Jacobi solver (std::thread) with
  a measured, diagnosed speedup study.
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
* [ ] Naive Jacobi kernel, validated against the CPU `jacobi_solve`.
* [ ] Benchmarking methodology, then a real GPU vs. CPU speedup study.

## Roadmap

**Heat solver (complete)**

* [x] Grid data structure for 2D scalar fields
* [x] Jacobi, Gauss-Seidel, and SOR iterative solvers
* [x] Closed-form auto-tuned SOR
* [x] 2D steady-state heat solver (`apps/heat`) with PPM heatmap output
* [x] Convergence and validation studies in `docs/`
* [x] Multithreaded Jacobi solver

**Possible extensions**

* [ ] Persistent worker pool (std::barrier) to eliminate thread respawn cost
* [ ] Red-black ordering for parallel Gauss-Seidel/SOR

## Building

Requires g++ and make. From the repo root:

    make

builds the test programs and the heat app. Run `./heat.exe` to solve
the demo problem and write `heat.ppm`.

The CUDA target needs `nvcc` and an NVIDIA GPU, which this repo
doesn't assume you have locally:

    make cuda-test

builds and runs the GPU grid-transfer test — do this on Colab, not
your own machine, unless you actually have an NVIDIA GPU. Every push
is checked by GitHub Actions: the CPU build is compiled and run for
real, the CUDA build is compile-checked only, since GitHub's runners
have no GPU.

## License

MIT
