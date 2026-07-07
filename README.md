# numkit

A numerical computing library written from scratch in modern C++, with no
external dependencies. numkit solves partial differential equations by finite
differences: it discretizes the steady-state heat equation (Laplace's
equation) on a 2D grid and solves the resulting sparse linear system with
iterative methods.

![Steady-state heat solution](docs/heatmap.png)

*Steady-state temperature field on a 300×200 plate with the top edge held at
100° and the remaining edges at 0°, solved with SOR and rendered to PPM by
the library.*

## What it does

* **Discretization** — the continuous PDE is sampled on a uniform grid, and
  second derivatives are replaced with central finite differences, reducing
  the problem to a sparse linear system (one equation per interior point).
* **Iterative solvers** — Jacobi, Gauss-Seidel, and SOR, implemented from
  scratch and validated against each other on identical problems.
* **Visualization** — solved fields are rendered to PPM images with a
  blue-to-red color map.

## Results so far

* On a 10×10 test problem (tolerance 1e-6), the three solvers converge in
  237 / 124 / 29 sweeps respectively — SOR at its measured optimal
  relaxation factor is 8× cheaper than Jacobi.
* The measured optimal SOR factor (ω = 1.50) matches the theoretical
  prediction (~1.53) for this grid size. Full study:
  [docs/omega_study.md](docs/omega_study.md).

## Roadmap

* [x] Grid data structure for 2D scalar fields
* [x] Jacobi, Gauss-Seidel, and SOR iterative solvers
* [x] 2D steady-state heat solver (`apps/heat`) with PPM heatmap output
* [ ] Convergence and validation studies in `docs/`
* [ ] Multithreaded Jacobi solver
* [ ] 2D truss finite element solver (fall 2026)

## Building

Requires g++ and make. From the repo root:

    make

builds the test programs and the heat app. Run `./heat.exe` to solve the
demo problem and write `heat.ppm`.

## License

MIT
