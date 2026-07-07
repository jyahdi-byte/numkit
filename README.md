# numkit

A numerical computing library written from scratch in modern C++. numkit provides
the core data structures and iterative solvers behind small physics applications,
starting with a 2D steady-state heat (Poisson) solver.

## Goals

* Build core numerical tools (grids, iterative linear solvers) with no external dependencies
* Validate everything: comparison against analytical solutions, residual tracking, grid refinement studies
* Grow the library through real applications rather than isolated exercises

## Roadmap

* [ ] Grid data structure for 2D scalar fields
* [ ] Jacobi, Gauss-Seidel, and SOR iterative solvers
* [ ] 2D steady-state heat solver (`apps/heat`) with PPM heatmap output
* [ ] Convergence and validation studies in `docs/`
* [ ] Multithreaded Jacobi solver (pthreads)
* [ ] 2D truss finite element solver (fall 2026)

## Status

In active development. The core library is working:

* `Grid` — bounds-checked 2D field container
* Three iterative solvers (Jacobi, Gauss-Seidel, SOR), validated against
  each other on identical problems
* SOR relaxation factor study ([docs/omega_study.md](docs/omega_study.md)) —
  measured optimum ω = 1.50 vs. theoretical ~1.53 for the test grid

Next up: the 2D steady-state heat application with PPM heatmap output.

## License

MIT
