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

Early development. The first milestone is the heat solver.

## License

MIT
