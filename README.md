# NumKit

**By Jonas Yahdi**  
Applied Mathematics / Computational Mathematics & Computer Science  
Santa Clara University

[![CI](https://github.com/jyahdi-byte/numkit/actions/workflows/build.yml/badge.svg)](https://github.com/jyahdi-byte/numkit/actions/workflows/build.yml)

**NumKit is a C++ numerical-computing project for solving and studying elliptic PDEs, with an emphasis on numerical verification, parallelism, and performance engineering.**

The project starts with a finite-difference model of steady-state heat flow and follows it all the way through multiple CPU solvers, irregular domains, multithreading, CUDA kernels, correctness validation, convergence studies, and performance measurements.

It is designed to demonstrate the kind of workflow used in scientific computing and computational engineering:

```text
Physical / Mathematical Model
            ↓
Finite-Difference Discretization
            ↓
Iterative Linear Solver
            ↓
Independent Verification
            ↓
Convergence / Parameter Studies
            ↓
CPU Parallelization
            ↓
CUDA Implementation
            ↓
Benchmarking and Optimization
```

The central engineering question is not only **“does it produce an answer?”**, but also **“why should I trust the answer, and what limits its performance?”**

---

## Why this project is useful as a portfolio project

NumKit combines several areas that are often shown separately in internship projects:

- **Numerical methods:** finite differences, elliptic PDEs, iterative methods, convergence and error analysis.
- **C++ systems work:** C++20, memory layout, reusable headers, assertions, multithreading, synchronization, and benchmarking.
- **Parallel computing:** Jacobi parallelization, persistent worker threads, `std::barrier`, red-black ordering, and race-free parallel updates.
- **GPU computing:** CUDA kernels, device/host transfers, shared-memory tiling, block-size experiments, and GPU validation.
- **Scientific software engineering:** unit-style tests, analytical validation, conservation checks, CI, reproducible studies, and documented experiments.

The repository also contains written studies explaining *why* measured behavior differs from the ideal theoretical model. That makes the project more than a collection of solver implementations.

---

## What is implemented

### Elliptic PDE solvers

NumKit solves steady-state elliptic problems such as Laplace's equation and its variable-conductivity form.

Implemented methods include:

- **Jacobi**
- **Gauss-Seidel**
- **Successive Over-Relaxation (SOR)**
- **Red-black Gauss-Seidel**
- **Red-black SOR**
- **Conjugate Gradient (CG)**
- **Preconditioned Conjugate Gradient (diagonal/Jacobi preconditioner)**

The grid is stored densely, but cells carry a type:

```text
INTERIOR  → unknown updated by the solver
FIXED     → prescribed Dirichlet value
HOLE      → excluded / insulated region
```

The same representation is reused across the CPU solvers and CUDA implementations.

### Variable conductivity

The `Grid` class also supports a spatially varying conductivity field. Face conductivities are computed with a harmonic average, allowing the solver to model materials with different thermal conductivities without changing the overall grid representation.

This is demonstrated by the `main_lens.cpp` application, which combines:

- low-conductivity background material,
- high-conductivity fin regions,
- fixed hot and cold regions,
- smoothly varying conductivity bumps,
- multithreaded red-black SOR.

### Irregular and embedded domains

The solver supports geometry that is not just a simple rectangular interior.

Examples include:

- circular insulated obstacles,
- rectangular insulated regions,
- internal fixed-temperature regions,
- narrow openings through an insulated wall,
- spatially varying material properties.

For insulated `HOLE` cells, the implementation uses a ghost-point style reflection: a missing neighbor contributes the center value to the stencil, preserving the four-point update structure while enforcing zero normal flux at steady state.

---

# Numerical verification

A numerical solver is only useful when its output can be checked independently. NumKit therefore treats **verification as a first-class part of the project**.

## Manufactured-solution convergence

The solver is checked against the analytical solution

\[
T(x,y)=\sin(\pi x)\sinh(\pi y)
\]

for a compatible boundary-value problem.

A grid-refinement study gives approximately second-order convergence:

| Grid change | Measured error ratio |
|---|---:|
| 25 → 50 | 4.17 |
| 50 → 100 | 4.11 |

For a second-order finite-difference discretization, halving the grid spacing should reduce the leading discretization error by roughly \(2^2=4\). The measured ratios are consistent with that expectation.

See [`docs/validation.md`](docs/validation.md) for the validation work and [`docs/omega_study.md`](docs/omega_study.md) for the SOR study.

## Solver cross-validation

Different algorithms should converge to the same physical solution even though they take different computational paths. The test suite compares solver outputs and checks pointwise agreement, symmetry, boundary behavior, and other invariants.

For one representative test, the measured sweep counts were:

| Method | Sweeps |
|---|---:|
| Jacobi | 298 |
| Gauss-Seidel | 154 |
| SOR (`ω = 1.5`) | 41 |
| SOR (auto/theoretical `ω`) | 40 |
| Conjugate Gradient | 28 |
| Preconditioned CG | 27 |

The important result is not that one method “wins,” but that the implementations can be compared on the **same problem while independently checking that they reach the same solution**.

## Conservation and boundary checks

The tests also verify properties such as:

- fixed cells remain fixed,
- masked cells behave as intended,
- symmetry is preserved when the geometry is symmetric,
- solver outputs agree across implementations,
- CUDA results match the CPU reference within tight numerical tolerances.

---

# Algorithmic studies

The repository includes focused studies rather than only code.

### SOR parameter selection

SOR depends on a relaxation parameter \(\omega\). NumKit sweeps \(\omega\) and compares measured convergence against the theoretical optimum for a rectangular grid.

This turns a textbook formula into an experimentally testable claim.

### Red-black ordering

Standard Gauss-Seidel and SOR perform in-place updates, which makes naive parallelization unsafe because neighboring cells can be read while they are being written.

NumKit uses **red-black ordering**:

```text
R B R B R
B R B R B
R B R B R
B R B R B
```

All cells of one color can be updated before advancing to the other color. This creates a synchronization-friendly update pattern suitable for both CPU threads and CUDA kernels.

A representative 80×80 comparison:

| Method | Sweeps |
|---|---:|
| Gauss-Seidel | 11,346 |
| Red-black Gauss-Seidel | 11,365 |
| SOR | 323 |
| Red-black SOR | 321 |

The near-identical convergence counts show why red-black ordering is useful: it changes the update schedule enough to enable parallel execution without materially changing convergence on this test.

See [`docs/red_black_study.md`](docs/red_black_study.md).

### Conjugate Gradient scaling

The CG implementation is **matrix-free**. Instead of explicitly assembling a large sparse matrix, it applies the discretized stencil directly to vectors.

On the study grid, doubling the problem dimension produced approximately:

| Grid size | Gauss-Seidel | SOR | CG |
|---:|---:|---:|---:|
| 10 | 179 | 39 | 24 |
| 20 | 756 | 80 | 60 |
| 40 | 2,977 | 161 | 124 |
| 80 | 11,346 | 323 | 252 |

The sweep-growth trends match the qualitative behavior expected from the methods, while also showing why CG can be attractive when an easily tuned SOR parameter is unavailable.

The repository also contains a study of diagonal preconditioning on irregular domains. In the included geometries, the simple diagonal preconditioner produced little or no improvement at larger sizes. That result is documented rather than hidden because negative or inconclusive experimental results are still useful numerical evidence.

See [`docs/conjugate_gradient_study.md`](docs/conjugate_gradient_study.md).

---

# CPU parallelization

Jacobi is naturally parallel because every update reads only the previous iteration. NumKit first implements a straightforward threaded version and then a **persistent worker model** using `std::barrier` to reduce repeated thread-management overhead.

The project also extends red-black Gauss-Seidel and SOR to threaded execution. The synchronization structure follows directly from the numerical dependency graph: one color is completed before the other begins.

A benchmark study showed an important practical lesson: increasing the thread count does not guarantee linear speedup.

| Threads | Runtime | Speedup |
|---:|---:|---:|
| 1 | 152 s | 1.00× |
| 2 | 93 s | 1.63× |
| 4 | 110 s | 1.38× |
| 8 | 107.8 s | 1.41× |

The study attributes the weak scaling to factors such as memory bandwidth, synchronization, thread-management overhead, and hardware limits rather than assuming that more threads must be faster.

See [`docs/threading_study.md`](docs/threading_study.md).

---

# CUDA acceleration

The GPU implementation was developed incrementally:

```text
CPU reference implementation
          ↓
Naive CUDA kernel
          ↓
CPU/GPU correctness validation
          ↓
Shared-memory tiled kernel
          ↓
Block-size sweep
          ↓
Benchmark and analysis
```

### Shared-memory Jacobi

The tiled Jacobi kernel loads a local grid tile and halo cells into CUDA shared memory so neighboring values can be reused without repeatedly reading the same data from global memory.

The project includes:

- naive Jacobi CUDA kernel,
- shared-memory tiled Jacobi kernel,
- CUDA implementations of red-black Gauss-Seidel and SOR,
- GPU/CPU correctness checks,
- block-size experiments,
- fixed-work throughput benchmarks.

### Representative GPU benchmark

A documented 100×100 fixed-work benchmark on a Tesla T4 produced:

| Implementation | Mean time | Speedup vs. CPU |
|---|---:|---:|
| CPU, single thread | 113.09 ms | 1.00× |
| CUDA, naive | 93.69 ms | 1.21× |
| CUDA, tiled | **85.48 ms** | **1.32×** |

The modest speedup is an intentional part of the analysis. Jacobi is memory-dependent, the problem is relatively small, and kernel-launch overhead matters. The study therefore focuses on understanding the bottleneck instead of reporting a GPU speedup without context.

See [`docs/gpu_speedup_study.md`](docs/gpu_speedup_study.md).

---

# Debugging and validation of parallel code

The repository includes CUDA validation programs that compare GPU results against the CPU reference solution with tight tolerances.

This caught a halo-loading bug in the tiled Jacobi kernel for ragged block sizes. The important part of the workflow was:

```text
Implement
   ↓
Validate against reference
   ↓
Find a failing / suspicious case
   ↓
Fix kernel
   ↓
Re-run validation
```

That workflow is especially important for numerical and parallel programs because a kernel can produce visually plausible output while still containing a race, indexing error, or boundary bug.

---

# Build and run

## Requirements

### CPU

- Linux/macOS/WSL or another environment with a modern GNU C++ toolchain
- **C++20**
- `make`

### GPU work

- NVIDIA GPU for execution
- CUDA toolkit / `nvcc`

The repository's CI uses CUDA Toolkit 12.4 to **compile-check** the CUDA targets. GitHub-hosted runners used by the workflow do not provide an NVIDIA GPU, so the CUDA CI job verifies compilation rather than running GPU kernels.

## Build the CPU project

```bash
make all
```

This builds the main solver tests, numerical studies, heat-flow demo, and multithreaded targets.

Run a few core correctness checks:

```bash
./test_solvers.exe
./test_grid.exe
./test_jacobi.exe
./test_omega_auto.exe
```

Some benchmark-oriented and multithreaded programs intentionally perform much larger workloads and can take substantially longer.

## Run the heat demo

Basic heat-flow example:

```bash
./heat.exe
```

With a circular insulated obstacle:

```bash
# after building heat.exe
./heat.exe --obstacle
```

The program writes a PPM image of the converged temperature field.

Other applications are in [`apps/heat/`](apps/heat/).

## CUDA targets

Examples:

```bash
make cuda-test
make jacobi-validate
make jacobi-tiled-validate
make gauss-seidel-rb-validate
make sor-rb-validate
make bench-gpu
make bench-tiled-sweep
```

These targets assume a working CUDA installation and NVIDIA GPU for execution.

---

# Repository structure

```text
numkit/
├── include/                  # Core numerical library headers
│   ├── grid.hpp              # Grid, geometry, cell types, conductivity field
│   ├── jacobi.hpp            # Jacobi solver
│   ├── gauss_seidel.hpp      # Gauss-Seidel solver
│   ├── sor.hpp               # SOR solver
│   ├── *_rb.hpp              # Red-black variants
│   ├── *_mt.hpp              # Multithreaded variants
│   ├── conjugate_gradient.hpp # CG / diagonal-PCG
│   ├── *_kernel.cuh          # CUDA kernels
│   └── ppm.hpp               # Output/visualization helpers
│
├── apps/
│   └── heat/                 # Demonstration applications
│
├── tests/                    # Correctness tests, sweeps, and benchmarks
│   └── cuda/                 # CUDA validation / benchmark programs
│
├── docs/                     # Numerical studies and analysis
│   ├── validation.md
│   ├── omega_study.md
│   ├── red_black_study.md
│   ├── threading_study.md
│   ├── conjugate_gradient_study.md
│   └── gpu_speedup_study.md
│
├── .github/workflows/
│   └── build.yml             # CPU test CI + CUDA compile checks
│
├── Makefile
└── README.md
```

---

# Technical stack

| Area | Technologies |
|---|---|
| Language | C++20 |
| GPU | CUDA / `nvcc` |
| Numerical methods | Finite differences, iterative solvers, CG, relaxation methods |
| PDEs | Elliptic / steady-state diffusion and heat conduction |
| Parallelism | `std::thread`, `std::barrier`, red-black ordering, CUDA |
| GPU optimization | Shared memory, tiling, block-size experiments |
| Validation | Assertions, manufactured solutions, cross-solver checks, CPU/GPU comparison |
| Performance | Fixed-work benchmarks, speedup studies, runtime variance analysis |
| Tooling | Make, GitHub Actions |

---

# What this project demonstrates

NumKit is best viewed as a **scientific-computing implementation and investigation**, not simply a PDE demo.

It demonstrates the ability to:

- translate a mathematical model into a working numerical algorithm,
- implement and compare multiple iterative solvers,
- reason about convergence and asymptotic behavior,
- validate numerical output against analytical or independent reference results,
- represent irregular domains and variable material properties,
- identify data dependencies that constrain parallelism,
- redesign algorithms for thread-safe CPU execution,
- move a stencil computation to CUDA,
- optimize memory access with shared-memory tiling,
- benchmark implementations and interpret non-ideal scaling,
- debug numerical and parallel implementations using targeted validation, and
- document experiments well enough for another engineer to reproduce and critique the results.

For internships in **scientific computing, computational engineering, numerical methods, HPC, GPU computing, or performance-oriented C++ development**, the strongest parts of the repository are the combination of numerical correctness, parallel implementation, and measured performance analysis.

---

# Current scope and next directions

NumKit intentionally focuses on **elliptic PDEs**. Time-dependent PDE solvers and a Black-Scholes implementation have been moved to a separate project.

Natural next steps for NumKit include:

- stronger automated coverage for edge cases and irregular geometries,
- broader CUDA execution tests on real GPUs in addition to compile-only CI,
- better benchmarking controls and repeated-run statistics,
- sparse-matrix representations for larger problem sizes,
- stronger preconditioning strategies for irregular or variable-coefficient systems,
- additional 3D data structures and solvers.

These are extensions of the existing architecture rather than a replacement for it.

---

## Documentation and studies

The most useful deep dives are:

- [`docs/validation.md`](docs/validation.md) — numerical correctness and error checks
- [`docs/omega_study.md`](docs/omega_study.md) — SOR parameter selection
- [`docs/red_black_study.md`](docs/red_black_study.md) — race-free red-black ordering
- [`docs/threading_study.md`](docs/threading_study.md) — CPU scaling and synchronization
- [`docs/conjugate_gradient_study.md`](docs/conjugate_gradient_study.md) — CG and preconditioning studies
- [`docs/gpu_speedup_study.md`](docs/gpu_speedup_study.md) — CPU vs. CUDA performance analysis

## Author

**Jonas Yahdi**  
Applied Mathematics / Computational Mathematics & Computer Science  
Santa Clara University

GitHub: [jyahdi-byte](https://github.com/jyahdi-byte)  
LinkedIn: [Jonas Yahdi](https://www.linkedin.com/in/jonas-yahdi-b00a1a226/)  
Email: [jyahdi@scu.edu](mailto:jyahdi@scu.edu) · [jyahdi@gmail.com](mailto:jyahdi@gmail.com)

## License

MIT — see [`LICENSE`](LICENSE).
