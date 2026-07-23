# GPU vs. CPU: Jacobi Speedup Study

Correctness for both CUDA kernels was already confirmed in earlier studies (`jacobi_validate.cu`, `jacobi_tiled_validate.cu`, `jacobi_convergence.cu`). What's left here is the comparison itself: how much faster the GPU actually is, and why that number came out the way it did.

## Setup

All runs solve a 100×100 grid with the top edge held at 100 degrees, for a fixed 10,000 sweeps rather than running to convergence — this measures per-sweep throughput, not solution accuracy. Every version uses the same methodology: full compiler optimization (`g++ -O2` for CPU, `nvcc -O3 -arch=sm_75` for GPU), 10 warm-up sweeps discarded, then independently timed trials reported as mean ± standard deviation. CPU is single-threaded Jacobi (`tests/bench_cpu.cpp`, sweep logic in `include/n_jacobi.hpp`), measured across two separate sessions, 10 trials each, after confirming the code itself was unchanged between sessions — 20 trials pooled into one number. GPU is a Tesla T4 on Colab, tested with both the naive one-thread-per-point kernel (`cuda/bench_gpu.cu`) and the shared-memory tiled kernel at its fastest block size, 8×8 (`cuda/bench_tiled_sweep.cu`), 10 trials each in a single session.

## Results

| Version         | Time (ms) | σ (ms) | Trials | Speedup vs. CPU |
|------------------|-----------|--------|--------|------------------|
| CPU (1 thread)   | 113.09    | 10.35  | 20     | 1.00×            |
| GPU naive        | 93.69     | 1.05   | 10     | 1.21×            |
| GPU tiled (8×8)  | 85.48     | 0.19   | 10     | 1.32×            |

## Analysis

The two CPU sessions disagreed by more than either session's own standard deviation would predict: 104.29 ms in the first, 121.90 ms in the second, against a within-session σ of about 5.5 ms in both. That gap only shows up across sessions, not within one — the pooled standard deviation (10.35 ms) is roughly double either individual session's own number. The code was checked and confirmed identical between the two runs, so this is real system noise: background processes, thermal state, or power mode differing between sessions, not the algorithm. The GPU numbers don't show this same effect here, but that's not a confirmed asymmetry yet — the GPU benchmarks were only run in a single session each, not repeated across separate sessions the way the CPU ended up being tested.

The GPU wins either way, but not by much. The more interesting number here might still be the standard deviation, not the mean. Both GPU kernels stay under 1.5% run-to-run variation within a session; the CPU sits at roughly 9% pooled across sessions. The T4 looks like a more consistent place to measure from than a general-purpose core sharing a machine with an OS and background processes, though a fair test of that claim would mean re-running the GPU benchmarks across separate sessions too.

The modest speedup has two real causes. First, every sweep is a separate kernel launch, and 10,000 launches means paying fixed launch overhead 10,000 times. At 100×100, that overhead eats a real fraction of each sweep's cost. A larger grid would do more actual work per launch while paying the same fixed cost per launch, which should push the ratio further toward the GPU — untested here, and the obvious next experiment. Second, Jacobi is memory-bound, which the CPU threading study already established. A T4's memory bandwidth is roughly 320 GB/s against a CPU core's 30-60 GB/s, and that gap is where a much larger speedup should eventually show up. At this grid size, the kernel isn't moving enough data per launch for that bandwidth gap to dominate over the fixed launch cost.

Tiling's edge over naive (1.09×) comes from real savings: each interior value gets read from global memory up to four times in the naive kernel and once in the tiled version. The same small-grid effect caps how much that shows up in total time — less traffic overall means less room for a memory-access optimization to matter.

The tiled kernel's first version, built in milestone 4, had a real bug. The halo-loading logic checked only local thread position within a block, which breaks for the ragged last block in any dimension that isn't an exact multiple of the block size — true for every grid/block combination tested, including the one supposedly validated at the time. It passed its correctness check by luck, not because the logic was right, and the bug only surfaced while building milestone 6's convergence study. The numbers above were re-measured after the fix and moved by under 0.1 ms from the broken version, confirming the bug hit correctness, not speed. Same memory traffic either way — whether the data being moved was right or garbage.

## Future work

Re-run this study at larger grid sizes (500×500, 1000×1000) to check whether the GPU's advantage grows as launch overhead becomes a smaller fraction of total runtime. Benchmark the existing multithreaded CPU Jacobi (`bench_mt.cpp`) under this same fixed-sweep, warm-up, multi-trial methodology — its 1.63× number came from a different setup (different grid size, run to convergence, no warm-up) and isn't directly comparable to the numbers above. Re-run the GPU benchmarks across separate sessions, the same way the CPU ended up being tested, to check whether the GPU's low variance holds up or was partly a single-session artifact. Beyond that, conjugate gradient on GPU is the original stretch goal, where convergence rate becomes the interesting question instead of raw per-sweep cost.