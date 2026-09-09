# GPU vs. CPU: 3D Jacobi Speedup Study

This extends `gpu_speedup_study.md` to three dimensions. The 2D study's own future-work section predicted the GPU's advantage should grow at a larger grid size, since launch overhead becomes a smaller fraction of the total work per launch. This is that test.

## Setup

All runs solve a 40×40×40 grid with the entire top face held at 100 degrees, for a fixed 5000 sweeps rather than running to convergence — same per-sweep-throughput methodology as the 2D study. CPU is single-threaded Jacobi (`tests/bench_cpu_3d.cpp`, sweep logic in `include/n_jacobi.hpp`), compiled with `g++ -O2`, 10 warm-up sweeps discarded, 10 timed trials. GPU is a T4 on Colab, compiled with `nvcc -O3 -arch=sm_75`, tested with both the naive kernel (`tests/cuda/bench_gpu_3d.cu`) and the shared-memory tiled kernel (`tests/cuda/bench_tiled_sweep_3d.cu`) across block sizes 4, 6, and 8, 10 trials each.

## Results

| Version              | Time (ms) | σ (ms) | Trials | Speedup vs. CPU |
|-----------------------|-----------|--------|--------|------------------|
| CPU (1 thread)        | 3216.47   | 210.47 | 10     | 1.00×            |
| GPU naive             | 111.37    | 15.35  | 10     | 28.9×            |
| GPU tiled (block 4)   | 100.09    | 16.77  | 10     | 32.1×            |

## Analysis

The 2D study's prediction held. 40×40×40 is 64,000 grid points against 100×100's 10,000, so each kernel launch here is doing over six times more work while paying the same fixed launch cost. The result is a 28-32× speedup instead of 2D's 1.32×, which lines up with the launch-overhead argument directly: more real work per launch means the fixed cost per launch matters less.

The CPU number carries more noise than the GPU numbers do — about 6.5% relative standard deviation against roughly 15% for the GPU runs in absolute terms, though the GPU's own run-to-run spread (99.5 ms to over 140 ms during warm-up-adjacent trials) is proportionally larger relative to its much smaller mean. Neither number is clean, but a 29-32× gap is well outside what noise this size could explain either way.

The block-size result is the part that doesn't have an obvious explanation yet: block 4 beat block 8, not the other way around. Bigger blocks usually amortize thread-launch and synchronization overhead better, so this looks backwards at first. The likely culprit is shared memory. The tiled kernel allocates `(blockSize + 2)³` doubles per block to hold the halo region — block 4 needs 216 doubles, block 8 needs 1000, almost five times as much. If that footprint limits how many blocks can sit resident on a single SM at once, block 8 could lose on occupancy even while doing more efficient halo reuse per block. That's a hypothesis, not something confirmed here — checking it properly would mean profiling actual occupancy, which Colab's default `nvidia-smi` doesn't expose without Nsight Compute.

## Future work

Profile actual SM occupancy for each block size (Nsight Compute or similar) to check whether the shared-memory-footprint explanation for block 4 beating block 8 actually holds. Run this same comparison at a second grid size — something between 40³ and 100³ — to see where the speedup curve starts to plateau as memory bandwidth rather than launch overhead becomes the bottleneck, the way the 2D study predicted it eventually would. Benchmark a multithreaded 3D CPU Jacobi under this same methodology, the same gap the 2D study left for its own threaded comparison. Re-run both CPU and GPU sides across separate sessions rather than one sitting, since the 2D study found CPU noise that only showed up across sessions and this study hasn't tested for that yet.