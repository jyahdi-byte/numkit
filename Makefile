all: test_solvers.exe omega_sweep.exe heat.exe validate.exe test_mt.exe bench_mt.exe test_omega_auto.exe

test_solvers.exe: include/grid.hpp include/jacobi.hpp include/gauss_seidel.hpp include/sor.hpp tests/test_solvers.cpp
	g++ -std=c++20 -Wall -I include tests/test_solvers.cpp -o test_solvers.exe

omega_sweep.exe: include/grid.hpp include/sor.hpp tests/omega_sweep.cpp
	g++ -std=c++20 -Wall -I include tests/omega_sweep.cpp -o omega_sweep.exe

heat.exe: include/grid.hpp include/sor.hpp include/ppm.hpp apps/heat/main.cpp
	g++ -std=c++20 -Wall -I include apps/heat/main.cpp -o heat.exe

validate.exe: include/grid.hpp include/sor.hpp tests/validate.cpp
	g++ -std=c++20 -Wall -I include tests/validate.cpp -o validate.exe

test_mt.exe: include/grid.hpp include/jacobi.hpp include/jacobi_mt.hpp tests/test_mt.cpp
	g++ -std=c++20 -Wall -pthread -I include tests/test_mt.cpp -o test_mt.exe

bench_mt.exe: include/grid.hpp include/jacobi.hpp include/jacobi_mt.hpp tests/bench_mt.cpp
	g++ -std=c++20 -Wall -pthread -I include tests/bench_mt.cpp -o bench_mt.exe

test_omega_auto.exe: include/grid.hpp include/sor.hpp tests/test_omega_auto.cpp
	g++ -std=c++20 -Wall -I include tests/test_omega_auto.cpp -o test_omega_auto.exe

cuda-test: include/grid.hpp cuda/grid_transfer_test.cu
	nvcc -I include cuda/grid_transfer_test.cu -o grid_transfer_test.exe
	./grid_transfer_test.exe

jacobi-validate: include/grid.hpp include/jacobi.hpp cuda/jacobi_validate.cu
	nvcc -I include cuda/jacobi_validate.cu -o jacobi_validate.exe
	./jacobi_validate.exe

bench-gpu: include/grid.hpp include/jacobi_kernel.cuh include/stats.hpp cuda/bench_gpu.cu
	nvcc -I include -O3 -arch=sm_75 cuda/bench_gpu.cu -o bench_gpu.exe
	./bench_gpu.exe

jacobi-tiled-validate: include/grid.hpp include/jacobi.hpp include/jacobi_tiled_kernel.cuh cuda/jacobi_tiled_validate.cu
	nvcc -I include cuda/jacobi_tiled_validate.cu -o jacobi_tiled_validate.exe
	./jacobi_tiled_validate.exe

bench-tiled-sweep: include/grid.hpp include/jacobi_tiled_kernel.cuh include/stats.hpp cuda/bench_tiled_sweep.cu
	nvcc -I include -O3 -arch=sm_75 cuda/bench_tiled_sweep.cu -o bench_tiled_sweep.exe
	./bench_tiled_sweep.exe
