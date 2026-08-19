all: test_solvers.exe omega_sweep.exe heat.exe validate.exe test_mt.exe bench_mt.exe test_omega_auto.exe test_advection.exe test_advection_exact.exe diffusion_study.exe test_grid.exe test_jacobi.exe test_diffusion.exe test_diffusion_exact.exe diffusion_convergence.exe test_black_scholes_exact.exe

test_solvers.exe: include/grid.hpp include/jacobi.hpp include/gauss_seidel.hpp include/sor.hpp include/jacobi_mt.hpp include/gauss_seidel_rb.hpp include/sor_rb.hpp tests/test_solvers.cpp
	g++ -std=c++20 -Wall -I include tests/test_solvers.cpp -o test_solvers.exe

test_grid.exe: include/grid.hpp tests/test_grid.cpp
	g++ -std=c++20 -Wall -I include tests/test_grid.cpp -o test_grid.exe

test_jacobi.exe: include/grid.hpp include/jacobi.hpp tests/test_jacobi.cpp
	g++ -std=c++20 -Wall -I include tests/test_jacobi.cpp -o test_jacobi.exe

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

jacobi-convergence: include/grid.hpp include/jacobi.hpp include/jacobi_tiled_kernel.cuh cuda/jacobi_convergence.cu
	nvcc -std=c++20 -I include cuda/jacobi_convergence.cu -o jacobi_convergence.exe
	./jacobi_convergence.exe

bench_cpu.exe: include/grid.hpp include/stats.hpp include/n_jacobi.hpp tests/bench_cpu.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/bench_cpu.cpp -o bench_cpu.exe
	./bench_cpu.exe

test_advection.exe: include/grid1d.hpp include/advection.hpp include/space_time_log.hpp tests/test_advection.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/test_advection.cpp -o test_advection.exe 

demo_cfl_violation.exe: include/grid1d.hpp include/advection.hpp tests/demo_cfl_violation.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/demo_cfl_violation.cpp -o demo_cfl_violation.exe 

test_advection_exact.exe: include/grid1d.hpp include/advection.hpp tests/test_advection_exact.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/test_advection_exact.cpp -o test_advection_exact.exe 

diffusion_study.exe: include/grid1d.hpp include/advection.hpp include/stats.hpp include/space_time_log.hpp include/ppm.hpp tests/diffusion_study.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/diffusion_study.cpp -o diffusion_study.exe

test_diffusion.exe: include/grid1d.hpp include/diffusion.hpp include/space_time_log.hpp include/ppm.hpp tests/test_diffusion.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/test_diffusion.cpp -o test_diffusion.exe

demo_r_violation.exe: include/grid1d.hpp include/diffusion.hpp include/space_time_log.hpp tests/demo_r_violation.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/demo_r_violation.cpp -o demo_r_violation.exe

test_diffusion_exact.exe: include/grid1d.hpp include/diffusion.hpp include/space_time_log.hpp include/ppm.hpp include/diffusion_exact.hpp tests/test_diffusion_exact.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/test_diffusion_exact.cpp -o test_diffusion_exact.exe

diffusion_convergence.exe: include/grid1d.hpp include/diffusion.hpp include/diffusion_exact.hpp include/space_time_log.hpp tests/diffusion_convergence.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/diffusion_convergence.cpp -o diffusion_convergence.exe

test_black_scholes_exact.exe: include/grid1d.hpp include/diffusion.hpp include/black_scholes.hpp include/black_scholes_exact.hpp include/space_time_log.hpp tests/test_black_scholes_exact.cpp
	g++ -std=c++20 -Wall -O2 -I include tests/test_black_scholes_exact.cpp -o test_black_scholes_exact.exe