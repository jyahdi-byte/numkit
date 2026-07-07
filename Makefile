all: test_solvers.exe omega_sweep.exe

test_solvers.exe: include/grid.hpp include/jacobi.hpp include/gauss_seidel.hpp include/sor.hpp tests/test_solvers.cpp
	g++ -std=c++17 -Wall -I include tests/test_solvers.cpp -o test_solvers.exe

omega_sweep.exe: include/grid.hpp include/sor.hpp tests/omega_sweep.cpp
	g++ -std=c++17 -Wall -I include tests/omega_sweep.cpp -o omega_sweep.exe