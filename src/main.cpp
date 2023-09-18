#include <iostream>
#include "Solvers/BasicSolver2D.hpp"
#include "Solvers/SolverMRTDimensional2D.hpp"
#include "Solvers/Solver3Phases2D.hpp"
#include "Solvers/SolverPVTsim2D.hpp"
#include "Solvers/SolverPVTsimDynamic2D.hpp"
using namespace Engine;

int main()
{
	int Time = 10000000;
	int Nx = 600;
	int	Ny = 100;
	int	numspec = 3;
	std::shared_ptr<BasicSolver2D> Solver;
	Solver = std::make_shared<SolverPVTsimDynamic2D>(Nx, Ny, numspec);
	for (int t = 0; t < Time; t++) 
	{
		std::cout << "time = " << t << std::endl;
		Solver->LBM_Step();
		if (t % 10 == 0)
		{
			Solver->SaveVTKFile(t);
		}
		if (t == 100)
		{
			Solver->set_stream(true);
		}
	}
	return 0;
}