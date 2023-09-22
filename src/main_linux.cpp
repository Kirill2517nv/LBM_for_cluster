#include <iostream>
#include <sstream>
#include <memory>
#include "Solvers/BasicSolver2D.hpp"
#include "Solvers/SolverMRTDimensional2D.hpp"
#include "Solvers/Solver3Phases2D.hpp"
#include "Solvers/SolverPVTsim2D.hpp"
#include "Solvers/SolverPVTsimDynamic2D.hpp"
using namespace Engine;

int main(int argc, char* argv[])
{
	const char* filenameini = "../../Init/ini.txt";
	std::stringstream folder;
	double kappa = ini_read<double>(filenameini, "kappa", 100);
	std::cout << kappa << std::endl;
	folder << "mkdir ../../VTK/kappa=" << kappa;
	std::system(folder.str().c_str());
	int Time = ini_read<int>(filenameini, "MAX_TIME", 1e7);
	int Nx = ini_read<int>(filenameini, "Nx", 100);
	int	Ny = ini_read<int>(filenameini, "Ny", 100);
	int	numspec = ini_read<int>(filenameini, "number_of_species", 2);
	int step_VTK = ini_read<int>(filenameini, "step_VTK", 1000);
	int start_pumping = ini_read<int>(filenameini, "start_pumping", 10000);
	std::shared_ptr<BasicSolver2D> Solver;
	Solver = std::make_shared<SolverPVTsimDynamic2D>(Nx, Ny, numspec);
	for (int t = 0; t < Time; t++) 
	{
		std::cout << "time = " << t << std::endl;
		Solver->LBM_Step();
		if (t % step_VTK == 0)
		{
			Solver->SaveVTKFile(t);
			std::cout << "time = " << t << std::endl;
		}
		if (t == start_pumping)
		{
			Solver->set_stream(true);
			std::cout << "time start pumping = " << t << std::endl;
		}
	}
	return 0;
}