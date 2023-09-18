#pragma once
#include "BasicSolver2D.hpp"

namespace Engine {
	class Solver3Phases2D : public Engine::BasicSolver2D {
	public:
		Solver3Phases2D() = default;
		Solver3Phases2D(int Nx, int Ny, int numspec);
		virtual void set_initial_conditions() override;
		virtual void set_border_conditions() override;
		virtual void movement_step() override;
		virtual void calculate_moments() override;
		virtual void calculate_pressure() override;
		virtual void set_Phi() override;
		virtual void calculate_force() override;
		virtual void collision_step() override;
	private:
		double B = - 0.0000;
		double RR = 0.000;
		bool stream = false;
		std::vector<double> s_i;
		std::vector<double> rho_in;
		std::vector<double> rho_out;
		std::vector<std::vector<double>> k_ij;
		std::vector<std::vector<double>> pressure_water;
		std::vector<std::vector<double>> effrho_water;
		std::vector<std::vector<double>> sqr_effrho_water;

	};
}