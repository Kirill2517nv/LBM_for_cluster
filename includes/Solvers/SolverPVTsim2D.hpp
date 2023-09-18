#pragma once
#include "BasicSolver2D.hpp"

namespace Engine {
	class SolverPVTsim2D : public Engine::BasicSolver2D {
	public:
		SolverPVTsim2D() = default;
		SolverPVTsim2D(int Nx, int Ny, int numspec);
		virtual void set_initial_conditions() override;
		virtual void set_border_conditions() override;
		virtual void movement_step() override;
		virtual void calculate_moments() override;
		virtual void calculate_pressure() override;
		virtual void set_Phi() override;
		virtual void calculate_force() override;
		virtual void collision_step() override;

		virtual void check_rho() override;
	private:
		std::vector<double> mixture(double per1, double per2);
		virtual void eq_func(double rho, double ux, double uy, double* f_eq) override;
		void eq_func_dimensional_moment(double rho, double ux, double uy, double* moment_eq);
		void set_surface_force(int i, int j, double* surface_force);

		std::vector<double> s_k = { 1.0, 0.8, 0.8, 1., 1.1, 1., 1.1, 1.72, 1.72 };
		//std::vector<double> s_k = { 0, 1 / tau, 1 / tau, 1., 1 / tau, 1., 1 / tau, 1 / tau, 1 / tau };
		std::vector<double> w_k = { 0, 1 / 3., 1 / 3., 1. / 3., 1 / 3., 1. / 12., 1 / 12., 1 / 12., 1 / 12. };

		std::vector<double> s_i;
		std::vector<std::vector<double>> k_ij;

		double M[9][9] =
		{
			{1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 },
			{-4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0},
			{4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0},
			{0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0},
			{0.0, -2.0,  0.0,  2.0,  0.0,  1.0, -1.0, -1.0,  1.0},
			{0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0},
			{0.0,  0.0, -2.0,  0.0,  2.0,  1.0,  1.0, -1.0, -1.0},
			{0.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0},
			{0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0},
		};

		double M_1[9][9] =
		{
			{1.0 / 9.0 , -1.0 / 9.0 ,  1.0 / 9.0 ,    0.0     ,     0.0    ,    0.0    ,     0.0    ,    0.0 ,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0,  1.0 / 6.0 , -1.0 / 6.0 ,    0.0    ,     0.0    ,    0.25,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0,    0.0     ,     0.0    , 1.0 / 6.0 , -1.0 / 6.0 ,   -0.25,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0, -1.0 / 6.0 ,  1.0 / 6.0 ,    0.0    ,     0.0    ,    0.25,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0,    0.0     ,     0.0    , -1.0 / 6.0,  1.0 / 6.0 ,   -0.25,  0.0 },
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0,  1.0 / 6.0 ,  1.0 / 12.0,  1.0 / 6.0,  1.0 / 12.0,    0.0 ,  0.25},
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0, -1.0 / 6.0 , -1.0 / 12.0,  1.0 / 6.0,  1.0 / 12.0,    0.0 , -0.25},
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0, -1.0 / 6.0 , -1.0 / 12.0, -1.0 / 6.0, -1.0 / 12.0,    0.0 ,  0.25},
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0,  1.0 / 6.0 ,  1.0 / 12.0, -1.0 / 6.0, -1.0 / 12.0,    0.0 , -0.25},
		};

	};
}