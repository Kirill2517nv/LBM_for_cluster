#pragma once
#include "BasicSolver2D.hpp"

namespace Engine {
	class SolverMRTDimensional2D : public Engine::BasicSolver2D{
	public:
		SolverMRTDimensional2D() = default;
		SolverMRTDimensional2D(int Nx, int Ny, int numspec);
		virtual void set_initial_conditions() override;
		virtual void set_border_conditions() override;
		virtual void movement_step() override;
		virtual void calculate_moments() override;
		virtual void calculate_pressure() override;
		virtual void set_Phi() override;
		virtual void calculate_force() override;
		virtual void collision_step() override;
	private:
		virtual void eq_func(double rho, double ux, double uy, double* f_eq) override;
		void eq_func_dimensional_moment(double rho, double ux, double uy, double* moment_eq);

		//std::vector<double> s_k = { 0, 1.63, 1.14, 1., 1.92, 1., 1.92, 1 / tau, 1 / tau };
		std::vector<double> s_k = { 0, 1 / tau, 1 / tau, 1., 1 / tau, 1., 1 / tau, 1 / tau, 1 / tau };

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