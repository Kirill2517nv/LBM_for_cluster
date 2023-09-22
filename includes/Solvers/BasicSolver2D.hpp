#pragma once
#include <vector>
#include "IniRead/ini.hpp"
#include <cmath>

namespace Engine {
	class BasicSolver2D {

	public:
		BasicSolver2D(int Nx, int Ny, int numspec);
		BasicSolver2D(const BasicSolver2D& other) = delete;
		BasicSolver2D(BasicSolver2D&& other) = delete;
		virtual void SaveVTKFile(int tStep);
		virtual ~BasicSolver2D();

		virtual void LBM_Step();
		virtual void set_initial_conditions();
		virtual void set_border_conditions() = 0;
		virtual void movement_step() = 0;
		virtual void calculate_moments() = 0;
		virtual void calculate_pressure() = 0;
		virtual void set_Phi() = 0;
		virtual void calculate_force() = 0;
		virtual void collision_step() = 0;

		virtual void check_rho();

		const int get_Nx() { return mNx; }
		const int get_Ny() { return mNy; }
		const double get_k() { return k; }
		const double get_A() { return A2; }
		const double get_gamma0() { return gamma[0]; }
		const double get_gamma1() { return gamma[1]; }
		const double get_gamma2() { return gamma[2]; }
		const double get_max_rho() { return max_rho; }
		const double get_min_rho() { return min_rho; }
		std::vector<std::vector<std::vector<double>>> get_rhomulticomponent() { return rhomulticomponent; }
		std::vector<std::vector<double>> get_rho() { return rho; }
		std::vector<std::vector<double>> get_pressure() { return pressure; }
		std::vector<double> get_mol_fraction() { return mol_fraction; };
		std::vector<double> get_fraction() { return fraction; };
		std::vector<double> get_fraction_liq() { return fraction_liq; };
		std::vector<double> get_fraction_vap() { return fraction_vap; };
		const bool get_leak() { return leak; }
		const bool get_saveVTK() { return saveVTK; }
		const bool get_stream() { return stream; }
		const double get_leak_coef() { return leak_coef; }
		const double get_ful_rho() { return ful_rho; }
		const double get_g() { return g; }
		const double get_kappa() { return kappa; }
		const double get_radius_of_droplet() { return radius_of_droplet; }

		void set_kappa(double k1) { kappa = k1; }
		void set_g(double k1) { g = k1; }
		void set_k(float k1) { k = k1; }
		void set_A(float A) { A2 = A; }
		void set_gamma0(double gamm) { gamma[0] = gamm; }
		void set_gamma1(double gamm) { gamma[1] = gamm; }
		void set_gamma2(double gamm) { gamma[2] = gamm; }
		void set_leak(bool leak1) { leak = leak1; }
		void set_saveVTK(bool leak1) { saveVTK = leak1; }
		void set_leak_coef(double leak1) { leak_coef = leak1; }
		void set_stream(bool stream1) { stream = stream1; }

	protected:
		virtual void eq_func(double rho, double ux, double uy, double* f_eq);
		double h = 1.0e-6;
		double delta_t = 1.0e-9;
		double wettability1 = 1.03; // walls
		double wettability2 = 1.03; //pores
		double b0 = 0.07780669;
		double a0 = 0.4572793;
		double rho_mixture = 300.;
		double ful_rho = 0.;
		double leak_coef = 0.995;
		double g = 0.0e1;
		double k = 1.;
		double R = 8.31446; // [J/(mol*K)]
		double max_rho = 0;
		double min_rho = 1000;
		double A2 = -0.580;  // -0.8080 -0.1381 Coefficient for calculating forces Peng-Robinson
		double tau = 1.0; // relaxation time
		double temperature; // Temperature of fluid [K]
		int mNx; // Number of nodes along the x-axis
		int mNy; // Number of nodes along the y-axis
		int number_of_species; // Number of species on fluid
		std::vector<double> gamma; // Coefficient ??????????
		std::vector<double> omega; // Acentric factor of each component
		std::vector<double> critical_temperatures; // Critical temperatures of each component [K]
		std::vector<double> critical_pressures; // Critical pressures of each component [MPa]
		std::vector<double> critical_rho; // Critical densities of each component [kg/m^3]
		std::vector<double> molarmass; // Molar masses of each component [kg/mol]
		std::vector<std::vector<double>> gamma_multiply_rho; // need for calculate multicomponent forces
		std::vector<std::vector<double>> rho; // Total density of fluid [kg/m^3]
		std::vector<std::vector<double>> effrho; // need for calculate forces
		std::vector<std::vector<double>> sqr_effrho;
		std::vector<std::vector<int>> mask;
		std::vector<std::vector<double>> ux; // Projection of the total fluid velocity on the x-axis at each node
		std::vector<std::vector<double>> uy; // Projection of the total fluid velocity on the y-axis at each node
		std::vector<std::vector<double>> dux_force; // Projection of the total force on the x-axis at each node
		std::vector<std::vector<double>> duy_force; // Projection of the total force on the y-axis at each node
		std::vector<std::vector<std::vector<double>>> ux_spec; // Projection of the component of fluid velocity on the x-axis at each node
		std::vector<std::vector<std::vector<double>>> uy_spec; // Projection of the component of fluid velocity on the y-axis at each node
		std::vector<std::vector<std::vector<double>>> dux_force_spec; // Projection of the component force on the x-axis at each node
		std::vector<std::vector<std::vector<double>>> duy_force_spec; // Projection of the component force on the y-axis at each node
		std::vector<std::vector<std::vector<double>>> rhomulticomponent; // Component density[]
		std::vector<std::vector<std::vector<std::vector<double>>>> fmulticomponent; // Function of distribution
		std::vector<std::vector<double>> pressure; // Pressure
		std::vector<double> Ai; // coefficient in PR EoS
		std::vector<double> Bi; // coefficient in PR EoS
		const int dx[9] = { 0,   1, 0, -1,  0,   1, -1, -1,  1 };
		const int dy[9] = { 0,   0, 1,  0, -1,   1,  1, -1, -1 };
		const int opposite_directions[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
		const int x_invert[9] = { 0, 3, 2, 1, 4, 6, 5, 8, 7 };
		const double G[9] = { 1, 1, 1, 1, 1,   0.25,  0.25, 0.25, 0.25 };
		std::vector<double> fraction_liq;
		std::vector<double> fraction_vap;
		std::vector<double> fraction;
		std::vector<double> mol_fraction;
		std::vector<double> rho_in;
		std::vector<double> rho_out;
		double kappa = 0.4;
		double radius_of_droplet = 0.;
		bool leak = false;
		bool stream = false;
		bool saveVTK = true;


	private:

	};
}