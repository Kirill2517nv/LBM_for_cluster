#include "Solvers/SolverPVTsimDynamic2D.hpp"

#include <iostream>

namespace Engine {


	SolverPVTsimDynamic2D::SolverPVTsimDynamic2D(int Nx, int Ny, int numspec) :
		BasicSolver2D(Nx, Ny, numspec),
		s_i(std::vector<double>(numspec))
	{

		set_initial_conditions();
		double f_eq[9];
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			for (int i = 1; i <= mNx; i++)
			{
				for (int j = 1; j <= mNy; j++)
				{
					eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j], uy_spec[numspec][i][j], f_eq);
					for (int k = 0; k < 9; k++)
					{
						fmulticomponent[numspec][k][i][j] = f_eq[k];
					}
				}
			}
		}
	}


	std::vector<double> SolverPVTsimDynamic2D::mixture(double per1, double per2)
	{
		double rho1, rho2, rho3;
		double a1 = molarmass[0] * molarmass[1] * per1 / (per1 * molarmass[0] + (1 - per1) * molarmass[1]);
		double a2 = (1 - per2) / molarmass[1] + per2 / molarmass[2] + a1 * ((1 - per2) / molarmass[1] / molarmass[2] +
			per2 / molarmass[0] / molarmass[2] - (1 - per2) / molarmass[1] / molarmass[1] - per2 / molarmass[0] / molarmass[1]);
		rho3 = 1 / a2 * rho_mixture * ((1 - per2) / molarmass[1] - a1 * ((1 - per2) / molarmass[1] / molarmass[1]
			+ per2 / molarmass[0] / molarmass[1]));
		rho1 = a1 * (rho_mixture / molarmass[1] + rho3 / molarmass[2] - rho3 / molarmass[1]);
		rho2 = rho_mixture - rho1 - rho3;
		return { rho1, rho2, rho3 };
	}

	void SolverPVTsimDynamic2D::set_initial_conditions()
	{
		//LOG_INFO("set_initial_conditions");
		temperature = 350.; /*376*/
		/*метан (omega = 0.01142), T cr = 190.564 K , Pcr = 4.5992 MPa , mol_mass = 0.016 kg/mol, rho cr = 162.66 kg/m^3  */
		/*бутан (omega = 0.201), T cr = 425.125 K , Pcr = 3.796 MPa , mol_mass = 0.05812 kg/mol, rho cr = 228 kg/m^3 */
		/*декан (omega = 0.4884), T cr = 617.7 K , Pcr = 2.103 MPa , mol_mass = 0.14229 kg/mol, rho cr = 233 kg/m^3  */
		/*пентан (omega = 0.2514), T cr = 469.7 K , Pcr = 3.3675 MPa , mol_mass = 0.07215 kg/mol, rho cr = 232 kg/m^3  */
		s_i = { -0.154, -0.06413, -0.082 };
		omega = { 0.01142, 0.201, 0.4884 };
		critical_pressures = { 4.5992e6, 3.796e6, 2.103e6 }; /* [Pa] */
		critical_temperatures = { 190.564, 425.125, 617.7 }; /* [k] */
		critical_rho = { 162.66, 228., 233. }; /*[kg/m^3]*/
		molarmass = { 0.016, 0.05812, 0.14229 }; /*[kg/mol]*/
		rho_mixture = 260.; // 231.3 400
		//ux_spec.swap(std::vector<std::vector<std::vector<double>>> (number_of_species, std::vector<std::vector<double>>(mNx + 2, std::vector<double>(mNy + 2, 0.05))));
		//uy_spec.swap(std::vector<std::vector<std::vector<double>>> (number_of_species, std::vector<std::vector<double>>(mNx + 2, std::vector<double>(mNy + 2, 0.05))));
		Bi = { b0 * R * critical_temperatures[0] / (critical_pressures[0] * molarmass[0]),
		   b0 * R * critical_temperatures[1] / (critical_pressures[1] * molarmass[1]) ,
		   b0 * R * critical_temperatures[2] / (critical_pressures[2] * molarmass[2]) };
		Ai = { R * critical_temperatures[0] / (sqrt(critical_pressures[0]) * molarmass[0]),
			   R * critical_temperatures[1] / (sqrt(critical_pressures[1]) * molarmass[1]),
			   R * critical_temperatures[2] / (sqrt(critical_pressures[2]) * molarmass[2]) };
		k_ij = { {0., 0.01, 0.045},
				 {0.01, 0., 0.005},
				 {0.045, 0.005, 0.} };

		rho_in = {139.3, 46.2, 5.9};
		rho_out = {97.6, 72.4, 253.5};
		gamma = { 0.28, 0.541 , 1.8 }; //{ 1., 0.482 }

		double rho0vap = 119.3; //  14.76;
		double rho1vap = 46.2; //	4.98 ;
		double rho2vap = 5.9; //	0.24 ;
		double rho0liq = 97.6; //	5.23 ;
		double rho1liq = 72.4; //	71.0 ;
		double rho2liq = 253.5; // 	514.5;


		//std::vector<double> temp(3, 0.0);
		//for (int i = 0; i <= mNx + 1; i++)
		//{
		//	for (int j = 0; j <= mNy + 1; j++)
		//	{
		//		/*temp = mixture(0.803, 0.13);
		//		for (int numspec = 0; numspec < number_of_species; numspec++)
		//		{
		//			rhomulticomponent[numspec][i][j] = temp[numspec] + 10 * (static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 0.5);
		//		}*/
		//		rhomulticomponent[0][i][j] = 123.;
		//		rhomulticomponent[1][i][j] = 90.;
		//		rhomulticomponent[2][i][j] = 70.;
		//	}
		//}

		int edge = 50;
		int delta = 10;
		for (int i = 1; i <= mNx; i++)
		{
			for (int j = 1; j <= mNy; j++)
			{
				if (i < edge - delta / 2)
				{
					rhomulticomponent[0][i][j] = rho0vap;
					rhomulticomponent[1][i][j] = rho1vap;
					rhomulticomponent[2][i][j] = rho2vap;
				}
				if ((i >= edge - delta / 2) && (i <= edge + delta / 2))
				{
					rhomulticomponent[0][i][j] = rho0vap + (rho0liq - rho0vap) * (i - edge + delta / 2) / delta;
					rhomulticomponent[1][i][j] = rho1vap + (rho1liq - rho1vap) * (i - edge + delta / 2) / delta;
					rhomulticomponent[2][i][j] = rho2vap + (rho2liq - rho2vap) * (i - edge + delta / 2) / delta;
				}
				if (i > edge + delta / 2)
				{
					rhomulticomponent[0][i][j] = rho0liq;
					rhomulticomponent[1][i][j] = rho1liq;
					rhomulticomponent[2][i][j] = rho2liq;
				}
			}
		}
	}


	void SolverPVTsimDynamic2D::set_border_conditions()
	{
		/*Periodic*/
		/*for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			fmulticomponent[numspec][5][0][0] = fmulticomponent[numspec][5][mNx][mNy];
			fmulticomponent[numspec][8][0][mNy + 1] = fmulticomponent[numspec][8][mNx][1];
			fmulticomponent[numspec][7][mNx + 1][mNy + 1] = fmulticomponent[numspec][7][1][1];
			fmulticomponent[numspec][6][mNx + 1][0] = fmulticomponent[numspec][6][1][mNy];
			for (int i = 1; i < mNy + 1; i++)
			{
				fmulticomponent[numspec][1][0][i] = fmulticomponent[numspec][1][mNx][i];
				fmulticomponent[numspec][3][mNx + 1][i] = fmulticomponent[numspec][3][1][i];
				fmulticomponent[numspec][5][0][i] = fmulticomponent[numspec][5][mNx][i];
				fmulticomponent[numspec][8][0][i] = fmulticomponent[numspec][8][mNx][i];
				fmulticomponent[numspec][6][mNx + 1][i] = fmulticomponent[numspec][6][1][i];
				fmulticomponent[numspec][7][mNx + 1][i] = fmulticomponent[numspec][7][1][i];
			}
			for (int i = 1; i < mNx + 1; i++)
			{
				fmulticomponent[numspec][2][i][0] = fmulticomponent[numspec][2][i][mNy];
				fmulticomponent[numspec][4][i][mNy + 1] = fmulticomponent[numspec][4][i][1];
				fmulticomponent[numspec][5][i][0] = fmulticomponent[numspec][5][i][mNy];
				fmulticomponent[numspec][6][i][0] = fmulticomponent[numspec][6][i][mNy];
				fmulticomponent[numspec][7][i][mNy + 1] = fmulticomponent[numspec][7][i][1];
				fmulticomponent[numspec][8][i][mNy + 1] = fmulticomponent[numspec][8][i][1];
			}
		}*/
		if (leak)
		{
			for (int numspec = 0; numspec < number_of_species; numspec++)
			{
				for (int i = 1; i <= mNx; i++)
				{
					for (int j = 1; j <= mNy; j++)
					{
						for (int k = 0; k < 9; k++)
						{
							fmulticomponent[numspec][k][i][j] = fmulticomponent[numspec][k][i][j] * leak_coef;
						}
					}
				}
			}
			leak = false;
		}


		/*Walls*/
		/*for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			for (int i = 1; i < mNy + 1; i++)
			{
				fmulticomponent[numspec][1][0][i] = fmulticomponent[numspec][3][1][i];
				fmulticomponent[numspec][3][mNx + 1][i] = fmulticomponent[numspec][1][mNx][i];
				fmulticomponent[numspec][5][0][i - 1] = fmulticomponent[numspec][7][1][i];
				fmulticomponent[numspec][8][0][i + 1] = fmulticomponent[numspec][6][1][i];
				fmulticomponent[numspec][6][mNx + 1][i - 1] = fmulticomponent[numspec][8][mNx][i];
				fmulticomponent[numspec][7][mNx + 1][i + 1] = fmulticomponent[numspec][5][mNx][i];
			}
			for (int i = 1; i < mNx + 1; i++)
			{
				fmulticomponent[numspec][2][i][0] = fmulticomponent[numspec][4][i][1];
				fmulticomponent[numspec][4][i][mNy + 1] = fmulticomponent[numspec][2][i][mNy];
				fmulticomponent[numspec][5][i - 1][0] = fmulticomponent[numspec][7][i][1];
				fmulticomponent[numspec][6][i + 1][0] = fmulticomponent[numspec][8][i][1];
				fmulticomponent[numspec][7][i + 1][mNy + 1] = fmulticomponent[numspec][5][i][mNy];
				fmulticomponent[numspec][8][i - 1][mNy + 1] = fmulticomponent[numspec][6][i][mNy];
			}
		}*/

		/*Stream*/
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			double f_eq_in[9];
			double f_eq_out[9];
			for (int i = 1; i < mNy + 1; i++)
			{
				if (stream)
				{
					if (mask[0][i] == 0)
					{
						eq_func(rho_in[numspec] - rhomulticomponent[numspec][1][i], 0.0, 0.0, f_eq_in);
						eq_func(rho_out[numspec] - rhomulticomponent[numspec][mNx][i], 0.0, 0.0, f_eq_out);
						for (int k = 0; k < 9; k++)
						{
							fmulticomponent[numspec][k][0][i] = fmulticomponent[numspec][x_invert[k]][1][i] + f_eq_in[k];
							fmulticomponent[numspec][k][mNx + 1][i] = fmulticomponent[numspec][x_invert[k]][mNx][i] + f_eq_out[k];
							//fmulticomponent[numspec][k][mNx + 1][i] = fmulticomponent[numspec][k][mNx][i];
						}
					}
					else
					{
						fmulticomponent[numspec][1][0][i] = fmulticomponent[numspec][3][1][i];
						fmulticomponent[numspec][3][mNx + 1][i] = fmulticomponent[numspec][1][mNx][i];
						fmulticomponent[numspec][5][0][i] = fmulticomponent[numspec][7][1][i + 1];
						fmulticomponent[numspec][8][0][i] = fmulticomponent[numspec][6][1][i - 1];
						fmulticomponent[numspec][6][mNx + 1][i] = fmulticomponent[numspec][8][mNx][i + 1];
						fmulticomponent[numspec][7][mNx + 1][i] = fmulticomponent[numspec][5][mNx][i - 1];
					}
				}
				else
				{
					fmulticomponent[numspec][1][0][i] = fmulticomponent[numspec][3][1][i];
					fmulticomponent[numspec][3][mNx + 1][i] = fmulticomponent[numspec][1][mNx][i];
					fmulticomponent[numspec][5][0][i - 1] = fmulticomponent[numspec][7][1][i];
					fmulticomponent[numspec][8][0][i + 1] = fmulticomponent[numspec][6][1][i];
					fmulticomponent[numspec][6][mNx + 1][i - 1] = fmulticomponent[numspec][8][mNx][i];
					fmulticomponent[numspec][7][mNx + 1][i + 1] = fmulticomponent[numspec][5][mNx][i];
				}
			}
			for (int i = 1; i < mNx + 1; i++)
			{
				fmulticomponent[numspec][2][i][0] = fmulticomponent[numspec][4][i][1];
				fmulticomponent[numspec][4][i][mNy + 1] = fmulticomponent[numspec][2][i][mNy];
				fmulticomponent[numspec][5][i - 1][0] = fmulticomponent[numspec][7][i][1];
				fmulticomponent[numspec][6][i + 1][0] = fmulticomponent[numspec][8][i][1];
				fmulticomponent[numspec][7][i + 1][mNy + 1] = fmulticomponent[numspec][5][i][mNy];
				fmulticomponent[numspec][8][i - 1][mNy + 1] = fmulticomponent[numspec][6][i][mNy];
				for (int j = 1; j < mNy + 1; j++)
				{
					if (mask[i][j] == 1)
					{
						for (int k = 0; k < 9; k++)
						{
							if (mask[i + dx[k]][j + dy[k]] == 0)
							{
								fmulticomponent[numspec][k][i][j] = fmulticomponent[numspec][opposite_directions[k]][i + dx[k]][j + dy[k]];
							}
						}
		/*				if (mask[i + 1][j] == 0)
						{
							fmulticomponent[numspec][1][i][j] = fmulticomponent[numspec][3][i + 1][j];
						}
						if (mask[i - 1][j] == 0)
						{
							fmulticomponent[numspec][3][i][j] = fmulticomponent[numspec][1][i - 1][j];
						}
						if (mask[i][j + 1] == 0)
						{
							fmulticomponent[numspec][2][i][j] = fmulticomponent[numspec][4][i][j + 1];
						}
						if (mask[i][j - 1] == 0)
						{
							fmulticomponent[numspec][4][i][j] = fmulticomponent[numspec][2][i][j - 1];
						}
						if (mask[i - 1][j + 1] == 0)
						{
							fmulticomponent[numspec][6][i][j] = fmulticomponent[numspec][8][i - 1][j + 1];
						}
						if (mask[i - 1][j - 1] == 0)
						{
							fmulticomponent[numspec][7][i][j] = fmulticomponent[numspec][5][i - 1][j - 1];
						}
						if (mask[i + 1][j - 1] == 0)
						{
							fmulticomponent[numspec][8][i][j] = fmulticomponent[numspec][6][i + 1][j - 1];
						}
						if (mask[i + 1][j + 1] == 0)
						{
							fmulticomponent[numspec][5][i][j] = fmulticomponent[numspec][7][i + 1][j + 1];
						}*/
					}
				}
			}
		}
	}


	void SolverPVTsimDynamic2D::movement_step()
	{
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			std::vector<std::vector<double>> f_temp(mNx + 2, std::vector<double>(mNy + 2));
			for (int k = 1; k < 9; k++)
			{
				for (int i = 1; i <= mNx; i++)
				{
					for (int j = 1; j <= mNy; j++)
					{
						if (mask[i][j] == 0)
						{
							f_temp[i][j] = fmulticomponent[numspec][k][i - dx[k]][j - dy[k]];
						}
						else
						{
							f_temp[i][j] = fmulticomponent[numspec][k][i][j];
						}

					}
				}
				fmulticomponent[numspec][k].swap(f_temp);
			}
		}
	}


	void SolverPVTsimDynamic2D::calculate_moments()
	{
		ful_rho = 0;
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				if (mask[i][j] == 0)
				{
					rho[i][j] = 0;
					gamma_multiply_rho[i][j] = 0;
					ux[i][j] = 0;
					uy[i][j] = 0;
					for (int numspec = 0; numspec < number_of_species; numspec++)
					{
						rhomulticomponent[numspec][i][j] = 0;
						ux_spec[numspec][i][j] = 0;
						uy_spec[numspec][i][j] = 0;
						for (int k = 0; k < 9; k++)
						{
							rhomulticomponent[numspec][i][j] += fmulticomponent[numspec][k][i][j];
							ux_spec[numspec][i][j] += M[3][k] * fmulticomponent[numspec][k][i][j];
							uy_spec[numspec][i][j] += M[5][k] * fmulticomponent[numspec][k][i][j];

						}
						gamma_multiply_rho[i][j] += gamma[numspec] * rhomulticomponent[numspec][i][j];
						rho[i][j] += rhomulticomponent[numspec][i][j];
						ux_spec[numspec][i][j] *= h / delta_t;
						uy_spec[numspec][i][j] *= h / delta_t;
						ux[i][j] += ux_spec[numspec][i][j];
						uy[i][j] += uy_spec[numspec][i][j];
						ux_spec[numspec][i][j] /= rhomulticomponent[numspec][i][j];
						uy_spec[numspec][i][j] /= rhomulticomponent[numspec][i][j];

					}
					ful_rho += rho[i][j];
					ux[i][j] /= rho[i][j];
					uy[i][j] /= rho[i][j];
				}
			}
		}
	}

	double a(const double& temperature, double omega) // need for PengRobinson EOS
	{
		double m = omega <= 0.49 ? 0.37464 + 1.54226 * omega - 0.269922 * omega * omega :
			0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega; //Classic
		/*double m = omega <= 0.49 ? 0.382144 + 1.476905 * omega - 0.134488 * omega * omega :
			0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;*/  // Kalashnikov
		double a = pow((1 + m * (1 - sqrt(temperature))), 2);
		return a;
	}


	void SolverPVTsimDynamic2D::calculate_pressure()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++)
			{
				double D = 0;
				double B = 0;
				double A = 0;
				double S = 0;
				for (int numspec = 0; numspec < number_of_species; numspec++)
				{
					D += rhomulticomponent[numspec][i][j] / molarmass[numspec];
					B += rhomulticomponent[numspec][i][j] * Bi[numspec];
					S += s_i[numspec] * Bi[numspec] * rhomulticomponent[numspec][i][j];
					for (int numspec2 = 0; numspec2 < number_of_species; numspec2++)
					{
						A += rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec2][i][j] * Ai[numspec] * Ai[numspec2] *
							sqrt(a(temperature / critical_temperatures[numspec], omega[numspec]) *
								a(temperature / critical_temperatures[numspec2], omega[numspec2])
							) * (1 - k_ij[numspec][numspec2]);
					}

				}
				D *= R;
				A *= a0;
				pressure[i][j] = D * temperature / (1 + S - B) - A / ((1 + S) * (1 + S) + 2 * B * (1 + S) - B * B);
			}
		}
	}


	void SolverPVTsimDynamic2D::set_Phi()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				sqr_effrho[i][j] = rho[i][j] * h * h / (3.0 * delta_t * delta_t) - k * pressure[i][j];
				effrho[i][j] = sqrt(sqr_effrho[i][j]);
			}
		}

		effrho[0][mNy + 1] = wettability1 * effrho[1][mNy];
		effrho[mNx + 1][mNy + 1] = wettability1 * effrho[mNx][mNy];
		effrho[0][0] = wettability1 * effrho[1][1];
		effrho[mNx + 1][0] = wettability1 * effrho[mNx][1];
		sqr_effrho[0][0] = effrho[0][0] * effrho[0][0];
		sqr_effrho[0][mNy + 1] = effrho[0][mNy + 1] * effrho[0][mNy + 1];
		sqr_effrho[mNx + 1][0] = effrho[mNx + 1][0] * effrho[mNx + 1][0];
		sqr_effrho[mNx + 1][mNy + 1] = effrho[mNx + 1][mNy + 1] * effrho[mNx + 1][mNy + 1];

		for (int i = 1; i <= mNx; i++)
		{
			for (int j = 1; j <= mNy; j++)
			{
				if (mask[i][j] == 1)
				{
					double sumG = 0.;
					effrho[i][j] = 0.;
					for (int k = 1; k < 9; k++)
					{
						if (mask[i - dx[k]][j - dy[k]] == 0)
						{
							effrho[i][j] += effrho[i - dx[k]][j - dy[k]] * G[k];
							sumG += G[k];
						}
					}
					if (sumG != 0)
					{
						effrho[i][j] *= wettability2 / sumG;
						sqr_effrho[i][j] = effrho[i][j] * effrho[i][j];
					}
				}
			}
		}
		for (int i = 1; i <= mNx; i++) {
			effrho[i][mNy + 1] = wettability1 * effrho[i][mNy];
			effrho[i][0] = wettability1 * effrho[i][1];
			sqr_effrho[i][mNy + 1] = effrho[i][mNy + 1] * effrho[i][mNy + 1];
			sqr_effrho[i][0] = effrho[i][0] * effrho[i][0];
		}
		for (int i = 1; i <= mNy; i++) {
			effrho[0][i] = effrho[1][i];
			effrho[mNx + 1][i] = effrho[mNx][i];
			sqr_effrho[0][i] = effrho[0][i] * effrho[0][i];
			sqr_effrho[mNx + 1][i] = effrho[mNx + 1][i] * effrho[mNx + 1][i];
		}
	}


	void SolverPVTsimDynamic2D::calculate_force()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				if (mask[i][j] == 0)
				{
					double force_x = 0;
					force_x = 2.0 / 3 / h * delta_t * ((1 - 2 * A2) * effrho[i][j] * (effrho[i + 1][j] - effrho[i - 1][j] +
						0.25 * (effrho[i + 1][j + 1] + effrho[i + 1][j - 1] - effrho[i - 1][j + 1] - effrho[i - 1][j - 1])) +
						A2 * (sqr_effrho[i + 1][j] - sqr_effrho[i - 1][j] + 0.25 *
							(sqr_effrho[i + 1][j + 1] + sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j + 1] - sqr_effrho[i - 1][j - 1])));
					dux_force[i][j] = force_x / rho[i][j];
					double force_y = 0;
					force_y = 2.0 / 3 / h * delta_t * ((1 - 2 * A2) * effrho[i][j] * (effrho[i][j + 1] - effrho[i][j - 1] +
						0.25 * (effrho[i + 1][j + 1] + effrho[i - 1][j + 1] - effrho[i + 1][j - 1] - effrho[i - 1][j - 1])) +
						A2 * (sqr_effrho[i][j + 1] - sqr_effrho[i][j - 1] + 0.25 *
							(sqr_effrho[i + 1][j + 1] + sqr_effrho[i - 1][j + 1] - sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j - 1])));
					duy_force[i][j] = force_y / rho[i][j];
					for (int numspec = 0; numspec < number_of_species; numspec++) {
						dux_force_spec[numspec][i][j] = force_x * gamma[numspec] / gamma_multiply_rho[i][j];
						duy_force_spec[numspec][i][j] = force_y * gamma[numspec] / gamma_multiply_rho[i][j];
					}
				}
			}
		}
	}


	void SolverPVTsimDynamic2D::collision_step()
	{
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			double f_eq_spec[9], f_eq_spec_n[9];
			double mom[9], mom_eq[9];
			double surface_force[9];
			for (int i = 1; i <= mNx; i++)
			{
				for (int j = 1; j <= mNy; j++)
				{
					if (mask[i][j] == 0)
					{
						for (int k = 0; k < 9; k++)
						{
							mom[k] = 0;
							for (int z = 0; z < 9; z++) {
								mom[k] += M[k][z] * fmulticomponent[numspec][z][i][j];
							}
						}
						eq_func_dimensional_moment(rhomulticomponent[numspec][i][j], ux[i][j], uy[i][j], mom_eq);
						set_surface_force(i, j, surface_force);
						for (int k = 0; k < 9; k++)
						{
							mom[k] -= s_k[k] * (mom[k] - mom_eq[k]) + delta_t / h * delta_t / h * surface_force[k];
						}
						eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j], uy_spec[numspec][i][j], f_eq_spec);
						eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j] + dux_force_spec[numspec][i][j],
							uy_spec[numspec][i][j] + duy_force_spec[numspec][i][j], f_eq_spec_n);
						for (int k = 0; k < 9; k++)
						{
							fmulticomponent[numspec][k][i][j] = 0;
							for (int z = 0; z < 9; z++) {
								fmulticomponent[numspec][k][i][j] += M_1[k][z] * mom[z];
							}
							fmulticomponent[numspec][k][i][j] += f_eq_spec_n[k] - f_eq_spec[k];
						}
					}
				}
			}
		}
	}

	void SolverPVTsimDynamic2D::check_rho()
	{
		double full_rho = 0;
		double full_rho_nC5 = 0;
		double full_rho_C1 = 0;
		double full_rho_C2 = 0;
		std::vector<double> rhol(number_of_species, 0.);
		std::vector<double> rhov(number_of_species, 0.);
		std::vector<double> ful_rhol(number_of_species, 0.);
		std::vector<double> ful_rhov(number_of_species, 0.);
		int S_rhol = 0;
		int S_rhov = 0;
		double max_rho_temp = 0;
		double min_rho_temp = 1000;
		int max_rho_x = 0;
		int max_rho_y = 0;
		int min_rho_x = 0;
		int min_rho_y = 0;
		double r_liquid = 0;
		for (int numspec = 0; numspec < number_of_species; numspec++) {
			for (int i = 1; i <= mNx; i++)
			{
				for (int j = 1; j <= mNy; j++)
				{
					if (numspec == 0) {
						full_rho += rho[i][j];
						full_rho_nC5 += rhomulticomponent[1][i][j];
						full_rho_C1 += rhomulticomponent[numspec][i][j];
						full_rho_C2 += rhomulticomponent[2][i][j];
						if (rho[i][j] < min_rho_temp)
						{
							min_rho = rho[i][j];
							min_rho_temp = rho[i][j];
						}
						if (rho[i][j] > max_rho_temp)
						{
							max_rho = rho[i][j];
							max_rho_temp = rho[i][j];
						}
						if (max_rho_temp == rho[i][j])
						{
							max_rho_x = i;
							max_rho_y = j;
						}
						if (min_rho_temp == rho[i][j])
						{
							min_rho_x = i;
							min_rho_y = j;
						}
					}
				}
			}
		}
		for (int i = 1; i <= mNx; i++)
		{
			for (int j = 1; j <= mNy; j++)
			{
				if (rho[i][j] > (max_rho + min_rho) / 2)
				{
					ful_rhol[0] += rhomulticomponent[0][i][j];
					ful_rhol[1] += rhomulticomponent[1][i][j];
					ful_rhol[2] += rhomulticomponent[2][i][j];
					S_rhol++;
				}
				if (rho[i][j] <= (max_rho + min_rho) / 2)
				{
					ful_rhov[0] += rhomulticomponent[0][i][j];
					ful_rhov[1] += rhomulticomponent[1][i][j];
					ful_rhov[2] += rhomulticomponent[2][i][j];
					S_rhov++;
				}
			}
		}
		/*mole fraction mixture in vapour*/
		mol_fraction[0] = (ful_rhov[0] * S_rhov / molarmass[0] + ful_rhov[1] * S_rhov / molarmass[1] +
			ful_rhov[2] * S_rhov / molarmass[2]) /
			(ful_rhov[0] * S_rhov / molarmass[0] + ful_rhov[1] * S_rhov / molarmass[1] + ful_rhov[2] * S_rhov / molarmass[2] +
				ful_rhol[0] * S_rhol / molarmass[0] + ful_rhol[1] * S_rhol / molarmass[1] + ful_rhol[2] * S_rhol / molarmass[2]);
		mol_fraction[1] = 1. - mol_fraction[0];
		/*components*/
		fraction[0] = full_rho_C1 / molarmass[0] /
			(full_rho_C1 / molarmass[0] + full_rho_nC5 / molarmass[1] + full_rho_C2 / molarmass[2]);
		fraction[1] = full_rho_nC5 / molarmass[1] /
			(full_rho_C1 / molarmass[0] + full_rho_nC5 / molarmass[1] + full_rho_C2 / molarmass[2]);
		fraction[2] = 1 - fraction[0] - fraction[1];

		fraction_vap[0] = ful_rhov[0] * S_rhov / molarmass[0] /
			(ful_rhov[0] * S_rhov / molarmass[0] + ful_rhov[1] * S_rhov / molarmass[1] + ful_rhov[2] * S_rhov / molarmass[2]);
		fraction_vap[1] = ful_rhov[1] * S_rhov / molarmass[1] /
			(ful_rhov[0] * S_rhov / molarmass[0] + ful_rhov[1] * S_rhov / molarmass[1] + ful_rhov[2] * S_rhov / molarmass[2]);
		fraction_vap[2] = 1. - fraction_vap[0] - fraction_vap[1];

		fraction_liq[0] = ful_rhol[0] * S_rhol / molarmass[0] /
			(ful_rhol[0] * S_rhol / molarmass[0] + ful_rhol[1] * S_rhol / molarmass[1] + ful_rhol[2] * S_rhol / molarmass[2]);
		fraction_liq[1] = ful_rhol[1] * S_rhol / molarmass[1] /
			(ful_rhol[0] * S_rhol / molarmass[0] + ful_rhol[1] * S_rhol / molarmass[1] + ful_rhol[2] * S_rhol / molarmass[2]);
		fraction_liq[2] = 1. - fraction_liq[0] - fraction_liq[1];

	}


	void SolverPVTsimDynamic2D::eq_func(double rho, double ux, double uy, double* f_eq)
	{
		double du2 = 1.0 - 1.5 * (ux * ux * delta_t / h * delta_t / h + uy * uy * delta_t / h * delta_t / h);
		f_eq[0] = 4.0 / 9.0 * rho * du2;
		f_eq[1] = rho / 9.0 * (du2 + ux * delta_t / h * (3.0 + 4.5 * ux * delta_t / h));
		f_eq[2] = rho / 9.0 * (du2 + uy * delta_t / h * (3.0 + 4.5 * uy * delta_t / h));
		f_eq[3] = rho / 9.0 * (du2 - ux * delta_t / h * (3.0 - 4.5 * ux * delta_t / h));
		f_eq[4] = rho / 9.0 * (du2 - uy * delta_t / h * (3.0 - 4.5 * uy * delta_t / h));
		f_eq[5] = rho / 36.0 * (du2 + (ux * delta_t / h + uy * delta_t / h) * (3.0 + 4.5 * (ux * delta_t / h + uy * delta_t / h)));
		f_eq[6] = rho / 36.0 * (du2 + (uy * delta_t / h - ux * delta_t / h) * (3.0 + 4.5 * (uy * delta_t / h - ux * delta_t / h)));
		f_eq[7] = rho / 36.0 * (du2 - (ux * delta_t / h + uy * delta_t / h) * (3.0 - 4.5 * (ux * delta_t / h + uy * delta_t / h)));
		f_eq[8] = rho / 36.0 * (du2 + (ux * delta_t / h - uy * delta_t / h) * (3.0 + 4.5 * (ux * delta_t / h - uy * delta_t / h)));
	}


	void SolverPVTsimDynamic2D::eq_func_dimensional_moment(double rho, double ux, double uy, double* moment_eq)
	{
		moment_eq[0] = rho;
		moment_eq[1] = -2. * rho + 3. * rho * (ux * ux + uy * uy) * delta_t * delta_t / h / h;
		moment_eq[2] = rho - 3. * rho * (ux * ux + uy * uy) * delta_t * delta_t / h / h;
		moment_eq[3] = rho * ux * delta_t / h;
		moment_eq[4] = -rho * ux * delta_t / h;
		moment_eq[5] = rho * uy * delta_t / h;
		moment_eq[6] = -rho * uy * delta_t / h;
		moment_eq[7] = rho * (ux * ux - uy * uy) * delta_t * delta_t / h / h;
		moment_eq[8] = rho * ux * uy * delta_t * delta_t / h / h;
	}
	void SolverPVTsimDynamic2D::set_surface_force(int i, int j, double* surface_force)
	{
		double Qxx = 0.;
		double Qxy = 0.;
		double Qyy = 0.;

		for (int k = 1; k < 9; k++)
		{
			Qxx += 2. / 3. * w_k[k] * (effrho[i + dx[k]][j + dy[k]] - effrho[i][j]) * dx[k] * dx[k];
			Qxy += 2. / 3. * w_k[k] * (effrho[i + dx[k]][j + dy[k]] - effrho[i][j]) * dx[k] * dy[k];
			Qyy += 2. / 3. * w_k[k] * (effrho[i + dx[k]][j + dy[k]] - effrho[i][j]) * dy[k] * dy[k];
		}

		Qxx *= kappa * effrho[i][j];
		Qxy *= kappa * effrho[i][j];
		Qyy *= kappa * effrho[i][j];

		surface_force[0] = 0.;
		surface_force[1] = 1.5 * s_k[1] * (Qxx + Qyy);
		surface_force[2] = -1.5 * s_k[2] * (Qxx + Qyy);
		surface_force[3] = 0.;
		surface_force[4] = 0.;
		surface_force[5] = 0.;
		surface_force[6] = 0.;
		surface_force[7] = -s_k[7] * (Qxx - Qyy);
		surface_force[8] = -s_k[8] * (Qxy);
	}
}