#include "Solvers/Solver3Phases2D.hpp"
#include <fstream>

namespace Engine {
    Solver3Phases2D::Solver3Phases2D(int Nx, int Ny, int numspec) :
        BasicSolver2D(Nx, Ny, numspec),
        s_i(std::vector<double>(numspec)),
		pressure_water(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
		effrho_water(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
		sqr_effrho_water(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2)))
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


    void Solver3Phases2D::set_initial_conditions()
	{
		k = 0.01;
		double Temperature = 400. / 190.564; /*376*/
		/*пентан (omega = 0.2514), T cr = 469.7 K , Pcr = 3.3675 MPa , mol_mass = 0.07215 kg/mol, rho cr = 232 kg/m^3  */
		/*метан (omega = 0.01142), T cr = 190.564 K , Pcr = 4.5992 MPa , mol_mass = 0.016 kg/mol, rho cr = 162.66 kg/m^3  */
		/*вода (omega = 0.3443, T cr = 647.096 K , Pcr = 22.0640 MPa , mol_mass = 0.018 kg/mol, rho cr = 322 kg/m^3  */
		s_i = { -0.154, -0.04183, 0.0 };
		omega = { 0.01142, 0.2514, 0.3443 };
		critical_pressures = { 4.5992e6, 3.3675e6, 22.0640e6 }; /* [Pa] */
		critical_temperatures = { 190.564, 469.7, 647.096 }; /* [k] */
		critical_rho = { 162.66, 232, 322 }; /*[kg/m^3]*/
		//vector<double> molar_mass = { 0.016 / 0.016, 0.07215 / 0.016, 0.018 / 0.016 }; /*[kg/mol]*/
		molarmass = { 0.016, 0.07215, 0.018 }; /*[kg/mol]*/
		double rho_mixture = 232. / 162.66; // 231.3 400
		k_ij = { {0., 0.03, 0.005},
			{0.03, 0., 0.010},
			{0.005, 0.01, 0.} };
		rho_in = { 0.08, 0.068, 0.1486 };
		rho_out = { 0.054, 0.417, 6.429 };

		int height = 0; //30
		int start = 40;
		int thickness = 6;
		double rho0liq = 0.017;
		double rho0vap = 0.05;
		double rho0water = 0.054;
		double rho1liq = 3.11;
		double rho1vap = 0.068;
		double rho1water = 0.417;
		double rho2liq = 0.15;
		double rho2vap = 0.1486;
		double rho2water = 6.429;

		for (int i = 1; i <= mNx; i++)
		{
			for (int j = 1; j <= mNy; j++)
			{
				/*if ((j < height - thickness) && (i >= edge1) && (i <= edge2))
				{
					rhomulticomponent[0][i][j] = rho0liq;
					rhomulticomponent[1][i][j] = rho1liq;
					rhomulticomponent[2][i][j] = rho2liq;
				}
				if ((j >= height - thickness) && (j <= height + thickness) && (i >= edge1) && (i <= edge2))
				{
					rhomulticomponent[0][i][j] = rho0liq + (rho0water - rho0liq) * (j - height + thickness) / thickness;
					rhomulticomponent[1][i][j] = rho1liq + (rho1water - rho1liq) * (j - height + thickness) / thickness;
					rhomulticomponent[2][i][j] = rho2liq + (rho2water - rho2liq) * (j - height + thickness) / thickness;
				}*/
				if ((j > height) && (i < start - thickness / 2))
				{
					rhomulticomponent[0][i][j] = rho0vap;
					rhomulticomponent[1][i][j] = rho1vap;
					rhomulticomponent[2][i][j] = rho2vap;
				}
				if ((j > height) && (i >= start - thickness / 2) && (i <= start + thickness / 2))
				{
					rhomulticomponent[0][i][j] = rho0vap + (rho0water - rho0vap) * (i - start + thickness / 2) / thickness;
					rhomulticomponent[1][i][j] = rho1vap + (rho1water - rho1vap) * (i - start + thickness / 2) / thickness;
					rhomulticomponent[2][i][j] = rho2vap + (rho2water - rho2vap) * (i - start + thickness / 2) / thickness;
				}
				if ((j > height) && (i > start + thickness / 2))
				{
					rhomulticomponent[0][i][j] = rho0water;
					rhomulticomponent[1][i][j] = rho1water;
					rhomulticomponent[2][i][j] = rho2water;
				}
			}
		}

		/* капля воды посреди углеводорода*/
		//for (int i = 1; i <= nx; i++)
		//{
		//	for (int j = 1; j <= ny; j++)
		//	{
		//		if (j < ny / 2 - thickness / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0liq;
		//			rhomulticomponent[1][i][j] = rho1liq;
		//			rhomulticomponent[2][i][j] = rho2liq;
		//		}
		//		if ((j >= ny / 2 - thickness / 2) && (j <= ny / 2 + thickness / 2))
		//			{
		//				rhomulticomponent[0][i][j] = rho0liq + (rho0vap - rho0liq) * (j - ny / 2 + thickness / 2) / thickness;
		//				rhomulticomponent[1][i][j] = rho1liq + (rho1vap - rho1liq) * (j - ny / 2 + thickness / 2) / thickness;
		//				rhomulticomponent[2][i][j] = rho2liq + (rho2vap - rho2liq) * (j - ny / 2 + thickness / 2) / thickness;
		//			}
		//		if ((j > ny / 2 + thickness / 2))
		//		{
		//			rhomulticomponent[0][i][j] = rho0vap;
		//			rhomulticomponent[1][i][j] = rho1vap;
		//			rhomulticomponent[2][i][j] = rho2vap;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= radius * radius)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water;
		//			rhomulticomponent[1][i][j] = rho1water;
		//			rhomulticomponent[2][i][j] = rho2water;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 2) * (radius + 2) && 
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > radius * radius && 
		//			j >= ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0vap - rho0water) * 0.2;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1vap - rho1water) * 0.2;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2vap - rho2water) * 0.2;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 2) * (radius + 2) &&
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > radius * radius &&
		//			j < ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0liq - rho0water) * 0.2;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1liq - rho1water) * 0.2;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2liq - rho2water) * 0.2;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 4) * (radius + 4) &&
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > (radius + 2) * (radius + 2) &&
		//			j >= ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0vap - rho0water) * 0.4;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1vap - rho1water) * 0.4;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2vap - rho2water) * 0.4;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 4) * (radius + 4) &&
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > (radius + 2) * (radius + 2) &&
		//			j < ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0liq - rho0water) * 0.4;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1liq - rho1water) * 0.4;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2liq - rho2water) * 0.4;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 6) * (radius + 6) && 
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > (radius + 4) * (radius + 4) &&
		//			j >= ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0vap - rho0water) * 0.6;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1vap - rho1water) * 0.6;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2vap - rho2water) * 0.6;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 6) * (radius + 6) &&
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > (radius + 4) * (radius + 4) &&
		//			j < ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0liq - rho0water) * 0.6;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1liq - rho1water) * 0.6;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2liq - rho2water) * 0.6;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 8) * (radius + 8) && 
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > (radius + 6) * (radius + 6) &&
		//			j >= ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0vap - rho0water) * 0.8;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1vap - rho1water) * 0.8;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2vap - rho2water) * 0.8;
		//		}
		//		if (((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) <= (radius + 8) * (radius + 8) &&
		//			((i - x_centre) * (i - x_centre) + (j - y_centre) * (j - y_centre)) > (radius + 6) * (radius + 6) &&
		//			j < ny / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0water + (rho0liq - rho0water) * 0.8;
		//			rhomulticomponent[1][i][j] = rho1water + (rho1liq - rho1water) * 0.8;
		//			rhomulticomponent[2][i][j] = rho2water + (rho2liq - rho2water) * 0.8;
		//		}
		//	}
		//}


		/*3 strips*/
		//for (int i = 1; i <= nx; i++)
		//{
		//	for (int j = 1; j <= ny; j++)
		//	{
		//		if (j < edge1 - thickness / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0up;
		//			rhomulticomponent[1][i][j] = rho1up;
		//			rhomulticomponent[2][i][j] = rho2down;
		//		}
		//		if ((j >= edge1 - thickness / 2) && (j <= edge1 + thickness / 2))
		//		{
		//			/*rhomulticomponent[0][i][j] = rho0down + (-rho0down + rho0up) * (1 + sin((i - edge1) * M_PI / thickness)) / 2;
		//			rhomulticomponent[1][i][j] = rho1down + (-rho1down + rho1up) * (1 + sin((i - edge1) * M_PI / thickness)) / 2;*/
		//			rhomulticomponent[0][i][j] = rho0up + (rho0down - rho0up) * (j - edge1 + thickness / 2) / thickness;
		//			rhomulticomponent[1][i][j] = rho1up + (rho1down - rho1up) * (j - edge1 + thickness / 2) / thickness;
		//			rhomulticomponent[2][i][j] = rho2down;
		//		}
		//		if ((j >= edge1 + thickness / 2) && (j <= edge1 + 3 * thickness / 2))
		//		{
		//			rhomulticomponent[2][i][j] = rho2down + (rho2up - rho2down) * (j - edge1 - thickness / 2) / thickness;
		//		}
		//		if ((j > edge1 + 3 * thickness / 2) && (j < edge2 - thickness / 2))
		//		{
		//			rhomulticomponent[2][i][j] = rho2up;
		//		}
		//		if ((j > edge1 + thickness / 2) && (j < edge2 - thickness / 2))
		//		{
		//			rhomulticomponent[0][i][j] = rho0down;
		//			rhomulticomponent[1][i][j] = rho1down;
		//		}
		//		if ((j >= edge2 - thickness / 2) && (j <= edge2 + thickness / 2))
		//		{
		//			rhomulticomponent[0][i][j] = rho0down;
		//			rhomulticomponent[1][i][j] = rho1down;
		//			rhomulticomponent[2][i][j] = rho2up + (rho2down - rho2up) * (j - edge2 + thickness / 2) / thickness;
		//		}
		//		if (j > edge2 + thickness / 2)
		//		{
		//			rhomulticomponent[0][i][j] = rho0down;
		//			rhomulticomponent[1][i][j] = rho1down;
		//			rhomulticomponent[2][i][j] = rho2down;
		//		}
		//	}
		//}
	}

	void Solver3Phases2D::set_border_conditions()
	{
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			double f_eq_in[9];
			double f_eq_out[9];
			for (int i = 1; i < mNy + 1; i++)
			{
				if (stream)
				{
					/*fmulticomponent[numspec][1][0][i] = leaks * fmulticomponent[numspec][3][1][i];
					fmulticomponent[numspec][3][mNx + 1][i] = fmulticomponent[numspec][3][mNx][i];
					fmulticomponent[numspec][5][0][i - 1] = leaks * fmulticomponent[numspec][7][1][i];
					fmulticomponent[numspec][8][0][i + 1] = leaks * fmulticomponent[numspec][6][1][i];
					fmulticomponent[numspec][6][mNx + 1][i - 1] = fmulticomponent[numspec][6][mNx][i];
					fmulticomponent[numspec][7][mNx + 1][i + 1] = fmulticomponent[numspec][7][mNx][i];*/
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
						if (mask[i + 1][j] == 0)
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
						}
					}
				}
			}
		}
	}


	void Solver3Phases2D::movement_step()
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
							int i_adj = i - dx[k];
							int j_adj = j - dy[k];
							f_temp[i][j] = fmulticomponent[numspec][k][i_adj][j_adj];
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


	void Solver3Phases2D::calculate_moments()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				if (mask[i][j] == 0) {
					rho[i][j] = 0;
					gamma_multiply_rho[i][j] = 0;
					ux[i][j] = 0;
					uy[i][j] = 0;
					for (int numspec = 0; numspec < number_of_species - 1; numspec++)
					{
						if (numspec == 0)
						{
							rhomulticomponent[number_of_species - 1][i][j] = fmulticomponent[number_of_species - 1][0][i][j];
							for (int k = 1; k < 9; k++)
							{
								rhomulticomponent[number_of_species - 1][i][j] += fmulticomponent[number_of_species - 1][k][i][j];
							}
							ux_spec[number_of_species - 1][i][j] = (fmulticomponent[number_of_species - 1][1][i][j] + fmulticomponent[number_of_species - 1][5][i][j] +
								fmulticomponent[number_of_species - 1][8][i][j] - fmulticomponent[number_of_species - 1][3][i][j] -
								fmulticomponent[number_of_species - 1][6][i][j] - fmulticomponent[number_of_species - 1][7][i][j]);
							uy_spec[number_of_species - 1][i][j] = (fmulticomponent[number_of_species - 1][2][i][j] + fmulticomponent[number_of_species - 1][5][i][j] +
								fmulticomponent[number_of_species - 1][6][i][j] - fmulticomponent[number_of_species - 1][4][i][j] -
								fmulticomponent[number_of_species - 1][7][i][j] - fmulticomponent[number_of_species - 1][8][i][j]);
							ux_spec[number_of_species - 1][i][j] /= rhomulticomponent[number_of_species - 1][i][j];
							uy_spec[number_of_species - 1][i][j] /= rhomulticomponent[number_of_species - 1][i][j];
						}
						rhomulticomponent[numspec][i][j] = fmulticomponent[numspec][0][i][j];
						for (int k = 1; k < 9; k++)
						{
							rhomulticomponent[numspec][i][j] += fmulticomponent[numspec][k][i][j];
						}
						gamma_multiply_rho[i][j] += gamma[numspec] * rhomulticomponent[numspec][i][j];
						rho[i][j] += rhomulticomponent[numspec][i][j];
						ux_spec[numspec][i][j] = (fmulticomponent[numspec][1][i][j] + fmulticomponent[numspec][5][i][j] +
							fmulticomponent[numspec][8][i][j] - fmulticomponent[numspec][3][i][j] -
							fmulticomponent[numspec][6][i][j] - fmulticomponent[numspec][7][i][j]);
						uy_spec[numspec][i][j] = (fmulticomponent[numspec][2][i][j] + fmulticomponent[numspec][5][i][j] +
							fmulticomponent[numspec][6][i][j] - fmulticomponent[numspec][4][i][j] -
							fmulticomponent[numspec][7][i][j] - fmulticomponent[numspec][8][i][j]);
						ux[i][j] += ux_spec[numspec][i][j];
						uy[i][j] += uy_spec[numspec][i][j];
						ux_spec[numspec][i][j] /= rhomulticomponent[numspec][i][j];
						uy_spec[numspec][i][j] /= rhomulticomponent[numspec][i][j];
					}
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

	void Solver3Phases2D::calculate_pressure()
	{
		double Z = 0.307;
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				double D = 0;
				double B = 0;
				double A = 0;
				double S = 0;
				for (int numspec = 0; numspec < number_of_species - 1; numspec++)
				{
					D += rhomulticomponent[numspec][i][j] * molarmass[0] / molarmass[numspec];
					B += rhomulticomponent[numspec][i][j] * critical_rho[0] / critical_rho[numspec];
					S += rhomulticomponent[numspec][i][j] * critical_rho[0] / critical_rho[numspec] * s_i[numspec];
					for (int numspec2 = 0; numspec2 < number_of_species - 1; numspec2++)
					{
						A += rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec2][i][j] *
							sqrt(
								(critical_temperatures[numspec] / critical_temperatures[0] *
									molarmass[0] / molarmass[numspec] *
									critical_rho[0] / critical_rho[numspec]) *
								(critical_temperatures[numspec2] / critical_temperatures[0] *
									molarmass[0] / molarmass[numspec2] *
									critical_rho[0] / critical_rho[numspec2]) *
								a(temperature * critical_temperatures[0] / critical_temperatures[numspec], omega[numspec]) *
								a(temperature * critical_temperatures[0] / critical_temperatures[numspec2], omega[numspec2])
							) * (1 - k_ij[numspec][numspec2]);
					}

				}
				A *= 1.487;
				B *= 0.253;
				pressure[i][j] = 1. / Z * (D * temperature / (1 + S - B) - A / ((1 + S) * (1 + S) + 2 * B * (1 + S) - B * B));
				D = 0;
				B = 0;
				A = 0;
				S = 0;
				for (int numspec = number_of_species - 1; numspec < number_of_species; numspec++)
				{
					D += rhomulticomponent[numspec][i][j] * molarmass[0] / molarmass[numspec];
					B += rhomulticomponent[numspec][i][j] * critical_rho[0] / critical_rho[numspec];
					S += rhomulticomponent[numspec][i][j] * critical_rho[0] / critical_rho[numspec] * s_i[numspec];
					for (int numspec2 = number_of_species - 1; numspec2 < number_of_species; numspec2++)
					{
						A += rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec2][i][j] *
							sqrt(
								(critical_temperatures[numspec] / critical_temperatures[0] *
									molarmass[0] / molarmass[numspec] *
									critical_rho[0] / critical_rho[numspec]) *
								(critical_temperatures[numspec2] / critical_temperatures[0] *
									molarmass[0] / molarmass[numspec2] *
									critical_rho[0] / critical_rho[numspec2]) *
								a(temperature * critical_temperatures[0] / critical_temperatures[numspec], omega[numspec]) *
								a(temperature * critical_temperatures[0] / critical_temperatures[numspec2], omega[numspec2])
							) * (1 - k_ij[numspec][numspec2]);
					}
				}
				A *= 1.487;
				B *= 0.253;
				pressure_water[i][j] = 1. / Z * (D * temperature / (1 + S - B) - A / ((1 + S) * (1 + S) + 2 * B * (1 + S) - B * B));
			}
		}
		pressure[0][0] = pressure[mNx][mNy];
		pressure[mNx + 1][0] = pressure[1][mNy];
		pressure[mNx + 1][mNy + 1] = pressure[1][1];
		pressure[0][mNy + 1] = pressure[mNx][1];
		pressure_water[0][0] = pressure_water[mNx][mNy];
		pressure_water[mNx + 1][0] = pressure_water[1][mNy];
		pressure_water[mNx + 1][mNy + 1] = pressure_water[1][1];
		pressure_water[0][mNy + 1] = pressure_water[mNx][1];


		for (int i = 1; i < mNy + 1; i++)
		{
			pressure_water[0][i] = pressure_water[mNx][i];
			pressure_water[mNx + 1][i] = pressure_water[1][i];

			pressure[0][i] = pressure[mNx][i];
			pressure[mNx + 1][i] = pressure[1][i];
		}
		for (int i = 1; i < mNx + 1; i++)
		{
			pressure_water[i][0] = pressure_water[i][mNy];
			pressure_water[i][mNy + 1] = pressure_water[i][1];

			pressure[i][0] = pressure[i][mNy];
			pressure[i][mNy + 1] = pressure[i][1];
		}
	}


	void Solver3Phases2D::set_Phi()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				sqr_effrho[i][j] = rho[i][j] / 3.0 - k * pressure[i][j];
				effrho[i][j] = sqrt(sqr_effrho[i][j]);
				sqr_effrho_water[i][j] = rhomulticomponent[number_of_species - 1][i][j] / 3.0 - k * pressure_water[i][j];
				effrho_water[i][j] = sqrt(sqr_effrho_water[i][j]);
			}
		}
		effrho_water[0][mNy + 1] = wettability1 * effrho_water[1][mNy];
		effrho_water[mNx + 1][mNy + 1] = wettability1 * effrho_water[mNx][mNy];
		effrho_water[0][0] = wettability1 * effrho_water[1][1];
		effrho_water[mNx + 1][0] = wettability1 * effrho_water[mNx][1];
		sqr_effrho_water[0][0] = effrho_water[0][0] * effrho_water[0][0];
		sqr_effrho_water[0][mNy + 1] = effrho_water[0][mNy + 1] * effrho_water[0][mNy + 1];
		sqr_effrho_water[mNx + 1][0] = effrho_water[mNx + 1][0] * effrho_water[mNx + 1][0];
		sqr_effrho_water[mNx + 1][mNy + 1] = effrho_water[mNx + 1][mNy + 1] * effrho_water[mNx + 1][mNy + 1];

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
					effrho_water[i][j] = 0.;
					for (int k = 1; k < 9; k++)
					{
						if (mask[i - dx[k]][j - dy[k]] == 0)
						{
							effrho[i][j] += effrho[i - dx[k]][j - dy[k]] * G[k];
							effrho_water[i][j] += effrho_water[i - dx[k]][j - dy[k]] * G[k];
							sumG += G[k];
						}
					}
					if (sumG != 0)
					{
						effrho[i][j] *= wettability2 / sumG;
						sqr_effrho[i][j] = effrho[i][j] * effrho[i][j];
						effrho_water[i][j] *= wettability2 / sumG;
						sqr_effrho_water[i][j] = effrho_water[i][j] * effrho_water[i][j];
					}
				}
			}
		}
		for (int i = 1; i <= mNx; i++) {
			effrho_water[i][mNy + 1] = wettability1 * effrho_water[i][mNy];
			effrho_water[i][0] = wettability1 * effrho_water[i][1];
			sqr_effrho_water[i][mNy + 1] = effrho_water[i][mNy + 1] * effrho_water[i][mNy + 1];
			sqr_effrho_water[i][0] = effrho_water[i][0] * effrho_water[i][0];

			effrho[i][mNy + 1] = wettability1 * effrho[i][mNy];
			effrho[i][0] = wettability1 * effrho[i][1];
			sqr_effrho[i][mNy + 1] = effrho[i][mNy + 1] * effrho[i][mNy + 1];
			sqr_effrho[i][0] = effrho[i][0] * effrho[i][0];

		}
		for (int i = 1; i <= mNy; i++) {
			effrho_water[0][i] = effrho_water[1][i];
			effrho_water[mNx + 1][i] = effrho_water[mNx][i];
			sqr_effrho_water[0][i] = effrho_water[0][i] * effrho_water[0][i];
			sqr_effrho_water[mNx + 1][i] = effrho_water[mNx + 1][i] * effrho_water[mNx + 1][i];

			effrho[0][i] = effrho[1][i];
			effrho[mNx + 1][i] = effrho[mNx][i];
			sqr_effrho[0][i] = effrho[0][i] * effrho[0][i];
			sqr_effrho[mNx + 1][i] = effrho[mNx + 1][i] * effrho[mNx + 1][i];
		}
	}


	void Solver3Phases2D::calculate_force()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				if (mask[i][j] == 0) {
					double force_x = 0, forceB_x = 0, forceRR_x = 0;
					forceRR_x = -RR * rho[i][j] * rhomulticomponent[number_of_species - 1][i][j] *
						(ux[i][j] - ux_spec[number_of_species - 1][i][j]);
					forceB_x = B * pressure[i][j] * (abs(pressure_water[i + 1][j]) - abs(pressure_water[i - 1][j]) +
						0.25 * (abs(pressure_water[i + 1][j + 1]) + abs(pressure_water[i + 1][j - 1]) -
							abs(pressure_water[i - 1][j + 1]) - abs(pressure_water[i - 1][j - 1])));
					force_x = 2.0 / 3 * ((1 - 2 * A2) * effrho[i][j] * (effrho[i + 1][j] - effrho[i - 1][j] +
						0.25 * (effrho[i + 1][j + 1] + effrho[i + 1][j - 1] - effrho[i - 1][j + 1] - effrho[i - 1][j - 1])) +
						A2 * (sqr_effrho[i + 1][j] - sqr_effrho[i - 1][j] + 0.25 *
							(sqr_effrho[i + 1][j + 1] + sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j + 1] - sqr_effrho[i - 1][j - 1])));
					dux_force[i][j] = (force_x + forceB_x + forceRR_x) / rho[i][j];
					double force_y = 0, forceB_y = 0, forceRR_y = 0;
					forceRR_y = -RR * rho[i][j] * rhomulticomponent[number_of_species - 1][i][j] *
						(uy[i][j] - uy_spec[number_of_species - 1][i][j]);
					forceB_y = B * pressure[i][j] * (abs(pressure_water[i][j + 1]) - abs(pressure_water[i][j - 1]) +
						0.25 * (abs(pressure_water[i + 1][j + 1]) + abs(pressure_water[i - 1][j + 1]) -
							abs(pressure_water[i + 1][j - 1]) - abs(pressure_water[i - 1][j - 1])));
					force_y = 2.0 / 3 * ((1 - 2 * A2) * effrho[i][j] * (effrho[i][j + 1] - effrho[i][j - 1] +
						0.25 * (effrho[i + 1][j + 1] + effrho[i - 1][j + 1] - effrho[i + 1][j - 1] - effrho[i - 1][j - 1])) +
						A2 * (sqr_effrho[i][j + 1] - sqr_effrho[i][j - 1] + 0.25 *
							(sqr_effrho[i + 1][j + 1] + sqr_effrho[i - 1][j + 1] - sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j - 1])));
					duy_force[i][j] = (force_y + forceB_y + forceRR_y) / rho[i][j];

					double force_x_water = 0, forceB_x_water = 0, forceRR_x_water = 0;
					forceRR_x_water = RR * rho[i][j] * rhomulticomponent[number_of_species - 1][i][j] *
						(ux[i][j] - ux_spec[number_of_species - 1][i][j]);
					forceB_x_water = B * abs(pressure_water[i][j]) * (pressure[i + 1][j] - pressure[i - 1][j] +
						0.25 * (pressure[i + 1][j + 1] + pressure[i + 1][j - 1] -
							pressure[i - 1][j + 1] - pressure[i - 1][j - 1]));
					force_x_water = 2.0 / 3 * ((1 - 2 * A2) * effrho_water[i][j] * (effrho_water[i + 1][j] - effrho_water[i - 1][j] +
						0.25 * (effrho_water[i + 1][j + 1] + effrho_water[i + 1][j - 1] - effrho_water[i - 1][j + 1] - effrho_water[i - 1][j - 1])) +
						A2 * (sqr_effrho_water[i + 1][j] - sqr_effrho_water[i - 1][j] + 0.25 *
							(sqr_effrho_water[i + 1][j + 1] + sqr_effrho_water[i + 1][j - 1] - sqr_effrho_water[i - 1][j + 1] - sqr_effrho_water[i - 1][j - 1])));
					dux_force_spec[number_of_species - 1][i][j] = (force_x_water + forceB_x_water + forceRR_x_water)
						/ rhomulticomponent[number_of_species - 1][i][j];
					double force_y_water = 0, forceB_y_water = 0, forceRR_y_water = 0;
					forceRR_y_water = RR * rho[i][j] * rhomulticomponent[number_of_species - 1][i][j] *
						(uy[i][j] - uy_spec[number_of_species - 1][i][j]);
					forceB_y_water = B * abs(pressure_water[i][j]) * (pressure[i][j + 1] - pressure[i][j - 1] +
						0.25 * (pressure[i + 1][j + 1] + pressure[i - 1][j + 1] -
							pressure[i + 1][j - 1] - pressure[i - 1][j - 1]));
					force_y_water = 2.0 / 3 * ((1 - 2 * A2) * effrho_water[i][j] * (effrho_water[i][j + 1] - effrho_water[i][j - 1] +
						0.25 * (effrho_water[i + 1][j + 1] + effrho_water[i - 1][j + 1] - effrho_water[i + 1][j - 1] - effrho_water[i - 1][j - 1])) +
						A2 * (sqr_effrho_water[i][j + 1] - sqr_effrho_water[i][j - 1] + 0.25 *
							(sqr_effrho_water[i + 1][j + 1] + sqr_effrho_water[i - 1][j + 1] - sqr_effrho_water[i + 1][j - 1] - sqr_effrho_water[i - 1][j - 1])));
					duy_force_spec[number_of_species - 1][i][j] = (force_y_water + forceB_y_water + forceRR_y_water)
						/ rhomulticomponent[number_of_species - 1][i][j];
					for (int numspec = 0; numspec < number_of_species - 1; numspec++) {
						dux_force_spec[numspec][i][j] = (force_x + forceB_x + forceRR_x) * gamma[numspec] / gamma_multiply_rho[i][j];
						duy_force_spec[numspec][i][j] = (force_y + forceB_y + forceRR_y) * gamma[numspec] / gamma_multiply_rho[i][j];

						/*dux_force_spec[numspec][i][j] = (force_x) * hi[numspec][i][j];
						duy_force_spec[numspec][i][j] = (force_y) * hi[numspec][i][j];*/
					}
				}
			}
		}
	}


	void Solver3Phases2D::collision_step()
	{
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			if (numspec == number_of_species - 1)
			{
				double f_eq_spec[9], f_eq_spec_n[9];
				for (int i = 1; i <= mNx; i++)
				{
					for (int j = 1; j <= mNy; j++)
					{
						if (mask[i][j] == 0) {
							eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j], uy_spec[numspec][i][j], f_eq_spec);
							eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j] + dux_force_spec[numspec][i][j],
								uy_spec[numspec][i][j] + duy_force_spec[numspec][i][j], f_eq_spec_n);
							for (int k = 0; k < 9; k++)
							{
								fmulticomponent[numspec][k][i][j] += f_eq_spec_n[k] - f_eq_spec[k] + (f_eq_spec[k] - fmulticomponent[numspec][k][i][j]) / tau;
								//fmulticomponent[numspec][k][i][j] += (f_eq[k] - fmulticomponent[numspec][k][i][j]) / tau;
							}
						}
					}
				}
			}
			else
			{
				double f_eq[9], f_eq_spec[9], f_eq_spec_n[9];
				for (int i = 1; i <= mNx; i++)
				{
					for (int j = 1; j <= mNy; j++)
					{
						if (mask[i][j] == 0) {
							eq_func(rhomulticomponent[numspec][i][j], ux[i][j], uy[i][j], f_eq);
							eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j], uy_spec[numspec][i][j], f_eq_spec);
							eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j] + dux_force_spec[numspec][i][j],
								uy_spec[numspec][i][j] + duy_force_spec[numspec][i][j], f_eq_spec_n);
							for (int k = 0; k < 9; k++)
							{
								fmulticomponent[numspec][k][i][j] += f_eq_spec_n[k] - f_eq_spec[k] + (f_eq[k] - fmulticomponent[numspec][k][i][j]) / tau;
								//fmulticomponent[numspec][k][i][j] += (f_eq[k] - fmulticomponent[numspec][k][i][j]) / tau;
							}
						}
					}
				}
			}
		}
	}
}