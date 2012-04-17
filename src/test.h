#ifndef TEST_H
#define TEST_H

#include <iostream>

#include "interpolation.h"

using namespace fasim;

// interpolation ==================================================================================
void test_interpolation()
{
	std::string input_path = "e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Er_upconversion_input\\";
	std::string output_path = "e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\";
	std::string i_filename = "interp_test.txt";
	std::string o_filename = "interp_test.txt";

	double *x = new double [sufficiently_big_memory_chunk],
	       *y = new double [sufficiently_big_memory_chunk];
	int n;
	f_read_table(x, y, &n, input_path + i_filename);

	SplineInterpolation spline(x, y, n);

	int n_out = 71;
	double x_min = 0.e-9,
		   x_max = 7.e-9;

	double dx = (x_max - x_min) / (n_out - 1);

	for(int i = 0; i < n_out; i++)
	{
		x[i] = x_min + i * dx;
		y[i] = spline(x[i]);
	}

	f_write_table(x, y, n_out, output_path + o_filename);

	delete [] x;
	delete [] y;
}


// ODE ============================================================================================

typedef Vector<double, 2> V2;
typedef Vector<V2, 6> V2int;

typedef Vector<double, 3> V3;
typedef Vector<V3, 6> V3int;

typedef Vector<double, 3> V4Base;
typedef Vector<V4Base, 4> V4;
typedef Vector<V4, 6> V4int;

class ODETest
{
public:

	V2 func(double z, const V2 &y)
	{
		 V2 dy(y[1], (1 - y[0] * y[0]) * y[1] - y[0]);
		// V2 dy(y[1], 100*(1 - y[0] * y[0]) * y[1] - y[0]);
		// V2 dy(y[1] * 0.1 * z, (1 - y[0] * y[0]) * y[1] - y[0]);
		return dy;
	};

	V3 dN3(double t, const V3 &N)
	{
		double W01 = 430.,
			   W10 = 0.03,
			   W12 = 107.,
			   A1 = 1. / 10e-3;
		V3 dN;
		dN[0] = W10 * N[1] - W01 * N[0] + A1 * N[2];
		dN[1] = W01 * N[0] - W10 * N[1] - A1 * N[1] - W12 * N[1];
		dN[2] = -1. * (dN[0] + dN[1]);
		return dN;
	}

	V4 dN(double z, const V4 &y)
	{
		V4 dy;
		
		double
			t1 = 1.e-3,
			t2 = 1.e-4,
			t3 = 0.5e-4,
			W01 = 1.e5,
			W12 = 1.e4,
			W13 = 1.e4,
			W23 = 1.e4;
		dy[1] = - 1. / t1 * y[1] + W01 * y[0] - (W12 + W13) * y[1];
		dy[2] = - 1. / t2 * y[2] + W12 * y[1] - W23 * y[2];
		dy[3] = - 1. / t3 * y[3] + W13 * y[1] + W23 * y[2];
		dy[0] = - 1. * (dy[1] + dy[2] + dy[3]);
		return dy;
	}
};


void test_VDP()
{
	std::string output_path = "e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\";
	std::string fname_all = "test_ODE_all.txt";
	std::string fname_dense = "test_ODE_dense.txt";
	
	V2 *y = new V2 [segment_size];
	double *z = new double [segment_size];
	V2int *dc = new V2int [segment_size];
	int grid_size = 0;

	ODETest A;
	//ODE_solve_RK45(A, &ODETest::func, 0, 400, V2(2., 0.), grid_size, z, y, 1.e-8, 1.e-8, 1.e-5, &dc);
	//ODE_solve_RK45(A, &ODETest::func, 20, 0, V2(-0.742528, 1.06929), grid_size, z, y, 1.e-9, 1.e-9, 1.e-5, &dc);
	ODE_solve_RK45(A, &ODETest::func, 20, 0, V2(2.00814976217545, -0.04250887527543), grid_size, z, y, 1.e-10, 1.e-10, 1.e-5, &dc);

	// dense output start
	V2 yi;
	double d_low = 2.5,
		d_high = 3.0;
	int d_n = 100;
	
	std::ofstream d_file_output(output_path + fname_dense);
	for(int i = 0; i < d_n; i++)
	{
		double x = d_low + i * (d_high - d_low) / (double)(d_n - 1);
		yi = ODE_interpolate(x, z, y, dc, grid_size);
		d_file_output << x << "\t" << yi;

		if (i != d_n - 1) d_file_output << "\n";
	}
	d_file_output.close();
	// dense output end

	std::ofstream file_output(output_path + fname_all);
	file_output.precision(12);
	for (int i = 0; i < grid_size; i++)
	{
		file_output << std::showpoint << z[i] <<"\t" << y[i];
		if (i != grid_size - 1) file_output << "\n";
	}	
	file_output.close();

	delete [] z;
	delete [] y;
	delete [] dc;
}


void test_dN()
{
	std::string output_path = "e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\";
	std::string fname_all = "test_ODE_all.txt";
	std::string fname_dense = "test_ODE_dense.txt";	
	
	V4 *y = new V4 [segment_size];
	double *z = new double [segment_size];
	V4int *dc = new V4int [segment_size];
	int grid_size = 0;

	double tMax = 5.e-4;
	double Ns =  1.e25;
	ODETest A;
	ODE_solve_RK45(A, &ODETest::dN, 0., infinity, V4(V4Base(Ns, Ns / 2., Ns / 4.), V4Base(0., 0., 0.) , V4Base(0., 0., 0.) , V4Base(0., 0., 0.)),
		grid_size, z, y, 1.e-9, 1.e-9, 1.e-10, &dc);

	std::cout << A.dN(0., y[grid_size - 1]) << "\n" ;
	
	// dense output start
	//V4 yi;
	//double d_low = 2.5,
	//	d_high = 3.0;
	//int d_n = 100;
	//
	//std::ofstream d_file_output(output_path + fname_dense);
	//for(int i = 0; i < d_n; i++)
	//{
	//	double x = d_low + i * (d_high - d_low) / (double)(d_n - 1);
	//	yi = ODE_interpolate(x, z, y, dc, grid_size);
	//	d_file_output << x << "\t" << yi;

	//	if (i != d_n - 1) d_file_output << "\n";
	//}
	//d_file_output.close();
	// dense output end

	std::ofstream file_output(output_path + fname_all);
	file_output.precision(12);
	for (int i = 0; i < grid_size; i++)
	{
		file_output << std::showpoint << z[i] <<"\t" << y[i];
		if (i != grid_size - 1) file_output << "\n";
	}	
	file_output.close();

	delete [] z;
	delete [] y;
	delete [] dc;
}

void test_dN3()
{
	std::string output_path = "e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\";
	std::string fname_all = "test_ODE_all.txt";
	std::string fname_dense = "test_ODE_dense.txt";	
	
	V3 *y = new V3 [segment_size];
	double *z = new double [segment_size];
	V3int *dc = new V3int [segment_size];
	int grid_size = 0;

	double tMax = 1.e-1;
	double Ns =  0.7e25;
	ODETest A;
	ODE_solve_RK45(A, &ODETest::dN3, 0., tMax, V3(0., 0., 0.),
		grid_size, z, y, 1.e-9, 1.e-9, 1.e-10, &dc);

	std::cout << A.dN3(0., y[grid_size - 1]) << "\n" ;

	V3 yi;
	double d_low = z[0],
		d_high = z[grid_size - 1];
	int d_n = 100;
	
	std::ofstream d_file_output(output_path + fname_dense);
	for(int i = 0; i < d_n; i++)
	{
		double x = d_low + i * (d_high - d_low) / (double)(d_n - 1);
		yi = ODE_interpolate(x, z, y, dc, grid_size);
		d_file_output << x << "\t" << yi;

		if (i != d_n - 1) d_file_output << "\n";
	}
	d_file_output.close();

	std::ofstream file_output(output_path + fname_all);
	file_output.precision(12);
	for (int i = 0; i < grid_size; i++)
	{
		file_output << std::showpoint << z[i] <<"\t" << y[i];
		if (i != grid_size - 1) file_output << "\n";
	}	
	file_output.close();

	delete [] z;
	delete [] y;
	delete [] dc;
}

#endif