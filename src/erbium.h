#ifndef ERBIUM_H
#define ERBIUM_H

#include "model_specification.h"

const int n_ASE = 120;
//const int n_signal = 1;
const int n_signal = 120;
const int max_iter = 100;

void net_cross_section(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	int n = 100;
	double *x = new double [n],
		   *y = new double [n];
	
	double wl_min = 1450.,
		   wl_max = 1650.;
	double R = 0.0;

	for(int i = 0; i < n; i++)
	{
		double wl = wl_min + i * (wl_max - wl_min) / (double) (n - 1);
		double cs = R * E.F.cs_spline[1][0]->operator()(wl) - (1. - R) * E.F.cs_spline[0][1]->operator()(wl);
		x[i] = wl;
		y[i] = cs;
	}

	f_write_table(x, y, n, E.output_path + "net_cross_section.txt");
	delete [] x;
	delete [] y;
}

void gain_vs_L(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	double P_goal = 10.e-3,
		   P_start = 60.e-3;
	E.set_boundary("pump", P_goal);
	E.set_boundary("signal", 1.e-7);
	E.reset_channel("pump", 980.);
	E.reset_channel("signal", 1550.);
	E.update();	
	
	double L_min = 15.,
		   L_max = 25.;

	int n = 100;
	double *x = new double [n],
		   *y = new double [n];

	for(int i = 0; i < n; i++)
	{
		double L = L_min + i * (L_max - L_min) / double(n - 1);
		std::cout <<  L << "\n" ;
		E.set_fiber_parameters(L, 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		E.update();
		
		if(i > 0)
			E.relaxation(1.e-8, 1.e-8, 100, E);
		else
			E.relaxation();
		//int n_P = 10;
		//ErbiumModel<1 + 1, 0, n_ASE>::ModelSolution sol = 0.;
		//for(int j = 0; j < n_P; j++)
		//{
		//	double P_curr = P_start - j * (P_start - P_goal) / (double) (n_P - 1);
		//	std::cout <<  P_curr << "\n" ;
		//	E.set_boundary("pump", P_curr);
		//	E.relaxation(1.e-8, 1.e-8, 100, sol);
		//	sol = E;
		//}
			
		x[i] = L;
		y[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]);
	}
	f_write_table(x, y, n, E.output_path + "gain_vs_L.txt");
}

void gain_vs_P(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	E.set_boundary("signal", 1.e-7);
	E.reset_channel("pump", 980.);
	E.reset_channel("signal", 1530.);
	
	double P_min = 1.e-3,
		   P_max = 30.e-3;
	
	int n = 20;
	double *x = new double [n],
		   *y = new double [n];
	
	for(int i = 0; i < n; i++)
	{
		double P = P_min + i * (P_max - P_min) / double(n - 1);
		std::cout << P << "\n";

		E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		E.set_boundary("pump", P);
		E.update();

		E.relaxation();
		x[i] = P;
		y[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]);
		std::cout << y[i] << "\n";
	}
	f_write_table(x, y, n, E.output_path + "gain_vs_P.txt");
}

void total_ASE_vs_position(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	E.set_boundary("pump", 10.e-3);
	E.set_boundary("signal", 0.e-7);
	E.reset_channel("pump", 980.);
	E.reset_channel("signal", 1530.);
	E.update();

	E.relaxation();

	//int n = E.F.n_z;
	//double *x = new double [n],
	//	   *y = new double [n];

	//for(int i = 0; i < n; i ++)
	//{
	//	x[i] = E.F.z[i];
	//	
	//	y[i] = 0.;
	//	for(int j = 0; j < E.F.n_ch; j++)
	//		if(E.F.type[j] == ASE)
	//			y[i] += E.F.P[i][j];
	//}

	int n = E.B.n_z;
	double *x = new double [n],
		   *y = new double [n];

	for(int i = 0; i < n; i ++)
	{
		x[i] = E.B.z[i];
		
		y[i] = 0.;
		for(int j = 0; j < E.B.n_ch; j++)
			if(E.B.type[j] == ASE)
				y[i] += E.B.P[i][j];
	}

	f_write_table(x, y, n, E.output_path + "total_ASE_vs_position.txt");
}

void total_ASE_vs_P(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	E.set_boundary("signal", 0.e-7);
	E.reset_channel("pump", 1480.);
	E.reset_channel("signal", 1530.);
	E.update();

	double P_min = 1.e-3,
		   P_max = 50.e-3;
	
	int n = 30;
	double *x = new double [n],
		   *y = new double [n];
	
	for(int i = 0; i < n; i++)
	{
		double P = P_min + i * (P_max - P_min) / double(n - 1);
		std::cout << "\n" << P ;

		E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		E.set_boundary("pump", P);
		E.update();

		E.relaxation();
		x[i] = P;
		y[i] = 0.;
		
		for(int j = 0; j < E.F.n_ch; j++)
			if(E.F.type[j] == ASE)
				y[i] += E.F.P[E.F.n_z - 1][j];
		
		for(int j = 0; j < E.B.n_ch; j++)
			if(E.B.type[j] == ASE)
				y[i] += E.B.P[E.B.n_z - 1][j];

	}
	f_write_table(x, y, n, E.output_path + "total_ASE_vs_P.txt");
}

void gain_vs_position(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	E.set_boundary("pump", 4.e-3);
	E.set_boundary("signal", 1.e-7);
	E.reset_channel("pump", 980.);
	E.reset_channel("signal", 1550.);
	E.update();

	E.relaxation();

	int n = E.F.n_z;
	double *x = new double [n],
		   *y = new double [n];

	for(int i = 0; i < n; i ++)
	{
		x[i] = E.F.z[i];
		y[i] = 10 * log10(E.F.P[i][1] / E.F.P[0][1]);
	}

	f_write_table(x, y, n, E.output_path + "gain_position.txt");
}

void P_vs_position(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	E.set_boundary("pump", 4.e-3);
	E.set_boundary("signal", 1.e-7);
	E.reset_channel("pump", 980.);
	E.reset_channel("signal", 1550.);
	E.update();

	E.relaxation();

	int n = E.F.n_z;
	double *x = new double [n],
		   *y = new double [n];

	for(int i = 0; i < n; i ++)
	{
		x[i] = E.F.z[i];
		y[i] = E.F.P[i][0];
	}

	f_write_table(x, y, n, E.output_path + "P_vs_position.txt");
}

void gain_vs_output_power(ErbiumModel<1 + 1, 0, n_ASE> &E)
{
	double P = 15.e-3;
	E.set_boundary("pump", P);
	E.reset_channel("pump", 1480.);
	E.reset_channel("signal", 1550.);

	double P_min = 1.e-10,
		   P_max = P / 2.;
	
	int n = 20;
	double *x = new double [n],
		   *y = new double [n];
	
	for(int i = 0; i < n; i++)
	{
		double Ps = P_min * pow(10, i * (log10(P_max) - log10(P_min)) / double(n - 1));
		std::cout << 10 * log10(Ps / 1.e-3) << "\n";

		E.set_fiber_parameters(20., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		E.set_boundary("signal", Ps);
		E.update();

		E.relaxation();
		x[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / 1.e-3);
		y[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]);
		std::cout << y[i] << "\n";
	}
	f_write_table(x, y, n, E.output_path + "gain_vs_output_power.txt");
}

void quantum_efficiency_vs_P()
{
	typedef ErbiumModel<1 + n_signal, 1, n_ASE> Erbium_model;
	Erbium_model E;

	E.set_paths("e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Er_input\\",
		"e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\");
	
	E.add_mono_channel(980., forward, "F_pump");
	E.add_mono_channel(980., backward, "B_pump");
	E.add_mono_channel(1550., forward, "signal");
	E.add_ASE_channel(1500., 1580., n_ASE);
	
	double scale = 1.e-25;
	E.add_transition(0, 1, "ErAbsorption.txt", scale);
	E.add_transition(1, 0, "ErEmission.txt", scale);

	E.set_h(0, "h_Fiber_A.txt");
	
	E.set_boundary("signal", 1.e-5);
	
	double P_min = 1.e-3,
		   P_max = 90.e-3;
	
	int n = 30;
	double *x = new double [n],
		   *y = new double [n];
	
	for(int i = 0; i < n; i++)
	{
		double P = P_min + i * (P_max - P_min) / double(n - 1);
		std::cout << P << "\n";

		E.set_fiber_parameters(14., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		
		double R = 0.0;
		E.set_boundary("F_pump", R * P );
		E.set_boundary("B_pump", (1 - R) * P );
		
		E.update();
		
		E.relaxation(1.e-8, 1.e-8, 100, E);
		
		x[i] = P;
		y[i] = ( (E.F.P[E.F.n_z - 1][1] - E.F.P[0][1]) * E.F.wl[1]) / (E.F.P[0][0] * E.F.wl[0] + E.B.P[0][0] * E.B.wl[0]);
		std::cout << y[i] << "\n";
	}
	f_write_table(x, y, n, E.output_path + "quantum_efficiency_vs_P.txt");
}

void run_erbium_ESA()
{
	ErbiumESAModel<1 + n_signal, 0, n_ASE> E;
	E.set_paths("e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Er_input\\",
		"e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\");
	
	E.add_mono_channel(980., forward, "pump");
	//E.add_mono_channel(1550., forward, "signal");
	E.add_mono_channel(1500., 1580., n_signal, forward, "signal");
	E.add_ASE_channel(1500., 1580., n_ASE);
	
	double scale = 1.e-25;
	E.add_transition(0, 1, "ErAbsorption.txt", scale);
	E.add_transition(1, 0, "ErEmission.txt", scale);
	E.add_transition(1, 2, "ErESA.txt", 0.0 * scale);

	E.set_fiber_parameters(8., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	E.set_boundary("signal", 1.e-7 / n_signal);
    E.set_boundary("pump", 0.e-3);
	
	//double P_min = 2.e-3,
	//	   P_max = 50.e-3;
	//
	//int n = 60;
	//double *x = new double [n],
	//	   *y = new double [n];
	//
	//for(int i = 0; i < n; i++)
	//{
	//	double P = P_min + i * (P_max - P_min) / double(n - 1);
	//	std::cout << P << "\n" ;

	//	E.set_fiber_parameters(8., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	//	E.set_boundary("pump", P);
	//	E.update();

	//	E.relaxation(1.e-8, 1.e-8, 100, E);
	//	x[i] = P;
	//	y[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]);
	//}
	//f_write_table(x, y, n, E.output_path + "gain_vs_P_ESA.txt");
	
	E.update();
	E.relaxation();
	E.output(50, dB);
}

void run_erbium()
{
	typedef ErbiumModel<1 + n_signal, 1, n_ASE> Erbium_model;
	Erbium_model E;

	E.set_paths("e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Er_input\\",
		"e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\");
	
	E.add_mono_channel(1480., forward, "pump");
	E.add_mono_channel(980., backward, "pump2");
	//E.add_mono_channel(1550., forward, "signal");
	E.add_mono_channel(1500., 1580., n_signal, forward, "signal");
	E.add_ASE_channel(1500., 1580., n_ASE);
	
	double scale = 1.e-25;
	E.add_transition(0, 1, "ErAbsorption.txt", scale);
	E.add_transition(1, 0, "ErEmission.txt", scale);

	E.set_h(0, "h_Fiber_A.txt");

	E.set_fiber_parameters(8., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	E.set_boundary("pump", 0.e-3);
	E.set_boundary("pump2", 4.e-3);
	E.set_boundary("signal", 1.e-7 / n_signal);

	//net_cross_section(E);
	//gain_vs_L(E);
	//gain_vs_P(E);
	//gain_vs_position(E);
	//total_ASE_vs_P(E);
	//total_ASE_vs_position(E);
	//P_vs_position(E);
	//gain_vs_output_power(E);
	//quantum_efficiency_vs_P();
	//run_erbium_ESA();

	E.update();
	E.relaxation(1.e-8, 1.e-8, 100, E);

	//Erbium_model::ModelSolution sol1 = E; 
	//Erbium_model::ModelSolution sol2 = sol1;

	//E.relaxation(1.e-8, 1.e-8, 100, E);

	E.output(50, dB);
}


#endif