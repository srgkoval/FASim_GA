#ifndef ERBIUM_UPCONVERSION
#define ERBIUM_UPCONVERSION

#include "model_specification.h"
#include "vector.h"
#include "ga.h"
#include "fiber.h"

const int n_signal_UC = 50; //50;
const int n_ASE_UC = 80;
const double N_Er = 6.4e26;

void gain_vs_L_and_P(ErbiumUpconversionModel<1 + 1 + 0, 0, n_ASE_UC> &E)
{
	E.set_boundary("signal", 1.e-7);
	E.reset_channel("pump", 1480.);
	E.reset_channel("signal", 551.5);
	E.update();	

	double P_min = 20.e-3,
		   P_max = 100.e-3,
		   L_min = 0.03,
		   L_max = 0.18;
	
	int P_points = 10;
	int L_points = 10;

	// prepare output
	std::string fname = "gain_vs_L_and_P.txt";
	
	std::ofstream g_out( (E.output_path + fname).c_str() );
	if (!g_out)
	{
		std::cout << "gain_vs_L_and_P: error opening " << fname <<"\n";
		wait_and_exit();
	}
	g_out.setf( std::ios::showpoint );
	g_out.precision(6);

	for(int i = 0; i < L_points; i++)
		g_out << L_min + i * (L_max - L_min) / (double) (L_points - 1) << (i != L_points - 1 ? "\t" : "\n");

	// main loop
	for(int p = 0; p < P_points; p++)
	{
		double P = P_min + p * (P_max - P_min) / (double) (P_points - 1);
		g_out << P;
		
		for(int l = 0; l < L_points; l++)
		{
			double L = L_min + l * (L_max - L_min) / (double) (L_points - 1);
			std::cout << "L = " << L << "\tP = " << P << "\nStep " << l + 1 + p * P_points << " of " << L_points * P_points << "\n";

			E.set_fiber_parameters(L, N_Er, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
			E.update();

			E.set_boundary("pump", P);

			E.relaxation(1.e-8, 1.e-8, 50);
			
			g_out << "\t" << 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]);
			std::cout << "\nG= " << 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]) << "\n";
		}
		if (p != P_points - 1 ) g_out << "\n";
	}

	g_out.close();
}

void gain_vs_P(ErbiumUpconversionModel<1 + 1, 0, n_ASE_UC> &E)
{
	E.set_boundary("signal", 1.e-7);
	E.reset_channel("pump", 1480.);
	E.reset_channel("signal", 551.5);
	
	double P_min = 20.e-3,
		   P_max = 200.e-3;
	
	int n = 20;
	double *x = new double [n],
		   *y = new double [n];
	
	for(int i = 0; i < n; i++)
	{
		double P = P_min + i * (P_max - P_min) / double(n - 1);
		std::cout << P << "\n";

		E.set_fiber_parameters(0.18, N_Er, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		E.set_boundary("pump", P);
		E.update();

		E.relaxation(1.e-8, 1.e-8, 50, E);
		x[i] = P;
		y[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]);
		std::cout << y[i] << "\n";
	}
	f_write_table(x, y, n, E.output_path + "gain_vs_P.txt");
}


void gain_vs_output_power(ErbiumUpconversionModel<1 + 1 + 0, 0, n_ASE_UC> &E)
{
	double P = 35.e-3;
	E.set_boundary("pump", P);
	E.reset_channel("pump", 980.);
	E.reset_channel("signal", 551.5);

	double P_min = 1.e-10,
		   P_max = P / 2.;
	
	int n = 20;
	double *x = new double [n],
		   *y = new double [n];
	
	for(int i = 0; i < n; i++)
	{
		double Ps = P_min * pow(10, i * (log10(P_max) - log10(P_min)) / double(n - 1));
		std::cout << 10 * log10(Ps / 1.e-3) << "\n";

		E.set_fiber_parameters(0.13, N_Er, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		E.set_boundary("signal", Ps);
		E.update();

		E.relaxation();
		x[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / 1.e-3);
		y[i] = 10 * log10(E.F.P[E.F.n_z - 1][1] / E.F.P[0][1]);
		std::cout << y[i] << "\n";
	}
	f_write_table(x, y, n, E.output_path + "gain_vs_output_power.txt");
}

void quantum_efficiency_vs_P_UC()
{
	typedef ErbiumUpconversionModel<1 + n_signal_UC, 1, n_ASE_UC> Erbium_model;
	Erbium_model E;

	E.set_paths("e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Er_upconversion_input_revised\\",
		"e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\");
	
	E.add_mono_channel(980., forward, "F_pump");
	E.add_mono_channel(980., backward, "B_pump");
	E.add_mono_channel(541., forward, "signal");
	E.add_ASE_channel(1500., 1620., 50);
	E.add_ASE_channel(515., 560., 50);
	
	double scale = 1.e-25;
	E.add_transition(0, 1, "cs01.txt", scale);
	E.add_transition(1, 0, "cs10.txt", scale);
	E.add_transition(0, 2, "cs02.txt", scale);
	E.add_transition(1, 3, "cs13.txt", scale);
	E.add_transition(2, 4, "cs24.txt", scale);

	double scale_upper = 10.e-25; //3.2e-25;
	//double scale_upper = 3.2e-25; 

	E.add_transition(0, 4, "cs04.txt", scale_upper);
	E.add_transition(4, 0, "cs40.txt", scale_upper);
	
	E.set_boundary("signal", 1.e-3);
	
	double P_min = 30.e-3,
		   P_max = 100.e-3;
	
	int n = 20;
	double *x = new double [n],
		   *y = new double [n];
	
	for(int i = 0; i < n; i++)
	{
		double P = P_min + i * (P_max - P_min) / double(n - 1);
		std::cout << P << "\n";

		E.set_fiber_parameters(0.06, N_Er, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
		
		double R = 0.5;
		E.set_boundary("F_pump", R * P );
		E.set_boundary("B_pump", (1 - R) * P );
		
		E.update();
		
		E.relaxation(1.e-8, 1.e-8, 50, E);
		
		x[i] = P;
		y[i] = ( (E.F.P[E.F.n_z - 1][1] - E.F.P[0][1]) * E.F.wl[1]) / (E.F.P[0][0] * E.F.wl[0] + E.B.P[0][0] * E.B.wl[0]);
		std::cout << y[i] << "\n";
	}
	f_write_table(x, y, n, E.output_path + "quantum_efficiency_vs_P.txt");
}


void run_erbium_upconversion_optim()
{
	ErbiumUpconversionModel<n_signal_UC + 1 + 0, 1, n_ASE_UC> E;
	
	E.set_paths("e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Er_upconversion_input_revised\\",
		"e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\");

	E.add_mono_channel(980., forward, "pump1");
	E.add_mono_channel(1480., backward, "pump2");

	//	E.add_mono_channel(537., 560., n_signal_UC, forward, "signal");
	E.add_mono_channel(541., 553., n_signal_UC, forward, "signal");

	E.add_ASE_channel(1500., 1620., 50);
	E.add_ASE_channel(515., 560., 30);

	double scale = 1.e-25;
	E.add_transition(0, 1, "cs01.txt", scale);
	E.add_transition(1, 0, "cs10.txt", scale);
	E.add_transition(0, 2, "cs02.txt", scale);
	E.add_transition(1, 3, "cs13.txt", scale);
	E.add_transition(2, 4, "cs24.txt", scale);

	double scale_upper = 10.e-25; //3.2e-25;
	E.add_transition(0, 4, "cs04.txt", scale_upper);
	E.add_transition(4, 0, "cs40.txt", scale_upper);
	
	E.set_boundary("pump1", 1 * 1.e-1);
	E.set_boundary("pump2", 0 * 1.e-1);
	E.set_boundary("signal", 1.e-7 / n_signal_UC);

	double P_total = 200.e-3;
	double L_min = 0.02;
	double L_max = 0.30;
	
	int r_points = 10;
	int l_points = 10;

	// memorize signal channel indexes
	int ch_num = 0;
	int *ch_list = new int [sufficiently_big_memory_chunk];
	for(int i = 0; i < E.F.n_ch; i++)
		if( E.F.type[i] == monochromatic && E.F.handle[i].find("signal") != string::npos)
			ch_list[ch_num++] = i;

	// prepare output
	double *g = new double [sufficiently_big_memory_chunk];
	
	std::string fname = "ripple_vs_L_and_R.txt";
	
	std::ofstream g_out( (E.output_path + fname).c_str() );
	if (!g_out)
	{
		std::cout << "run_erbium_upconversion_optim: error opening " << fname <<"\n";
		wait_and_exit();
	}
	g_out.setf( std::ios::showpoint );
	g_out.precision(6);

	for(int i = 0; i < l_points; i++)
		g_out << L_min + i * (L_max - L_min) / (double) (l_points - 1) << (i != l_points - 1 ? "\t" : "\n");

	// main loop
	for(int r = 0; r < r_points; r++)
	{
		double R = r * 1.0 / (double) (r_points - 1);
		g_out << R;
		
		for(int l = 0; l < l_points; l++)
		{
			double L = L_min + l * (L_max - L_min) / (double) (l_points - 1) - 1.e-5;
			std::cout << "L = " << L << "\tP 980 = " << R * P_total << "\nStep " << l + 1 + r * r_points << " of " << l_points * r_points << "\n";

			E.set_fiber_parameters(L, N_Er, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
			E.update();

			E.set_boundary("signal", 1.e-7 / n_signal_UC);
			E.set_boundary("pump1", R * P_total);
			E.set_boundary("pump2", (1 - R) * P_total);

			E.relaxation(1.e-8, 1.e-8, 50, E);
			
			// gain and maximal gain --------------------------------------------------------------
			double g_max = -1000., 
				   g_min = +1000., 
				   g_average = 0.;

			for(int i = 0; i < ch_num; i++)
			{
				g[i] = 10 * log10( E.F.P[E.F.n_z - 1][ ch_list[i] ] / E.F.P[0][ ch_list[i] ] );
				if (g[i] > g_max && E.F.wl[ch_list[i]] > 0.)	g_max = g[i];
				if (g[i] < g_min)	g_min = g[i];
				g_average += g[i];
			}
			g_average /= (double) ch_num;
			
			double g_ripple = g_average / (g_max - g_min); 
			
			// spectrum width ---------------------------------------------------------------------
			// level for spectrum width calculation
			
			//double g_level = 1. / 2.;
			//double g_goal = g_max * g_level;

			////calculate spectrum width
			//double rise_wl = 0.,
			//	   fall_wl = 0.;

			//for(int i = 0; i < ch_num; i++)
			//	if( g[i] > g_goal )
			//	{
			//		double k = (g_goal - g[i - 1]) / (g[i] - g_goal);
			//		rise_wl = ( E.F.wl[ ch_list[i] ] * k + E.F.wl[ch_list[i - 1]] ) / (k + 1);
			//		break;
			//	}

			//for(int i = ch_num - 1; i >= 0; i--)
			//	if( g[i] > g_goal )
			//	{
			//		double k = (g_goal - g[i - 1]) / (g[i] - g_goal);
			//		fall_wl = ( E.F.wl[ ch_list[i] ] * k + E.F.wl[ch_list[i - 1]] ) / (k + 1);
			//		break;
			//	}
			//	
			//double g_width = fall_wl - rise_wl;

			// -------------------------------------------------------------------------------------------------------------

			g_out << "\t" << g_max;
			std::cout << "\nG_max = " << g_max << "\n";
		}
		if (r != r_points - 1 ) g_out << "\n";
	}

	delete [] ch_list;
	delete [] g;
	g_out.close();
}



void run_erbium_upconversion()
{
	ErbiumUpconversionModel<n_signal_UC + 1 + 0, 1, n_ASE_UC> E;
	
	E.set_paths("e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Er_upconversion_input_revised\\",
		"e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\");

	E.add_mono_channel(980., forward, "pump");
	E.add_mono_channel(980., backward, "pump2");
	E.add_mono_channel(537., 560., n_signal_UC, forward, "signal");
	//E.add_mono_channel(541., forward, "signal");
	E.add_ASE_channel(1500., 1620., 50);
	E.add_ASE_channel(515., 560., 30);

	double scale = 1.e-25;
	E.add_transition(0, 1, "cs01.txt", scale);
	E.add_transition(1, 0, "cs10.txt", scale);
	E.add_transition(0, 2, "cs02.txt", scale);
	E.add_transition(1, 3, "cs13.txt", scale);
	E.add_transition(2, 4, "cs24.txt", scale);

	double scale_upper = 10.e-25; //3.2e-25;
	//double scale_upper = 3.2e-25; 

	E.add_transition(0, 4, "cs04.txt", scale_upper);
	E.add_transition(4, 0, "cs40.txt", scale_upper);
	
	E.set_boundary("pump", 0.e-9);
	E.set_boundary("pump2", 100 * 1.e-3);
	E.set_boundary("signal", 1.e-7 / n_signal_UC);

	E.set_fiber_parameters(0.13, N_Er, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	E.update();

	//gain_vs_L_and_P(E);
	//gain_vs_P(E);
	//gain_vs_output_power(E);
	//quantum_efficiency_vs_P_UC();

	//run_erbium_upconversion_optim();
	
	E.relaxation(1.e-8, 1.e-8, 100);
	E.output(50, dB);
}


#endif