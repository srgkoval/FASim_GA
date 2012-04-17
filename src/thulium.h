#ifndef THULIUM_H
#define THULIUM_H

#include "model_specification.h"

void run_thulium()
{
	const int n_signal = 16;
	typedef ThuliumModel<1 + n_signal , 1, 30 + 30> thulium_engine;
	thulium_engine T;

	T.add_mono_channel(1064., forward, "pump1");
	T.add_mono_channel(1064., backward, "pump2");
	//T.add_mono_channel(1470., forward, "signal");
	T.add_mono_channel(1450., 1525., n_signal, forward, "signal");
	T.add_ASE_channel(785., 820., 30);
	T.add_ASE_channel(1400., 2100., 30);
	
	T.set_paths("e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Thulium_input_new\\",
		"e:\\Programming - code and etc\\Fiber Amplifier Simulation\\!Spectra\\Output\\");
	double scale = 1.e-25;
	T.add_transition(0, 1, "cs01.txt", scale);
	T.add_transition(0, 2, "cs02.txt", scale);
	T.add_transition(0, 3, "cs03.txt", scale);
	T.add_transition(1, 0, "cs10.txt", scale);
	T.add_transition(1, 3, "cs13.txt", scale);
	T.add_transition(1, 4, "cs14.txt", scale);
	T.add_transition(3, 0, "cs30.txt", scale);
	T.add_transition(3, 1, "cs31.txt", scale);

	// T.set_G("G_factor_Peterka.txt");
	T.set_boundary("pump1", 0.0);
	T.set_boundary("pump2", 1.0);

	T.set_boundary("signal", 1.e-7 / n_signal);
	
	T.set_fiber_parameters(6, 1.56e25, 0.30, 1.3e-6, 1.3e-6 / 1.3e-6);
	T.update();

	T.relaxation();
	T.output(0, dB);
	
    //int n = T.F.n_z;
	//double *x = new double [n],
	//	   *y = new double [n];

	//for(int i = 0; i < n; i ++)
	//{
	//	x[i] = T.F.z[i];
	//	y[i] = 10 * log10(T.F.P[i][2] / T.F.P[0][2]);
	//}

	//f_write_table(x, y, n, T.output_path + "gain_vs_position.txt");

	// gain vs P additional for Kasamatsu scheme --------------------------------------------------
	//double P_min = 1.e-3,
	//	   P_max = 1000.e-3;
	//
	//int n = 20;
	//double *x = new double [n],
	//	   *y = new double [n];
	//
	//for(int i = 0; i < n; i++)
	//{
	//	double Pa = P_min * pow(10, i * (log10(P_max) - log10(P_min)) / double(n - 1));
	//	std::cout << Pa << "\n";

	//	T.set_fiber_parameters(15., 0.7e25, 0.28, 1.4e-6, 1.05e-6 / 1.4e-6);
	//	T.set_boundary("pump2", Pa);
	//	T.update();

	//	T.relaxation(1.e-8, 1.e-8, 100, T);
	//	x[i] = 10*log10(Pa / 1.e-3);
	//	y[i] = 10 * log10(T.F.P[T.F.n_z - 1][1] / T.F.P[0][1]);
	//	std::cout << y[i] << "\n";
	//}
	//f_write_table(x, y, n, T.output_path + "gain_vs_P.txt");
	// --------------------------------------------------------------------------------------------
}


#endif