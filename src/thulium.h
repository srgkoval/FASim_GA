#ifndef THULIUM_H
#define THULIUM_H

#include "model_specification.h"
#include "vector.h"
#include "ga.h"
#include "fiber.h"

const int n_ASE_800 = 30,
          n_ASE_2000 = 30;

const int n_signal_thulium = 16;

const double N_thulium = 1.56e25;
const double P_total_thulium = 1.;

template <int F_mono_num, int B_mono_num, int ASE_num, int N_param = 4, int N_obj = 2> class Thulium_ga_wrapper
{
public:
    typedef ThuliumModel<F_mono_num, B_mono_num, ASE_num> TM;
    TM *tm;

    Thulium_ga_wrapper(TM  * _tm) : tm(_tm) {}

    Vector<double, N_obj> operator()(const typename GA<N_param, N_obj>::Individual &x)
    {
        Vector<double, N_obj> res;

        tm->set_fiber_parameters(x[0], N_thulium, fiber_NA, fiber_r, fiber_r_doping / fiber_r);
        tm->set_boundary("pump1", x[1] * P_total_thulium);
        tm->set_boundary("pump2", (1. - x[1]) * P_total_thulium);
        tm->reset_channel("pump1", x[2]);
        tm->reset_channel("pump2", x[3]);
        
        tm->update();

        tm->relaxation();

	    int ch_num = 0;
	    int *ch_list = new int [F_mono_num];
	    for(int i = 0; i < tm->F.n_ch; i++)
		    if( tm->F.type[i] == monochromatic && tm->F.handle[i].find("signal") != string::npos)
			    ch_list[ch_num++] = i;

		double g_max = -1000., 
			   g_min = +1000., 
			   g_average = 0.,
               g [F_mono_num];

		for(int i = 0; i < ch_num; i++)
		{
			g[i] = 10 * log10( tm->F.P[tm->F.n_z - 1][ ch_list[i] ] / tm->F.P[0][ ch_list[i] ] );
			if (g[i] > g_max)	g_max = g[i];
			if (g[i] < g_min)	g_min = g[i];
			g_average += g[i];
		}
		g_average /= (double) ch_num;
			
		double g_ripple = g_average / (g_max - g_min);
        
        delete [] ch_list;

        // multiobjective: max gain + gain ripple -------------------------------------------------
        if(g_average >= 0.)
        {
            res[0] = - g_max;
            res[1] = g_max - g_min;
        }
        else
        {
            res[0] = 0.;
            res[1] = 100.;
        }
        
        // single-objective: max gain -------------------------------------------------------------
        //res = - g_max;
        
        
        return res;
    }
};


void run_thulium()
{
	typedef ThuliumModel<1 + 0 + n_signal_thulium , 1, n_ASE_800 + n_ASE_2000> thulium_engine;
	thulium_engine T;

	T.add_mono_channel(798., forward, "pump1");
	T.add_mono_channel(1020., backward, "pump2");
	//T.add_mono_channel(1470., forward, "signal");
	T.add_mono_channel(1450., 1490., n_signal_thulium, forward, "signal");
	T.add_ASE_channel(785., 820., n_ASE_800);
	T.add_ASE_channel(1400., 2100., n_ASE_2000);
	
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
	
    double R = 0.25;
    T.set_boundary("pump1", R * 1.);
	T.set_boundary("pump2", (1. - R) * 1.);

	T.set_boundary("signal", 1.e-4 / n_signal_thulium);
	
	T.set_fiber_parameters(10., N_thulium, fiber_NA, fiber_r, fiber_r_doping / fiber_r);
	T.update();

	T.relaxation();
	T.output(0, dB);

    // run GA multiobjective-----------------------------------------------------------------------
    //const int n_o = 2,
    //          n_p = 4;

    //Thulium_ga_wrapper<1 + 1 + n_signal_thulium, 0, n_ASE_800 + n_ASE_2000, n_p, n_o> thulium_ga_wrapper(&T);

    //Vector<double, n_p> lower(2., 0., 720., 1000.), upper(20., 1., 850., 1100.);
    //
    //GA<n_p, n_o> ga;
    //ga.run_multiobjective(thulium_ga_wrapper, lower, upper);

    // end GA -------------------------------------------------------------------------------------


    // run GA single-objective---------------------------------------------------------------------
    //const int n_o = 2,
    //          n_p = 4;

    //Thulium_ga_wrapper<1 + 0 + n_signal_thulium, 1, n_ASE_800 + n_ASE_2000, n_p, n_o> thulium_ga_wrapper(&T);
    //GA<n_p, n_o> ga;

    //Vector<double, n_p> lower(2., 0., 1000., 1350.), upper(20., 1., 1100., 1450.);  // scheme a, b
    ////Vector<double, n_p> lower(2., 0., 1330., 1450.), upper(20., 1., 1450., 1650.);  // scheme c 
    ////Vector<double, n_p> lower(2., 0., 720., 1330.), upper(20., 1., 850., 1450.);  // scheme d
    ////Vector<double, n_p> lower(2., 0., 720., 1000.), upper(20., 1., 850., 1100.);  // scheme e
    //
    //ga.options.output_directory = "e:\\Programming - code and etc\\Genetic Algorithm\\Output\\1\\";
    //ga.run_multiobjective(thulium_ga_wrapper, lower, upper);

    // end GA -------------------------------------------------------------------------------------



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