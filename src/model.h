#ifndef MODEL_ENGINE_H
#define MODEL_ENGINE_H

//#define NEGATIVE_CS_WARNING

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
//#include <map>
//#include <algorithm>

#include "vector.h"
#include "support.h"
#include "ode.h"
#include "constants.h"
#include "interpolation.h"

#include "boost/math/special_functions/bessel.hpp"

using namespace fasim;
using std::string;
using std::cout;
using std::ofstream;


enum SignalOutput {dB, linear};											// dB - signal gain is output in decibels, else just signal power in Watts is returned
enum ChannelType {monochromatic, ASE};
enum TransitionDirection {emission = 1, absorption = -1};
enum PropagationDirection {forward = 1, backward = -1, forward_and_backward = 0};


// Grid ===================================================================================================================================

template <int r_nodes, int phi_nodes> struct Grid
{
	typedef Vector<double, r_nodes * phi_nodes> Values;
	
	static inline double integrate_uniform(const Values &g, double r_core)			// row-major order R x Phi
	{
		double dr = r_core / (double)(r_nodes - 1),
			   dphi = 2 * pi / (double)(phi_nodes);		// and no -1 here
		
		// angular Simpson integration
		double p[r_nodes];
		for(int i = 0; i < r_nodes; i++)
		{
			int shift = i * phi_nodes;
			p[i] = 2. * g[shift];
			
			for(int j = 1; j < phi_nodes; j += 2)
				p[i] += 4. * g[shift + j];
			for(int j = 2; j < phi_nodes; j += 2)
				p[i] += 2. * g[shift + j];
			p[i] *= (dphi / 3.);
		}
		
		// radial Simpson integration
		double res = 0. + p[r_nodes - 1] * r_core;
		
		for(int i = 1; i < r_nodes - 1; i += 2)
			res += 4. * p[i] * (dr * i);
		for(int i = 2; i < r_nodes - 1; i += 2)
			res += 2. * p[i] * (dr * i);
		
		return res * dr / 3.;
	}
};


template <int r_nodes> struct Grid<r_nodes, 1>
{
	typedef Vector<double, r_nodes> Values;
	
	static inline double integrate_uniform(const Values &g, double r_core)			// row-major order R x Phi
	{
		double dr = r_core / (double)(r_nodes - 1);
		double res = g[r_nodes - 1] * r_core;
		
		for(int i = 1; i < r_nodes - 1; i += 2)
			res += 4. * g[i] * (dr * i);
		for(int i = 2; i < r_nodes - 1; i += 2)
			res += 2. * g[i] * (dr * i);
		
		return 2. * pi * res * dr / 3.;
	}
};


template <> struct Grid<1, 1>
{
	typedef Vector<double, 1> Values;

	static inline double integrate_uniform(const Values &g, double r_core)
	{
		return g[0];
	}
};


// Set_envelopes_wrapper ==================================================================================================================

template <int n_level, int n_channel, int r_nodes, int phi_nodes> class Subsystem;


template <int n_level, int n_channel, int r_nodes, int phi_nodes> struct Set_envelopes_wrapper
{
	static void set_envelopes(Subsystem< n_level, n_channel, r_nodes, phi_nodes > & s, double r_core, double DRF, double NA, double n1);
};


template <int n_level, int n_channel> struct Set_envelopes_wrapper< n_level, n_channel, 1, 1 >
{
	static void set_envelopes(Subsystem< n_level, n_channel, 1, 1> & s, double r_core, double DRF, double NA, double n1);
};


// Subsystem ==============================================================================================================================

template <int n_level, int n_channel, int r_nodes, int phi_nodes> class Subsystem
{
public:
	static const int n_ch = n_channel;
	typedef Vector<double, n_channel> Channel;
	typedef Vector<Channel, ODE_dense_num> Dense_coefficients;

	Channel	*P,
			P_boundary,
			fr,
			wl,
			fr_width,
			wl_width;
	
	Dense_coefficients *dense_coeff;
	double *z;
	int n_z;

	SplineInterpolation *cs_spline[n_level][n_level];

	double cs[n_level][n_level][n_channel];
	TransitionDirection direction[n_level][n_level];

	int n_nonzero_cs[n_channel];
	int nonzero_cs_from[n_level * n_level][n_channel];
	int nonzero_cs_to[n_level * n_level][n_channel];

	int n_nonzero_cs_channel[n_level][n_level];
	int nonzero_cs_channel[n_level][n_level][n_channel];

	ChannelType	type[n_channel];
	
	SplineInterpolation *G_spline;
	SplineInterpolation *h_spline[max_mode_index];

	Vector< typename Grid<r_nodes, phi_nodes>::Values, n_channel> mode_envelope;
	int	mode_index[n_channel];

	friend struct Set_envelopes_wrapper<n_level, n_channel, r_nodes, phi_nodes>;
		void set_envelopes(double r_core, double DRF, double NA, double n1)
		{
			Set_envelopes_wrapper<n_level, n_channel, r_nodes, phi_nodes>::set_envelopes(*this, r_core, DRF, NA, n1);
		}

	string handle[n_channel];

	void free()
	{
		delete [] P;
		delete [] dense_coeff;
		delete [] z;
	}

	// ctor, copy ctor, assignment operator, dtor -------------------------------------------------
	void clear()
	{
		delete [] P;
		delete [] dense_coeff;
		delete [] z;

		for(int i = 0; i < n_level; i++)
			for(int j = 0; j < n_level; j++) 
				delete cs_spline[i][j];

		delete G_spline;
		for(int i = 0; i < max_mode_index; i++)
			delete h_spline[i];
	}


	Subsystem(): n_ch_added(0) 
	{
		P_boundary = 0.0;
		
		for(int i = 0; i < n_level; i++)
			for(int j = 0; j < n_level; j++)
				cs_spline[i][j] = NULL;
		P = NULL;
		dense_coeff = NULL;
		z = NULL;

		G_spline = NULL;
		for(int i = 0; i < max_mode_index; i++)
			h_spline[i] = NULL;
	}


	Subsystem & operator=(const Subsystem & s)
	{
		clear();

		P_boundary = s.P_boundary;
		fr = s.fr;
		wl = s.wl;
		fr_width = s.fr_width;
		wl_width = s.wl_width;

		n_z = s.n_z;		
		if(n_z == 0)
		{
			P = NULL;
			dense_coeff = NULL;
			z = NULL;
		}
		else
		{
			P = new Channel [n_z];
			dense_coeff = new Dense_coefficients [n_z];
			z = new double [n_z];
			
			for(int i = 0; i < n_z; i++)
			{
				P[i] = s.P[i];
				dense_coeff[i] = s.dense_coeff[i];
				z[i] = s.z[i];
			}
		}

		for(int i = 0; i < n_level; i++)
			for(int j = 0; j < n_level; j++)
			{
				if(s.cs_spline[i][j] == NULL)
					cs_spline[i][j] = NULL;
				else
					cs_spline[i][j] = new SplineInterpolation( *(s.cs_spline[i][j]) );
				
				direction[i][j] = s.direction[i][j];
				n_nonzero_cs_channel[i][j] = s.n_nonzero_cs_channel[i][j];
				for(int k = 0; k < n_channel; k++)
				{
					cs[i][j][k] = s.cs[i][j][k];
					nonzero_cs_channel[i][j][k] = s.nonzero_cs_channel[i][j][k];
				}
			}

		for(int k = 0; k < n_channel; k++)
		{
			n_nonzero_cs[k] = s.n_nonzero_cs[k];
			type[k] = s.type[k];
			for(int i = 0; i < n_level * n_level; i++)
			{
				nonzero_cs_from[i][k] = s.nonzero_cs_from[i][k];
				nonzero_cs_to[i][k] = s.nonzero_cs_to[i][k];
			}
			mode_index[k] = s.mode_index[k];
			handle[k] = s.handle[k];
		}
			
		if (s.G_spline == NULL)
			G_spline = NULL;
		else
			G_spline = new SplineInterpolation(*s.G_spline);

		for(int i = 0; i < max_mode_index; i++)
			if (s.h_spline[i] == NULL)
				h_spline[i] = NULL;
			else
				h_spline[i] = new SplineInterpolation(*s.h_spline[i]);
	
		return *this;
	}

	Subsystem(Subsystem & s)
	{
		for(int i = 0; i < n_level; i++)
			for(int j = 0; j < n_level; j++)
				cs_spline[i][j] = NULL;
		P = NULL;
		dense_coeff = NULL;
		z = NULL;

		G_spline = NULL;
		for(int i = 0; i < max_mode_index; i++)
			h_spline[i] = NULL;
		
		*this = s;
	}

	~Subsystem()
	{
		clear();
	}

	// i/o members --------------------------------------------------------------------------------

	int n_ch_added;
	void add_channel(double _wl, ChannelType ch_type, double ch_width, string _handle = "none");

	void update();					// reset and update all cross-section info

	bool reset_channel(string _handle, double wl_new);

	void output(int refined_points, SignalOutput signal_output, 
		ofstream &pump_out, ofstream &signal_out, ofstream &ASE_out, ofstream &ASE_total_out, ofstream &ASE_spectrum_out); 
};


// solution container =====================================================================================================================

template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes> class ModelEngine;


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes> struct SolutionContainer
{
	typedef ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes> ME;

	typedef Subsystem< n_level, n_F_mono + n_ASE, r_nodes, phi_nodes > F_subsystem;
	typedef Subsystem< n_level, n_B_mono + n_ASE, r_nodes, phi_nodes > B_subsystem;
	
	F_subsystem F;
	B_subsystem B;

	SolutionContainer(double L = magic_L)							// see constants.h for magic_L definition
	{
		F.n_z = 10;
		F.z = new double [F.n_z];
		F.P = new F_subsystem::Channel [F.n_z];
		F.dense_coeff = new F_subsystem::Dense_coefficients [F.n_z];

		for(int i = 0; i < F.n_z; i++)
		{
			double z = i * L  / (F.n_z - 1);			
			F.z[i] = z;												
			F.P[i] = 0.;
			F.dense_coeff[i] = 0.;
		}

		B.n_z = 10;
		B.z = new double [B.n_z];
		B.P = new B_subsystem::Channel [B.n_z];
		B.dense_coeff = new B_subsystem::Dense_coefficients [B.n_z];
	
		for(int i = 0; i < B.n_z; i++)
		{
			double z = L - i * L / (double)(B.n_z - 1);		// L / 100. is added to avoid error accumulation issue
			B.z[i] = z;												// when L >> 1 z[B.n_z - 1] may be not zero
			B.P[i] = 0.;
			B.dense_coeff[i] = 0.;
		}		
	}
	
	SolutionContainer(const SolutionContainer &s)
	{
		F = s.F;
		B = s.B;
	}

	SolutionContainer(const ME &m)
	{
		if( m.F.dense_coeff != NULL && m.B.dense_coeff != NULL)
		{
			F = m.F;
			B = m.B;
		}
		else
		{
			SolutionContainer tmp;
			*this = tmp;
		}
	}

	SolutionContainer & operator=(const ME & m)
	{
		F = m.F;
		B = m.B;
		return *this;
	}

	SolutionContainer & operator=(const SolutionContainer & s)
	{
		F = s.F;
		B = s.B;
		return *this;
	}
};


// calculate_W_val_wrapper ================================================================================================================


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes> struct Calculate_W_val_wrapper
{
	typedef ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes> ME;
	
	static void calculate_W_val(ME & m, const typename ME::F_subsystem::Channel & _PF, const typename ME::B_subsystem::Channel & _PB )
		{
			for(int i = 0; i < n_level; i++)
				for(int j = 0; j < n_level; j++)
				{
					m.W_val[i][j] = 0.0;
					for(int q = 0; q < m.F.n_nonzero_cs_channel[i][j]; q++)
					{
						int k = m.F.nonzero_cs_channel[i][j][q];
						m.W_val[i][j] += ( m.F.cs[i][j][k] * _PF[k] / ( h_planck * m.F.fr[k]) ) * m.F.mode_envelope[k];
					}

					for(int q = 0; q < m.B.n_nonzero_cs_channel[i][j]; q++)
					{
						int k = m.B.nonzero_cs_channel[i][j][q];
						m.W_val[i][j] += ( m.B.cs[i][j][k] * _PB[k] / ( h_planck * m.B.fr[k]) ) * m.B.mode_envelope[k];
					}
				}
		}
};


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE> struct Calculate_W_val_wrapper<n_level, n_F_mono, n_B_mono, n_ASE, 1, 1>
{
	typedef ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, 1, 1> ME;

	static void calculate_W_val(ME &m, const typename ME::F_subsystem::Channel & _PF, const typename ME::B_subsystem::Channel & _PB )
		{
			for(int i = 0; i < n_level; i++)
				for(int j = 0; j < n_level; j++)
				{
					m.W_val[i][j] = 0.0;
					for(int q = 0; q < m.F.n_nonzero_cs_channel[i][j]; q++)
					{
						int k = m.F.nonzero_cs_channel[i][j][q];
						m.W_val[i][j] += ( m.F.cs[i][j][k] * _PF[k] / ( h_planck * m.F.fr[k]) ) * m.F.mode_envelope[k];
					}

					for(int q = 0; q < m.B.n_nonzero_cs_channel[i][j]; q++)
					{
						int k = m.B.nonzero_cs_channel[i][j][q];
						m.W_val[i][j] += ( m.B.cs[i][j][k] * _PB[k] / ( h_planck * m.B.fr[k]) ) * m.B.mode_envelope[k];
					}
					
					m.W_val[i][j] /= (pi * (m.r_core * m.DRF) * (m.r_core * m.DRF));
				}
		}
};

// MODEL_ENGINE ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes> class ModelEngine
{
public:
	// physics ------------------------------------------------------------------------------------

	double L, N_activator, NA, r_core, DRF, n1;
	void set_fiber_parameters(double _L, double _N_activator, double _NA, double _r_core, double _DRF, double _n1 = 1.52069) 
		{ L = _L; N_activator = _N_activator; NA = _NA; r_core = _r_core; DRF = _DRF; n1 = _n1;
		N_previous_step = 0.;
		N_previous_step[0] = _N_activator;
	}

	typedef Subsystem< n_level, n_F_mono + n_ASE, r_nodes, phi_nodes > F_subsystem;
	typedef Subsystem< n_level, n_B_mono + n_ASE, r_nodes, phi_nodes > B_subsystem;
	F_subsystem F;
	B_subsystem B;
	
	typename F_subsystem::Channel dPF(double _z, const typename F_subsystem::Channel & PF);
	typename B_subsystem::Channel dPB(double _z, const typename B_subsystem::Channel & PB);

	typedef Vector<typename Grid<r_nodes, phi_nodes>::Values, n_level> N_vector;
	N_vector N_previous_step;
	N_vector *N;
	double *t;
	int n_t;
	
	virtual N_vector dN(double _t, const N_vector & _N) = 0;

	typename Grid<r_nodes, phi_nodes>::Values W_val[n_level][n_level];
	const typename Grid<r_nodes, phi_nodes>::Values &W(int i, int j) const	{ return W_val[i][j]; }

	friend struct Calculate_W_val_wrapper<n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>;
	void calculate_W_val(const typename F_subsystem::Channel & _PF, const typename B_subsystem::Channel & _PB)
	{
		Calculate_W_val_wrapper<n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::calculate_W_val(*this, _PF, _PB);
	}

	typedef SolutionContainer< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes> ModelSolution;

	bool relaxation(double abs_tolerance = 1.e-8, double rel_tolerance = 1.e-8, int max_iterations = 100, 
		ModelSolution sol = magic_L);		// see constants.h for magic_L definition	do NEVER use SOL with BIDIRECTIONAL pump

	ModelEngine();
	~ModelEngine();

	// counters -----------------------------------------------------------------------------------
	
	int dPF_called, dPB_called, dN_called, dN_average_called, dN_min_called, dN_max_called;

	// input/output -------------------------------------------------------------------------------

	string input_path, output_path;
	void set_paths(string _input_path, string _output_path);

	void add_transition(int i, int j, string filename, double cs_scale = 1., double wl_scale = 1.);

	void add_channel(double wl, ChannelType ch_type, PropagationDirection direction, double ch_width, string handle = "none");
	void add_mono_channel(double wl, PropagationDirection direction, string handle = "none");
	void add_mono_channel(double wl0, double wl1, int n_ch, PropagationDirection direction, string handle = "none");
	void add_ASE_channel(double wl0, double wl1, int n_ch, string handle = "none");

	void set_h(int _mode_index, string filename, double h_scale = 1., double wl_scale = 1.);
	void set_G(string filename, double G_scale = 1., double wl_scale = 1.);

	void set_boundary(string handle, double P0);

	void update() 				// reset and update all cross-section info
	{
		F.update();
		B.update();
		F.set_envelopes(r_core, DRF, NA, n1);
		B.set_envelopes(r_core, DRF, NA, n1);
	}

	void reset_channel(string handle, double wl_new);

	std::ofstream set_output_stream(string &filename);

	void output(int refined_points = 0, SignalOutput signal_output = linear, 
		string f_pump_name = "f_pump.txt", string b_pump_name = "b_pump.txt",
		string f_signal_name = "f_signal.txt", string b_signal_name = "b_signal.txt",
		string f_ASE_name = "f_ASE.txt", string b_ASE_name = "b_ASE.txt",
		string f_ASE_total_name = "f_ASE_total.txt", string b_ASE_total_name = "b_ASE_total.txt",
		string f_ASE_spectrum_name = "f_ASE_spectrum.txt", string b_ASE_spectrum_name = "b_ASE_spectrum.txt",
		string N_name = "N.txt", 
		string NF_name = "NF.txt");
};	


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
	ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::ModelEngine()
		: input_path(""), output_path("")
{
	if ( r_nodes % 2 == 0 || (phi_nodes % 2 == 1 && phi_nodes != 1) )	// evenness check, required for Simpson's rule integration	
		wait_and_exit("ModelEngine: incorrect grid dimensions");		// r_nodes must be odd, and phi_nodes must be even
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::~ModelEngine()
{
	// no need to delete N and t, since they are created and deleted immedialely after use in ModelEngine::dPF and ModelEngine::dPB
}


#include "model.cpp"
#include "model_io.cpp"

#endif