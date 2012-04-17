#ifndef MODEL_IO_CPP
#define MODEL_IO_CPP

#include "model.h"
#include "boost/math/special_functions/bessel.hpp"

// INPUT //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// subsystem methods ======================================================================================================================

template<int n_level, int n_channel, int r_nodes, int phi_nodes>
void Subsystem<n_level, n_channel, r_nodes, phi_nodes>::
		add_channel(double _wl, ChannelType ch_type, double ch_width, string _handle)
{
	if(n_ch_added >= n_channel)
		wait_and_exit("Subsystem::add_channel: attempt to add too many channels");
	
	int i = n_ch_added;
	fr[i] = nm_to_Hz(_wl);
	wl[i] = _wl;
	fr_width[i] = nm_to_Hz(_wl - ch_width / 2.) - nm_to_Hz(_wl + ch_width / 2.);
	wl_width[i] = ch_width;
	type[i] = ch_type;
	handle[i] = _handle;
	
	n_ch_added++;	
}


template<int n_level, int n_channel, int r_nodes, int phi_nodes>
void Subsystem<n_level, n_channel, r_nodes, phi_nodes>::
		update()
{
	if(n_ch_added != n_channel)
	{
		std::cout << "Subsystem::update: added " << n_ch_added << " channels instead of " << n_channel << " required";
		wait_and_exit();
	}

	for(int k = 0; k < n_channel; k++)
		n_nonzero_cs[k] = 0;
	
	for(int i = 0; i < n_level; i++)
		for(int j = 0; j < n_level; j++)
		{
			n_nonzero_cs_channel[i][j] = 0;
			
			if(i > j) 
				direction[i][j] = emission;
			else
				direction[i][j] = absorption;
			
			if(cs_spline[i][j] != NULL)
			{
				double cs_max_ij = 0.,
					   tmp = 0.;

				for(int k = 0; k < n_channel; k++)
				{
					tmp = (*cs_spline[i][j])(wl[k]);
					if( tmp > cs_max_ij )
						cs_max_ij = tmp;
				}

				for(int k = 0; k < n_channel; k++)
				{
					double tmp = (*cs_spline[i][j])(wl[k]);
					if( tmp / cs_max_ij >= cs_sensitivity_threshold)
					{
						cs[i][j][k] = tmp;

						nonzero_cs_from[n_nonzero_cs[k]][k] = i;
						nonzero_cs_to[n_nonzero_cs[k]][k] = j;
						n_nonzero_cs[k]++;

						nonzero_cs_channel[i][j][n_nonzero_cs_channel[i][j]] = k;
						n_nonzero_cs_channel[i][j]++;
					}
					else
						cs[i][j][k] = 0.;
				}
			}
		}
}


template<int n_level, int n_channel, int r_nodes, int phi_nodes>
bool Subsystem<n_level, n_channel, r_nodes, phi_nodes>::
	reset_channel(string _handle, double wl_new)
{
	bool found = false;
	for(int i = 0; i < n_channel; i++)
		if(handle[i] == _handle)
		{
			found = true;
			
			if(type[i] == ASE)
				wait_and_exit("Subsystem::reset_channel: only monochromatic channels can be reset");

			for(int j = i + 1; j < n_channel; j++)
				if(handle[j] == _handle)
				{
					std::cout << "Subsystem::reset_channel: multiple channels \"" << handle << "\"";
					wait_and_exit();
				}
			
			wl[i] = wl_new;
			fr[i] = nm_to_Hz(wl_new);
		}
	return found;
}

// Set_envelopes_wrapper ==================================================================================================================

template<int n_level, int n_channel, int r_nodes, int phi_nodes>
void Set_envelopes_wrapper<n_level, n_channel, r_nodes, phi_nodes>::
	set_envelopes(Subsystem< n_level, n_channel, r_nodes, phi_nodes > & s, double r_core, double DRF, double NA, double n1)
{
	using boost::math::cyl_bessel_j;
	using boost::math::cyl_bessel_k;

	double U, V, W, h;

	for(int k = 0; k < n_channel; k++)
	{
		s.mode_index[k] = 0;

		V = 2 * pi * r_core * NA / (1.e-9 * s.wl[k]);
			
		if(s.h_spline[s.mode_index[k]] == NULL)
			U = V * (1 + sqrt2) / ( 1 + pow(4 + pow(V, 4), 0.25) );
		else
		{
			h = (*s.h_spline[s.mode_index[k]])(s.wl[k]);
			U = r_core * sqrt( pow(2 * pi * n1 /(1.e-9 * s.wl[k]), 2) - h * h );
		}
			
		W = sqrt(V * V - U * U);

		double psi_int_S = pi * r_core * r_core * ( pow(cyl_bessel_j(0, U), 2) * (pow(cyl_bessel_k(1, W), 2) / pow(cyl_bessel_k(0, W), 2)) + pow(cyl_bessel_j(1, U), 2));			
		
		double dr = r_core * DRF / (double) (r_nodes - 1);
		
		for(int i = 0; i < r_nodes; i++)
			for(int j = 0; j < phi_nodes; j++)
				s.mode_envelope[k][i * phi_nodes + j] = pow( cyl_bessel_j(0, U * (dr * i) / r_core), 2);
		
		s.mode_envelope[k] *= (1. / psi_int_S);
	}
}


template<int n_level, int n_channel>
void Set_envelopes_wrapper<n_level, n_channel, 1, 1>::
	set_envelopes(Subsystem< n_level, n_channel, 1, 1 > & s, double r_core, double DRF, double NA, double n1)
{
	using boost::math::cyl_bessel_j;
	using boost::math::cyl_bessel_k;

	double U, V, W, h;
	double b = r_core * DRF;

	for(int k = 0; k < n_channel; k++)
	{
		s.mode_index[k] = 0;

		if (s.G_spline == NULL)
		{
			V = 2 * pi * r_core * NA / (1.e-9 * s.wl[k]);
			if(s.h_spline[s.mode_index[k]] == NULL)
				U = V * (1 + sqrt2) / ( 1 + pow(4 + pow(V, 4), 0.25) );
			else
			{
				h = (*s.h_spline[s.mode_index[k]])(s.wl[k]);
				U = r_core * sqrt( pow(2 * pi * n1 /(1.e-9 * s.wl[k]), 2) - h * h );
			}
			W = sqrt(V * V - U * U);

			double psi_int_S = pi * r_core * r_core * ( pow(cyl_bessel_j(0, U), 2) * (pow(cyl_bessel_k(1, W), 2) / pow(cyl_bessel_k(0, W), 2)) + pow(cyl_bessel_j(1, U), 2));
			s.mode_envelope[k][0] = pi * b * b * ( pow( cyl_bessel_j(0, U * b / r_core), 2) + pow( cyl_bessel_j(1, U * b / r_core), 2) ) / psi_int_S;
		}
		else
			s.mode_envelope[k][0] = (*s.G_spline)(s.wl[k]);
	}
};


// model engine methods ===================================================================================================================

template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	set_paths(string _input_path, string _output_path)
{
	input_path = _input_path;
	output_path = _output_path;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	add_transition(int i, int j, string filename, double cs_scale, double wl_scale)
{	
	if(F.cs_spline[i][j] != NULL || B.cs_spline[i][j] != NULL)
	{
		std::cout << "ModelEngine::add_transition: transition " << i << "->" << j << "already exists\n";
		wait_and_exit();
	}

	filename = input_path + filename;
	int n;
	double *x = new double [sufficiently_big_memory_chunk];
	double *y = new double [sufficiently_big_memory_chunk];
	f_read_table(x, y, &n, filename);

	for(int ii = 0; ii < n; ii++)
	{
		x[ii] *= wl_scale;
		y[ii] *= cs_scale;
	}

	F.cs_spline[i][j] = new SplineInterpolation(x, y, n);
	B.cs_spline[i][j] = new SplineInterpolation(x, y, n);

	delete [] x;
	delete [] y;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	add_channel(double wl, ChannelType ch_type, PropagationDirection direction, double ch_width, string handle)
{
	switch(direction)
	{
		case forward:
			F.add_channel(wl, ch_type, ch_width, handle);
			break;
		case backward:
			B.add_channel(wl, ch_type, ch_width, handle);
			break;
		case forward_and_backward:
			F.add_channel(wl, ch_type, ch_width, handle);
			B.add_channel(wl, ch_type, ch_width, handle);
			break;
	}
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	add_mono_channel(double wl, PropagationDirection direction, string handle)
{
	add_channel(wl, monochromatic, direction, 0.0, handle); 
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	add_mono_channel(double wl0, double wl1, int n_ch, PropagationDirection direction, string handle)
{
	assert(wl0 < wl1 && n_ch > 0);
	
	double ch_separation = (wl1 - wl0) / (double)(n_ch - 1);
	for(int i = 0; i < n_ch; i++)
		add_mono_channel(wl0 + i * ch_separation, direction, handle);
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	add_ASE_channel(double wl0, double wl1, int n_ch, string handle = "none")
{
	assert(wl0 < wl1 && n_ch > 0);

	double ch_width = (wl1 - wl0) / (double)n_ch;
	for(int i = 0; i < n_ch; i++)
		add_channel(wl0 + ch_width / 2. + i * ch_width, ASE, forward_and_backward, ch_width, handle);
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	set_h(int _mode_index, string filename, double h_scale, double wl_scale)
{
	if(F.h_spline[_mode_index] != NULL || B.h_spline[_mode_index] != NULL)
	{
		std::cout << "ModelEngine::set_h: h_spline for mode index = " << _mode_index << " already exists\n";
		wait_and_exit();
	}

	filename = input_path + filename;
	int n;
	double *x = new double [sufficiently_big_memory_chunk];
	double *y = new double [sufficiently_big_memory_chunk];
	f_read_table(x, y, &n, filename);

	for(int i = 0; i < n; i++)
	{
		x[i] *= wl_scale;
		y[i] *= h_scale;
	}

	F.h_spline[_mode_index] = new SplineInterpolation(x, y, n, false);
	B.h_spline[_mode_index] = new SplineInterpolation(x, y, n, false);

	delete [] x;
	delete [] y;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	set_G(string filename, double G_scale, double wl_scale)
{
	if(F.G_spline != NULL || B.G_spline != NULL)
	{
		std::cout << "ModelEngine::set_G: G_spline already exists\n";
		wait_and_exit();
	}

	filename = input_path + filename;
	int n;
	double *x = new double [sufficiently_big_memory_chunk];
	double *y = new double [sufficiently_big_memory_chunk];
	f_read_table(x, y, &n, filename);

	for(int i = 0; i < n; i++)
	{
		x[i] *= wl_scale;
		y[i] *= G_scale;
	}

	F.G_spline = new SplineInterpolation(x, y, n, false);
	B.G_spline = new SplineInterpolation(x, y, n, false);

	delete [] x;
	delete [] y;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	set_boundary(string handle, double P0)
{
	bool f_found = false,
		 b_found = false;
	
	for(int i = 0; i < F.n_ch; i++)
		if(F.handle[i] == handle)
		{
			f_found = true;
			F.P_boundary[i] = P0;
		}

	for(int i = 0; i < B.n_ch; i++)
		if(B.handle[i] == handle)
		{
			b_found = true;
			B.P_boundary[i] = P0;
		}

	if(!f_found && !b_found)
	{
		std::cout << "ModelEngine::set_boundary: handle \"" << handle <<" \" not found\n";
		wait_and_exit();
	}
	if(f_found && b_found)
	{
		std::cout << "ModelEngine::set_boundary: handle \"" << handle <<" \" found in both forward and backward subsystems\n";
		wait_and_exit();
	}
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	reset_channel(string handle, double wl_new)
{
	bool f = F.reset_channel(handle, wl_new);
	bool b = B.reset_channel(handle, wl_new);
	if( f && b)
	{	
		std::cout << "ModelEngine::reset_channel: both forward and backward systems have channel \"" << handle << "\"";
		wait_and_exit();
	}
	if(!f && !b)
	{
		std::cout << "ModelEngine::reset_channel: no channel \"" << handle << "\"found";
		wait_and_exit();
	}
	
	update();
}

// OUTPUT /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<int n_level, int n_channel, int r_nodes, int phi_nodes>
void Subsystem<n_level, n_channel, r_nodes, phi_nodes>::
	output(int refined_points, SignalOutput signal_output, 
		ofstream &pump_out, ofstream &signal_out, ofstream &ASE_out, ofstream &ASE_total_out, ofstream &ASE_spectrum_out)
{
	std::string separator = "\t";
	double dz = fabs(z[0] - z[n_z - 1]) / (double)(refined_points - 1);
	
	int ch_num = 0;
	int *ch_list = new int [sufficiently_big_memory_chunk];
	
	// pump
	for(int i = 0; i < n_channel; i++)
		if(type[i] == monochromatic && handle[i].find("pump") != string::npos)
			ch_list[ch_num++] = i;

	for(int i = 0; i < ch_num; i++)
		pump_out << wl[ch_list[i]] << (i != ch_num - 1 ? separator : "\n");
	
	if(refined_points == 0)
		for(int i = 0; i < n_z; i++)
		{
			pump_out << z[i];
			for(int j = 0; j < ch_num; j++)
				pump_out << separator << P[i][ch_list[j]];
			if(i != n_z - 1) pump_out << "\n";
		}
	else
		for(int i = 0; i < refined_points; i++)
		{
			double zi = i * dz;
			Channel Pi = ODE_interpolate(zi, z, P, dense_coeff, n_z);
			
			pump_out << zi;
			for(int j = 0; j < ch_num; j++)
				pump_out << separator << Pi[ch_list[j]];
			if(i != refined_points - 1) pump_out << "\n";
		};

	// signal
	ch_num = 0;
	for(int i = 0; i < n_channel; i++)
		if(type[i] == monochromatic && handle[i].find("signal") != string::npos)
			ch_list[ch_num++] = i;

	for(int i = 0; i < ch_num; i++)
		signal_out << wl[ch_list[i]] << (i != ch_num - 1 ? separator : "\n");
	
	if(refined_points == 0)
		for(int i = 0; i < n_z; i++)
		{
			signal_out << z[i];
			for(int j = 0; j < ch_num; j++)
				if(signal_output == linear)
					signal_out << separator << P[i][ch_list[j]];
				else
					signal_out << separator << 10 * log10(P[i][ch_list[j]] / P[0][ch_list[j]]);
			if(i != n_z - 1) signal_out << "\n";
		}
	else
		for(int i = 0; i < refined_points; i++)
		{
			double zi = i * dz;
			Channel Pi = ODE_interpolate(zi, z, P, dense_coeff, n_z);
			
			signal_out << zi;
			for(int j = 0; j < ch_num; j++)
				if(signal_output == linear)
					signal_out << separator << Pi[ch_list[j]];
				else
					signal_out << separator << 10 * log10(Pi[ch_list[j]] / P[0][ch_list[j]]);
			if(i != refined_points - 1) signal_out << "\n";
		};

	// ASE and ASE total
	ch_num = 0;
	for(int i = 0; i < n_channel; i++)
		if(type[i] == ASE)
			ch_list[ch_num++] = i;
	
	double average_wl = 0.;
	for(int i = 0; i < ch_num; i++)
	{
			ASE_out << wl[ch_list[i]] << (i != ch_num - 1 ? separator : "\n");
			average_wl += wl[ch_list[i]];
	}
	ASE_total_out << average_wl / ch_num << "\n";	

	if(refined_points == 0)
		for(int i = 0; i < n_z; i++)
		{
			ASE_out << z[i];
			ASE_total_out << z[i];
			
			double pow_total = 0.;

			for(int j = 0; j < ch_num; j++)
			{
				ASE_out << separator << P[i][ch_list[j]];
				pow_total += P[i][ch_list[j]];
			}
			ASE_total_out << separator << pow_total;

			if(i != n_z - 1)
			{
				ASE_out << "\n";
				ASE_total_out << "\n";
			}
		}
	else
		for(int i = 0; i < refined_points; i++)
		{
			double zi = i * dz;
			Channel Pi = ODE_interpolate(zi, z, P, dense_coeff, n_z);
			
			ASE_out << zi;
			ASE_total_out << zi;
			
			double pow_total = 0.;

			for(int j = 0; j < ch_num; j++)
			{
				ASE_out << separator << Pi[ch_list[j]];
				pow_total += Pi[ch_list[j]];
			}
			ASE_total_out << separator << pow_total;

			if(i != refined_points - 1)
			{
				ASE_out << "\n";
				ASE_total_out << "\n";
			}
		};

	//ASE spectrum
	ASE_spectrum_out << ch_num << "\n";
	for(int i = 0; i < ch_num; i++)
	{
		int j = ch_list[i];

		ASE_spectrum_out << wl[j];
		if(signal_output == linear)
			ASE_spectrum_out << separator << P[n_z - 1][j] / wl_width[j];
		else
			ASE_spectrum_out << separator << 10 * log10(P[n_z - 1][j] / wl_width[j]);
		
		if(i != ch_num - 1) ASE_spectrum_out << "\n";
	}

	delete [] ch_list;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
std::ofstream ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	 set_output_stream(string &filename)
{
	filename = output_path + filename;
	std::ofstream fstream(filename.c_str());
	
	if (!fstream)
	{
		std::cout << "ModelEngine::set_output_stream: error opening " << filename <<"\n";
		wait_and_exit();
	}
	fstream.setf( std::ios::showpoint );
	fstream.precision( 5 );

	return fstream;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
void ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	output(int refined_points, SignalOutput signal_output,
		string f_pump_name, string b_pump_name,
		string f_signal_name, string b_signal_name,
		string f_ASE_name, string b_ASE_name, 
		string f_ASE_total_name, string b_ASE_total_name,
		string f_ASE_spectrum_name, string b_ASE_spectrum_name,
		string N_name,
		string NF_name)
	{
		if (refined_points != 0 && refined_points <= 1)
		{
			std::cout << "ModelEngine::output: incorrect number of refined points = " << refined_points << "\n";
			wait_and_exit();
		}

		std::ofstream f_pump_out = set_output_stream(f_pump_name);
		std::ofstream b_pump_out = set_output_stream(b_pump_name);
		std::ofstream f_signal_out = set_output_stream(f_signal_name);
		std::ofstream b_signal_out = set_output_stream(b_signal_name);
		std::ofstream f_ASE_out = set_output_stream(f_ASE_name);
		std::ofstream b_ASE_out = set_output_stream(b_ASE_name);
		std::ofstream f_ASE_total_out = set_output_stream(f_ASE_total_name);
		std::ofstream b_ASE_total_out = set_output_stream(b_ASE_total_name);
		std::ofstream f_ASE_spectrum_out = set_output_stream(f_ASE_spectrum_name);
		std::ofstream b_ASE_spectrum_out = set_output_stream(b_ASE_spectrum_name);
		std::ofstream N_out = set_output_stream(N_name);
		std::ofstream NF_out = set_output_stream(NF_name);

		F.output(refined_points, signal_output, f_pump_out, f_signal_out, f_ASE_out, f_ASE_total_out, f_ASE_spectrum_out);
		B.output(refined_points, signal_output, b_pump_out, b_signal_out, b_ASE_out, b_ASE_total_out, b_ASE_spectrum_out);

		// N output ===============================================================================
		double dz = L / (double)(refined_points - 1);

		for(int i = 0; i <= n_level; i++)
			N_out << i << (i != n_level ? "\t" : "\n");

		if (refined_points == 0)
		{
			for(int i = 0; i < F.n_z; i++)
			{
				double _z = F.z[i];
				N_out << _z;
				calculate_W_val(F.P[i], ODE_interpolate(_z, B.z, B.P, B.dense_coeff, B.n_z));

				n_t = 0;
				t = new double[segment_size];
				N = new N_vector[segment_size];
				ODE_solve_RK45(*this, & ModelEngine::dN, 0, infinity, N_previous_step, n_t, t, N);

				for (int  j = 0; j < n_level; j++)
				{
					double N_a = Grid<r_nodes, phi_nodes>::integrate_uniform(N[n_t - 1][j], r_core * DRF) / (pi * (r_core * DRF) * (r_core * DRF));
					N_out << "\t" << N_a;
				}
				N_out << "\t" << N_activator;
				if (i != F.n_z - 1) N_out << "\n";

				delete [] t;
				delete [] N;
			}
		}
		else
		{
			for(int i = 0; i < refined_points; i++)
			{
				double _z = i * dz;
				N_out << _z;
				calculate_W_val(ODE_interpolate(_z, F.z, F.P, F.dense_coeff, F.n_z), ODE_interpolate(_z, B.z, B.P, B.dense_coeff, B.n_z));

				n_t = 0;
				t = new double[segment_size];
				N = new N_vector[segment_size];
				ODE_solve_RK45(*this, & ModelEngine::dN, 0, infinity, N_previous_step, n_t, t, N);

				for (int  j = 0; j < n_level; j++)
				{
					double N_a = Grid<r_nodes, phi_nodes>::integrate_uniform(N[n_t - 1][j], r_core * DRF) / (pi * (r_core * DRF) * (r_core * DRF));
					N_out << "\t" << N_a;
				}
				N_out << "\t" << N_activator;
				if (i != F.n_z - 1) N_out << "\n";

				delete [] t;
				delete [] N;
			}
		}

		// Noise figure (NF) output
		int ch_num = 0;
		int *ch_list = new int [sufficiently_big_memory_chunk];

		double wl_min = infinity,
			   wl_max = 0.;
		for(int i = 0; i < F.n_ch; i++)
			if( F.type[i] == monochromatic && F.handle[i].find("signal") != string::npos)
			{
				ch_list[ch_num++] = i;

				if(F.wl[i] > wl_max) wl_max = F.wl[i];
				if(F.wl[i] < wl_min) wl_min = F.wl[i];
			};

		if (ch_num <= 1)
			std::cout << "\nWith " << ch_num << " signal channel(s) no NF can be calculated\n";
		else
		{
			// first we interpolate gain values to fit ASE slots
			double *x = new double [sufficiently_big_memory_chunk],
				   *y = new double [sufficiently_big_memory_chunk];

			for(int i = 0; i < ch_num; i++)
			{
				x[i] = F.wl[ch_list[i]];
				y[i] = F.P[F.n_z - 1][ch_list[i]] / F.P[0][ch_list[i]];
			}
			
			SplineInterpolation gain_interp(x, y, ch_num);

			delete [] x;
			delete [] y;

			NF_out << 0.;

			for(int i = 0; i < F.n_ch; i++)
				if (F.type[i] == ASE && F.wl[i] >= wl_min && F.wl[i] <= wl_max)
				{
					double G = gain_interp.interpolate(F.wl[i]);
					// std::cout << G << "\n";
					NF_out << "\n" << F.wl[i] << "\t"
						<< 10 * log10( (1 +  /* 2 * */ F.P[F.n_z - 1][i] / (h_planck * F.fr[i] * F.fr_width[i]) ) / G);		// Becker, eq.7.76, p.222, also p. 443
				}
		}	
			
		delete [] ch_list;

		f_pump_out.close();
		b_pump_out.close();
		f_signal_out.close();
		b_signal_out.close();
		f_ASE_out.close();
		b_ASE_out.close();
		f_ASE_total_out.close();
		b_ASE_total_out.close();
		N_out.close();
		NF_out.close();
	}	// ModelEngine::output


#endif