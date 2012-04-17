#ifndef MODEL_CPP
#define MODEL_CPP

//#define SHOW_CALL_INFO

#include "model.h"


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
	typename ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::F_subsystem::Channel 
	ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::dPF(double _z, const typename F_subsystem::Channel & PF)
{
	//std::cout << "z = " << _z << "\r";
	dPF_called++;

	B_subsystem::Channel _PB = ODE_interpolate(_z, B.z, B.P, B.dense_coeff, B.n_z);
	calculate_W_val(PF, _PB);

    if(exact_stationary_solution_possible)
    {
        N_previous_step = N_stationary();
    }
    else
    {
	    t = new double[segment_size];
	    N = new N_vector[segment_size];
	    n_t = 0;

	    dN_called = 0;
	
	    ODE_solve_RK45(*this, & ModelEngine::dN, 0, infinity, N_previous_step, n_t,	t, N);
	
	    dN_average_called += dN_called;
	    if(dN_max_called < dN_called) dN_max_called = dN_called;
	    if(dN_min_called > dN_called) dN_min_called = dN_called;
	
	    N_previous_step = N[n_t - 1];

	    delete [] t;
	    delete [] N;
    }

	F_subsystem::Channel dPF;
	for(int k = 0; k < F_subsystem::n_ch; k++)
	{
		typename Grid<r_nodes, phi_nodes>::Values G;
		G = 0.;
		for(int q = 0; q < F.n_nonzero_cs[k]; q++)
		{
			int i = F.nonzero_cs_from[q][k];
			int j = F.nonzero_cs_to[q][k];
			
			if( F.direction[i][j] == emission )
			{
				if(F.type[k] == ASE)
					G += ((PF[k] + 2 * h_planck * F.fr[k] * F.fr_width[k]) * F.cs[i][j][k]) * N_previous_step[i];
				else
					G += (PF[k] * F.cs[i][j][k]) * N_previous_step[i];
			}
			else
				G -= (PF[k] * F.cs[i][j][k]) * N_previous_step[i];
		}
		
		//std::cout << Grid<r_nodes, phi_nodes>::integrate_uniform(F.mode_envelope[k], r_core * DRF) << "\n";

		dPF[k] = Grid<r_nodes, phi_nodes>::integrate_uniform(G * F.mode_envelope[k], r_core * DRF);
	}

	return dPF;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
	typename ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::B_subsystem::Channel 
	ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::dPB(double _z, const typename B_subsystem::Channel & PB)
{
	//std::cout << "z = " << _z << "\r";

	dPB_called++;

	F_subsystem::Channel _PF = ODE_interpolate(_z, F.z, F.P, F.dense_coeff, F.n_z);
	calculate_W_val(_PF, PB);

    if(exact_stationary_solution_possible)
    {
        N_previous_step = N_stationary();
    }
    else
    {
	    t = new double[segment_size];
	    N = new N_vector[segment_size];
	    n_t = 0;

	    dN_called = 0;

	    ODE_solve_RK45(*this, & ModelEngine::dN, 0, infinity, N_previous_step, n_t,	t, N);
	    
	    dN_average_called += dN_called;
	    if(dN_max_called < dN_called) dN_max_called = dN_called;
	    if(dN_min_called > dN_called) dN_min_called = dN_called;	
	    
	    N_previous_step = N[n_t - 1];

	    delete [] t;
	    delete [] N;
    }

    N_previous_step = N_stationary();

	B_subsystem::Channel dPB;
	for(int k = 0; k < B_subsystem::n_ch; k++)
	{
		typename Grid<r_nodes, phi_nodes>::Values G;
		G = 0.;

		for(int q = 0; q < B.n_nonzero_cs[k]; q++)
		{
			int i = B.nonzero_cs_from[q][k];
			int j = B.nonzero_cs_to[q][k];
			
			if( B.direction[i][j] == emission )
			{
				if(B.type[k] == ASE)
					G += ((PB[k] + 2 * h_planck * B.fr[k] * B.fr_width[k]) * B.cs[i][j][k]) * N_previous_step[i];
				else
					G += (PB[k] * B.cs[i][j][k]) * N_previous_step[i];
			}
			else
				G -= (PB[k] * B.cs[i][j][k]) * N_previous_step[i];
		}
		
		dPB[k] = -1. * Grid<r_nodes, phi_nodes>::integrate_uniform(G * B.mode_envelope[k], r_core * DRF);
	}
	return dPB;
}


template <int n_level, int n_F_mono, int n_B_mono, int n_ASE, int r_nodes, int phi_nodes>
bool ModelEngine< n_level, n_F_mono, n_B_mono, n_ASE, r_nodes, phi_nodes>::
	relaxation(double abs_tolerance, double rel_tolerance, int max_iterations, ModelSolution sol)
{
	cout.precision(8);
	//cout.setf(std::ios::showpoint);

	F.free();
	B.free();
	N_previous_step = 0.;
	N_previous_step[0] = N_activator;

	bool F_convergence = false,
		 B_convergence = false,
		 FB_convergence =  false,
		 converged = true;
	int iter = 0;
	
	// preset B according to zero backward powers assumption
	B.n_z = sol.B.n_z;
	B.z = new double [B.n_z];
	B.P = new B_subsystem::Channel [B.n_z];
	B.dense_coeff = new B_subsystem::Dense_coefficients [B.n_z];
	
	double K_scale = L / fabs(sol.B.z[sol.B.n_z - 1] - sol.B.z[0]);

	for(int i = 0; i < B.n_z; i++)
	{
		B.z[i] = K_scale * sol.B.z[i];					
														
		B.P[i] = sol.B.P[i];
		B.dense_coeff[i] = (1. / K_scale) * sol.B.dense_coeff[i];
	}
	B.z[0] += L * 1.e-8;					// z interval is widened to avoid error accumulation issue
	B.z[B.n_z - 1] -= L * 1.e-8;			// ... so that (L, 0) belongs to (B.z[0], B.z[B.z_n - 1])

	//B.n_z = 10;
	//B.z = new double [B.n_z];
	//B.P = new B_subsystem::Channel [B.n_z];
	//B.dense_coeff = new B_subsystem::Dense_coefficients [B.n_z];
	//
	//for(int i = 0; i < B.n_z; i++)
	//{
	//	B.z[i] = L - i * (L + 1.e-9 * L) / (double) (B.n_z - 1);					
	//	B.P[i] = 0.;
	//	B.dense_coeff[i] = 0.;
	//}


	F_subsystem::Channel F_previous_step, F_err, F_threshold, F_error_scale;
	F_threshold = abs_tolerance / rel_tolerance;
	
	B_subsystem::Channel B_previous_step, B_err, B_threshold, B_error_scale;
	B_threshold = abs_tolerance / rel_tolerance;
	double error;

	do
	{
		if(iter > 0)
		{
			delete [] F.z;
			delete [] F.P;
			delete [] F.dense_coeff;
		}
		
		F.z = new double [segment_size];
		F.P = new F_subsystem::Channel [segment_size];
		F.dense_coeff = new F_subsystem::Dense_coefficients [segment_size];
		F.n_z = 0;
		
		dPF_called = dN_average_called = dN_max_called = 0;
		dN_min_called = 10000000;

		ODE_solve_RK45(*this, &ModelEngine::dPF, 0., L, F.P_boundary, F.n_z, F.z, F.P, 1.e-9, 1.e-9, 1.e-5, &F.dense_coeff);

		if(max_iterations == 0) break;			// max_iterations == 0 corresponds to 1 forward integration

		if (iter > 0)
		{
			F_err = F.P[F.n_z - 1] - F_previous_step;
			F_error_scale = vector_norm(F_err) / max_abs(F_threshold, max_abs(F.P[F.n_z - 1], F_previous_step));
			error = scalar_norm_inf(F_error_scale);

			if(error <= rel_tolerance)
				F_convergence = true;
			else
				F_convergence = false;

            // output
//			std::cout << "\trelaxation iteration " << iter + 1 << "\tF rel error = " << error / rel_tolerance;
//#ifdef SHOW_CALL_INFO
//            std::cout << "\t\t\tcalls: dPF " << dPF_called << "\tdN_min " << dN_min_called << "\tdN_max " << dN_max_called << "\tdN_average " << dN_average_called / dPF_called << "\n";
//#else
//			std::cout << "\n";
//#endif
		}
		
		F_previous_step = F.P[F.n_z - 1];

		delete [] B.z;
		delete [] B.P;
		delete [] B.dense_coeff;

		B.z = new double [segment_size];
		B.P = new B_subsystem::Channel [segment_size];
		B.dense_coeff = new B_subsystem::Dense_coefficients [segment_size];
		B.n_z = 0;

		dPB_called = dN_average_called = dN_max_called = 0;
		dN_min_called = 10000000;
		
		ODE_solve_RK45(*this, &ModelEngine::dPB, L, 0., B.P_boundary, B.n_z, B.z, B.P, 1.e-9, 1.e-9, 1.e-5, &B.dense_coeff);

		if (iter > 0)
		{
			B_err = B.P[B.n_z - 1] - B_previous_step;
			B_error_scale = vector_norm(B_err) / max_abs(B_threshold, max_abs(B.P[B.n_z - 1], B_previous_step));
			error = scalar_norm_inf(B_error_scale);

			if(error <= rel_tolerance)
				B_convergence = true;
			else
				B_convergence = false;

        // output
//			std::cout << "\t\t\t\tB rel error = " << error / rel_tolerance;
//#ifdef SHOW_CALL_INFO
//            std::cout << "\t\t\tcalls: dPB " << dPB_called << "\tdN_min " << dN_min_called << "\tdN_max " << dN_max_called << "\tdN_average " << dN_average_called / dPB_called << "\n";
//#else
//			std::cout << "\n";
//#endif
		}

		B_previous_step = B.P[B.n_z - 1];

		FB_convergence = F_convergence && B_convergence;

		iter++;
		if (iter == max_iterations - 1) 
			converged = false;
	} while (iter < max_iterations && !FB_convergence);
	
	return converged;
}


#endif