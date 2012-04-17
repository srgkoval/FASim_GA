#ifndef MODEL_SPECIFICATION_H
#define MODEL_SPECIFICATION_H

#include "boost\math\special_functions\fpclassify.hpp"

#include "model.h"

// GRID PARAMETERS
// r_nodes = phi_nodes = 1 corresponds to overlap factor approximation
// phi_nodes = 1 - axial symmetry

const int r_nodes = 1; //21;
const int phi_nodes = 1;


template <int F_mono_num, int B_mono_num, int ASE_num> class ThuliumModel: public ModelEngine<5, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>
{
public:
	ThuliumModel(): 
	  ModelEngine<5, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>::ModelEngine()
      {
          exact_stationary_solution_possible = true;
      }	
	
    N_vector N_stationary()
    {
		N_vector N;

        double
			E3 = 1 / 650.e-6,
			An1 = 1 / 430.e-6,
			An3 = 1 / (45.e-6),
//			An5 = 1 / 784.e-6,
			A31 = E3 * 0.09,
			A32 = E3 * 0.03,
			A10 = 1 / 3500.e-6;

        double C10, C11, C13,
               C30, C31, C33,
               common_denominator;

		for(int i = 0; i < r_nodes * phi_nodes; i++)
		{
            C10 = W(0,1)[i] + W(0,2)[i];        C11 = - (W(1,0)[i] + W(1,3)[i] + W(1,4)[i]) - (An1 + A10);    C13 = W(3,1)[i] + An3 + A31 + A32;
            C30 = W(0,3)[i];                    C31 = W(1,3)[i] + W(1,4)[i];                              C33 = - (An3 + E3) - (W(3,1)[i] + W(3,0)[i]);
            common_denominator = 1. / ( C13 * (-C30 + C31) + C11 * (C30 - C33) + C10 * (-C31 + C33) );

            N[1][i] = -1. * (C13 * C30 - C10 * C33) * N_activator * common_denominator;
            N[2][i] = 0.;
            N[3][i] = (C11 * C30 - C10 * C31) * N_activator * common_denominator;
            N[4][i] = 0.;
            N[0][i] = N_activator - N[1][i] - N[3][i];
		}

        return N;
    }


	N_vector dN(double z, const N_vector & N)
	{
		dN_called++;
		
		N_vector dN;
	
		double
			E3 = 1 / 650.e-6,
			An1 = 1 / 430.e-6,
			An3 = 1 / (45.e-6),
//			An5 = 1 / 784.e-6,
			A31 = E3 * 0.09,
			A32 = E3 * 0.03,
			A10 = 1 / 3500.e-6;

		for(int i = 0; i < r_nodes * phi_nodes; i++)
		{
			dN[1][i] = (W(0,1)[i] + W(0,2)[i]) * N[0][i] - (W(1,0)[i] + W(1,3)[i] + W(1,4)[i]) * N[1][i] 
				- (An1 + A10) * N[1][i] + (W(3,1)[i]) * N[3][i] + (An3 + A31 + A32) * N[3][i]; 
			dN[2][i] = 0.;
			dN[3][i] = W(0,3)[i] * N[0][i] + (W(1,3)[i] + W(1,4)[i]) * N[1][i] - (An3 + E3) * N[3][i] - (W(3,1)[i] + W(3,0)[i]) * N[3][i];
			dN[4][i] = 0.;
			dN[0][i] = (-1.) * (dN[1][i] + dN[3][i]);

            if( boost::math::isnan<double>(dN[1][i]))
                std::cin.get();
		}
		return dN;
	}
};


template <int F_mono_num, int B_mono_num, int ASE_num> class ErbiumModel: public ModelEngine<2, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>
{
public:
	ErbiumModel(): 
	  ModelEngine<2, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>::
		  ModelEngine()	{};	
	
    N_vector N_stationary()
    {
		N_vector N;
        N = 0.;
        return N;
    }


	N_vector dN(double z, const N_vector & N)
	{
		dN_called++;

		N_vector dN;
		double A1 = 1. / 10.e-3;
		
		//double t_uc = 5.e-3;
		//double C_uc = 1. / (t_uc * N_activator);
	
		dN[0] = W(1,0) * N[1] - W(0,1) * N[0] + A1 * N[1]; // + C_uc * N[1] * N[1];
		dN[1] = (-1.) * dN[0];	

		return dN;
	}
};


template <int F_mono_num, int B_mono_num, int ASE_num> class ErbiumESAModel: public ModelEngine<3, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>
{
public:
	ErbiumESAModel(): 
	  ModelEngine<3, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>::
		  ModelEngine()	{};	

    N_vector N_stationary()
    {
		N_vector N;
        N = 0.;
        return N;
    }

	N_vector dN(double _t, const N_vector & N)
	{
		//std::cout << "t = " << _t << "\r";
		
		dN_called++;
		
		N_vector dN;
		double A1 = 1. / 10.e-3;
		double A2 = 1. / 0.58e-5;
		for(int i = 0; i < r_nodes * phi_nodes; i++)
		{
			dN[0][i] = W(1,0)[i] * N[1][i] - W(0,1)[i] * N[0][i] + A1 * N[1][i];
			//dN[1] = (-1.) * dN[0];	
			//dN[1] = W(0,1) * N[0] - W(1,0) * N[1] - A1 * N[1] - W(1,2) * N[1] + A2 * N[2];
			dN[2][i] = W(1,2)[i] * N[1][i] - A2 * N[2][i];
			//dN[2] = 0.; //(-1.) * (dN[0] + dN[1]);
			dN[1][i] = (-1.) * (dN[0][i] + dN[2][i]);
		}
		return dN;
	}
};


template <int F_mono_num, int B_mono_num, int ASE_num> class ErbiumUpconversionModel: public ModelEngine<5, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>
{
public:
	ErbiumUpconversionModel(): 
	  ModelEngine<5, F_mono_num, B_mono_num, ASE_num, r_nodes, phi_nodes>::
		  ModelEngine()	{};	

    N_vector N_stationary()
    {
		N_vector N;
        N = 0.;
        return N;
    }


	N_vector dN(double z, const N_vector & N)
	{
		dN_called++;

		N_vector dN;
		double
			A1 = 1. / 9.0e-3,
			A2 = 1. / 6.9e-3,
			A3 = 1. / 5.0e-6,
			A4 = 1. / 0.58e-3,
			C_ETU1 = 6.7e-23,
			C_ETU2 = 1.9e-23,
			C_CR = 3.3e-23;
		
		for(int i = 0; i < r_nodes * phi_nodes; i++)
		{	
			dN[1][i] = W(0,1)[i] * N[0][i] - (W(1,3)[i] + W(1,0)[i]) * N[1][i] - A1 * N[1][i] + A2 * N[2][i]	- 2 * C_ETU1 * N[1][i] * N[1][i] + C_CR * N[0][i] * N[4][i];
			dN[2][i] = W(0,2)[i] * N[0][i] - W(2,4)[i] * N[2][i] - A2 * N[2][i] + A3 * N[3][i]					- 2 * C_ETU2 * N[2][i] * N[2][i];
			dN[3][i] = W(1,3)[i] * N[1][i] - A3 * N[3][i] + A4 * N[4][i]										+ C_ETU1 * N[1][i] * N[1][i] + C_CR * N[0][i] * N[4][i];
			dN[4][i] = W(0,4)[i] * N[0][i] + W(2, 4)[i] * N[2][i] - W(4,0)[i] * N[4][i] - A4 * N[4][i]			+ C_ETU2 * N[2][i] * N[2][i] - C_CR * N[0][i] * N[4][i];
			dN[0][i] = (-1.) * (dN[1][i] + dN[2][i] + dN[3][i] + dN[4][i]);
		}
		
		return dN;
	}
};


#endif