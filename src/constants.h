#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

namespace fasim
{
	const int segment_size = 1000;						// used for memory allocation in ode routines
	const int sufficiently_big_memory_chunk = 10000;	// used for data input/output memory allocation

	const double infinity = std::numeric_limits<double>::max();			// infinity integration flag
	const double neg_infinity = std::numeric_limits<double>::min();		// negative infinity integration flag

	const double default_zero_tolerance = 1.e-7;		// ODE_solve_RK45 infinity integration parameter
	
	const int max_integration_steps = 10000;			// used for infinite integrations
	const int ODE_dense_num = 6;						// dense output interpolation coefficients number

	const double cs_sensitivity_threshold = 1.e-2;		// cross-section is considered to be zero, if it is lower than <peak_cs> * cs_s_threshold

	const double C0 = 299792458;						// lightspeed in vacuum
	const double h_planck = 6.62606957e-34;				// Planck constant
	const double pi = 3.14159265358979323846;			// guess what?
	const double sqrt2 = 1.41421356237310;

	const int max_mode_index = 5;						// maximal amount of modes

	const double magic_L = 10.;							// default amplifier length for SolutionContainer generation
}

#endif