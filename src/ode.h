#ifndef ODE_H
#define ODE_H

//#define ODE_DEBUG
//#define ODE_DEBUG_INFINITE
//#define ODE_DEBUG_LAST_STEP
//#define NON_NEGATIVE

#include <cmath>

#include "vector.h"
#include "support.h"
#include "constants.h"

using namespace fasim;


template <class C, class X, int x_size> void ODE_solve_RK45(C &obj, Vector<X, x_size> (C::*ODE_fun)(double , const Vector<X, x_size> &),
	double z0, double z1, const Vector<X, x_size> &y0, int &grid_points, double *&z, Vector<X, x_size> *&y,
	double abs_tolerance = 1.e-9, double rel_tolerance = 1.e-9, double zero_tolerance = default_zero_tolerance,
	Vector< Vector<X, x_size>, 6> **dense_coeff = NULL);

template <typename X, int x_size> Vector<X, x_size>
	ODE_interpolate(double x, double const *z, Vector<X, x_size> * const & y, Vector< Vector<X, x_size>, 6> * const &dense_coeff, int grid_points);

template <class X> void reallocate(X *& ,int , int );

template <class X, int x_size> void reallocate(Vector<X, x_size> *& ,int , int );


template <class C, class X, int x_size> void ODE_solve_RK45(C &obj, Vector<X, x_size> (C::*ODE_fun)(double , const Vector<X, x_size> &),
	double z0, double z1, const Vector<X, x_size> &y0, int &grid_points, double *&z, Vector<X, x_size> *&y,
	double abs_tolerance, double rel_tolerance, double zero_tolerance,
	Vector< Vector<X, x_size>, 6> **dense_coeff)
{
	double
	a21 = 1./5.,
	a31 = 3./40.,		a32 = 9./40.,
	a41 = 44./45.,		a42 =- 56./15.,			a43 = 32./9.,
	a51 = 19372./6561.,	a52 =- 25360./2187.,	a53 = 64448./6561.,		a54 =- 212./729.,
	a61 = 9017./3168.,	a62 =- 355./33.,		a63 = 46732./5247.,		a64 = 49./176.,		a65 =- 5103./18656.,
	a71 = 35./384.,		a72 = 0.,				a73 = 500./1113.,		a74 = 125./192.,	a75 =- 2187./6784.,		a76 = 11./84.,
		
	e1 = 71./57600.,	e2 = 0.,				e3 =- 71./16695.,		e4 = 71./1920.,		e5 =- 17253./339200.,	e6 = 22./525.,	e7 =- 1./40.,

	c2 = 1./5.,			c3 = 3./10.,			c4 = 4./5.,				c5 = 8./9.,			c6 = 1.,				c7 = 1.;
	
	double
		method_order = 5.,
		method_pow = 1.0 / method_order,
//		PI_alpha = 0.7 / method_order,
//		PI_beta = 0.4 / method_order,
		safety_factor = 0.8;
	
	if (abs_tolerance < 1.e-12 || rel_tolerance < 1.e-12 || abs_tolerance > 0.1 || rel_tolerance > 0.1)
		wait_and_exit("ODE_solve_RK45: all integration tolerances must lie in (0.1, 1.0e-12) range");

	bool infinite_integration = false;
	if ( z1 == infinity || z1 == neg_infinity)
		infinite_integration = true;
	
	double 
		threshold  = abs_tolerance / rel_tolerance,
		error,
		error_previous;

	Vector<double, x_size> 
		abs_tolerance_vector, 
		threshold_vector,
		error_scale;

	abs_tolerance_vector = abs_tolerance;
	threshold_vector = threshold;

	int	y_size = segment_size;
//	int steps_rejected = 0;

	Vector<X, x_size>	k1, k2, k3, k4, k5, k6, k7, ynew, yerr;

	y[0] = y0;
	z[0] = z0;

	if (z0 == z1)
	{
		grid_points = 1;
		if (dense_coeff)
			(*dense_coeff)[0][0] = (obj.*ODE_fun)( z[0], y[0] );
		return;
	}

	double z_current = z0;
	int i = 0;							// index of current step
	k7 = (obj.*ODE_fun)( z_current, y[i] );	
	
	int dir = (z1 > z0 ? 1: -1);		// integration direction: poisitve = 1, negative = -1 
	
	double hmax;
	if(!infinite_integration)
		hmax = abs(z1 - z0) / 10.;
	else
		hmax = 1.0e-1;					// hmax for INFINITE case
	double hmin = 16 * eps(z_current);

	double absh = hmax;
	double rh = scalar_norm( vector_norm(k7) / max_abs(threshold_vector, vector_norm(y0)) ) / (safety_factor * pow(rel_tolerance, method_pow));
	if (absh * rh > 1)
		absh = 1 / rh;
	absh = std::max(absh, hmin);
	
	// here we limit absh for dN calculations
	absh = std::max(absh, 1.e-6);
	
	double h = dir * absh;
	double h_scale;

	bool done = false;
//	bool previous_step_successful = true;

	// MAIN LOOPS <==================================================================================================================================
	
	// FINITE INTEGRATION ---------------------------------------------------------------------------------------------------------------------------
	if (!infinite_integration)
	{
		//integration loop
		while ( !done )
		{
			if (dir == 1)
			{
				if (z_current + h >= z1)	
				{	
					h = z1 - z_current;
					done = true;
				}
			}
			else
			{
				if (z_current + h <= z1)
				{
					h = z1 - z_current;
					done = true;
				}
			}

#ifdef ODE_DEBUG
			std::cout << i << "\th:" << h << "\tz:" << z_current << "\t" << y[i] << "\n";
#endif			
			k1 = k7;
			while (true) // attempt to advance a step
			{
				k2 = (obj.*ODE_fun)( z_current + c2 * h, linear_combination2(y[i], h * a21, k1) );
				k3 = (obj.*ODE_fun)( z_current + c3 * h, linear_combination3(y[i], h * a31, k1, h * a32, k2) );
				k4 = (obj.*ODE_fun)( z_current + c4 * h, linear_combination4(y[i], h * a41, k1, h * a42, k2, h * a43, k3) );
				k5 = (obj.*ODE_fun)( z_current + c5 * h, linear_combination5(y[i], h * a51, k1, h * a52, k2, h * a53, k3, h * a54, k4) );
				k6 = (obj.*ODE_fun)( z_current + c6 * h, linear_combination6(y[i], h * a61, k1, h * a62, k2, h * a63, k3, h * a64, k4, h * a65, k5) );
		
				ynew = linear_combination7(y[i], h * a71, k1, h * a72, k2, h * a73, k3, h * a74, k4, h * a75, k5, h * a76, k6);
				k7 = (obj.*ODE_fun)( z_current + c7 * h, ynew );

				yerr = linear_combination7(e1, k1, e2, k2, e3, k3, e4, k4, e5, k5, e6, k6, e7, k7);
				error_scale = vector_norm(yerr) / max_abs(threshold_vector, max_abs(ynew, y[i]));
				error = abs(h) * scalar_norm_inf(error_scale);

				if (error > rel_tolerance) // failed step
				{
					absh = abs(h);
					if (absh <= hmin)
					{					
						std::cout << "ODE_solve_RK45: unable to meet integration tolerances without reducing the step size \n below the smallest value allowed "
							<< hmin << " at z =" << z_current << std::endl;
						wait_and_exit();
					}

					h_scale = std::max(0.2, safety_factor * pow((rel_tolerance / error), method_pow));
					absh = std::max(hmin, h_scale * absh);
					h = dir * absh;
					done = false;

#ifdef ODE_DEBUG
				std::cout <<"\t\terror/r_tol = " << error / rel_tolerance 
					<< "\th decreased to " << h << ",\t steps rejected: " << ++steps_rejected << "\tk7 = " << k7 << std::endl;
#endif

				}
				else		// successful step
				{
					break;
				}
			} // while true

#ifdef ODE_DEBUG
				std::cout <<"\terror/r_tol = " << error / rel_tolerance 
					<< "\tk7 = " << k7 << "\n";
#endif

			// check if there is not enough memory to save the data
			if ( i == y_size - 1)
			{
				reallocate(y, y_size, segment_size);
				reallocate(z, y_size, segment_size);
				if (dense_coeff) 
					reallocate(*dense_coeff, y_size, segment_size);
				y_size += segment_size;
			}

			if (dense_coeff)
			{
				(*dense_coeff)[i][0] = k1;
				(*dense_coeff)[i][1] = k2;
				(*dense_coeff)[i][2] = k3;
				(*dense_coeff)[i][3] = k4;
				(*dense_coeff)[i][4] = k5;
				(*dense_coeff)[i][5] = k6;
			};

			i++;	// advance a step
			z_current += h;
			z[i] = z_current;

            y[i] = ynew;
#ifdef NON_NEGATIVE
            y[i].nullify_negative();
#endif
			error_previous = error;

			h_scale = (1 / safety_factor) * pow((error / rel_tolerance), method_pow);
			h_scale = std::max( 1. / 5., h_scale);
			h /= h_scale;
			absh = std::max(hmin, abs(h));
			h = dir * absh;
		}			// main while loop, finite integration case
	}
	else
	// INFINITE integration -------------------------------------------------------------------------------------------------------------------------
	{
		while ( !done )
		{

#ifdef ODE_DEBUG_INFINITE
			std::cout << i << "\th:" << h << "\tz:" << z_current << "\t" << y[i] << "\n";
#endif			
			k1 = k7;
			while (true) // attempt to advance a step
			{
				k2 = (obj.*ODE_fun)( z_current + c2 * h, linear_combination2(y[i], h * a21, k1) );
				k3 = (obj.*ODE_fun)( z_current + c3 * h, linear_combination3(y[i], h * a31, k1, h * a32, k2) );
				k4 = (obj.*ODE_fun)( z_current + c4 * h, linear_combination4(y[i], h * a41, k1, h * a42, k2, h * a43, k3) );
				k5 = (obj.*ODE_fun)( z_current + c5 * h, linear_combination5(y[i], h * a51, k1, h * a52, k2, h * a53, k3, h * a54, k4) );
				k6 = (obj.*ODE_fun)( z_current + c6 * h, linear_combination6(y[i], h * a61, k1, h * a62, k2, h * a63, k3, h * a64, k4, h * a65, k5) );
		
				ynew = linear_combination7(y[i], h * a71, k1, h * a72, k2, h * a73, k3, h * a74, k4, h * a75, k5, h * a76, k6);
				k7 = (obj.*ODE_fun)( z_current + c7 * h, ynew );

				yerr = linear_combination7(e1, k1, e2, k2, e3, k3, e4, k4, e5, k5, e6, k6, e7, k7);

				error_scale = vector_norm(yerr) / max_abs(threshold_vector, max_abs(ynew, y[i]));
				error = abs(h) * scalar_norm_inf(error_scale); 

				if (error > rel_tolerance) // failed step
				{
					absh = abs(h);
					if (absh <= hmin)
					{					
						std::cout << "ODE_solve_RK45: unable to meet integration tolerances without reducing the step size \n below the smallest value allowed "
							<< hmin << " at z =" << z_current << std::endl;
						wait_and_exit();
					}

					h_scale = std::max(0.2, safety_factor * pow((rel_tolerance / error), method_pow));
					absh = std::max(hmin, h_scale * absh);
					h = dir * absh;
					done = false;

#ifdef ODE_DEBUG_INFINTE
				std::cout <<"\t\terror/r_tol = " << error / rel_tolerance 
					<< "\th decreased to " << h << ",\t steps rejected: " << ++steps_rejected << "\tk7 = " << k7 << std::endl;
#endif

				}
				else		// successful step
				{
					break;
				}
			} // while true

#ifdef ODE_DEBUG_INFINITE
				std::cout <<"\t\terror/r_tol = " << error / rel_tolerance 
					<< "\tk7 = " << k7 << "\n";
#endif

			// check if there is not enough memory to save the data
			if ( i == y_size - 1)
			{
				reallocate(y, y_size, segment_size);
				reallocate(z, y_size, segment_size);
				if (dense_coeff) 
					reallocate(*dense_coeff, y_size, segment_size);
				y_size += segment_size;
			}

			if (dense_coeff)
			{
				(*dense_coeff)[i][0] = k1;
				(*dense_coeff)[i][1] = k2;
				(*dense_coeff)[i][2] = k3;
				(*dense_coeff)[i][3] = k4;
				(*dense_coeff)[i][4] = k5;
				(*dense_coeff)[i][5] = k6;
			};

			i++;	// advance a step
			z_current += h;
			z[i] = z_current;
			y[i] = ynew;
			error_previous = error;

			h_scale = (1 / safety_factor) * pow((error / rel_tolerance), method_pow);
			h_scale = std::max( 1. / 5., h_scale);
			h /= h_scale;
			absh = std::max(hmin, abs(h));
			h = dir * absh;

			int hindsight = 4;
				if (i > hindsight)
				{
					Vector<double, x_size> norm_current = vector_norm(y[i]);
					double rel_delta = std::min( scalar_norm_inf(norm_current), scalar_norm_inf( vector_norm(y[i] - y[i-hindsight]) / norm_current));
#ifdef ODE_DEBUG_INFINITE
					std::cout << "\t\tdelta = " << rel_delta << "\n";
#endif
					if ( rel_delta <= zero_tolerance)
						done = true;
				}
		
		};			// main while loop, infinite integration case
	}; // MAIN LOOP'S END ============================================================================================================================

#if defined ODE_DEBUG || defined ODE_DEBUG_INFINITE
#if defined ODE_DEBUG_LAST_STEP
			std::cout << i << "\th:" << h << "\tz:" << z_current << "\t" << y[i] << "\n";
#endif
#endif

	if (dense_coeff)
		(*dense_coeff)[i][0] = (obj.*ODE_fun)( z_current, y[i] );
	grid_points = i + 1;
}


// support functions ================================================================================================================================
// ==================================================================================================================================================

template <typename X, int x_size> Vector<X, x_size>
	ODE_interpolate(double x, const double *z, Vector<X, x_size> * const & y, Vector< Vector<X, x_size>, 6> * const &k, int grid_points)
{
	int// lower = 0,
		upper = grid_points - 1,
		mid;
	if(z[0] < z[grid_points - 1])
	{
		if( x < z[0] || x > z[upper] )
		{
			std::cout << "\nODE_interpolate: x = " << x << " is out of (" << z[0] << ", " << z[upper] << ")";
			wait_and_exit();
		}
		for(mid = 0; mid < grid_points; mid++)
			if( z[mid] > x) break;
	}
	else
	{
		if( x > z[0] || x < z[upper] )
		{
			std::cout << "\nODE_interpolate: x = " << x << " is out of (" << z[0] << ", " << z[upper] << ")";
			wait_and_exit();
		}
		for(mid = 0; mid < grid_points; mid++)
			if( z[mid] < x) break;
	}
	mid--;

	if(mid == grid_points - 1)
		return y[mid];

	// dense output (Hairer, Norsett, Wanner)
	// actually taken from Matlab ode45 dense output
	double absh = z[mid + 1] - z[mid];
	double s = (x - z[mid]) / absh,
		s2 = s * s,
		s3 = s2 * s,
		s4 = s3 * s;
	
	double b1, b2, b3, b4, b5, b6, b7;

	b1 = 1. * s	- 183./64. * s2	+ 37./12. * s3 - 145./128 * s4;
	b2 = 0.;
	b3 = 1500./371. * s2 - 1000./159 * s3 + 1000./371. * s4;
	b4 = -125./32. * s2 + 125./12. * s3 - 375./64. * s4;
	b5 = 9477./3392. * s2 - 729./106. * s3 + 25515./6784. * s4;
	b6 = -11./7. * s2 + 11./3. * s3 - 55./28. * s4;
	b7 = 3./2. *s2 - 4. * s3 + 5./2. * s4;

	return y[mid] + absh * 
		( b1 * k[mid][0] + /* b2 = 0*/ b3 * k[mid][2] + b4 * k[mid][3] + b5 * k[mid][4] + b6 * k[mid][5] + b7 * k[mid+1][0] );
}


template <class X> void reallocate(X *&y, int y_size, int y_additional_size)
{
	X *temp;
	temp = new X [y_size + y_additional_size];
	for ( int i = 0; i < y_size; i++)
		temp[i] = y[i];
	delete [] y;
	y = temp;
}

template <class X, int x_size> void reallocate(Vector<X, x_size> *&y, int y_size, int y_additional_size)
{
	Vector<X, x_size> *temp;
	temp = new Vector<X, x_size> [y_size + y_additional_size];
//	int dimension = y[0].size(); // и нафига это вообще тут было? странно
	for ( int i = 0; i < y_size; i++)
		temp[i] = y[i];
	delete [] y;
	y = temp;
}


#endif