#include "interpolation.h"
#include <iostream>
#include "support.h"

using namespace fasim;


SplineInterpolation::SplineInterpolation(double *_x, double *_y, int _n, bool _zero_if_out_of_range):
	n(_n), omit_warnings(false), zero_if_out_of_range(_zero_if_out_of_range)
	{
		x = new double [n];
		y = new double [n];
		y2 = new double [n];
		u = new double [n];

		for(int i = 0; i < n; i++)
		{	
			x[i] = _x[i];
			y[i] = _y[i];
		}

		y2[0] = u[0] = 0.0;

		for(int i = 0; i < n - 1; i++)
		{
			if(x[i] == x[i+1])
			{
				std::cout << "SplineInterpolation::SplineInterpolation: two equal {x}=" << x[i]
					<< " entries in the function table\n";
				wait_and_exit();
			}
			if(x[i] >= x[i+1])
			{
				std::cout << "SplineInterpolation::SplineInterpolation: unsorted x array at i = " 
					<< i << ", x = " << x[i] << "\n";
				wait_and_exit();
			}
		}

		for (int  i = 1; i < n - 1;  i++)
		{
			double sgn = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
			double  p = sgn * y2[i-1] + 2.0;
			y2[i] = (sgn - 1.0) / p;
			u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
			u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sgn * u[i-1]) / p;
		}
		y2[n-1] = 0;	

		for (int i = n - 2; i >= 0 ; i--)
			y2[i] = y2[i] * y2[i+1] + u[i];
		
		// debug output
		// for (int i = 0; i < n; i++)
		// std::cout << i + 1 << "\t" << x[i] << "\t" << y[i] << "\t" <<
		//	y2[i] << "\t" << u[i] << "\n";

	};


SplineInterpolation::~SplineInterpolation()
{
	delete [] x;
	delete [] y;
	delete [] y2;
	delete [] u;
}


SplineInterpolation::SplineInterpolation(const SplineInterpolation &s)
{
	n = s.n;
	omit_warnings = false;
	zero_if_out_of_range = s.zero_if_out_of_range;

	x = new double [n];
	y = new double [n];
	y2 = new double [n];
	u = new double [n];

	for(int i = 0; i < n; i++)
	{	
		x[i] = s.x[i];
		y[i] = s.y[i];
		y2[i] = s.y2[i];
		u[i] = s.u[i];
	}
}


const SplineInterpolation & SplineInterpolation::operator=(const SplineInterpolation &s)
{
	if( this != &s)
	{
		delete [] x;
		delete [] y;
		delete [] y2;
		delete [] u;

		n = s.n;
		omit_warnings = false;
		zero_if_out_of_range = s.zero_if_out_of_range;

		x = new double [n];
		y = new double [n];
		y2 = new double [n];
		u = new double [n];

		for(int i = 0; i < n; i++)
		{	
			x[i] = s.x[i];
			y[i] = s.y[i];
			y2[i] = s.y2[i];
			u[i] = s.u[i];
		}
	}
	return *this;
}


double SplineInterpolation::interpolate(double xx)
{
	if( xx > x[n-1] || xx < x[0])
	{
		if (zero_if_out_of_range)
		{
#ifdef VERBOSE_INTERPOLATION
			if (!omit_warnings)
			{
				std::cout << "Warning! SplineInterpolation::interpolate: x = " << xx << " is out of range (" << x[0] << ", " << x[n-1] << "), n = " << n << ". Zero is assumed\n";
				std::cout << "\t...further warnings are omitted\n";
				omit_warnings = true;
			};
#endif
			return 0;
		}
		else
		{
			std::cout << "SplineInterpolation::interpolate: x = " << xx << " is out of range (" << x[0] << ", " << x[n-1] << "), n = " << n << "\n";
			wait_and_exit();
		}
	}

	int jl, jm , ju;
	jl = 0;
	ju = n - 1;

	// some binary search
	while (ju - jl > 1)
	{
		jm = (jl + ju) / 2;
		if (xx > x[jm])
			jl = jm;
		else
			ju = jm;
	}

	int kl = jl, kh = jl + 1;
	double h = x[kh] - x[kl];
	
	if (h == 0.0)
		wait_and_exit("SplineInterpolation::interpolate: two equal {x} entries in the function table");

	double A = (x[kh] - xx) / h;
	double B = (xx - x[kl]) / h;
	double res = A * y[kl] + B * y[kh] + ( A * (A * A - 1) * y2[kl] + B * (B * B - 1) * y2[kh]) * h * h / 6.0;

	////debug code
	//std::cout << "bin_search i = " << jl << " for x = " << xx << ", res" << res 
	//	<< " in {(" << x[kl] << ", " << y[kl] << "), (" << x[kh] << ", " << y[kh] << ")}\n";
	
	return res;
}
