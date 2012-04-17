#ifndef INTERPOLATION_H
#define INTERPOLATION_H

//#define VERBOSE_INTERPOLATION

namespace fasim{

class SplineInterpolation
{
public:
	double *x, *y;			// x and y array that represent the function to be interpolated
	int n;					// x and y size
	double *y2, *u;			// auxiliary vectors
	bool omit_warnings;		// prevents warning flood
	bool zero_if_out_of_range;	// if set to false, will report an error if x is out of range, if true will return zero
	
	SplineInterpolation(double *_x, double *_y, int _n, bool _zero_if_out_of_range = true);
	~SplineInterpolation();

	SplineInterpolation(const SplineInterpolation &s);
	const SplineInterpolation & operator=(const SplineInterpolation &s);

	double interpolate(double xx);
	
	double operator()(double xx)
	{
		return interpolate(xx);
	}
};

}

#endif