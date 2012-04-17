#include "vector.h"

namespace fasim
{

double scalar_norm(double op)
{
	return fabs(op);
}

double linear_combination2(double v0, double k1, double v1)
{
	return v0 + k1 * v1;
}

double linear_combination3(double v0, double k1, double v1, double k2, double v2)
{
	return v0 + k1 * v1 + k2 * v2;
}

double linear_combination4(double v0, double k1, double v1, double k2, double v2, double k3, double v3)
{
	return v0 + k1 * v1 + k2 * v2 + k3 * v3;
}

double linear_combination5(double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4)
{
	return v0 + k1 * v1 + k2 * v2 + k3 * v3 + k4 * v4;
}

double linear_combination6(double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4, double k5, double v5)
{
	return v0 + k1 * v1 + k2 * v2 + k3 * v3 + k4 * v4 + k5 * v5;
}

double linear_combination7(double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4, double k5, double v5, double k6, double v6)
{
	return v0 + k1 * v1 + k2 * v2 + k3 * v3 + k4 * v4 + k5 * v5 + k6 * v6;
}

double linear_combination7(double k0, double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4, double k5, double v5, double k6, double v6)
{
	return k0 * v0 + k1 * v1 + k2 * v2 + k3 * v3 + k4 * v4 + k5 * v5 + k6 * v6;
}

}