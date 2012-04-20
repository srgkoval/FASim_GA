#ifndef TEST_GA_H
#define TEST_GA_H

#include <boost/math/constants/constants.hpp>
#include "ga.h"

using namespace fasim;

// single-objective problems

template<int N> double rastrigin(const typename GA<N>::Individual &x)
{
    double res = 0.;
    for(int i = 0; i < N; i++)
        res += 10. + x[i] * x[i] - 10. * pow(cos(2 * boost::math::constants::pi<double>() * x[i]), 2);

    return res;
}

template<int N> double schwefel(const typename GA<N>::Individual &x)
{
    double res = 0.;
    for(int i = 0; i < N; i++)
        res += -1. * x[i] * sin(sqrt(fabs(x[i])));
    
    return res;
}

template<int N> double rosenbrock(const typename GA<N>::Individual &x)
{
    double res = 0.;
    for(int i = 0; i < N - 1; i++)
        res += 100. * pow((x[i + 1] - x[i] * x[i]), 2) + pow(1 - x[i], 2);
    
    return res;
}


// multiobjective problems

Vector<double, 2> mymulti1(const GA<2, 2>::Individual &x)
{
    Vector<double, 2> res;
    res[0] = pow(x[0], 4) - 10.*x[0]*x[0] + x[0]*x[1] + pow(x[1], 4) - x[0]*x[0]*x[1]*x[1];
    res[1] = pow(x[1], 4) - x[0]*x[0]*x[1]*x[1] + pow(x[0], 4) + x[0]*x[1];

    return res;
}

template <int N> Vector<double, 2> KUR(const typename GA<N, 2>::Individual &x)
{
    Vector<double, 2> res;
    res = 0.;
    for(int i = 0; i < N; i++)
    {
        if(i != N - 1)    
            res[0] += -10. * exp(-0.2 * sqrt(x[i]*x[i] + x[i+1]*x[i+1]));
        res[1] += pow(fabs(x[i]), 0.8) + 5. *  sin(pow(x[i], 3));
    }

    return res;
}


const double KUR_scale_dim1 = 1.e-3,
             KUR_scale_dim2 = 1.e-1,
             KUR_scale_dim3 = 1.e3;

template <int N> Vector<double, 2> KUR_scaled(const typename GA<N, 2>::Individual & _x)
{
    Vector<double, 2> res;
    res = 0.;

    GA<N, 2>::Individual x(_x);
    x[0] *= KUR_scale_dim1;
    x[1] *= KUR_scale_dim2;
    x[2] *= KUR_scale_dim3;

    for(int i = 0; i < N; i++)
    {
        if(i != N - 1)    
            res[0] += -10. * exp(-0.2 * sqrt(x[i]*x[i] + x[i+1]*x[i+1]));
        res[1] += pow(fabs(x[i]), 0.8) + 5. *  sin(pow(x[i], 3));
    }

    return res;
}


template <int N> Vector<double, 2> ZDT1(const typename GA<N, 2>::Individual &x)
{
    Vector<double, 2> res;
    double g, h;
    res = g = h = 0.;

    for(int i = 1; i < N; i++)
        g += 9. / (N - 1) * x[i];
    g += 1.;
    h = 1. - sqrt(x[0] / g);

    res[0] = x[0];
    res[1] = g * h;

    return res;
}

template <int N> Vector<double, 2> ZDT2(const typename GA<N, 2>::Individual &x)
{
    Vector<double, 2> res;
    double g, h;
    res = g = h = 0.;

    for(int i = 1; i < N; i++)
        g += 9. / (N - 1) * x[i];
    g += 1.;
    h = 1. - pow((x[0] / g), 2);

    res[0] = x[0];
    res[1] = g * h;

    return res;
}

template <int N> Vector<double, 2> ZDT3(const typename GA<N, 2>::Individual &x)
{
    Vector<double, 2> res;
    double g, h;
    res = g = h = 0.;

    for(int i = 1; i < N; i++)
        g += 9. / (N - 1) * x[i];
    g += 1.;
    h = 1. - sqrt(x[0] / g) - x[0] / g * sin(10. * boost::math::constants::pi<double>() * x[0]);

    res[0] = x[0];
    res[1] = g * h;

    return res;
}

template <int N> Vector<double, 2> ZDT4(const typename GA<N, 2>::Individual &x)
{
    Vector<double, 2> res;
    double g, h;
    res = g = h = 0.;

    for(int i = 1; i < N; i++)
        g += x[i] * x[i] - 10. * cos(4. * boost::math::constants::pi<double>() * x[i]);
    g += 1. + 10. * (N - 1);
    h = 1. - sqrt(x[0] / g);

    res[0] = x[0];
    res[1] = g * h;

    return res;
}


template <int N> Vector<double, 2> ZDT6(const typename GA<N, 2>::Individual &x)
{
    Vector<double, 2> res;
    double g, h;
    res = g = h = 0.;

    for(int i = 1; i < N; i++)
        g += x[i] / double (N - 1);
    g += 1. + 9. * pow(g, 0.25);

    res[0] = 1. - exp(-4. * x[0]) * pow( sin(4. * boost::math::constants::pi<double>() * x[0]), 6);
    h = 1. - pow((res[0] / g), 2);    
    res[1] = g * h;

    return res;
}


//  test run routines =============================================================================

void run_MATLAB_example()
{
    const int n = 2;
    Vector<double, n> lower, upper;
    lower = -5.;
    upper = 5.;

    GA<n,2> ga;
    ga.run_multiobjective(mymulti1, lower, upper);
    ga.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt"); 
}


void run_KUR()
{
    const int n = 3;
    Vector<double, n> lower, upper;
    lower = -5.;
    upper = 5.;
    
    GA<n,2> kur;
    kur.run_multiobjective(KUR<n>, lower, upper);
    kur.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt"); 
}


void run_KUR_scaled()
{
    const int n = 3;
    Vector<double, n> lower, upper;
    lower = -5.;    lower[0] /= KUR_scale_dim1; lower[1] /= KUR_scale_dim2; lower[2] /= KUR_scale_dim3;
    upper = 5.;     upper[0] /= KUR_scale_dim1; upper[1] /= KUR_scale_dim2; upper[2] /= KUR_scale_dim3;
    
    GA<n,2> kur;
    kur.run_multiobjective(KUR_scaled<n>, lower, upper);
    kur.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt"); 
}


void run_ZDT1()
{
    const int n_ZDT1 = 30;
    GA<n_ZDT1,2> zdt1;

    Vector<double, n_ZDT1> lower_ZDT1, upper_ZDT1;
    lower_ZDT1 = 0.;
    upper_ZDT1 = 1.;

    zdt1.run_multiobjective(ZDT1<n_ZDT1>, lower_ZDT1, upper_ZDT1);

    zdt1.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt");
}


void run_ZDT2()
{
    const int n = 30;
    GA<n,2> zdt2;

    Vector<double, n> lower, upper;
    lower = 0.;
    upper = 1.;

    zdt2.run_multiobjective(ZDT2<n>, lower, upper);

    zdt2.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt");
}


void run_ZDT3()
{
    const int n = 30;
    GA<n,2> zdt3;

    Vector<double, n> lower, upper;
    lower = 0.;
    upper = 1.;

    zdt3.run_multiobjective(ZDT3<n>, lower, upper);

    zdt3.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt");
}


void run_ZDT4()
{
    const int n = 10;
    GA<n,2> zdt4;

    Vector<double, n> lower, upper;
    lower = -5.;    lower[0] = 0.;
    upper = 5.;     upper[0] = 1.;

    zdt4.run_multiobjective(ZDT4<n>, lower, upper);

    zdt4.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt");
}


void run_ZDT6()
{
    const int n = 10;
    GA<n,2> zdt6;

    Vector<double, n> lower, upper;
    lower = 0.;
    upper = 1.;

    zdt6.run_multiobjective(ZDT6<n>, lower, upper);

    zdt6.output("e:\\Programming - code and etc\\Genetic Algorithm\\Output\\ga.txt");
}


#endif