#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <fstream>
#include <stdarg.h>
#include "support.h"

#define RANGE_CHECK


namespace fasim{

double scalar_norm(double op);

template <typename X, int N> class Vector{
	int n;
	X data[N];

public:
	explicit Vector(): n(N) 
		{ 
			//for(int i = 0; i < N; i++)
			//	data[i] = 0;
		}

	explicit Vector(X v0, ...);				
	Vector(const Vector<X, N> &op);

	int size() const {return n;}

	X & operator[](int i); 
	const X & operator[](int i) const;

	Vector<X, N> & operator=(const Vector<X, N> & op);
	Vector<X, N> & operator=(double k);
	Vector<X, N> & operator=(int k);
	
	Vector<X, N> & operator+=(const Vector<X, N> &op);
	Vector<X, N> & operator-=(const Vector<X, N> &op);
	Vector<X, N> & operator*=(const Vector<X, N> &op);
	Vector<X, N> & operator/=(const Vector<X, N> &op);

	Vector<X, N> & operator*=(double k);
	Vector<X, N> & operator/=(double k);
	
	const Vector<X, N> operator+(const Vector<X, N> &op) const;
	const Vector<X, N> operator-(const Vector<X, N> &op) const;
	const Vector<X, N> operator*(const Vector<X, N> &op) const;
	const Vector<X, N> operator/(const Vector<X, N> &op) const;

	template <class Y, int y_size> friend const Vector<Y, y_size> operator*(double k, const Vector<Y, y_size> &op);
	template <class Y, int y_size> friend const Vector<Y, y_size> operator*(const Vector<Y, y_size> &op, double k);

	template <class Y, int y_size> friend double scalar_norm(const Vector<Y, y_size> &op);
	template <class Y, int y_size> friend double scalar_norm_inf(const Vector<Y, y_size> &op);				// returns the largest absolute value of all vector elements
	template <class Y, int y_size> friend Vector<double, y_size> vector_norm(const Vector<Y, y_size> &op);

	template <class Y, int y_size> friend Vector<double, y_size> max_abs(const Vector<Y, y_size> &op1, const Vector<Y, y_size> &op2);

	template <class Y, int y_size> friend Y sum(const Vector<Y, y_size> & op);

    bool operator==(const Vector<X, N> &op);
    bool operator!=(const Vector<X, N> &op);
    bool operator<(const Vector<X, N> &op);
    bool operator>(const Vector<X, N> &op);
    bool operator<=(const Vector<X, N> &op);
    bool operator>=(const Vector<X, N> &op);

    bool operator==(double op);
    bool operator!=(double op);
    bool operator<(double op);
    bool operator>(double op);
    bool operator<=(double op);
    bool operator>=(double op);

    void nullify_negative();
};


double linear_combination2(double v0, double k1, double v1);
double linear_combination3(double v0, double k1, double v1, double k2, double v2);
double linear_combination4(double v0, double k1, double v1, double k2, double v2, double k3, double v3);
double linear_combination5(double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4);
double linear_combination6(double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4, double k5, double v5);
double linear_combination7(double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4, double k5, double v5, double k6, double v6);
double linear_combination7(double k0, double v0, double k1, double v1, double k2, double v2, double k3, double v3, double k4, double v4, double k5, double v5, double k6, double v6);

template <typename T> inline T linear_combination2(T v0, double k1, T v1);
template <typename T, int n> inline Vector<T,n> linear_combination3(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2);
template <typename T, int n> inline Vector<T,n> linear_combination4(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3);
template <typename T, int n> inline Vector<T,n> linear_combination5(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4);
template <typename T, int n> inline Vector<T,n> linear_combination6(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4, double k5, Vector<T,n> v5);
template <typename T, int n> inline Vector<T,n> linear_combination7(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4, double k5, Vector<T,n> v5, double k6, Vector<T,n> v6);
template <typename T, int n> inline Vector<T,n> linear_combination8(double k0, Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4, double k5, Vector<T,n> v5, double k6, Vector<T,n> v6);

// VECTOR =================================================================================================================================

template <typename X, int N> std::ostream & operator<<(std::ostream &os, const Vector<X, N> &op);


template <typename X, int N> std::ostream & operator<<(std::ostream &os, const Vector<X, N> &op)
{
//	os << " (";
	if (op.size() > 0)
	{
		for (int i = 0; i < op.size() - 1; i++)
			os << op[i] << "\t";
		os << op[op.size() - 1];
	}
	else
		os << "empty vector" << std::endl;
//	os << ") ";
	return os;
}


template <class X, int N> Vector<X, N>::Vector(X v0, ...)
	:n(N)
// to make sure this constructor functions properly, pass doubles via (double) X or use decimal point
{
	va_list components;
	X v(v0);
	va_start(components, v0);
	data[0] = v;
	for (int i = 1; i < N ; i++)
	{
		v = va_arg(components, X);
		data[i] = v;
		double v_norm = scalar_norm(v);
		if (v_norm > 1.e39 || (v_norm > 0. && v_norm < 1.e-39))
		{
			std::cout << "Vector<X, N>::Vector(X v0, ...): argument number " << i << " is " << v 
				<< "\n\tmaybe this is not intended and constructor is called with incorrect number of arguments?\n";
			wait_and_exit();
		}
	}
	va_end(components);
}

template <class X, int N> Vector<X, N>::Vector(const Vector<X, N> &op)
// copy constructor
	:n(N)
{
	for (int i = 0; i < N; i++)
		data[i] = op.data[i];
}

template <class X, int N> inline X & Vector<X, N>::operator[](int i)
{
#ifdef RANGE_CHECK	
	if ( i >= 0 && i < n)
		return data[i];
	else
	{
		std::cout << "& Vector<X>::operator[]: index " << i << " out of range(0," << n << ")\n";
		wait_and_exit();

		return data[0];	// bogus string to get rid of compiler warning "not all control paths return a value"
	}
#elif
	return data[i];
#endif
}

template <class X, int N> inline const X & Vector<X, N>::operator[](int i) const
{
#ifdef RANGE_CHECK	
	if ( i >= 0 && i < n)
		return data[i];
	else
	{
		std::cout << "& Vector<X>::operator[]: index " << i << " out of range(0," << n << ")\n";
		wait_and_exit();

		return data[0];	// bogus string to get rid of compiler warning "not all control paths return a value"
	}
#elif
	return data[i];
#endif
}

template <class X, int N> Vector<X, N> & Vector<X, N>::operator=(const Vector<X, N> & op)
{
	if (&op != this)
	{
		for (int i = 0; i < n; i++)
			data[i] = op.data[i];
	}
	return *this;
}

template <class X, int N> Vector<X, N> & Vector<X,N>::operator=(double k)
{
	for(int i = 0; i < n; i++)
		data[i] = k;
	return *this;
}

template <class X, int N> Vector<X, N> & Vector<X,N>::operator=(int k)
{
	for(int i = 0; i < n; i++)
		data[i] = k;
	return *this;
}

template <typename X, int N> Vector<X, N> & Vector<X, N>::operator+=(const Vector<X, N> &op)
{
	for (int i = 0; i < n; i++)
		data[i] += op.data[i];
	return *this;
}

template <typename X, int N> Vector<X, N> & Vector<X, N>::operator-=(const Vector<X, N> &op)
{
	for (int i = 0; i < n; i++)
		data[i] -= op.data[i];
	return *this;
}

template <typename X, int N> Vector<X, N> & Vector<X, N>::operator*=(const Vector<X, N> &op)
{
	for (int i = 0; i < n; i++)
		data[i] *= op.data[i];
	return *this;
}

template <typename X, int N> Vector<X, N> & Vector<X, N>::operator/=(const Vector<X, N> &op)
{
	for (int i = 0; i < n; i++)
		data[i] /= op.data[i];
	return *this;
}

template <typename X, int N> Vector<X, N> & Vector<X, N>::operator*=(double k)
{
	for (int i = 0; i < n; i++)
		data[i] *= k;
	return *this;
}

template <typename X, int N> Vector<X, N> & Vector<X, N>::operator/=(double k)
{
	for (int i = 0; i < n; i++)
		data[i] /= k;
	return *this;
}


template <class X, int N> const Vector<X, N> Vector<X, N>::operator+(const Vector<X, N> &op) const
{
	Vector <X, N> temp(*this);
	return temp += op;
}

template <class X, int N> const Vector<X, N> Vector<X, N>::operator-(const Vector<X, N> &op) const
{
	Vector <X, N> temp(*this);
	return temp -= op;
}

template <class X, int N> const Vector<X, N> Vector<X, N>::operator*(const Vector<X, N> &op) const
{
	Vector <X, N> temp(*this);
	return temp *= op;
}

template <class X, int N> const Vector<X, N> Vector<X, N>::operator/(const Vector<X, N> &op) const
{
	Vector <X, N> temp(*this);
	return temp /= op;
}

template <class Y, int y_size> const Vector<Y, y_size> operator*(double k, const Vector<Y, y_size> &op)
{
	Vector <Y, y_size> temp(op);
	return temp *= k;
}

template <class Y, int y_size> const Vector<Y, y_size> operator*(const Vector<Y, y_size> &op, double k)
{
	Vector <Y, y_size> temp(op);
	return temp *= k;
}

template <class Y, int y_size> Vector<double, y_size> vector_norm(const Vector<Y, y_size> &op)
{
	Vector <double, y_size> temp;
	for(int i = 0; i < y_size; i++)
		temp[i] = scalar_norm(op[i]);
	return temp;
}


template <class Y, int y_size> double scalar_norm(const Vector<Y, y_size> &op)
{
	double res = 0.;
	Vector<double, y_size> temp = vector_norm(op);

	for(int i = 0; i < y_size; i++)
		res += (temp.data[i] * temp.data[i]);
	return sqrt(res);	// / ((double) y_size);
}

template <class Y, int y_size> double scalar_norm_inf(const Vector<Y, y_size> &op)
{
	double res = 0.;
	Vector<double, y_size> temp = vector_norm(op);

	for(int i = 0; i < y_size; i++)
		if (fabs(temp[i]) > res) res = fabs(temp[i]);
	return res;
}

template <class Y, int y_size> Vector<double, y_size> max_abs(const Vector<Y, y_size> &op1, const Vector<Y, y_size> &op2)
{
	Vector<double, y_size> tmp;
	for(int i = 0; i < y_size; i++)
		tmp[i] = std::max(scalar_norm(op1[i]), scalar_norm(op2[i]));
	return tmp;
}

template <class Y, int y_size> Y sum(const Vector<Y, y_size> & op)
{
	Y tmp = 0;
	for(int i = 0; i < y_size; i++)
		tmp += op.data[i];
	return tmp;
}


// LINEAR COMBINATION =====================================================================================================================

template <typename T> T linear_combination2(T v0, double k1, T v1)
{
	T V;
	for(int i = 0; i < V.size(); i++)
		V[i] = linear_combination2(v0[i], k1, v1[i]);
	return V;
}

template <typename T, int n> Vector<T,n> linear_combination3(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2)
{
	Vector<T,n> V;
	for(int i = 0; i < V.size(); i++)
		V[i] = linear_combination3(v0[i], k1, v1[i], k2, v2[i]);
	return V;
}

template <typename T, int n> Vector<T,n> linear_combination4(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3)
{
	Vector<T,n> V;
	for(int i = 0; i < V.size(); i++)
		V[i] = linear_combination4(v0[i], k1, v1[i], k2, v2[i], k3, v3[i]);
	return V;
}

template <typename T, int n> Vector<T,n> linear_combination5(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4)
{
	Vector<T,n> V;
	for(int i = 0; i < V.size(); i++)
		V[i] = linear_combination5(v0[i], k1, v1[i], k2, v2[i], k3, v3[i], k4, v4[i]);
	return V;
}

template <typename T, int n> Vector<T,n> linear_combination6(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4, double k5, Vector<T,n> v5)
{
	Vector<T,n> V;
	for(int i = 0; i < V.size(); i++)
		V[i] = linear_combination6(v0[i], k1, v1[i], k2, v2[i], k3, v3[i], k4, v4[i], k5, v5[i]);
	return V;
}

template <typename T, int n> Vector<T,n> linear_combination7(Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4, double k5, Vector<T,n> v5, double k6, Vector<T,n> v6)
{
	Vector<T,n> V;
	for(int i = 0; i < V.size(); i++)
		V[i] = linear_combination7(v0[i], k1, v1[i], k2, v2[i], k3, v3[i], k4, v4[i], k5, v5[i], k6, v6[i]);
	return V;
}

template <typename T, int n> Vector<T,n> linear_combination7(double k0, Vector<T,n> v0, double k1, Vector<T,n> v1, double k2, Vector<T,n> v2, double k3, Vector<T,n> v3, double k4, Vector<T,n> v4, double k5, Vector<T,n> v5, double k6, Vector<T,n> v6)
{
	Vector<T,n> V;
	for(int i = 0; i < V.size(); i++)
		V[i] = linear_combination7(k0, v0[i], k1, v1[i], k2, v2[i], k3, v3[i], k4, v4[i], k5, v5[i], k6, v6[i]);
	return V;
}

// comparison operators ===================================================================================================================

template<typename X, int N>  bool Vector<X, N>::operator==(const Vector<X, N> &op)
{
    for(int i = 0; i < N; i++)
        if(data[i] != op.data[i])
            return false;
    return true;
}

template<typename X, int N>  bool Vector<X, N>::operator==(double op)
{
    for(int i = 0; i < N; i++)
        if(data[i] != op)
            return false;
    return true;
}


template<typename X, int N>  bool Vector<X, N>::operator!=(const Vector<X, N> &op)
{
    for(int i = 0; i < N; i++)
        if(data[i] != op.data[i])
            return true;
    return false;
}

template<typename X, int N>  bool Vector<X, N>::operator!=(double op)
{
    for(int i = 0; i < N; i++)
        if(data[i] != op)
            return true;
    return false;
}


template<typename X, int N>  bool Vector<X, N>::operator<(const Vector<X, N> &op)
{
    for(int i = 0; i < N; i++)
        if(data[i] >= op.data[i])
            return false;
    return true;
}

template<typename X, int N>  bool Vector<X, N>::operator<(double op)
{
    for(int i = 0; i < N; i++)
        if(data[i] >= op)
            return false;
    return true;
}


template<typename X, int N>  bool Vector<X, N>::operator>(const Vector<X, N> &op)
{
    for(int i = 0; i < N; i++)
        if(data[i] <= op.data[i])
            return false;
    return true;
}

template<typename X, int N>  bool Vector<X, N>::operator>(double op)
{
    for(int i = 0; i < N; i++)
        if(data[i] <= op)
            return false;
    return true;
}

template<typename X, int N>  bool Vector<X, N>::operator<=(const Vector<X, N> &op)
{
    for(int i = 0; i < N; i++)
        if(data[i] > op.data[i])
            return false;
    return true;
}

template<typename X, int N>  bool Vector<X, N>::operator<=(double op)
{
    for(int i = 0; i < N; i++)
        if(data[i] > op)
            return false;
    return true;
}

template<typename X, int N>  bool Vector<X, N>::operator>=(const Vector<X, N> &op)
{
    for(int i = 0; i < N; i++)
        if(data[i] < op.data[i])
            return false;
    return true;
}

template<typename X, int N>  bool Vector<X, N>::operator>=(double op)
{
    for(int i = 0; i < N; i++)
        if(data[i] < op)
            return false;
    return true;
}


template<typename X, int N>  void Vector<X, N>::nullify_negative()
{
    for(int i = 0; i < N; i++)
        if(data[i] < 0.)
            data[i] = 0.;

}


};	// namespace fasim

#endif