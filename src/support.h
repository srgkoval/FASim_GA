#ifndef SUPPORT_H
#define SUPPORT_H

#include <string>

void wait_and_exit();
void wait_and_exit(std::string error_message);

//	<nm_to_Hz>
//	converts nanometer wavelength (say, 1540, without any exponent) into corresponding frequency
double nm_to_Hz(double wavelength);
double Hz_to_nm(double frequency);

double round(double x);
double eps(double x);

void f_read_table(double *x, double *y, int *n, const std::string &filename);
void f_write_table(double *x, double *y, int n, const std::string &filename);

#endif