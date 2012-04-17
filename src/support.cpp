#include "support.h"
#include <iostream>
#include <fstream>
#include <conio.h>

#include "constants.h"

void wait()
{
	std::cout << "\npress Enter to continue...";
	std::cin.get();
}

void wait_and_exit()
{
	wait();
	exit(EXIT_FAILURE);
}

void wait_and_exit(std::string error_message)
{
	std::cout << error_message << "\n";
	wait_and_exit();
}

//	<nm_to_Hz>
//	converts nanometer wavelength (say, 1540, without any exponent) into corresponding frequency
	double nm_to_Hz(double wavelength)
	{
		if (wavelength > 2500. || wavelength < 100.)
		{
			std::cout << "nm_to_Hz: wavelength = " << wavelength << " provided, maybe that's not intended?\n";
			wait_and_exit();
		}
		return fasim::C0 / (wavelength * 1.e-9);
	}

	double Hz_to_nm(double frequency)
	{
		if (frequency > 3.e15 || frequency < 1.2e14)
		{
			std::cout << "Hz_to_nm: frequency = " << frequency << " provided, maybe that's not intended?\n";
			wait_and_exit();
		}		
		return fasim::C0 / (frequency * 1.e-9);
	}


double round(double x)
{
    return x < 0.0 ? ceil(x - 0.5) : floor(x + 0.5);
}


// eps returns the minimal floating point number that can be added to a given double <x> so that the result is different from <x>
// implementation exploits the fact that IEEE double mantissa is exactly 53 bits (hidden bit included)

double eps(double x)
{
	if (x == 0. || x == -1.*0.)
	{
		return pow(2., -1022);
	}
	else
	{
		int e = (int) floor( log(abs(x)) / log(2.) );
		return pow(2., e - 52);
	};
}


// file import/export =============================================================================

void f_read_table(double *x, double *y, int *n, const std::string &filename)
{
	std::ifstream file_input(filename.c_str());
		
	if (!file_input) 
	{
		std::cout << "f_read_table: error reading " << filename.c_str() << "\n";
		wait_and_exit();
	}
		
	int input_length = 0;
	for( ; !file_input.eof(); input_length++)
	{
		file_input >> x[input_length] >> y[input_length];
		//std::cout << x[input_length] << "\t" << y[input_length] << "\n";
		//std::cin.get();
	}


	*n = input_length;

	file_input.close();
}


void f_write_table(double *x, double *y, int n, const std::string &filename)
{
	std::ofstream file_output(filename.c_str());
		
	if (!file_output) 
	{
		std::cout << "f_write_table: error reading " << filename.c_str() << "\n";
		wait_and_exit();
	}
		
	for(int i = 0; i < n; i++)
	{
		file_output << x[i] << "\t" << y[i];
		if (i != n - 1)
			file_output << "\n";
	}

	file_output.close();
}