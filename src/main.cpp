#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <conio.h>
#include "boost/lexical_cast.hpp"

#include "ode.h"
#include "model_specification.h"
#include "test.h"
#include "test_ga.h"

#include "thulium.h"
#include "erbium_upconversion.h"
#include "erbium.h"

void main()
{
	clock_t start, finish;
	start = clock();
	std::cout.setf( std::ios::showpoint );
	std::cout.precision(12);

	//run_erbium_ESA();
	//run_erbium();
	//test_dN3();
	run_thulium();
	//run_erbium_upconversion();
	//run_erbium_upconversion_optim();
	
    //run_KUR();

	finish = clock();
	std::cout << "\ntimer:\t" << (double (finish - start) / CLOCKS_PER_SEC) << std::endl;
	std::cin.get();
}