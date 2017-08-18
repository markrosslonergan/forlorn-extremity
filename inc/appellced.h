#ifndef APPELLCED_H_
#define APPELLCED_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <complex>
#include <gsl/gsl_sf_hyperg.h>

#include "hyp_2F1.h"

//The Appell F1 hypergeometric function, defined for all real arguments by analytic cont.
double AppellF1(double, double, double, double, double, double, double*); 

//Wrapper for the MS hyp_2F1 function. 
double G2F1(double, double, double, double);

//Old implementation of the analytically continued hyp_2F1 function. Redundant as it has numerical issues. 
double G2F1_Forrey(double, double, double, double);

//Old implementation of 2F1 implemented by its series. Used to replace GSL implementation, which is odd.
double G2F1_sum(double, double, double, double);

//Pochhammer symbol.
double poch(double, int);

//Two currently redundant functions. Might be used later.

//Not really an Appell function: the G2 function which is part of the analytic
//continuation of the Appell F1.
double AppellG2(double, double, double, double, double, double);

//Appell F2 hypergeometric defined by sum.
double AppellF2(double, double, double, double, double, double, double);

#endif
