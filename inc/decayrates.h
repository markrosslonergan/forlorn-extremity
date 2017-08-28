#ifndef DECAYRATES_H_
#define DECAYRATES_H_

#include <cmath>
#include <cstdio>
#include <iostream>
#include "cuba.h"
#include "appellced.h"

struct decay_params
{
	double chi;
	double mZprime;
	double mS;
	double Ue4;
	double Um4;
	double Ut4;
	double numeval;
};

struct master_integrals 
{ 
	double I00; 
	double I01; 
	double I02; 
	double I10; 
	double I11; 
	double I12; 
	double I20; 
	double I21; 
	double I22; 
	double J01_0; 
	double J02_0; 
	double J11_0; 
	double J01_1; 
	double J02_1; 
	double J11_1; 
};
 
class decay_obj {

public:
	decay_obj(decay_params*);

	void update(decay_params*);
	void compute_master_integrals();
	void compute_gamma_total();
		
	decay_params params;
	master_integrals master_ints;

	double Gamma_nununu;
	double Gamma_nugamma;
	double Gamma_ee;
	double Gamma_epi;
	double Gamma_nupi0;
	double Gamma_numue;
	double Gamma_mupi;
	double Gamma_numumu;
	double Gamma_total;
};

double Gamma_EE(decay_params*, master_integrals*);
double Gamma_EPI(decay_params*);
double Gamma_MUPI(decay_params*);
double Gamma_NUPI0(decay_params*);
double Gamma_NUGAMMA(decay_params*);
double Gamma_NUMUE(decay_params*);
double Gamma_NUMUMU(decay_params*);
double Gamma_NUNUNU(decay_params*);

//Integrand and integrator for Gamma_ee and wrapper for the ee decay integration.
//Note that "matrix_element" includes all prefactors. Integrate first two arguments over 0,1 to get the full decay rate.
double Gamma_EE_integral(decay_params*);
double matrix_element_ee(double, double, decay_params*);
static int Integrand_ee(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);

//for testing and consistency purposes.
double Gamma_NUPI0_old(double, double, double, double, double, double); 
double Gamma_EE_old(double, double, double, double, double, double);
double Gamma_EE_oldest(double, double, double, double);
double Gamma_NUMUMU_old(double, double, double, double);
double Gamma_NUMUE_old(double, double, double, double);

//These declarations are for integral approximations given at the bottom of this file.
double I1(double, double);
double I2(double, double);
double I3(double, double);
double I1_C5(double,double);
double I2_C7(double);
double I1_C6(double);

//Old integral approximations
double I1_3arg_e(double);
double I1_3arg_mue(double);
double I1_3arg_mu(double);
double I2_3arg_e(double);
double I2_3arg_mu(double);

//Here are some nasty mathematica expressions
double sbeta_mathematica(double, double);
double cbeta_mathematica(double, double);
//This odd combination cancels analytically in a numerically nasty way.
inline double sinbeta_plus_cosbeta_sW_tanchi(double, double); 

double branch_plotter(decay_params*);

#define GF  1.16636e-5 /* Fermi's constant in GeV^-2 */
#define fPi  0.1307 /* pion decay constant in GeV */
#define Vud  0.97425 /* CKM element (for charged pion production) */
#define pi_mass   0.13957
#define	e_mass    0.000511
#define	Z_mass    91.1876 //GeV
#define	W_mass    80.385 //GeV
#define pi0_mass   0.135
#define mu_mass   0.105
#define xW  0.231 /* sin^2\theta_W -- the Weinberg angle */
#define v_vev 254.0 // Higgs vev.


//Tedious variables for the integration algorithm
#define NDIM 		2
#define NCOMP 		1
#define USERDATA 	NULL
#define NVEC 		1
#define EPSREL 		1e-3
#define EPSABS 		1e-12
#define VERBOSE 	0
#define LAST 		4
#define SEED 		0
#define MINEVAL 	5000
#define MAXEVAL 	500000000
#define NSTART 		5000
#define NINCREASE 	500
#define NBATCH 		1000
#define GRIDNO 		0
#define STATEFILE 	NULL
#define SPIN 		NULL


#endif
