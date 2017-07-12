#ifndef DECAYRATES_H_
#define DECAYRATES_H_

#include <cmath>
#include <cstdio>
#include <iostream>

double Gamma_total(double, double, double, double, double, double);
double Gamma_ZP_EE(double, double, double, double, double, double);
double Gamma_EE(double, double, double, double);
double Gamma_EPI(double, double, double, double);
double Gamma_MUPI(double, double, double, double);
double Gamma_NUPI0(double, double, double, double);
double Gamma_NUGAMMA(double,double , double , double );
double Gamma_NUMUMU(double, double, double, double);
double Gamma_NUNUNU(double, double, double, double);

//These declarations are for integral approximations given at the bottom of this file.
double I1(double, double);
double I2(double, double);
double I1_3arg_e(double);
double I1_3arg_mue(double);
double I1_3arg_mu(double);
double I2_3arg_e(double);
double I2_3arg_mu(double);

double branch_plotter(double, double);

// Verrified by mark
#define GF  1.16636e-5 /* Fermi's constant in GeV^-2 */
#define fPi  0.1307 /* pion decay constant in GeV */
#define Vud  0.97425 /* CKM element (for charged pion production) */
#define pi_mass   0.13957
#define	e_mass    0.000511
#define	Z_mass    91.1876 //GeV
#define pi0_mass   0.135
#define xW  0.231 /* sin^2\theta_W -- the Weinberg angle */
#define mu_mass   0.105

#endif
