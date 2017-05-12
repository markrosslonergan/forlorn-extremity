#ifndef MIN_H_
#define MIN_H_


#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TError.h"

#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

#include "minInstance.h"


class minInstance {
	
	public:
	double norm_nu;
	double norm_nubar;

	double sigma_zeta_nu;
	double sigma_zeta_nubar;

	bool use_bounds;

	std::string f_minimizer_mode;
	std::string f_minimizer_algo;
	
	std::vector<double> ebins;
	std::vector<double> cbins;

	std::vector<double > sig_E_nu ;
	std::vector<double > sig_C_nu  ;
	std::vector<double > sig_E_nubar ;
	std::vector<double > sig_C_nubar  ;

	std::vector<double > bkg_E_nu;
	std::vector<double > bkg_C_nu;
	std::vector<double > bkg_E_nubar;
	std::vector<double > bkg_C_nubar;

	std::vector<double > obs_E_nu;
	std::vector<double > obs_C_nu;
	std::vector<double > obs_E_nubar;
	std::vector<double > obs_C_nubar;

	TH1D * h_obs_E_nu;
	TH1D * h_obs_C_nu;
	TH1D * h_obs_E_nubar;
	TH1D * h_obs_C_nubar;

	TH1D * h_excess_E_nu;
	TH1D * h_excess_C_nu;
	TH1D * h_excess_E_nubar;
	TH1D * h_excess_C_nubar;



	TH1D * h_bf_E_nu;
	TH1D * h_bf_C_nu;
	TH1D * h_bf_E_nubar;
	TH1D * h_bf_C_nubar;

	TH1D * h_bkg_E_nu;
	TH1D * h_bkg_C_nu;
	TH1D * h_bkg_E_nubar;
	TH1D * h_bkg_C_nubar;

	double Nbkg_E_nu;
	double Nobs_E_nu;
	double Nbkg_E_nubar;
	double Nobs_E_nubar;

	double Nbkg_C_nu;
	double Nobs_C_nu;
	double Nbkg_C_nubar;
	double Nobs_C_nubar;


	double X_E_nu;
	double X_C_nu;
	double X_E_nubar;
	double X_C_nubar;


	double bf_chi;
	double bf_up;
	double bf_ud;
	double bf_zeta_b_nu;
	double bf_zeta_b_nubar;


	double Current_Chi;
	//ROOT::Math::Minimizer* min ; 
	
   /* 
   min.SetMaxFunctionCalls(100000);
   min.SetMaxIterations(10000);
   min.SetTolerance(0.01);
 
   ROOT::Math::Functor f(&RosenBrock,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = { -1.,1.2};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   //    min.SetVariable(0,"x",variable[0], step[0]);
   //       min.SetVariable(1,"y",variable[1], step[1]);
   //        
   //           min.Minimize(); 
   //            
   //               const double *xs = min.X();
   //                  cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
   //                          << RosenBrock(xs) << endl;
   //                           
   //                              return 0;
*/

	minInstance(double,double,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>);


	std::vector<double> calc_chi(double,double,double,double,double);// returns a vector of 4 chis Enu Cnu Enubar Cnubar

	double minim_calc_chi(const double * xx);
	
	int init_minim();
	double minimize();
	int reset_minim();

	int clear_all();
};
#endif
