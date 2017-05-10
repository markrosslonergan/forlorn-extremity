#include "minInstance.h"


minInstance::minInstance(double normn, double normb, std::vector<double> sen, std::vector<double> scn, std::vector<double> senb, std::vector<double> scnb) :
       		sig_E_nu(sen),
		sig_E_nubar(senb),
		sig_C_nu(scn),
		sig_C_nubar(scnb),
		norm_nu(normn),
		norm_nubar(normb)

		{	

		X_E_nu=0;
		X_C_nu=0;
		X_E_nubar=0;
		X_C_nubar=0;


	f_minimizer_mode ="GSLMultiMin"; //"GSLSimAn"
	f_minimizer_algo= "BFGS2";


		 obs_E_nu = {204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3};
		 bkg_E_nu = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9};
		 obs_E_nubar ={93,130,85,68,45,40,14,18,11,14,12,12,12,2,4,7,3,2,4};
		 bkg_E_nubar ={ 74.2,107.5,73.5,49.3,36.7,27.8,25.1,20.4,18.6,13.9,13.5,9.8,8.9,7.8,5.3,5,3.9,3.8,1.9};
			
		 obs_C_nu = {22,34,43,41,60,87,90,139,237,429};
		 bkg_C_nu = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390};
		 obs_C_nubar = {10,13,16,20,24,36,41,70,94,263};
		 bkg_C_nubar = {9.2,11.2,13.5,16,18.7,24.2,36,52.1,94.9,237.1};


		 Nbkg_E_nu= std::accumulate(bkg_E_nu.begin(), bkg_E_nu.end(), 0);
		 Nobs_E_nu= std::accumulate(obs_E_nu.begin(), obs_E_nu.end(), 0);
		 Nbkg_E_nubar= std::accumulate(bkg_E_nubar.begin(), bkg_E_nubar.end(), 0);
		 Nobs_E_nubar= std::accumulate(obs_E_nubar.begin(), obs_E_nubar.end(), 0);
	
		 Nbkg_C_nu= std::accumulate(bkg_C_nu.begin(), bkg_C_nu.end(), 0);
		 Nobs_C_nu= std::accumulate(obs_C_nu.begin(), obs_C_nu.end(), 0);
		 Nbkg_C_nubar= std::accumulate(bkg_C_nubar.begin(), bkg_C_nubar.end(), 0);
		 Nobs_C_nubar= std::accumulate(obs_C_nubar.begin(), obs_C_nubar.end(), 0);
	

		std::cout<<"Total Events Expected in Neutrino Mode, Angle: "<<Nbkg_C_nu<<" observed: "<<Nobs_C_nu<<" Excess: "<<-Nbkg_C_nu+Nobs_C_nu<<std::endl;
		std::cout<<"Total Events Expected in Neutrino Mode, Energy: "<<Nbkg_E_nu<<" observed: "<<Nobs_E_nu<<" Excess: "<<-Nbkg_E_nu+Nobs_E_nu<<std::endl;
		std::cout<<"Total Events Expected in AntiNeutrino Mode, Angle: "<<Nbkg_C_nubar<<" observed: "<<Nobs_C_nubar<<" Excess: "<<-Nbkg_C_nubar+Nobs_C_nubar<<std::endl;
		std::cout<<"Total Events Expected in AntiNeutrino Mode, Energy: "<<Nbkg_E_nubar<<" observed: "<<Nobs_E_nubar<<" Excess: "<<-Nbkg_E_nubar+Nobs_E_nubar<<std::endl;


}//end minInstance constructor;


int minInstance::clear_all(){

return 1;
}

int minInstance::reset_minim(){
	//min->Clear();
	//min->~Minimizer();
return 1;

}


int minInstance::init_minim(){

//	min = new ROOT::Math::GSLMinimizer(ROOT::Math::kConjugateFR);	
//	min = new ROOT::Math::GSLMinimizer(ROOT::Math::kVectorBFGS2);	
  // 	min->SetMaxIterations(250);     // for GSL
//	min->SetTolerance(0.001); 	//times 4 for normal
//	min->SetPrintLevel(0);
//	min->SetPrecision(0.0001);	//times 4 for normal

return 1;
}



double minInstance::minim_calc_chi(const double * x){
	double ans = 99999;

	double v_chi = x[0];
	double v_up = x[1];
	double v_ud = x[2];
	double v_zeta_b_nu = x[3];
	double v_zeta_b_nubar = x[4];


	std::vector<double> vchi;
	vchi = this->calc_chi(v_chi, v_ud, v_up, v_zeta_b_nu, v_zeta_b_nubar );	

	double chi = vchi[0]+vchi[1]+vchi[2]+vchi[3];
	


	/*if(ans<0){std::cout<<"WARNING: chi^2 is less than 0 in minimizer: "<<ans<<std::endl;
			std::cout<<"mass: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
			std::cout<<"ue: "<<Iue[0]<<" "<<Iue[1]<<" "<<Iue[2]<<std::endl;
			std::cout<<"um: "<<Ium[0]<<" "<<Ium[1]<<" "<<Ium[2]<<std::endl;
			std::cout<<"phi: "<<Iphi[0]<<" "<<Iphi[1]<<" "<<Iphi[2]<<std::endl;
		exit(EXIT_FAILURE);
	}*/

	return chi;

}



double minInstance::minimize(){

	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(f_minimizer_mode, f_minimizer_algo);
	min->SetMaxIterations(250);     // for GSL
	min->SetTolerance(0.001); 	//times 4 for normal
	min->SetPrintLevel(0);
	min->SetPrecision(0.0001);	//times 4 for normal




        ROOT::Math::Functor f( this, &minInstance::minim_calc_chi,5); 
	TRandom3 *rangen    = new TRandom3(0);


	double variable[5] ={-4,-4,-4,0.01,0.01};


	double step[5] = {0.01,0.01,0.01, 0.005,0.005};
	double lower[5] = {-5,-5,-5,-4,4};
	double upper[5] = {0,0,0,4,4 };	

	
	std::string name[5] ={"chi\0","Up\0","Ud\0","zeta_nu\0","zeta_nubar\0"};
	int isfixed[15]={0,0,0,0,0};

		
 	min->SetFunction(f);

   for(int i=0;i<5;i++){
	if(isfixed[i]){
	//	std::cout<<"Setting Fixed Variable: "<<i<<" "<<name[i]<<" value: "<<variable[i]<<std::endl;
	   	min->SetFixedVariable(i,name[i],variable[i]);
	} else {
	//	std::cout<<"Setting Variable: "<<i<<" "<<name[i]<<" value: "<<variable[i]<<" lower: "<<lower[i]<<" upper: "<<upper[i]<<" step: "<<step[i]<<std::endl;
   		min->SetLimitedVariable(i,name[i],variable[i], step[i], lower[i],upper[i]);
	}

   }
  min->Minimize(); 
   
  const double *xs = min->X();

  double valAns = minim_calc_chi(xs);
  
  return valAns;

}



std::vector<double> minInstance::calc_chi(double inchi, double inUp, double inUd, double zeta_b_nu, double zeta_b_nubar){
	X_E_nu=0;
	X_C_nu=0;
	X_E_nubar=0;
	X_C_nubar=0;


        for(int b=0;b<sig_E_nu.size();b++)
                        {
                                double temp_sig = sig_E_nu[b];
                                double temp_bg = (1.0+zeta_b_nu)*bkg_E_nu[b];
                                double lambda = temp_sig + temp_bg;
                                X_E_nu+= 2.0*(lambda-obs_E_nu[b]) + 2.0*obs_E_nu[b]*log(obs_E_nu[b]/lambda);

                        }
        X_E_nu+= pow((zeta_b_nu/sigma_zeta_nu),2.0);
	
        for(int b=0;b<sig_C_nu.size();b++)
                        {
                                double temp_sig = sig_C_nu[b];
                                double temp_bg = (1.0+zeta_b_nu)*bkg_C_nu[b];
                                double lambda = temp_sig + temp_bg;
                                X_C_nu+= 2.0*(lambda-obs_C_nu[b]) + 2.0*obs_C_nu[b]*log(obs_C_nu[b]/lambda);

                        }
        X_C_nu+= pow((zeta_b_nu/sigma_zeta_nu),2.0);

        for(int b=0;b<sig_E_nubar.size();b++)
                        {
                                double temp_sig = sig_E_nubar[b];
                                double temp_bg = (1.0+zeta_b_nubar)*bkg_E_nubar[b];
                                double lambda = temp_sig + temp_bg;
                                X_E_nubar+= 2.0*(lambda-obs_E_nubar[b]) + 2.0*obs_E_nubar[b]*log(obs_E_nubar[b]/lambda);

                        }
        X_E_nubar+= pow((zeta_b_nubar/sigma_zeta_nubar),2.0);

        for(int b=0;b<sig_C_nubar.size();b++)
                        {
                                double temp_sig = sig_C_nubar[b];
                                double temp_bg = (1.0+zeta_b_nu)*bkg_C_nubar[b];
                                double lambda = temp_sig + temp_bg;
                                X_C_nubar+= 2.0*(lambda-obs_C_nubar[b]) + 2.0*obs_C_nubar[b]*log(obs_C_nubar[b]/lambda);

                        }
        X_C_nubar+= pow((zeta_b_nu/sigma_zeta_nu),2.0);

	std::vector<double> ans = {X_E_nu,X_C_nu,X_E_nubar, X_C_nubar};

return ans;
}


