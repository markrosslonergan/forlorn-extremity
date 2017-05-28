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


//	f_minimizer_mode ="GSLMultiMin"; //"GSLSimAn"
//	f_minimizer_algo= "BFGS2";
	f_minimizer_mode ="GSLSimAn";
	f_minimizer_algo= "";



		sigma_zeta_nu = 0.025;
		sigma_zeta_nubar = 0.025;

		 ebins=	{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.};
		 cbins={-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}; 


		 obs_E_nu = {204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3};
		 bkg_E_nu = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9};
		 obs_E_nubar ={93,130,85,68,45,40,14,18,11,14,12,12,12,2,4,7,3,2,4};
		 bkg_E_nubar ={ 74.2,107.5,73.5,49.3,36.7,27.8,25.1,20.4,18.6,13.9,13.5,9.8,8.9,7.8,5.3,5,3.9,3.8,1.9};

		 h_bkg_E_nu = new TH1D("Bkg_Evis_Nu","",19,&ebins[0]); 
		 h_bkg_E_nubar = new TH1D("Bkg_Evis_Nubar","",19,&ebins[0]); 

		 h_bkg_C_nu = new TH1D("Bkg_Cos_Nu","",10,&cbins[0]); 
		 h_bkg_C_nubar = new TH1D("Bkg_Cos_Nubar","",10,&cbins[0]); 


		 h_bf_E_nu = new TH1D("Bf_Evis_Nu","",19,&ebins[0]); 
		 h_bf_E_nubar = new TH1D("Bf_Evis_Nubar","",19,&ebins[0]); 

		 h_bf_C_nu = new TH1D("Bf_Cos_Nu","",10,&cbins[0]); 
		 h_bf_C_nubar = new TH1D("Bf_Cos_Nubar","",10,&cbins[0]); 


		 h_obs_E_nu = new TH1D("Obs_Evis_Nu","",19,&ebins[0]); 
		 h_obs_E_nubar = new TH1D("Obs_Evis_Nubar","",19,&ebins[0]); 

		 h_obs_C_nu = new TH1D("Obs_Cos_Nu","",10,&cbins[0]); 
		 h_obs_C_nubar = new TH1D("Obs_Cos_Nubar","",10,&cbins[0]); 

		h_excess_E_nu = new TH1D("Excess_Evis_Nu","",19,&ebins[0]); 
		 h_excess_E_nubar = new TH1D("Excess_Evis_Nubar","",19,&ebins[0]); 

		 h_excess_C_nu = new TH1D("Excess_Cos_Nu","",10,&cbins[0]); 
		 h_excess_C_nubar = new TH1D("Excess_Cos_Nubar","",10,&cbins[0]); 


		 obs_C_nu = {22,34,43,41,60,87,90,139,237,429};
		 bkg_C_nu = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390};
		 obs_C_nubar = {10,13,16,20,24,36,41,70,94,263};
		 bkg_C_nubar = {9.2,11.2,13.5,16,18.7,24.2,36,52.1,94.9,237.1};

		for(int i=1; i<=19; i++){
			h_bkg_E_nu->SetBinContent(i,bkg_E_nu[i-1]);
			h_bkg_E_nubar->SetBinContent(i,bkg_E_nubar[i-1]);
			h_obs_E_nu->SetBinContent(i,obs_E_nu[i-1]);
			h_obs_E_nubar->SetBinContent(i,obs_E_nubar[i-1]);
		
			h_excess_E_nu->SetBinContent(i,obs_E_nu[i-1]- bkg_E_nu[i-1]  );
			h_excess_E_nubar->SetBinContent(i,obs_E_nubar[i-1]-bkg_E_nubar[i-1]);

		}

		for(int i=1; i<=10; i++){
			h_bkg_C_nu->SetBinContent(i,bkg_C_nu[i-1]);
			h_bkg_C_nubar->SetBinContent(i,bkg_C_nubar[i-1]);
			h_obs_C_nu->SetBinContent(i,obs_C_nu[i-1]);
			h_obs_C_nubar->SetBinContent(i,obs_C_nubar[i-1]);
	
			h_excess_C_nu->SetBinContent(i,obs_C_nu[i-1]- bkg_C_nu[i-1]  );
			h_excess_C_nubar->SetBinContent(i,obs_C_nubar[i-1]-bkg_C_nubar[i-1]);


		}



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

int minInstance::setMass(double mz, double ms){

	mass_z = mz;
	mass_s=ms;

	return 0;
}


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
	vchi = this->calc_chi(v_chi, v_up, v_ud, v_zeta_b_nu, v_zeta_b_nubar );	

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
	min->SetMaxIterations(500000);     // for GSL
	min->SetTolerance(0.000001); 	//times 4 for normal
	min->SetPrecision(0.000001);	//times 4 for normal
	min->SetPrintLevel(2);



        ROOT::Math::Functor f( this, &minInstance::minim_calc_chi,5); 
	TRandom3 *rangen    = new TRandom3(0);


	double variable[5] ={-3,-4,-4,0.01,0.01};


	double step[5] = {0.005,0.005,0.005, 0.001,0.001};
	double lower[5] = {-10,-10,-10,-1,-1};
	double upper[5] = {0,0,0,1,1 };	

	
	std::string name[5] ={"chi\0","Up\0","Ud\0","zeta_nu\0","zeta_nubar\0"};
	int isfixed[5]={0,0,1,0,0};

		
 	min->SetFunction(f);

   for(int i=0;i<5;i++){
	if(isfixed[i]){
		std::cout<<"Setting Fixed Variable: "<<i<<" "<<name[i]<<" value: "<<variable[i]<<std::endl;
	   	min->SetFixedVariable(i,name[i],variable[i]);
	} else {
		std::cout<<"Setting Variable: "<<i<<" "<<name[i]<<" value: "<<variable[i]<<" lower: "<<lower[i]<<" upper: "<<upper[i]<<" step: "<<step[i]<<std::endl;
   		min->SetLimitedVariable(i,name[i],variable[i], step[i], lower[i],upper[i]);
	}

   }


  min->Minimize(); 
   
  const double *xs = min->X();

	 bf_chi = xs[0];
	 bf_up = xs[1];
	 bf_ud = xs[2];
         bf_zeta_b_nu=xs[3];
	 bf_zeta_b_nubar=xs[4];
	 double diam_miniboone = 10;
 
	 double sc =  pow(pow(10,bf_up),2)* pow(pow(10,bf_chi),4)*100*pow(1.973,-1)*pow(10.0,5.0)*pow(10.0,9.0)*diam_miniboone*bound_vector[0].myRate(mass_s, mass_z); 
;


  for(int i=1; i<=19; i++){
	h_bf_E_nu->SetBinContent(i, sig_E_nu[i-1]*norm_nu*sc );
	h_bf_E_nubar->SetBinContent(i,sig_E_nubar[i-1]*norm_nubar*sc );

  }  
  for(int i=1; i<=11; i++){
	h_bf_C_nu->SetBinContent(i,sig_C_nu[i-1]*norm_nu*sc );
	h_bf_C_nubar->SetBinContent(i,sig_C_nubar[i-1]*norm_nubar*sc );
  } 


  double valAns = minim_calc_chi(xs);
  
  return valAns;






}



std::vector<double> minInstance::calc_chi(double inchi, double inUp, double inUd, double zeta_b_nu, double zeta_b_nubar){
	X_E_nu=0;
	X_C_nu=0;
	X_E_nubar=0;
	X_C_nubar=0;

	//this is just guessork, have to actuall see what the scaling is
	double ch2 = pow(pow(10,inchi),2.0);
	double u2 = pow(pow(10,inUp),2.0);

	double diam_miniboone =10;
	double pdec = 100*pow(1.973,-1)*pow(10.0,5.0)*pow(10.0,9.0)*diam_miniboone*bound_vector[0].myRate(mass_s, mass_z)*ch2*ch2; 

	//std::cout<<"PDEC: "<<pdec<<std::endl;
	double UX=pdec*u2;


        for(int b=0;b<sig_E_nu.size();b++)
                        {
                                double temp_sig = norm_nu*UX*sig_E_nu[b];
                                double temp_bg = (1.0+zeta_b_nu)*bkg_E_nu[b];
                                double lambda = temp_sig + temp_bg;
                                X_E_nu+= 2.0*(lambda-obs_E_nu[b]) + 2.0*obs_E_nu[b]*log(obs_E_nu[b]/lambda);

                        }
        X_E_nu+= pow((zeta_b_nu/sigma_zeta_nu),2.0);
	
        for(int b=0;b<sig_C_nu.size();b++)
                        {
                                double temp_sig = norm_nu*UX*sig_C_nu[b];
                                double temp_bg = (1.0+zeta_b_nu)*bkg_C_nu[b];
                                double lambda = temp_sig + temp_bg;
                                X_C_nu+= 2.0*(lambda-obs_C_nu[b]) + 2.0*obs_C_nu[b]*log(obs_C_nu[b]/lambda);

                        }
        X_C_nu+= pow((zeta_b_nu/sigma_zeta_nu),2.0);

        for(int b=0;b<sig_E_nubar.size();b++)
                        {
                                double temp_sig = norm_nubar*UX*sig_E_nubar[b];
                                double temp_bg = (1.0+zeta_b_nubar)*bkg_E_nubar[b];
                                double lambda = temp_sig + temp_bg;
                                X_E_nubar+= 2.0*(lambda-obs_E_nubar[b]) + 2.0*obs_E_nubar[b]*log(obs_E_nubar[b]/lambda);

                        }
        X_E_nubar+= pow((zeta_b_nubar/sigma_zeta_nubar),2.0);

        for(int b=0;b<sig_C_nubar.size();b++)
                        {
                                double temp_sig = norm_nubar*UX*sig_C_nubar[b];
                                double temp_bg = (1.0+zeta_b_nu)*bkg_C_nubar[b];
                                double lambda = temp_sig + temp_bg;
                                X_C_nubar+= 2.0*(lambda-obs_C_nubar[b]) + 2.0*obs_C_nubar[b]*log(obs_C_nubar[b]/lambda);

                        }
        X_C_nubar+= pow((zeta_b_nu/sigma_zeta_nu),2.0);


	if(use_bounds){

		double penalty = 0 ;
		
			bool is_ok = bound_vector.at(0).ps191(mass_s,mass_z,ch2,u2);
			if(!is_ok){
				penalty=1e10;
			}
			is_ok = bound_vector.at(1).asIs(mass_s,u2);
			if(!is_ok){
				penalty=1e10;
			}
			is_ok = bound_vector.at(2).asIs(mass_z,ch2);
			if(!is_ok){
				penalty=1e10;
			}





		X_E_nu += penalty;
		X_C_nu += penalty ;
		X_E_nubar += penalty;
		X_C_nubar += penalty;


	}

	
	std::vector<double> ans = {X_E_nu,X_C_nu,X_E_nubar, X_C_nubar};

return ans;
}


