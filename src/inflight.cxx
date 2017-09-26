/*
 *_________ _        _______  _       _________ _______          _________
 *\__   __/( (    /|(  ____ \( \      \__   __/(  ____ \|\     /|\__   __/
 *   ) (   |  \  ( || (    \/| (         ) (   | (    \/| )   ( |   ) (   
 *   | |   |   \ | || (__    | |         | |   | |      | (___) |   | |   
 *   | |   | (\ \) ||  __)   | |         | |   | | ____ |  ___  |   | |   
 *   | |   | | \   || (      | |         | |   | | \_  )| (   ) |   | |   
 *___) (___| )  \  || )      | (____/\___) (___| (___) || )   ( |   | |   
 *\_______/|/    )_)|/       (_______/\_______/(_______)|/     \|   )_(   
 *
 *		Inflight, Event generator for sterile decays at SBL facilities
 *
 *		If you have any questions, queries or comments please contact the authors;
 *			 mark.ross-lonergan@durham.ac.uk 
 *			 or 
 *			 peter.ballett@durham.ac.uk
 *
 *		The authors make no guarrentee of the behaviour, stability or bug-free-ness of this code.
 *		Use is at own risk.
 *
 * Copyright (2016) 
 * Mark Ross-Lonergan mark.ross-lonergan@durham.ac.uk 
 * Peter Ballett peter.ballett@durham.ac.uk
 */
#include <iostream>
#include <numeric>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define no_argument 0
#define required_argument 1
#define optional_argument 2

#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1.h"
#include "TCanvas.h"
#include "THStack.h"


#include "fourmomentum.h" // Defines a class for fourmomenta; makes getting directions and 3-momenta (and
// all that) easier.

#include "sterile_flux.h" // Defines getEvents, and describes a class whose objects denote a single 
// incoming sterile (mass, fourmomentum).

#include "channel.h"	  // Includes the mother class for two body final state
// decays, and the derived classes for the two channels we've studied so far.

#include "detector.h"	  // This includes the detector-specific cut functions.


#include "zprimecrosssection.h"


#include "minInstance.h"
#include "bounds.h"


#define MPION  0.13957
#define MPI0   0.13497
#define MKAON  0.49367
#define MMU    0.10566
#define	ME     0.00051




/* ########## Main function ############### */
int main(int argc, char* argv[])
{

	const gsl_rng_type * T; // Standard invocation prayer to the RNG
	gsl_rng * r;
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	int channel_flag = CHAN_ELECPI; 	// This one too.
	int detector_flag = DET_NOCUTS;	 	// This one too.
	int matrix_flag = 0; 	 		// This one too.
	int eff_flag = 0; 			// This one too.
	bool hepevt = false;
	std::string flux_name = "flux.dat";

	int NUMEVENTS = 100;

	int c;
	opterr = 0;

	double mS=0.1;
	double mZ=1.0;

	bool genMode = true;
	bool anaMode = false;
	bool statMode =false;
	bool boundMode =false;


	int index; 
	int iarg = 0;
	opterr=1;
	const struct option longopts[] = 
	{
		{"mass",	 	required_argument, 	0, 'm'},
		{"zprime",		required_argument,	0, 'z'},
		{"num",			required_argument,	0, 'n'},
		{"mode",		required_argument,	0, 'M'},
		{"help",		no_argument,		0, 'h'},
		{0,			no_argument, 		0,  0}
	};
	const char * tin;

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "M:hm:n:z:", longopts, &index);

		switch(iarg)
		{
			case 'z':
				mZ =strtof(optarg,NULL);
				break;
			case 'm':
				mS =strtof(optarg,NULL);
				break;
			case 'n':
				NUMEVENTS = strtod(optarg,NULL);
				break;
			case 'M':
				tin = optarg;//way to read in null terminated c strings and compare to known list
				if(strcmp(tin, "gen")==0){ genMode = true; }
				if(strcmp(tin, "ana")==0) { genMode = false; anaMode=true;}
				if(strcmp(tin, "stat")==0 ){ genMode = false; statMode=true;}
				if(strcmp(tin, "bound")==0 ){ genMode = false; boundMode=true;}
				break;
			case '?':
				std::cout<<"Abandon hope all ye who enter this value. "<<std::endl;
			case 'h':
				std::cout<<"\nInflight, an event generator for heavy sterile decays at SBL facilities"<<std::endl;
				std::cout<<"If you have any questions please contact the authors mark.ross-lonergan@durham.ac.uk or peter.ballett@durham.ac.uk"<<std::endl;
				std::cout<<"The authors make no guarrentee of this code"<<std::endl;
				std::cout<<"******************************************"<<std::endl;
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-m, --mass\t\t\tSets the parent sterile mass. [default = 0.1000]"<<std::endl;
				std::cout<<"\t-z, --zprime\t\t\tSets the scattering Zprime mass. [default = 1.000]"<<std::endl;
				std::cout<<"\t-n, --number\t\t\tHow many events will we generate? [Default 100]"<<std::endl;
				std::cout<<"\t-M, --mode\t\t\tMode is either ``gen'', ``ana'' or ``stat''  [Default gen]"<<std::endl;
				std::cout<<"******************************************"<<std::endl;

				return 0;
		}

	}
	int EBINS = 19;
	int COSBINS = 10;

	double EMIN = 0.1;
	double EMAX = 2.0;
	double COSMIN = -1.0;
	double COSMAX = 1.0;

	double EBINWIDTH=fabs(EMIN-EMAX)/((double)EBINS);
	double COSBINWIDTH=fabs(COSMIN-COSMAX)/((double)COSBINS);


	std::string core ="/home/mark/work/durham/zprime_lee_signal/forlorn-extremity/";

	std::string filename = "zmass_" + std::to_string(mZ) + "_smass_" + std::to_string(mS)+"_sim.root";

	if(false){//testing

		bound bound_ps191("../../data/bounds/PS191_UM4_EE_BOTH.dat",0.01,100);


		bound_ps191.setTypicalEnergy(5.0);

		bound_ps191.ps191(0.1, 91, 0.01, 1);


		return 0;
	}
	/************************************************************************
	 *
	 *
	 *		GENMODE genMode GeNMode
	 *
	 *
	 *
	 * **********************************************************************/



	if(genMode){

		std::cout<<"Running genMode with mZ: "<<mZ<<" GeV, and mS: "<<mS<<" GeV."<<std::endl;
		twoIP_channel *CHAN;
		std::vector<double> model_params; //This should include any theoretical parameters which the model needs to know.

		model_params.push_back(mZ); //mediator mass
		model_params.push_back((double) CHAN_ELECPOSI); // the pion. 
		model_params.push_back((double) CHAN_ELECPOSI); // the pion. 
		CHAN = new threebody(r,model_params); 

		detector * DETECTOR;

		DETECTOR = new nocuts(); 	// pseudo-detector that just allows every event.


		static OBSERVABLES Obs; //This struct is contained in "decay.h"; it specifically gives variables for a two body event (e+,e-)
		Obs.chan_identifier= CHAN->chan_identifier; // needs this so it knows what resolutions to use in smearing


		TFile *f1 = new TFile(filename.c_str(),"RECREATE");
		TTree *t1 = new TTree("nuTree","nuTree");
		TTree *t3 = new TTree("nubarTree","nubarTree");
		TTree *t2 = new TTree("summeryTree","summeryTree");

		int split = 0;
		int bsize = 64000;

		double sEn;
		double sPhi;
		double sCos;
		double vEn;
		double e1En;
		double e2En;
		double e1Cos;
		double e2Cos;
		double csTmandel;
		double dAngSep;
		double dFsAngSep;
		double dEnSum;
		double dThSum;
		double sPx;
		double sPz;
		double sPy;
		double e1Px;
		double e1Pz;
		double e1Py;
		double e2Px;
		double e2Pz;
		double e2Py;


		double sEn_B;
		double sPhi_B;
		double sCos_B;
		double vEn_B;
		double e1En_B;
		double e2En_B;
		double e1Cos_B;
		double e2Cos_B;
		double csTmandel_B;
		double dAngSep_B;
		double dFsAngSep_B;
		double dEnSum_B;
		double dThSum_B;
		double sPx_B;
		double sPz_B;
		double sPy_B;
		double e1Px_B;
		double e1Pz_B;
		double e1Py_B;
		double e2Px_B;
		double e2Pz_B;
		double e2Py_B;

		double sM;
		double zM;
		double csTotalNum;
		double csTotalNum_B;

		t2->Branch("sM", &sM);	
		t2->Branch("zM", &zM);	
		t2->Branch("csTotalNumNu", &csTotalNum);	
		t2->Branch("csTotalNumNuBar", &csTotalNum_B);	

		t1->Branch("csTmandel", &csTmandel);	
		t1->Branch("sEn", &sEn);	
		t1->Branch("sPhi",&sPhi);
		t1->Branch("sCos", &sCos);	
		t1->Branch("vEn", &vEn);	

		t1->Branch("e1En", &e1En);	
		t1->Branch("e1Cos", &e1Cos);	
		t1->Branch("e2En", &e2En);	
		t1->Branch("e2Cos", &e2Cos);	

		t1->Branch("dAngSep",&dAngSep);
		t1->Branch("dFsAngSep",&dFsAngSep);
		t1->Branch("dThSum",&dThSum);
		t1->Branch("dEnSum",&dEnSum);

		t1->Branch("sPx",&sPx);
		t1->Branch("sPy",&sPy);
		t1->Branch("sPz",&sPz);
		t1->Branch("e1Px",&e1Px);
		t1->Branch("e1Py",&e1Py);
		t1->Branch("e1Pz",&e1Pz);
		t1->Branch("e2Px",&e2Px);
		t1->Branch("e2Py",&e2Py);
		t1->Branch("e2Pz",&e2Pz);

		// branches
		t3->Branch("csTmandel", &csTmandel_B);	
		t3->Branch("sEn", &sEn_B);	
		t3->Branch("sPhi",&sPhi_B);
		t3->Branch("sCos", &sCos_B);	
		t3->Branch("vEn", &vEn_B);	

		t3->Branch("e1En", &e1En_B);	
		t3->Branch("e1Cos", &e1Cos_B);	
		t3->Branch("e2En", &e2En_B);	
		t3->Branch("e2Cos", &e2Cos_B);	

		t3->Branch("dAngSep",&dAngSep_B);
		t3->Branch("dFsAngSep",&dFsAngSep_B);
		t3->Branch("dThSum",&dThSum_B);
		t3->Branch("dEnSum",&dEnSum_B);

		t3->Branch("sPx",&sPx_B);
		t3->Branch("sPy",&sPy_B);
		t3->Branch("sPz",&sPz_B);
		t3->Branch("e1Px",&e1Px_B);
		t3->Branch("e1Py",&e1Py_B);
		t3->Branch("e1Pz",&e1Pz_B);
		t3->Branch("e2Px",&e2Px_B);
		t3->Branch("e2Py",&e2Py_B);
		t3->Branch("e2Pz",&e2Pz_B);


		//Bit from zprime
		double norm =totalevents2(1,1,mZ,mS, POSNU);
		double normBar =totalevents2(1,1,mZ,mS, NEGNU);

		//Maximum Total CS for normalisation sampling. 
		std::vector<double > CSmax = MaxDiffCs(mZ,mS,1,1);


		csTotalNum=norm;
		csTotalNum_B=normBar;
		zM = mZ;
		sM = mS;

		for(int i=0; i<NUMEVENTS; i++){

			if(i%500==0) std::cout<<"#: "<<i<<std::endl;
			double  phiS = 2.0*M_PI*gsl_rng_uniform(r);

			std::vector<double> event = MCkin(mZ,mS,r,CSmax[1],POSNU); //Get Event
			std::vector<double> eventBar = MCkin(mZ,mS,r,CSmax[1],NEGNU); //Get Event

			if(eventBar[2]<mS){std::cout<<"ERROR: Event with energy less than mass, nuBarMode, En "<<eventBar[2]<<std::endl; ;return 0;}
			if(event[2]<mS){std::cout<<"ERROR: Event with energy less than mass, nuMODE, En "<<event[2]<<std::endl; ;return 0;}


			// for note,  event = {Evin,tout,Esout,costh,Edecay,Esout/ms, howmanyrngs}; in case your wondering
			vEn = event[0];
			csTmandel= event[1];
			sEn = event[2];
			sCos = event[3];
			sPhi = phiS; 

			vEn_B = eventBar[0];
			csTmandel_B= eventBar[1];
			sEn_B = eventBar[2];
			sCos_B = eventBar[3];
			sPhi_B = phiS; 


			initial_sterile nus(mS, sEn, sCos, sPhi);
			initial_sterile nusBar(mS, sEn_B, sCos_B, sPhi_B);

			sPx = nus.labframeP.p[0] ;
			sPy = nus.labframeP.p[1] ;
			sPz = nus.labframeP.p[2] ;
			sPx_B = nusBar.labframeP.p[0] ;
			sPy_B = nusBar.labframeP.p[1] ;
			sPz_B = nusBar.labframeP.p[2] ;



			//	CHAN->decayfunction(nus);			//oldest approximate one
			//CHAN->decayfunctionMassive(nus,ME,ME,0.0); // Marks' old one
			decay_params params;
			params.mS = mS;
			params.mZprime = mZ;
			
			params.chi = 0.01;
			params.Ue4 = 0.1;
			params.Um4 = 0.1;
			params.Ut4 = 0.1;
			CHAN->decayfunction(nus,&params); // a new one (4 Aug 2017)
			CHAN->observables(&Obs, r);


			//printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5f %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", Obs.E_sum, Obs.Th_sum, Obs.AngSep, Obs.E_sterile, Obs.Th_sterile, Obs.E_high, Obs.Th_high, Obs.E_low, Obs.Th_low, Obs.FS_AngSep,Obs.P_high[0],Obs.P_high[1],Obs.P_high[2],Obs.P_low[0],Obs.P_low[1],Obs.P_low[2],  Obs.E_sum_smear, Obs.Th_sum_smear, Obs.AngSep_smear, Obs.E_sterile_smear, Obs.Th_sterile_smear, Obs.E_high_smear, Obs.Th_high_smear, Obs.E_low_smear, Obs.Th_low_smear, Obs.FS_AngSep_smear, Obs.Minvar_smear);

			//quick check how its calculting angles



			if(Obs.E_low <= Obs.E_high){
				e1En = Obs.E_low;
				e1Cos = cos(Obs.Th_low*M_PI/180.0);
				e2En = Obs.E_high;
				e2Cos = cos(Obs.Th_high*M_PI/180.0);
				e1Px = Obs.P_low[0];
				e1Py = Obs.P_low[1];
				e1Pz = Obs.P_low[2];
				e2Px = Obs.P_high[0];
				e2Py = Obs.P_high[1];
				e2Pz = Obs.P_high[2];



			}else{
				e2En = Obs.E_low;
				e2Cos = cos(Obs.Th_low*M_PI/180.0);
				e1En = Obs.E_high;
				e1Cos = cos(Obs.Th_high*M_PI/180.0);
				e2Px = Obs.P_low[0];
				e2Py = Obs.P_low[1];
				e2Pz = Obs.P_low[2];
				e1Px = Obs.P_high[0];
				e1Py = Obs.P_high[1];
				e1Pz = Obs.P_high[2];


			}

			dEnSum  = Obs.E_sum;
			dAngSep = Obs.AngSep;
			dFsAngSep = Obs.FS_AngSep;
			dThSum = Obs.Th_sum;			

			t1->Fill();

			CHAN->decayfunction(nusBar,&params); // a new one (4 Aug 2017)
			CHAN->observables(&Obs, r);

			//CHAN->decayfunctionMassive(nusBar,ME,ME,0.0); //my new one, tobetested
			//		CHAN->decayfunction(nus);			//old approximate one
			//CHAN->observables(&Obs, r);


			if(Obs.E_low <= Obs.E_high){
				e1En_B = Obs.E_low;
				e1Cos_B = cos(Obs.Th_low*M_PI/180.0);
				e2En_B = Obs.E_high;
				e2Cos_B= cos(Obs.Th_high*M_PI/180.0);
				e1Px_B = Obs.P_low[0];
				e1Py_B = Obs.P_low[1];
				e1Pz_B = Obs.P_low[2];
				e2Px_B = Obs.P_high[0];
				e2Py_B = Obs.P_high[1];
				e2Pz_B = Obs.P_high[2];


			}else{
				e2En_B = Obs.E_low;
				e2Cos_B = cos(Obs.Th_low*M_PI/180.0);
				e1En_B = Obs.E_high;
				e1Cos_B = cos(Obs.Th_high*M_PI/180.0);
				e2Px_B  = Obs.P_low[0];
				e2Py_B  = Obs.P_low[1];
				e2Pz_B  = Obs.P_low[2];
				e1Px_B  = Obs.P_high[0];
				e1Py_B  = Obs.P_high[1];
				e1Pz_B  = Obs.P_high[2];


			}
			dEnSum_B  = Obs.E_sum;
			dAngSep_B = Obs.AngSep;
			dFsAngSep_B = Obs.FS_AngSep;
			dThSum_B = Obs.Th_sum;			

			t3->Fill();


		}


		t2->Fill();

		// Make some historgams





		//f1->Write();
		t1->Write("",TObject::kWriteDelete);
		t2->Write("",TObject::kWriteDelete);
		t3->Write("",TObject::kWriteDelete);
		//Cleaning up.
		f1->Close();
		gsl_rng_free(r);
		delete CHAN;
		delete DETECTOR;

		return 0;
	}
	/************************************************************************
	 *
	 *
	 *		ANAMODE anaMode AnAMode
	 *
	 *
	 *
	 * **********************************************************************/



	if(anaMode){

		double type1_low_thresh = 0.02;
		double type1_high_thresh = 0.1;

		double type2_ang_thresh = 5;	


		std::cout<<"Running anaMode with mZ: "<<mZ<<" GeV, and mS: "<<mS<<" GeV."<<std::endl;

		std::cout<<"Loading up correct file"<<std::endl;
		TFile *fana= new TFile(filename.c_str(),"UPDATE");
		std::cout<<"Accessing SummeryTree to check if file is correct"<<std::endl;
		TTree *tsum = (TTree *)fana->Get("summeryTree");
		std::cout<<"SummeryTree loaded"<<std::endl;	
		//some simple checking of the summeryTree
		double tsum_mS, tsum_mZ, tsum_TotalNumNu, tsum_TotalNumNuBar;
		tsum->SetBranchAddress("zM",&tsum_mZ);
		tsum->SetBranchAddress("sM",&tsum_mS);
		tsum->SetBranchAddress("csTotalNumNu",&tsum_TotalNumNu);
		tsum->SetBranchAddress("csTotalNumNuBar",&tsum_TotalNumNuBar);

		for(int i=0;i<tsum->GetEntries();i++){
			tsum->GetEntry(i);
			std::cout<<"From CommandLine, mZ: "<<mZ<<", mS: "<<mS<<std::endl;
			std::cout<<"From Summery tree, mZ: "<<tsum_mZ<<", mS: "<<tsum_mS<<std::endl;
			if(tsum_mZ!=mZ || tsum_mS != mS){
				std::cout<<"ERROR: Suposed to be, mZ: "<<mZ<<", mS: "<<mS<<std::endl;

			}

			if(i>0) std::cout<<"ERROR: anaMode, more that 1 entry in summery tree"<<std::endl;
		}

		std::cout<<"Accessing nuTree and nubarTree"<<std::endl;
		TTree *tnu = (TTree *)fana->Get("nuTree");
		int Nnu = tnu->GetEntries();
		std::cout<<"nuTree loaded, "<<Nnu<<" events."<<std::endl;	
		TTree *tnubar = (TTree *)fana->Get("nubarTree");
		int Nnubar = tnubar->GetEntries();
		std::cout<<"nubarTree loaded, "<<Nnubar<<" events."<<std::endl;	

		//And creating a new ana subdirectory, and histograms that will go in said directory
		//If directory already exists, then simpley wipe it and work from there!
		TDirectory *anaDir = fana->GetDirectory("anaDir");
		if (anaDir) {
			std::cout<<"WARNING: anaDir exists, wiping previous run"<<std::endl;
			anaDir->cd();
			anaDir->Delete("*;*");
		} else {
			anaDir = fana->mkdir("anaDir");	
		}

		//make new TTree's for neutrino and antineutrino mode
		TTree *t_nu = new TTree("anatree_nu","anatree_nu");
		TTree *t_nubar = new TTree("anatree_nubar","anatree_nubar");
		//and their associated things
		double evis_nu, evis_nubar;
		double enu_nu, enu_nubar;
		double cos_nu, cos_nubar;

		t_nu->Branch("Es",&enu_nu);
		t_nu->Branch("Evis",&evis_nu);
		t_nu->Branch("Cos",&cos_nu);
		t_nubar->Branch("Es",&enu_nubar);
		t_nubar->Branch("Evis",&evis_nubar);
		t_nubar->Branch("Cos",&cos_nubar);


		TH1D  *hType1_Evis_nu = new TH1D("hType1Signal_Evis_nu","Type1Signal Visible En Nu",EBINS,EMIN,EMAX);
		TH1D  *hType1_Evis_nubar = new TH1D("hType1Signal_Evis_nubar","Type1Signal Visible En NuBar",EBINS,EMIN,EMAX);
		TH1D  *hType1_Cos_nu = new TH1D("hType1Signal_Cos_nu","Type1Signal Cos Nu",COSBINS,COSMIN,COSMAX);
		TH1D  *hType1_Cos_nubar = new TH1D("hType1Signal_Cos_nubar","Type1Signal Cos NuBar",COSBINS,COSMIN,COSMAX);

		TH1D  *hType2_Evis_nu = new TH1D("hType2Signal_Evis_nu","Type2Signal Visible En Nu",EBINS,EMIN,EMAX);
		TH1D  *hType2_Evis_nubar = new TH1D("hType2Signal_Evis_nubar","Type2Signal Visible En NuBar",EBINS,EMIN,EMAX);
		TH1D  *hType2_Cos_nu = new TH1D("hType2Signal_Cos_nu","Type2Signal Cos Nu",COSBINS,COSMIN,COSMAX);
		TH1D  *hType2_Cos_nubar = new TH1D("hType2Signal_Cos_nubar","Type2Signal Cos NuBar",COSBINS,COSMIN,COSMAX);

		TH1D  *hSignal_Evis_nu = new TH1D("hSignalSignal_Evis_nu","SignalSignal Visible En Nu",EBINS,EMIN,EMAX);
		TH1D  *hSignal_Evis_nubar = new TH1D("hSignalSignal_Evis_nubar","SignalSignal Visible En NuBar",EBINS,EMIN,EMAX);
		TH1D  *hSignal_Cos_nu = new TH1D("hSignalSignal_Cos_nu","SignalSignal Cos Nu",COSBINS,COSMIN,COSMAX);
		TH1D  *hSignal_Cos_nubar = new TH1D("hSignalSignal_Cos_nubar","SignalSignal Cos NuBar",COSBINS,COSMIN,COSMAX);

		TH1D *hType1_Ester_nu= new TH1D("hType1_Ester_nu","Type1 Parent Sterile En Nu",25,0,3);
		TH1D *hType2_Ester_nu= new TH1D("hType2_Ester_nu","Type2 Parent Sterile En Nu",25,0,3);
		TH1D *hType1_Ester_nubar= new TH1D("hType1_Ester_nubar","Type1 Parent Sterile En NuBar",25,0,3);
		TH1D *hType2_Ester_nubar= new TH1D("hType2_Ester_nubar","Type2 Parent Sterile En NuBar",25,0,3);


		double sEn;
		double sPhi;
		double sCos;
		double vEn;
		double e1En;
		double e2En;
		double e1Cos;
		double e2Cos;
		double csTmandel;
		double dAngSep;
		double dFsAngSep;
		double dEnSum;
		double dThSum;

		double sEn_B;
		double sPhi_B;
		double sCos_B;
		double vEn_B;
		double e1En_B;
		double e2En_B;
		double e1Cos_B;
		double e2Cos_B;
		double csTmandel_B;
		double dAngSep_B;
		double dFsAngSep_B;
		double dEnSum_B;
		double dThSum_B;

		tnu->SetBranchAddress("csTmandel", &csTmandel);	
		tnu->SetBranchAddress("sEn", &sEn);	
		tnu->SetBranchAddress("sPhi",&sPhi);
		tnu->SetBranchAddress("sCos", &sCos);	
		tnu->SetBranchAddress("vEn", &vEn);	

		tnu->SetBranchAddress("e1En", &e1En);	
		tnu->SetBranchAddress("e1Cos", &e1Cos);	
		tnu->SetBranchAddress("e2En", &e2En);	
		tnu->SetBranchAddress("e2Cos", &e2Cos);	

		tnu->SetBranchAddress("dAngSep",&dAngSep);
		tnu->SetBranchAddress("dFsAngSep",&dFsAngSep);
		tnu->SetBranchAddress("dThSum",&dThSum);
		tnu->SetBranchAddress("dEnSum",&dEnSum);

		// branches
		tnubar->SetBranchAddress("csTmandel", &csTmandel_B);	
		tnubar->SetBranchAddress("sEn", &sEn_B);	
		tnubar->SetBranchAddress("sPhi",&sPhi_B);
		tnubar->SetBranchAddress("sCos", &sCos_B);	
		tnubar->SetBranchAddress("vEn", &vEn_B);	

		tnubar->SetBranchAddress("e1En", &e1En_B);	
		tnubar->SetBranchAddress("e1Cos", &e1Cos_B);	
		tnubar->SetBranchAddress("e2En", &e2En_B);	
		tnubar->SetBranchAddress("e2Cos", &e2Cos_B);	

		tnubar->SetBranchAddress("dAngSep",&dAngSep_B);
		tnubar->SetBranchAddress("dFsAngSep",&dFsAngSep_B);
		tnubar->SetBranchAddress("dThSum",&dThSum_B);
		tnubar->SetBranchAddress("dEnSum",&dEnSum_B);

		if(Nnu!=Nnubar){ std::cout<<"ERROR: anaMode, nu and Nubar trees have different numbers of entries"<<std::endl;}

		for(int j=0; j< Nnu; j++){
			if(j%10000==0) std::cout<<"#: "<<j<<std::endl;
			tnu->GetEntry(j);
			tnubar->GetEntry(j);

			bool is_this_an_event_nu = false;
			bool is_this_an_event_nubar = false;

			double e1En_smr = e1En;
			double e2En_smr = e2En;
			double e1Cos_smr = e1Cos;
			double e2Cos_smr = e2Cos;
			double dAngSep_smr = dAngSep;
			double dThSum_smr =dThSum;

			double e1En_B_smr = e1En_B;
			double e2En_B_smr = e2En_B;
			double e1Cos_B_smr = e1Cos_B;
			double e2Cos_B_smr = e2Cos_B;
			double dAngSep_B_smr = dAngSep_B;
			double dThSum_B_smr =dThSum_B;




			double weight_nu = 1;//mS/sqrt(sEn*sEn-mS*mS);  
			double weight_nubar = 1;//mS/sqrt(sEn_B*sEn_B-mS*mS);  
			if(weight_nu != weight_nu || weight_nubar != weight_nubar){std::cout<<"ERROR:Infinite weight, mS: "<<mS<<" En_nu "<<sEn<<" En_nubar " <<sEn_B<<" , run "<<j<<std::endl;}


			// Neutrino run
			if( e1En_smr < type1_low_thresh && e2En_smr >=type1_high_thresh ) // So Called Type 1 event
			{
				hType1_Evis_nu->Fill(e2En_smr, weight_nu);	
				hType1_Cos_nu->Fill(e2Cos_smr, weight_nu);
				hSignal_Evis_nu->Fill(e2En_smr, weight_nu);//and add to signal too	
				hSignal_Cos_nu->Fill(e2Cos_smr, weight_nu);	

				hType1_Ester_nu->Fill(sEn);
				is_this_an_event_nu = true;
				enu_nu = sEn;
				evis_nu = e2En_smr;
				cos_nu = e2Cos_smr;

			}else if(fabs(dAngSep_smr) < type2_ang_thresh) // So called Type 2 event
			{
				hType2_Evis_nu->Fill(e2En_smr+e1En_smr, weight_nu);	
				hType2_Cos_nu->Fill(cos(dThSum_smr*M_PI/180.0), weight_nu);	
				hSignal_Evis_nu->Fill(e2En_smr+e1En_smr, weight_nu);	
				hSignal_Cos_nu->Fill(cos(dThSum_smr*M_PI/180.0), weight_nu);	

				hType2_Ester_nu->Fill(sEn);
				is_this_an_event_nu = true;
				enu_nu = sEn;
				evis_nu =e2En_smr+e1En_smr;
				cos_nu = cos(dThSum_smr*M_PI/180.0);
			}


			// Anti-Neutrino run
			if( e1En_B_smr < type1_low_thresh && e2En_B_smr >=type1_high_thresh ) // So Called Type 1 event
			{
				hType1_Evis_nubar->Fill(e2En_B_smr, weight_nubar);	
				hType1_Cos_nubar->Fill(e2Cos_B_smr, weight_nubar);
				hSignal_Evis_nubar->Fill(e2En_B_smr, weight_nubar);//and add to signal too	
				hSignal_Cos_nubar->Fill(e2Cos_B_smr, weight_nubar);

				hType1_Ester_nubar->Fill(sEn_B);
				is_this_an_event_nubar = true;
				enu_nubar = sEn_B;
				evis_nubar = e2En_B_smr;
				cos_nubar = e2Cos_B_smr;

			}else if(fabs(dAngSep_B_smr) < type2_ang_thresh) // So called Type 2 event
			{
				hType2_Evis_nubar->Fill(e2En_B_smr+e1En_B_smr, weight_nubar);	
				hType2_Cos_nubar->Fill(cos(dThSum_B_smr*M_PI/180.0), weight_nubar);	
				hSignal_Evis_nubar->Fill(e2En_B_smr+e1En_B_smr, weight_nubar);	
				hSignal_Cos_nubar->Fill(cos(dThSum_B_smr*M_PI/180.0), weight_nubar);	

				hType2_Ester_nubar->Fill(sEn_B);
				is_this_an_event_nubar = true;
				enu_nubar = sEn_B;
				evis_nubar =e2En_B_smr+e1En_B_smr;
				cos_nubar = cos(dThSum_B_smr*M_PI/180.0);

			}

			if(is_this_an_event_nu){
				t_nu->Fill();
			}

			if(is_this_an_event_nubar){
				t_nubar->Fill();
			}



		}


		fana->cd("anaDir");
		hType1_Evis_nu->Write();
		hType1_Cos_nu->Write();
		hSignal_Evis_nu->Write();
		hSignal_Cos_nu->Write();
		hType2_Evis_nu->Write();
		hType2_Cos_nu->Write();

		hType1_Evis_nubar->Write();
		hType1_Cos_nubar->Write();
		hSignal_Evis_nubar->Write();
		hSignal_Cos_nubar->Write();
		hType2_Evis_nubar->Write();
		hType2_Cos_nubar->Write();

		hType1_Ester_nu->Write();
		hType2_Ester_nu->Write();
		hType1_Ester_nubar->Write();
		hType2_Ester_nubar->Write();

		t_nu->Write();
		t_nubar->Write();
		

		fana->Close();


		return 0;
	}


	if(statMode){

		std::cout<<"Running statMode with mZ: "<<mZ<<" GeV, and mS: "<<mS<<" GeV."<<std::endl;
		std::cout<<"Loading up correct file"<<std::endl;
		TFile *fstat= new TFile(filename.c_str(),"UPDATE");
		std::cout<<"Accessing SummeryTree to check if file is correct"<<std::endl;
		TTree *tsum = (TTree *)fstat->Get("summeryTree");
		std::cout<<"SummeryTree loaded"<<std::endl;	
		//some simple checking of the summeryTree
		double Norm_nu=0;
		double Norm_nubar =0;

		double tsum_mS, tsum_mZ, tsum_TotalNumNu, tsum_TotalNumNuBar;
		tsum->SetBranchAddress("zM",&tsum_mZ);
		tsum->SetBranchAddress("sM",&tsum_mS);
		tsum->SetBranchAddress("csTotalNumNu",&tsum_TotalNumNu);
		tsum->SetBranchAddress("csTotalNumNuBar",&tsum_TotalNumNuBar);

		for(int i=0;i<tsum->GetEntries();i++){
			tsum->GetEntry(i);
			std::cout<<"From CommandLine, mZ: "<<mZ<<", mS: "<<mS<<std::endl;
			std::cout<<"From Summery tree, mZ: "<<tsum_mZ<<", mS: "<<tsum_mS<<std::endl;
			if(tsum_mZ!=mZ || tsum_mS != mS){
				std::cout<<"ERROR: Suposed to be, mZ: "<<mZ<<", mS: "<<mS<<std::endl;

			}


			Norm_nu =tsum_TotalNumNu;
			Norm_nubar =tsum_TotalNumNuBar;

			if(i>0) std::cout<<"ERROR: statMode, more that 1 entry in summery tree"<<std::endl;
		}

		//And creating a new ana subdirectory, and histograms that will go in said directory
		//If directory already exists, then simpley wipe it and work from there!
		TDirectory *statDir = fstat->GetDirectory("statDir");
		if (statDir) {
			std::cout<<"WARNING: statDir exists, wiping previous run"<<std::endl;
			statDir->cd();
			statDir->Delete("*;*");
		} else {
			statDir = fstat->mkdir("statDir");	
		}



		//Ok from here on its going to be pure stats baby
		//Enter the anaDir subdirectory

		TDirectory *anaDir = fstat->GetDirectory("anaDir");
		if (anaDir) {
			anaDir->cd();
		} else {
			std::cout<<"ERROR: anaDir does not exit, cannot run stat mode, please run --mode ana first"<<std::endl;
			exit(EXIT_FAILURE);
		}

		//Load up the ana_tree
		std::cout<<"Accessing ana_nuTree and ana_nubarTree"<<std::endl;
		fstat->ls();
		TTree *tAna_nu = (TTree *)anaDir->Get("anatree_nu");
		TTree *tAna_nubar = (TTree *)anaDir->Get("anatree_nubar");
		std::cout<<"Done. Gotten pointers to TTrees of interest"<<std::endl;
		tAna_nu->Print();
		tAna_nubar->Print();


		//create sum good signals
		TH1D * hSignal_Evis_nu = (TH1D*)anaDir->Get("hSignalSignal_Evis_nu");	
		TH1D * hSignal_Cos_nu = (TH1D*)anaDir->Get("hSignalSignal_Cos_nu");	
		TH1D * hSignal_Evis_nubar = (TH1D*)anaDir->Get("hSignalSignal_Evis_nubar");	
		TH1D * hSignal_Cos_nubar = (TH1D*)anaDir->Get("hSignalSignal_Cos_nubar");	

		bound bound_ps191(core+"data/bounds/PS191_UM4_EE_BOTH.dat",0.01,128);
		bound_ps191.setTypicalEnergy(5.0);
		bound bound_peak(core+"data/bounds/peak_um4.dat",0.00,128);
		bound bound_babar(core+"data/bounds/b1_babar2014.csv",0.00,128);

		bound bound_gm2(core+"data/bounds/b1_mg2.csv",0.00,128);
		bound bound_nutev(core+"data/bounds/nutev_muon.csv",0.213,1450);
		bound_nutev.setTypicalEnergy(100);



		//Right so returning to statDir, we basicall have all ingredients, 4 vSignal vectors, 4 Observed Vectors and 4 Expected background vector and 2 norms
		// Norm_nu and Norm_nubar
		std::cout<<"Going to initilize statInatance now"<<std::endl;
		minInstance statInstance(Norm_nu, Norm_nubar, tAna_nu, tAna_nubar  , (double)NUMEVENTS);
		std::cout<<"Done. initilized statInatance"<<std::endl;

		statInstance.setMass(mZ,mS);
		statInstance.bound_vector.push_back(bound_ps191);
		statInstance.bound_vector.push_back(bound_peak);
		statInstance.bound_vector.push_back(bound_babar);
		statInstance.bound_vector.push_back(bound_gm2);
		statInstance.bound_vector.push_back(bound_nutev);

		statInstance.use_bounds=true;
		std::cout<<"Set Bounds"<<std::endl;


		std::cout<<"Norm nu: "<< statInstance.norm_nu<<std::endl;
		std::cout<<"Norm nu_bar: "<< statInstance.norm_nubar<<std::endl;
		double ans_chi = statInstance.minimize();

		//Todo, get number correct using histogrammer from old LR.h file
		//get bounds implemented..
		//tweak minimizer!

		std::cout<<"*********************************************************"<<std::endl;
		std::cout<<"#: mZ = "<<mZ<<std::endl;
		std::cout<<"#: mS = "<<mS<<std::endl;
		std::cout<<"#: Um4^2 = "<<2*(statInstance.bf_up)<<std::endl;
		std::cout<<"#: chi^2 = "<<2*(statInstance.bf_chi)<<std::endl;
		std::cout<<"#: ANS "<<mZ<<" "<<mS<<" "<<2*(statInstance.bf_up)<<" "<<2*(statInstance.bf_chi)<<" "<<ans_chi<<std::endl;
//Peter changed ths 11/July
//		double thisrate =bound_ps191.myRate( mS,mZ)*pow(pow(10,statInstance.bf_chi),2.0);
	//	double thisrate =bound_ps191.myRate(pow(10,statInstance.bf_chi),mS,mZ);
		
	//	std::cout<<"#: DecayRate = "<<thisrate<<std::endl;
	//	std::cout<<"#: DecayLength = "<<bound_ps191.gev2meters(1.0/thisrate)<<std::endl;
		std::cout<<"#: NumScatters (neutrino mode) = "<<statInstance.norm_nu*pow(pow(10,statInstance.bf_chi),2.0)*pow(pow(10,statInstance.bf_up),2.0)<<std::endl;
		std::cout<<"#: NumScatters (antineutrino mode) = "<<statInstance.norm_nubar*pow(pow(10,statInstance.bf_chi),2.0)*pow(pow(10,statInstance.bf_up),2.0)<<std::endl;
		std::cout<<"#: NumDecays (neutrino mode) = "<<statInstance.N_bf<<std::endl;
		std::cout<<"#: NumDecays (antineutrino mode) = "<<statInstance.N_bf_bar<<std::endl;



		fstat->cd("statDir");
		statInstance.h_bkg_E_nu->Write();
		statInstance.h_bkg_E_nubar->Write();
		statInstance.h_bkg_C_nu->Write();
		statInstance.h_bkg_C_nubar->Write();

		statInstance.h_bf_E_nu->Write();
		statInstance.h_bf_E_nubar->Write();
		statInstance.h_bf_C_nu->Write();
		statInstance.h_bf_C_nubar->Write();

		statInstance.h_obs_E_nu->Write();
		statInstance.h_obs_E_nubar->Write();
		statInstance.h_obs_C_nu->Write();
		statInstance.h_obs_C_nubar->Write();


		double bkg_bf_zeta_nu = (statInstance.bf_zeta_b_nu);
		double bkg_bf_zeta_nubar = (statInstance.bf_zeta_b_nubar);

		//###################################################################################
		//			Plot some useful things for checks
		//###################################################################################
		TCanvas *c1 = new TCanvas("MiniBooNE","c1",1200,800);
		c1->Divide(2,2);
		c1->cd(1);
		THStack * h1 	   = new THStack("Evis Neutrino Mode","Evis Neutrino Mode");	
		statInstance.h_bkg_E_nu->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_E_nu->SetFillColor(kBlue-7);
		statInstance.h_bkg_E_nu->SetLineColor(kBlack);
		statInstance.h_bkg_E_nu->SetTitle("Evis Neutrino Mode");

		statInstance.h_bf_E_nu->SetMarkerColor(kOrange+6);
		statInstance.h_bf_E_nu->SetFillColor(kOrange+6);
		statInstance.h_bf_E_nu->SetLineColor(kBlack);


		h1->Add(statInstance.h_bkg_E_nu);
		h1->Add(statInstance.h_bf_E_nu);

		h1->Draw();
		statInstance.h_obs_E_nu->Scale(1);
		statInstance.h_obs_E_nu->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_obs_E_nu->GetXaxis()->SetTitle("E_{vis} (GeV)");
		statInstance.h_obs_E_nu->SetMarkerColor(kBlack);
		statInstance.h_obs_E_nu->Draw("same");			

		c1->cd(2);
		THStack * h2 	   = new THStack("cos #theta Neutrino Mode","cos #theta Neutrino Mode");	
		statInstance.h_bkg_C_nu->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_C_nu->SetFillColor(kBlue-7);
		statInstance.h_bkg_C_nu->SetLineColor(kBlack);
		statInstance.h_bkg_C_nu->SetTitle("cos #theta Neutrino Mode");

		statInstance.h_bf_C_nu->SetMarkerColor(kOrange+6);
		statInstance.h_bf_C_nu->SetFillColor(kOrange+6);
		statInstance.h_bf_C_nu->SetLineColor(kBlack);


		h2->Add(statInstance.h_bkg_C_nu);
		h2->Add(statInstance.h_bf_C_nu);

		h2->Draw();
		//	statInstance.h_obs_C_nu->Scale(1,"width");
		statInstance.h_obs_C_nu->Scale(1);
		statInstance.h_obs_C_nu->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_obs_C_nu->GetXaxis()->SetTitle("cos #theta (GeV)");
		statInstance.h_obs_C_nu->SetMarkerColor(kBlack);
		statInstance.h_obs_C_nu->Draw("same");			
		c1->cd(3);
		THStack * h3 	   = new THStack("Evis Anti Neutrino Mode","Evis Anti Neutrino Mode");	
		statInstance.h_bkg_E_nubar->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_E_nubar->SetFillColor(kBlue-7);
		statInstance.h_bkg_E_nubar->SetLineColor(kBlack);
		statInstance.h_bkg_E_nubar->SetTitle("Evis Anti Neutrino Mode");

		statInstance.h_bf_E_nubar->SetMarkerColor(kOrange+6);
		statInstance.h_bf_E_nubar->SetFillColor(kOrange+6);
		statInstance.h_bf_E_nubar->SetLineColor(kBlack);


		h3->Add(statInstance.h_bkg_E_nubar);
		h3->Add(statInstance.h_bf_E_nubar);

		h3->Draw();
		statInstance.h_obs_E_nubar->Scale(1);
		statInstance.h_obs_E_nubar->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_obs_E_nubar->GetXaxis()->SetTitle("E_{vis} (GeV)");
		statInstance.h_obs_E_nubar->SetMarkerColor(kBlack);
		statInstance.h_obs_E_nubar->Draw("same");			
		c1->cd(4);
		THStack * h4 	   = new THStack("cos #theta Anti Neutrino Mode","cos #theta Anti Neutrino Mode");	
		statInstance.h_bkg_C_nubar->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_C_nubar->SetFillColor(kBlue-7);
		statInstance.h_bkg_C_nubar->SetLineColor(kBlack);
		statInstance.h_bkg_C_nubar->SetTitle("cos #theta Anti  Neutrino Mode");

		statInstance.h_bf_C_nubar->SetMarkerColor(kOrange+6);
		statInstance.h_bf_C_nubar->SetFillColor(kOrange+6);
		statInstance.h_bf_C_nubar->SetLineColor(kBlack);


		h4->Add(statInstance.h_bkg_C_nubar);
		h4->Add(statInstance.h_bf_C_nubar);

		h4->Draw();
		statInstance.h_obs_C_nubar->Scale(1);
		statInstance.h_obs_C_nubar->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_obs_C_nubar->GetXaxis()->SetTitle("cos #theta (GeV)");
		statInstance.h_obs_C_nubar->SetMarkerColor(kBlack);
		statInstance.h_obs_C_nubar->Draw("same");			
		c1->Write();	


		//###################################################################################
		//			Plot some useful things for checks
		//###################################################################################

		TCanvas *c2 = new TCanvas("MiniBooNE Excess","c2",1200,800);
		c2->Divide(2,2);
		c2->cd(1);
		THStack * hh1 	   = new THStack("Evis Neutrino Mode","Evis Neutrino Mode");	
		statInstance.h_excess_E_nu->Scale(1);
		statInstance.h_excess_E_nu->Draw();	

		statInstance.h_bkg_E_nu->Scale( bkg_bf_zeta_nu ,"nosw2");
		statInstance.h_bkg_E_nu->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_E_nu->SetFillColor(kBlue-7);
		statInstance.h_bkg_E_nu->SetLineColor(kBlack);
		statInstance.h_bkg_E_nu->SetTitle("Evis Neutrino Mode");

		statInstance.h_bf_E_nu->SetMarkerColor(kOrange+6);
		statInstance.h_bf_E_nu->SetFillColor(kOrange+6);
		statInstance.h_bf_E_nu->SetLineColor(kBlack);


		hh1->Add(statInstance.h_bkg_E_nu);
		hh1->Add(statInstance.h_bf_E_nu);

		hh1->Draw("same");

		statInstance.h_excess_E_nu->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_excess_E_nu->GetXaxis()->SetTitle("E_{vis} (GeV)");
		statInstance.h_excess_E_nu->SetMarkerColor(kBlack);
		statInstance.h_excess_E_nu->SetMarkerStyle(20);	

		statInstance.h_excess_E_nu->Draw("same");	


		c2->cd(2);
		THStack * hh2 	   = new THStack("cos #theta Neutrino Mode","cos #theta Neutrino Mode");	
		statInstance.h_excess_C_nu->Scale(1);
		statInstance.h_excess_C_nu->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_excess_C_nu->GetXaxis()->SetTitle("cos #theta (GeV)");
		statInstance.h_excess_C_nu->Draw();

		statInstance.h_bkg_C_nu->Scale( bkg_bf_zeta_nu ,"nosw2");
		statInstance.h_bkg_C_nu->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_C_nu->SetFillColor(kBlue-7);
		statInstance.h_bkg_C_nu->SetLineColor(kBlack);
		statInstance.h_bkg_C_nu->SetTitle("cos #theta Neutrino Mode");

		statInstance.h_bf_C_nu->SetMarkerColor(kOrange+6);
		statInstance.h_bf_C_nu->SetFillColor(kOrange+6);
		statInstance.h_bf_C_nu->SetLineColor(kBlack);


		hh2->Add(statInstance.h_bkg_C_nu);
		hh2->Add(statInstance.h_bf_C_nu);

		hh2->Draw("same");
		//	statInstance.h_excess_C_nu->Scale(1,"width");

		statInstance.h_excess_C_nu->SetMarkerStyle(20);	
		statInstance.h_excess_C_nu->SetMarkerColor(kBlack);
		statInstance.h_excess_C_nu->Draw("same");


		c2->cd(3);
		THStack * hh3 	   = new THStack("Evis Anti Neutrino Mode","Evis Anti Neutrino Mode");	
		statInstance.h_excess_E_nubar->Scale(1);
		statInstance.h_excess_E_nubar->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_excess_E_nubar->GetXaxis()->SetTitle("E_{vis} (GeV)");
		statInstance.h_excess_E_nubar->Draw();			

		statInstance.h_bkg_E_nubar->Scale( bkg_bf_zeta_nubar ,"nosw2");
		statInstance.h_bkg_E_nubar->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_E_nubar->SetFillColor(kBlue-7);
		statInstance.h_bkg_E_nubar->SetLineColor(kBlack);
		statInstance.h_bkg_E_nubar->SetTitle("Evis Anti Neutrino Mode");

		statInstance.h_bf_E_nubar->SetMarkerColor(kOrange+6);
		statInstance.h_bf_E_nubar->SetFillColor(kOrange+6);
		statInstance.h_bf_E_nubar->SetLineColor(kBlack);


		hh3->Add(statInstance.h_bkg_E_nubar);
		hh3->Add(statInstance.h_bf_E_nubar);

		hh3->Draw("same");


		statInstance.h_excess_E_nubar->SetMarkerColor(kBlack);
		statInstance.h_excess_E_nubar->SetMarkerStyle(20);
		statInstance.h_excess_E_nubar->Draw("same");			



		c2->cd(4);
		THStack * hh4 	   = new THStack("cos #theta Anti Neutrino Mode","cos #theta Anti Neutrino Mode");	

		statInstance.h_excess_C_nubar->Scale(1);
		statInstance.h_excess_C_nubar->GetYaxis()->SetTitle("Events/GeV");
		statInstance.h_excess_C_nubar->GetXaxis()->SetTitle("cos #theta (GeV)");
		statInstance.h_excess_C_nubar->SetMarkerColor(kBlack);
		statInstance.h_excess_C_nubar->Draw();			


		statInstance.h_bkg_C_nubar->Scale( bkg_bf_zeta_nubar ,"nosw2");

		statInstance.h_bkg_C_nubar->SetMarkerColor(kBlue-7);
		statInstance.h_bkg_C_nubar->SetFillColor(kBlue-7);
		statInstance.h_bkg_C_nubar->SetLineColor(kBlack);
		statInstance.h_bkg_C_nubar->SetTitle("cos #theta Anti  Neutrino Mode");

		statInstance.h_bf_C_nubar->SetMarkerColor(kOrange+6);
		statInstance.h_bf_C_nubar->SetFillColor(kOrange+6);
		statInstance.h_bf_C_nubar->SetLineColor(kBlack);


		hh4->Add(statInstance.h_bkg_C_nubar);
		hh4->Add(statInstance.h_bf_C_nubar);

		hh4->Draw("same");

		statInstance.h_excess_C_nubar->SetMarkerStyle(20);
		statInstance.h_excess_C_nubar->SetMarkerColor(kBlack);
		statInstance.h_excess_C_nubar->Draw("same");			


		c2->Write();	


		fstat->Close();
		return 0;
	}


	if(boundMode){
		bound bound_ps191("../../data/bounds/PS191_UM4_EE_BOTH.dat",0.01,128);
		bound_ps191.setTypicalEnergy(5.0);
		bound bound_peak("../../data/bounds/peak_um4.dat",0.00,128);
		bound bound_babar("../../data/babar_data.dat",0.00,128);
		bound bound_gm2("../../data/bounds/b1_mg2.csv",0.00,128);
		bound bound_nutev("../../data/bounds/nutev_muon.csv",0.213,1450);
		bound_nutev.setTypicalEnergy(100);

		if(false){ 
			for(double ims=log10(0.3); ims<log10(0.5); ims+=0.01){
				double ms= pow(10,ims);
				//	std::cout<<ms<<" "<<bound_ps191.bound_file.getFlux(ms)<<" "<<bound_nutev.bound_file.getFlux(ms)<<" "<<bound_peak.bound_file.getFlux(ms)<<std::endl;
			}

			double mZ=1.0;
			double bf_ms_um4=-20;
			double bf_ms_chi2=-20;
			double bf_ms_mz=0;

			for(double imS = log10(0.0101); imS<log10(0.5); imS+=0.05){
				double mS=pow(10,imS);


				double bf_ch2=-10;
				double bf_um2=-10;

				double cm = std::max( 0.5*log10(bound_babar.bound_file.getFlux(mZ)), log10(bound_gm2.bound_file.getFlux(mZ)) );
				for(double ich2=cm; ich2>-10; ich2-=0.05){
					for(double im2 = 0.5*log10(bound_peak.bound_file.getFlux(mS)); im2 >-10   ; im2-=0.05){

						bool isPS191 = bound_ps191.ps191(mS,mZ, pow(pow(10,im2),2), pow(pow(10,ich2),2));
						bool isNUTEV = bound_nutev.ps191(mS,mZ, pow(pow(10,im2),2), pow(pow(10,ich2),2));

						bool isPEAK = bound_peak.asIs(mS, pow(pow(10,im2),2)   );
						bool isBABAR = bound_babar.asIs(mZ, pow(pow(10,ich2),2) );
						bool isMGMIN2 = bound_gm2.asIs(mZ, pow(10,ich2));//not on square for now

						if(isPS191 && isPEAK && isBABAR && isNUTEV && isMGMIN2){
							if( pow(10,ich2)*pow(10,im2) >  pow(10, bf_ch2)*pow(10,bf_um2)){

								bf_ch2 = ich2;
								bf_um2  = im2;
							}

							if( im2 > bf_ms_um4){

								bf_ms_um4 = im2;
								bf_ms_mz = mZ;
								bf_ms_chi2 = ich2;
							}

						}

					}
				}


				std::cout<<mS<<" "<<mZ<<" "<<bf_ch2<<" "<<bf_um2<<" "<<pow(pow(10, bf_ch2)*pow(10,bf_um2),2.0)<<std::endl;
			}//end mz





			return 0;
		}


		if(false){
			for(double mZ = 0.5; mZ<=2; mZ=mZ+0.02){
				double bf_ms_um4=-20;
				double bf_ms_chi2=-20;
				double bf_ms_mz=0;

				for(double imS = log10(0.0101); imS<log10(0.5); imS+=0.05){
					double mS=pow(10,imS);


					double bf_ch2=-10;
					double bf_um2=-10;

					double cm = std::max( 0.5*log10(bound_babar.bound_file.getFlux(mZ)), log10(bound_gm2.bound_file.getFlux(mZ)) );
					for(double ich2=cm; ich2>-10; ich2-=0.025){
						for(double im2 = 0.5*log10(bound_peak.bound_file.getFlux(mS)); im2 >-10   ; im2-=0.025){

							bool isPS191 = bound_ps191.ps191(mS,mZ, pow(pow(10,im2),2), pow(pow(10,ich2),2));
							bool isNUTEV = bound_nutev.ps191(mS,mZ, pow(pow(10,im2),2), pow(pow(10,ich2),2));

							bool isPEAK = bound_peak.asIs(mS, pow(pow(10,im2),2)   );
							bool isBABAR = bound_babar.asIs(mZ, pow(pow(10,ich2),2) );
							bool isMGMIN2 = bound_gm2.asIs(mZ, pow(10,ich2));//not on square for now

							if(isPS191 && isPEAK && isBABAR && isNUTEV && isMGMIN2){
								if( pow(10,ich2)*pow(10,im2) >  pow(10, bf_ch2)*pow(10,bf_um2)){

									bf_ch2 = ich2;
									bf_um2  = im2;
								}

								if( im2 > bf_ms_um4){

									bf_ms_um4 = im2;
									bf_ms_mz = mZ;
									bf_ms_chi2 = ich2;
								}

							}

						}
					}


					std::cout<<mS<<" "<<mZ<<" "<<bf_ch2<<" "<<bf_um2<<" "<<pow(pow(10, bf_ch2)*pow(10,bf_um2),2.0)<<std::endl;
				}//end mz


				std::cout<<"#mz: "<<mZ<<" "<<bf_ms_um4<<" "<<bf_ms_mz<<" "<<bf_ms_chi2 <<" "<<pow(pow(10,bf_ms_chi2),2.0)<<std::endl;
			}//end mS





			return 0;
		}

		for(double imS = log10(0.5); imS > log10(0.0101); imS-=0.05){
			double mS=pow(10,imS);

			double bf_ms_um4=-20;
			double bf_ms_chi2=-20;
			double bf_ms_mz=0;

			for(double mZ = 0.6; mZ<=2; mZ=mZ+0.05){
				double bf_ch2=-10;
				double bf_um2=-10;

				double cm = std::max( 0.5*log10(bound_babar.bound_file.getFlux(mZ)), log10(bound_gm2.bound_file.getFlux(mZ)) );
				for(double ich2=cm; ich2>-10; ich2-=0.033){
					for(double im2 = 0.5*log10(bound_peak.bound_file.getFlux(mS)); im2 >-10   ; im2-=0.033){

						bool isPS191 = bound_ps191.ps191(mS,mZ, pow(pow(10,im2),2), pow(pow(10,ich2),2));
						bool isNUTEV = bound_nutev.ps191(mS,mZ, pow(pow(10,im2),2), pow(pow(10,ich2),2));
						bool isPEAK = bound_peak.asIs(mS, pow(pow(10,im2),2)   );
						bool isBABAR = bound_babar.asIs(mZ, pow(pow(10,ich2),2) );
						bool isMGMIN2 = bound_gm2.asIs(mZ, pow(10,ich2));//not on square for now

						if(isPS191 && isPEAK && isBABAR && isNUTEV && isMGMIN2){
							if( pow(10,ich2)*pow(10,im2) >  pow(10, bf_ch2)*pow(10,bf_um2)){

								bf_ch2 = ich2;
								bf_um2  = im2;
							}

							if( im2 > bf_ms_um4){

								bf_ms_um4 = im2;
								bf_ms_mz = mZ;
								bf_ms_chi2 = ich2;
							}

						}

					}
				}


				std::cout<<mS<<" "<<mZ<<" "<<bf_ch2<<" "<<bf_um2<<" "<<pow(pow(10, bf_ch2)*pow(10,bf_um2),2.0)<<std::endl;
			}//end mz


			std::cout<<"#ms: "<<mS<<" "<<bf_ms_um4<<" "<<bf_ms_mz<<" "<<bf_ms_chi2 <<" "<<pow(pow(10,bf_ms_um4),2.0)<<std::endl;
		}//end mS




	}//end boundmode



}//end main











/*****************************************************************************
 *
 *		What follows is the main bit of the code of actual producing events
 *
 *
 * ***************************************************************************/


