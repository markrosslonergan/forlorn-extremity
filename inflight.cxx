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

#include "fourmomentum.h" // Defines a class for fourmomenta; makes getting directions and 3-momenta (and
			  // all that) easier.

#include "sterile_flux.h" // Defines getEvents, and describes a class whose objects denote a single 
			  // incoming sterile (mass, fourmomentum).

#include "channel.h"	  // Includes the mother class for two body final state
			  // decays, and the derived classes for the two channels we've studied so far.

#include "detector.h"	  // This includes the detector-specific cut functions.


#include "zprimecrosssection.h"


#include "minInstance.h"


#define MPION  0.13957
#define MPI0   0.13497
#define MKAON  0.49367
#define MMU    0.10566
#define	ME     0.00051



#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1.h"




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
	

	std::string filename = "zmass_" + std::to_string(mZ) + "_smass_" + std::to_string(mS)+"_sim.root";

	
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
	
			if(i%5000==0) std::cout<<"#: "<<i<<std::endl;
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



		//	CHAN->decayfunction(nus);			//old approximate one
			CHAN->decayfunctionMassive(nus,ME,ME,0.0);
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
			
			CHAN->decayfunctionMassive(nusBar,ME,ME,0.0); //my new one, tobetested
	//		CHAN->decayfunction(nus);			//old approximate one
			CHAN->observables(&Obs, r);

						
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




			double weight_nu = mS/sqrt(sEn*sEn-mS*mS);  
			double weight_nubar = mS/sqrt(sEn_B*sEn_B-mS*mS);  
			if(weight_nu != weight_nu || weight_nubar != weight_nubar){std::cout<<"ERROR:Infinite weight, mS: "<<mS<<" En_nu "<<sEn<<" En_nubar " <<sEn_B<<" , run "<<j<<std::endl;}


			// Neutrino run
			if( e1En_smr < type1_low_thresh && e2En_smr >=type1_high_thresh ) // So Called Type 1 event
			{
				hType1_Evis_nu->Fill(e2En_smr, weight_nu);	
				hType1_Cos_nu->Fill(e2Cos_smr, weight_nu);
				hSignal_Evis_nu->Fill(e2En_smr, weight_nu);//and add to signal too	
				hSignal_Cos_nu->Fill(e2Cos_smr, weight_nu);	

				hType1_Ester_nu->Fill(sEn);
	
			}else if(fabs(dAngSep_smr) < type2_ang_thresh) // So called Type 2 event
			{
				hType2_Evis_nu->Fill(e2En_smr+e1En_smr, weight_nu);	
				hType2_Cos_nu->Fill(cos(dThSum_smr*M_PI/180.0), weight_nu);	
				hSignal_Evis_nu->Fill(e2En_smr+e1En_smr, weight_nu);	
				hSignal_Cos_nu->Fill(cos(dThSum_smr*M_PI/180.0), weight_nu);	


				hType2_Ester_nu->Fill(sEn);
			}


			// Anti-Neutrino run
			if( e1En_B_smr < type1_low_thresh && e2En_B_smr >=type1_high_thresh ) // So Called Type 1 event
			{
				hType1_Evis_nubar->Fill(e2En_B_smr, weight_nubar);	
				hType1_Cos_nubar->Fill(e2Cos_B_smr, weight_nubar);
				hSignal_Evis_nubar->Fill(e2En_B_smr, weight_nubar);//and add to signal too	
				hSignal_Cos_nubar->Fill(e2Cos_B_smr, weight_nubar);

				hType1_Ester_nubar->Fill(sEn_B);
	
			}else if(fabs(dAngSep_B_smr) < type2_ang_thresh) // So called Type 2 event
			{
				hType2_Evis_nubar->Fill(e2En_B_smr+e1En_B_smr, weight_nubar);	
				hType2_Cos_nubar->Fill(cos(dThSum_B_smr*M_PI/180.0), weight_nubar);	
				hSignal_Evis_nubar->Fill(e2En_B_smr+e1En_B_smr, weight_nubar);	
				hSignal_Cos_nubar->Fill(cos(dThSum_B_smr*M_PI/180.0), weight_nubar);	

				hType2_Ester_nubar->Fill(sEn_B);
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

		//Enter the anaDir subdirectly

		TDirectory *anaDir = fstat->GetDirectory("anaDir");
		if (anaDir) {
			anaDir->cd();
		} else {
			std::cout<<"ERROR: anaDir does not exit, cannot run stat mode, please run --mode ana first"<<std::endl;
			exit(EXIT_FAILURE);
		}


			TH1D * hSignal_Evis_nu = (TH1D*)anaDir->Get("hSignalSignal_Evis_nu");	
			TH1D * hSignal_Cos_nu = (TH1D*)anaDir->Get("hSignalSignal_Evis_nu");	
			TH1D * hSignal_Evis_nubar = (TH1D*)anaDir->Get("hSignalSignal_Cos_nubar");	
			TH1D * hSignal_Cos_nubar = (TH1D*)anaDir->Get("hSignalSignal_Cos_nubar");	


			std::vector<double> vSignal_Evis_nu;
			std::vector<double> vSignal_Cos_nu;
			std::vector<double> vSignal_Evis_nubar;
			std::vector<double> vSignal_Cos_nubar;

			for(int i=1;i<=EBINS;i++){ //Dont forget, root hists count from 1 as 0 is underflow.
				vSignal_Evis_nu.push_back( hSignal_Evis_nu->GetBinContent(i));				
				vSignal_Evis_nubar.push_back( hSignal_Evis_nubar->GetBinContent(i));				
			}
			for(int i=1;i<=COSBINS;i++){
				vSignal_Cos_nu.push_back( hSignal_Cos_nu->GetBinContent(i));				
				vSignal_Cos_nubar.push_back( hSignal_Cos_nubar->GetBinContent(i));				
			}

		

			//Right so returning to statDir, we basicall have all ingredients, 4 vSignal vectors, 4 Observed Vectors and 4 Expected background vector and 2 norms
			// Norm_nu and Norm_nubar
			
			minInstance statInstance(Norm_nu,Norm_nubar, vSignal_Evis_nu, vSignal_Cos_nu, vSignal_Evis_nubar, vSignal_Cos_nubar);

			std::cout<<statInstance.sig_C_nubar[2]<<" "<<statInstance.bkg_C_nu[3]<<std::endl;


return 0;
}
}//end main











/*****************************************************************************
 *
 *		What follows is the main bit of the code of actual producing events
 *
 *
 * ***************************************************************************/


