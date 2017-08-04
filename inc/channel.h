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
 */

#ifndef CHANNEL_H_
#define CHANNEL_H_

#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "fourmomentum.h" // defines class fourmomentum
#include "sterile_flux.h" // defines class initial_sterile
#include "decayrates.h" //needs decay_params

#define CHAN_ELECPOSI 0
#define CHAN_ELECPI 1
#define CHAN_MUONPI 2
#define CHAN_NUPI0 3
#define CHAN_GAMMA 4
#define CHAN_MUMU 5
#define CHAN_MUE 6

class initial_sterile;

typedef struct OBSERVABLES { // this is a struct of relevant observables for two (visible) body final states
	double E_sum; 	
	double Th_sum; 
	double AngSep; 
	double E_sterile; 
	double Th_sterile; 
	double E_high; 
	double Th_high; 
	double E_low; 
	double Th_low;
        std::vector<double > P_high;
	std::vector<double > P_low;
	double Minvar;
	double FS_AngSep; //The foreshortened angular separation.

	int chan_identifier;

	double E_high_smear;
	double E_low_smear;
	double Th_high_smear;
	double Th_low_smear;
	double AngSep_smear;
	double FS_AngSep_smear;
	double E_sum_smear;
	double Th_sum_smear;
	double E_sterile_smear;
	double Th_sterile_smear;
	double Minvar_smear;

	double Th_high_2;
	double Th_high_3;


	double Th_low_2;
	double Th_low_3;

	} OBSERVABLES;

class twoIP_channel { //This is the mother class for all decay channels (into two Ionising Particles)

public:
	twoIP_channel(gsl_rng * g, std::vector<double> input);

	fourmomentum IP1;	//first outgoing particle 4 momentum.
	fourmomentum IP2;	//second outgoing particle 4 momentum.
	int chan_identifier;
	gsl_rng * r;
	std::vector<double> model_params;

	int observables(OBSERVABLES * output, gsl_rng *g);
	virtual int decayfunction(initial_sterile nuS, decay_params *);

};

/* ###############################
   
   Below here we have a derived class for each channel

   ############################### */

/* ########################################################################

	This is the nu_s \to \nu e+ e- channel (off-shell Zprime).

   ######################################################################## */

class threebody : public twoIP_channel {

public:
	threebody(gsl_rng * g, std::vector<double> input);
	int decayfunction(initial_sterile nuS, decay_params* params);
	int computeLabFrameVariables(initial_sterile nuS, double p0[4], double p1[4]);
	int drawRestFrameDist(gsl_rng * r, decay_params* params, double out0[4], double out1[4]);
private:	
	std::vector<double > generic_boost(double Ep, double px, double py, double pz, double Er, double rx, double ry, double rz);

}; 

#endif
