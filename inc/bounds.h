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

#ifndef BOUNDS_H_
#define BOUNDS_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>

#include <gsl/gsl_sf_lambert.h>
#include "sterile_flux.h"
#include "decayrates.h"

class bound {

public:
	fluxfile bound_file;

	double mymin;	
	double length;
	double typical_E;

	double meters2gev(double l);
	bound(std::string,double min,double L);
	int setTypicalEnergy(double in);

	std::vector<double> lambertBounds(double A, double B,double,double);
	double probDecay(double mass,double gam_c, double gam_t, double lambda);
	double gev2meters(double gev);

	double assumedRate(double mS);
	double old_assumedRate(double mS);
	bool ps191(double mS, double mZ, double Um, double chi);
	double myRate(double chi, double mS, double mZprime);
	double old_myRate(double chi, double mS, double mZprime);
	bool asIs(double mass, double Us);


};
#endif
