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

#include "bounds.h"

bound::bound(std::string namein, double min, double L) :  bound_file(namein,min), length(L), mymin(min){
typical_E=1.0;

}

int bound::setTypicalEnergy(double in){

typical_E = in;
return 0;
}

double bound::meters2gev(double len){

return 100*pow(1.973,-1)*pow(10.0,5.0)*pow(10.0,9.0)*len;
}


double bound::gev2meters(double gev){

return gev/meters2gev(1);

}
std::vector<double> bound::lambertBounds(double A, double B, double mass, double gamma_total_prime){

	double UtildeSq = bound_file.getFlux(mass);	

	double gamma = typical_E/mass;
	double gamma_beta = gamma*sqrt(1.0-1.0/(gamma*gamma));

	double kap =  -UtildeSq*gamma_total_prime*meters2gev(length)/(2*gamma_beta);

	double upp = UtildeSq/(B*kap)*gsl_sf_lambert_W0(exp(kap)*B/sqrt(A)*kap  );
	double low = UtildeSq/(B*kap)*gsl_sf_lambert_Wm1(exp(kap)*B/sqrt(A)*kap );
	double simple = UtildeSq/sqrt(B);

	std::vector<double> temp;
	temp.push_back(low);
	temp.push_back(upp);
	temp.push_back(UtildeSq);
	temp.push_back(simple);

	return temp;

}

double bound::probDecay(double mass, double gam_c, double gam_t, double lambda){

	double gamma = typical_E/mass;
	double gamma_beta = gamma*sqrt(1.0-1.0/(gamma*gamma));
	double kap = meters2gev(length)/gamma_beta;
	double epsilon = lambda/length;

	return exp(-gam_t *kap)*(1-exp(-gam_t*kap*epsilon))*gam_c/gam_t;

}

bool bound::ps191(double mS, double mZ, double Um, double chi){
	double UtildeSq = 1;
	if(mS > mymin){
		UtildeSq = bound_file.getFlux(mS);	
	}

//	double mr = myRate(chi,mS,mZ) + Um*Um*assumedRate(mS);
//	double tr = UtildeSq*assumedRate(mS);
	double my_rate = Gamma_ZP_EE(chi,mS,mZ,1/sqrt(3),1/sqrt(3),1/sqrt(3)) + Gamma_EE(mS,0.0,Um,0.0);
	double their_rate = UtildeSq*Gamma_EE(mS,1/sqrt(3),1/sqrt(3),1/sqrt(3)); 
	double my_fullrate = Gamma_total(chi,mZ,mS,0.0,Um,0.0);
	double their_fullrate = Gamma_total(chi,-1.0,mS,0.0,sqrt(UtildeSq),0.0);

	//Going to multiply by U^2 each side as I think we lost "flux folded bit"
	double my_flux_folded_prob = Um*Um*probDecay(mS, my_rate, my_fullrate,16);
	double their_flux_folded_prob = UtildeSq*probDecay(mS, their_rate, their_fullrate ,16);

	bool ans = false;

	if(my_flux_folded_prob < their_flux_folded_prob){
		ans = true;		
	}

	return ans;

}

bool bound::asIs(double mass, double Us){
	double UtildeSq = 1;
	if(mass > mymin){
	 UtildeSq = bound_file.getFlux(mass);	
	}
	bool ans = false;

	//std::cout<<"PEAK: "<<UtildeSq<<" "<<Us<<std::endl;
	if(Us <= UtildeSq){ ans = true;}

	return ans;
}


double bound::myRate(double chi, double mS, double mZprime)
{
//This should be the Zprime decay into e e.

	double ret = Gamma_ZP_EE(chi, mZprime, mS, 1/sqrt(3),  1/sqrt(3),  1/sqrt(3)); 

return ret; 
}

double bound::assumedRate(double mS)
{
//This should be the total decay rate of the particle.

	return Gamma_total(1.0, -1.0, mS, 1/sqrt(3), 1/sqrt(3), 1/sqrt(3));
}

double bound::old_assumedRate(double mS)
{
//This should be the total decay rate of the particle.

	double gf = 0.00001166;
	return gf*gf*pow(mS,5.0)/(96.0*pow(3.14159,3.0) )*(0.126469);	

}

double bound::old_myRate(double chi, double mS, double mZprime)
{
//This should be the Zprime decay into e e.

	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1.0-mu_s*mu_s);
	double g1 =0.0874287 ;
	double fine_struct=1.0/137.0;

	double invnorm;

	if(alpha < 0.01)
	{
		double invnorm_perturb_0 = (1.0/(1.0+alpha))*(7.0/4.0 + 41.0/60.0*alpha); 

		double invnorm_perturb_rest = -0.18333*pow(alpha,2.0)+0.22857*pow(alpha,3.0)-0.23274*pow(alpha,4.0)+0.22421*pow(alpha,5.0)-0.21190*pow(alpha,6.0)+0.19899*pow(alpha,7.0)-0.18662*pow(alpha,8.0)+0.17517*pow(alpha,9.0)-0.16475*pow(alpha,10.0)+0.15531*pow(alpha,11.0);
		invnorm = invnorm_perturb_0+invnorm_perturb_rest;
	}
	else 
	{		
		invnorm = (3.0/(2.0*pow(alpha,3.0)))*(2+alpha-3*pow(alpha,2.0))/(1.0+alpha) + (4.0*alpha*alpha-3.0)*(log(1+alpha)/log(exp(1.0)))/pow(alpha,4.0);
	}

		double ret = chi*chi*pow(fine_struct,2.0)*pow(mS,5.0)*invnorm/(96*pow(3.14159,3.0)*g1*g1*pow(mZprime*mZprime-mS*mS,2.0));

return ret; 
}


