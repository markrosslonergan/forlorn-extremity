#include "decayrates.h"

double Gamma_EE(decay_params * params)
{
	std::cout<<"Gamma_EE: starting"<<std::endl;
	int comp, nregions, neval, fail;
	cubareal integral, error, prob;

	double function = 1;
	params->numeval=0;
	Vegas(NDIM, NCOMP, Integrand_ee, params, NVEC,
			EPSREL, EPSABS, VERBOSE, SEED,MINEVAL, MAXEVAL, 					   NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE, SPIN,
			&neval, &fail, &integral, &error, &prob);

	//Note this factor of 1e-21 is removing the one I add to the integrand in Integrand_ee to help the integrator. 
	return (1e-21)*(double)integral;
}

double matrix_element_ee(double x, double y, decay_params * params)
{

	double mS = (*params).mS;

	double f = 0.0;

	if(mS>2*e_mass)
	{

		double mZprime = (*params).mZprime;
		double chi = (*params).chi;
		double Ue4 = (*params).Ue4;
		double Um4 = (*params).Um4;
		double Ut4 = (*params).Ut4;

		double alpha = mS*mS/(mZprime*mZprime-mS*mS);
		double beta = mS*mS/(Z_mass*Z_mass-mS*mS);
		double gamma = mS*mS/(W_mass*W_mass-mS*mS);

		//Cuba integrals go over 0,1.
		double u_min = e_mass*e_mass;
		double u_max = (mS-e_mass)*(mS-e_mass);
		double u = u_min + (u_max - u_min)*x; 

		double E_2 = (u + e_mass*e_mass)/(2.0*sqrt(u));
		double E_3 = (mS*mS - u - e_mass*e_mass)/(2.0*sqrt(u));
		double p2 = sqrt(E_2*E_2 - e_mass*e_mass);
		double p3 = sqrt(E_3*E_3 - e_mass*e_mass);

		double t_min = pow(E_2 + E_3,2.0) - (p2+p3)*(p2+p3);  
		double t_max = pow(E_2 + E_3,2.0) - (p2-p3)*(p2-p3);  
		double t = t_min + (t_max - t_min)*y;

		double s = mS*mS + 2.0*e_mass*e_mass - t - u;

		double jacobian = fabs((u_max-u_min)*(t_max-t_min));

		double sW = sqrt(xW); //I can't pretend the letter s doesn't exist.

		double tanchi = chi + pow(chi,3.0)/3.0 + (2.0/15.0)*pow(chi,5.0) + (17.0/315.0)*pow(chi,7.0) + + (62.0/2835.0)*pow(chi,9.0); //Just the tan(chi) but chi is small... works fine up to \chi=1.	
		double coschi = 1 - 0.5*pow(chi,2.0) + (1.0/24.0)*pow(chi,4.0) - (1.0/720.0)*pow(chi,6.0) + (1.0/40320.0)*pow(chi,8.0) - (1.0/3628800.0)*pow(chi,10.0); 	
		double sinbeta = sbeta_mathematica(mZprime/v_vev,chi);	
		double cosbeta = cbeta_mathematica(mZprime/v_vev,chi);	
		double Qs=1.0; //This is an assumption.
		double g=sqrt(4*M_PI/137.0)/sW; //The weak force coupling: g=e/sW.	
		double gX=g; //See QX.

		// Note that the d4i c4i and CC coefficients below do not have PMNS matrix elements included. I am doing that by hand to exploit unitarity relations.
		double d4i = -0.5*(cosbeta*sW*tanchi+sinbeta) + Qs*(gX/g)*sqrt(1-sW*sW)*cosbeta/coschi; 
		double c4i = 0.5*(-sinbeta*sW*tanchi+cosbeta) - Qs*(gX/g)*sqrt(1-sW*sW)*sinbeta/coschi; 

		double deV = 1.5*cosbeta*sW*tanchi+sinbeta*(-0.5+2.0*sW*sW);
		double deA = -0.5*(cosbeta*sW*tanchi+sinbeta);
		double ceV = -1.5*sinbeta*sW*tanchi+cosbeta*(-0.5+2.0*sW*sW);
		double ceA = -0.5*(-sinbeta*sW*tanchi+cosbeta);

		double A = (mZprime*mZprime-mS*mS)*sqrt(2)*pow(Z_mass/mS,2.0)*GF*2.0*alpha*d4i*(deV-deA); 
		double B = (Z_mass*Z_mass-mS*mS)*sqrt(2)*pow(Z_mass/mS,2.0)*GF*2.0*beta*c4i*(ceV-ceA); 
		double D = (mZprime*mZprime-mS*mS)*sqrt(2)*pow(Z_mass/mS,2.0)*GF*2.0*alpha*d4i*(deV+deA); 
		double E = (Z_mass*Z_mass-mS*mS)*sqrt(2)*pow(Z_mass/mS,2.0)*GF*2.0*beta*c4i*(ceV+ceA); 
		double F = -(W_mass*W_mass-mS*mS)*2.0*sqrt(2)*pow(Z_mass/mS,2.0)*GF*gamma*(1-sW*sW); 

		double func_f = 4.0*(t - e_mass*e_mass)*(mS*mS + e_mass*e_mass - t); 
		double func_g = 4.0*(u - e_mass*e_mass)*(mS*mS + e_mass*e_mass - u); 
		double func_h = 4.0*e_mass*e_mass*(t+u-2*e_mass*e_mass);

		double f1 = ( A/(t - mZprime*mZprime) + B/(t - Z_mass*Z_mass)  )*((A*func_f + D*func_h)/(t - mZprime*mZprime) + (B*func_f + E*func_h)/(t - Z_mass*Z_mass));
		double f2 = pow(F/(u - W_mass*W_mass),2.0)*func_h;
		double f3 = ((A*func_h + 2.0*D*func_g)/(t - mZprime*mZprime) + (B*func_h + 2.0*E*func_g)/(t - Z_mass*Z_mass))*F/(u-W_mass*W_mass);

		f = (1-Ue4*Ue4-Um4*Um4-Ut4*Ut4)*(Ue4*Ue4*(f1+f3)+Um4*Um4*f1+Ut4*Ut4*f1) + (1-Ue4*Ue4)*Ue4*Ue4*f2;
		f *= jacobian/(32.0*pow(2*M_PI*mS,3.0));

	}

	return f;
}

static int Integrand_ee(const int *ndim, const cubareal xx[], 
		const int *ncomp, cubareal ff[], void *userdata) {


	double x=xx[0];
	double y=xx[1];
	double f = 0.0;

	decay_params * params = (decay_params *)userdata;
	params->numeval++;
	//Note this factor of 1e21 should be removed again in Gamma_EE. Make sure it is!
	ff[0] = 1e21*matrix_element_ee(x,y,params);

	//std::cout<<"Integrand_ee: val "<<ff[0]<<" x "<<x<<" y "<<y<<" numeval "<<params->numeval<<std::endl;
	return 0;
}

double Gamma_NUPI0(decay_params * params)
{
	double mS = (*params).mS;
	double mZprime = (*params).mZprime;
	double chi = params->chi;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

	double G_nu_pi0 = 0.0; 

	if(mS>pi0_mass+1e-8)
	{ 
		double Us4_sq = 1.0 - Ue4*Ue4 - Um4*Um4 - Ut4*Ut4;
		double tanchi = chi + pow(chi,3.0)/3.0 + (2.0/15.0)*pow(chi,5.0) + (17.0/315.0)*pow(chi,7.0) + + (62.0/2835.0)*pow(chi,9.0); //Just the tan(chi) but chi is small... works fine up to \chi=1.	
		double coschi = 1 - 0.5*pow(chi,2.0) + (1.0/24.0)*pow(chi,4.0) - (1.0/720.0)*pow(chi,6.0) + (1.0/40320.0)*pow(chi,8.0) - (1.0/3628800.0)*pow(chi,10.0); 	
		double sinbeta = sbeta_mathematica(mZprime/v_vev,chi);	
		double cosbeta = cbeta_mathematica(mZprime/v_vev,chi);	
		double QX=1.0; //This is an assumption.
		double gX=1.0; //See QX.
		double sW = sqrt(xW); //Renamed for my own sanity. This is sin(thetaW).

		//what follows are the coupling constants for the Zprime and then the Z (_std).

		double cq_Z = sqrt(sqrt(2)*Z_mass*Z_mass*GF)*(cosbeta-sW*sinbeta*tanchi);
		double cq_Zp = sqrt(sqrt(2)*Z_mass*Z_mass*GF)*(sinbeta+sW*cosbeta*tanchi);

		double cnu_Z = -sqrt(sqrt(2)*Z_mass*Z_mass*GF)*(cosbeta - sW*sinbeta*tanchi) - QX*gX*sinbeta/coschi; 
		double cnu_Zp = -sqrt(sqrt(2)*Z_mass*Z_mass*GF)*(sinbeta + sW*cosbeta*tanchi) + QX*gX*cosbeta/coschi; 

		double F_sq = Us4_sq*(1.0-Us4_sq)*pow( cq_Zp*cnu_Zp/(2.0*(pi0_mass*pi0_mass-mZprime*mZprime)) + cq_Z*cnu_Z/(2.0*(pi0_mass*pi0_mass-Z_mass*Z_mass)),2.0);

		G_nu_pi0 =  (F_sq*fPi*fPi/(32.0*M_PI))*pow(mS,3.0)*pow(1.0 - pi0_mass*pi0_mass/(mS*mS), 2.0); 
	}

	return G_nu_pi0;
}


double Gamma_NUMUMU(decay_params * params)
{
	double mS = (*params).mS;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

	double G_nu_mu_mu = 0.0;

	if(mS>2.0*mu_mass)
	{
		double gRe = xW;
		double gRm = xW;
		double gLe = -0.5 + xW;
		double gLm = -0.5 + xW;

		G_nu_mu_mu = (pow(GF,2)*pow(mS,5)/(96.0*pow(M_PI,3)))*(Um4*Um4*((gLm*gRm +gRm)*I2_C7(mu_mass/mS) + (gLm*gLm + gRm*gRm + 1.0 + 2.0*gLm)*I1_C6(mu_mass/mS)) + (Ue4*Ue4 + Ut4*Ut4)*( gLm*gRm*I2_C7(mu_mass/mS) + (gLm*gLm + gRm*gRm)*I1_C6(mu_mass/mS) ));
	}

	return G_nu_mu_mu;
}


//Channels below here are unaltered (pure CC).
double Gamma_NUMUE(decay_params * params)
{
	//This channel is the sum of mu^+e^- and mu^-e^+

	double mS = (*params).mS;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;

	double G = 0.0;

	if(mS>=mu_mass+e_mass+1e-6)
	{
		G = 2.0*(GF*GF*pow(mS,5.0)/(192.0*pow(M_PI,3.0)))*(Ue4*Ue4*I1_C5(e_mass/mS,mu_mass/mS) + Um4*Um4*I1_C5(mu_mass/mS,e_mass/mS));
	}

	return G;
}

double Gamma_EPI(decay_params * params)
{
	double mS = (*params).mS;
	double Ue4 = (*params).Ue4;

	double G_e_pi = 0.0; 

	if(mS>pi_mass+e_mass+1e-8)
	{
		G_e_pi = 2.0*I1(e_mass*e_mass/(mS*mS), pi_mass*pi_mass/(mS*mS))*Ue4*Ue4*pow(Vud,2)*pow(GF,2)*pow(fPi,2)*pow(mS,3)/(16.0*M_PI);
	}

	return G_e_pi;
}

double Gamma_MUPI(decay_params * params)
{
	double mS = (*params).mS;
	double Um4 = (*params).Um4;

	double G_mu_pi = 0.0; 

	if(mS>mu_mass+pi_mass+1e-8)
	{ 
		G_mu_pi = 2.0*I1(mu_mass*mu_mass/(mS*mS), pi_mass*pi_mass/(mS*mS))*Um4*Um4*pow(Vud,2)*pow(GF,2)*pow(fPi,2)*pow(mS,3)/(16.0*M_PI);
	}

	return G_mu_pi;
}


double Gamma_NUGAMMA(decay_params * params)
{
	double mS = (*params).mS;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

	double G_nu_gamma = 0.0;

	G_nu_gamma = GF*GF/(96*pow(M_PI,3))*(Ue4*Ue4+Um4*Um4+Ut4*Ut4)*pow(mS,5)*(27/(32*M_PI*137));

	return G_nu_gamma;
}

double Gamma_NUNUNU(decay_params * params)
{
	double mS = (*params).mS;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

	/* taken from 0901.3589 */
	double G_3nu = (Ue4*Ue4+Um4*Um4+Ut4*Ut4)*pow(GF,2)*pow(mS,5)/(96.0*pow(M_PI,3)); 

	return G_3nu;
}

//This should report the total decay width in GeV.
double Gamma_total(decay_params * params)
{
	/* Total decay width */
	double Gamma =0.0;

	/* Now we add the relevant contributions according to the mass of nu_s */
	Gamma = Gamma_NUNUNU(params);

	//Add also radiative channel, why not, 
	Gamma += Gamma_NUGAMMA(params);

	double mS = (*params).mS;

	if(mS>=2.0*e_mass+1e-8)
	{
		Gamma += Gamma_EE(params); 

		if(mS>=mu_mass+e_mass)
		{
			Gamma += Gamma_NUMUE(params);

			if(mS>=pi0_mass)
			{
				Gamma += Gamma_NUPI0(params);

				if(mS>=pi_mass+e_mass)
				{
					Gamma += Gamma_EPI(params);

					if(mS>=2.0*mu_mass)
					{ 
						Gamma += Gamma_NUMUMU(params);

						if(mS>mu_mass + pi_mass)
						{
							Gamma += Gamma_MUPI(params);
						}
					}
				}
			}

		}
	}

	return Gamma;
}


//This can plot the branching ratios
double branch_plotter(decay_params * params)
{
	double GT = 0.0;
	double G1 = 0.0;
	double G2 = 0.0;
	double G3 = 0.0;
	double G4 = 0.0;
	double G5 = 0.0;
	double G6 = 0.0;
	double G7 = 0.0;
	double G8 = 0.0;

	double log_mS = 0.0;

	for(log_mS=-5.0; log_mS<=log(0.490)/log(10)+1e-5; log_mS+=0.01)
	{
		params->mS=pow(10.0,log_mS);

		GT = Gamma_total(params);
		G1 = Gamma_NUNUNU(params); 
		G2 = Gamma_NUGAMMA(params);
		G3 = Gamma_EE(params);
		G4 = Gamma_NUMUE(params);
		G5 = Gamma_NUPI0(params);
		G6 = Gamma_EPI(params);
		G7 = Gamma_NUMUMU(params);
		G8 = Gamma_MUPI(params);

		std::cout<<params->mS<<" "<<GT<<" "<<G1/GT<<" "<<G2/GT<<" "<<G3/GT<<" "<<G4/GT<<" "<<G5/GT<<" "<<G6/GT<<" "<<G7/GT<<" "<<G8/GT<<" "<<std::endl;
	}

	return 0;
}


//A bunch of very tedious approximations to make the total decay rate function.
//

double I1(double x, double y)
{
	return ((1+x-y)*(1+x) - 4.0*x)*sqrt(1.0 + x*x + y*y - 2.0*x - 2.0*y -2.0*x*y);
}

double I2(double x, double y)
{
	return ((1+x-y)*(1+x+2.0*y) - 4.0*x)*sqrt(1.0 + x*x + y*y - 2.0*x - 2.0*y -2.0*x*y);
}

double I3(double x, double y)
{
	return (1+2.0*y)*(1-y)*sqrt(1.0 + x*x + y*y - 2.0*x - 2.0*y -2.0*x*y);
}

double I1_C5(double x, double z)
{

	double X = 1 - 10*pow(z,2) + (33*pow(z,3))/64. + 
		pow(x,6)*(8 + 4*pow(z,2) + 8*pow(z,4) + 8*pow(z,6) + 8*pow(z,8) + 
				8*pow(z,10) + 8*pow(z,12) + 8*pow(z,14)) + 
		pow(x,8)*(-1 + 4*pow(z,4) + 10*pow(z,6) + 18*pow(z,8) + 28*pow(z,10) + 
				40*pow(z,12) + 54*pow(z,14)) + 
		pow(x,10)*((12*pow(z,4))/5. + (56*pow(z,6))/5. + 32*pow(z,8) + 
				72*pow(z,10) + 140*pow(z,12) + (1232*pow(z,14))/5.) + 
		pow(x,12)*((8*pow(z,4))/5. + 12*pow(z,6) + 50*pow(z,8) + 154*pow(z,10) + 
				392*pow(z,12) + (4368*pow(z,14))/5.) + 
		pow(x,14)*((8*pow(z,4))/7. + (88*pow(z,6))/7. + 72*pow(z,8) + 
				(1456*pow(z,10))/5. + (4704*pow(z,12))/5. + 2592*pow(z,14)) + 
		pow(x,4)*((33*pow(z,3))/64. + 18*pow(z,4) + (39*pow(z,5))/640. + 
				(2081*pow(z,6))/512. + (2211*pow(z,7))/35840. + (527*pow(z,8))/256. + 
				(1201*pow(z,9))/21504. + (10269*pow(z,10))/8192. + 
				(13955*pow(z,11))/270336. + (13927*pow(z,12))/16384. + 
				(114189*pow(z,13))/2.342912e6 + (284021*pow(z,14))/458752. + 
				(79847*pow(z,15))/1.703936e6 - 24*log(x) + pow(z,2)*(-30 + 24*log(x))) + 
		(pow(z,4)*(5029.8828125 + (-9361968 - 3440640*log(2) - 3440640*log(z))/2048.))/
		35. + (pow(z,5)*(16082.7734375 + 
					(11*(-9361968 - 3440640*log(2) - 3440640*log(z)))/4096. + 
					(18803624 + 18923520*log(2) + 18923520*log(z))/2048.))/35. + 
		(pow(z,6)*(40002.67578125 + (-28004480 - 49029120*log(2) - 49029120*log(z))/2048. + 
			   (33*(-9361968 - 3440640*log(2) - 3440640*log(z)))/4096. + 
			   (11*(18803624 + 18923520*log(2) + 18923520*log(z)))/4096.))/35. + 
		(pow(z,7)*(69947.8076171875 + (11*(-28004480 - 49029120*log(2) - 49029120*log(z)))/
			   4096. + (143*(-9361968 - 3440640*log(2) - 3440640*log(z)))/8192. + 
			   (33*(18803624 + 18923520*log(2) + 18923520*log(z)))/4096. + 
			   (35262640 + 80424960*log(2) + 80424960*log(z))/2048.))/35. + 
		(pow(z,8)*(107826.57470703125 + 
			   (-39153380 - 94617600*log(2) - 94617600*log(z))/2048. + 
			   (33*(-28004480 - 49029120*log(2) - 49029120*log(z)))/4096. + 
			   (1001*(-9361968 - 3440640*log(2) - 3440640*log(z)))/32768. + 
			   (143*(18803624 + 18923520*log(2) + 18923520*log(z)))/8192. + 
			   (11*(35262640 + 80424960*log(2) + 80424960*log(z)))/4096.))/35. + 
		(pow(z,9)*(146151.865234375 + (11*(-39153380 - 94617600*log(2) - 94617600*log(z)))/
			   4096. + (143*(-28004480 - 49029120*log(2) - 49029120*log(z)))/8192. + 
			   (3003*(-9361968 - 3440640*log(2) - 3440640*log(z)))/65536. + 
			   (1001*(18803624 + 18923520*log(2) + 18923520*log(z)))/32768. + 
			   (33*(35262640 + 80424960*log(2) + 80424960*log(z)))/4096. + 
			   (37511110 + 85155840*log(2) + 85155840*log(z))/2048.))/35. + 
		(pow(z,10)*(179559.2626953125 + 
			    (33*(-39153380 - 94617600*log(2) - 94617600*log(z)))/4096. + 
			    (-29813210 - 60318720*log(2) - 60318720*log(z))/2048. + 
			    (1001*(-28004480 - 49029120*log(2) - 49029120*log(z)))/32768. + 
			    (1001*(-9361968 - 3440640*log(2) - 3440640*log(z)))/16384. + 
			    (3003*(18803624 + 18923520*log(2) + 18923520*log(z)))/65536. + 
			    (143*(35262640 + 80424960*log(2) + 80424960*log(z)))/8192. + 
			    (11*(37511110 + 85155840*log(2) + 85155840*log(z)))/4096.))/35. + 
		(pow(z,11)*(203711.30920410156 + 
			    (143*(-39153380 - 94617600*log(2) - 94617600*log(z)))/8192. + 
			    (11*(-29813210 - 60318720*log(2) - 60318720*log(z)))/4096. + 
			    (3003*(-28004480 - 49029120*log(2) - 49029120*log(z)))/65536. + 
			    (2431*(-9361968 - 3440640*log(2) - 3440640*log(z)))/32768. + 
			    (1001*(18803624 + 18923520*log(2) + 18923520*log(z)))/16384. + 
			    (19128967 + 33707520*log(2) + 33707520*log(z))/2048. + 
			    (1001*(35262640 + 80424960*log(2) + 80424960*log(z)))/32768. + 
			    (33*(37511110 + 85155840*log(2) + 85155840*log(z)))/4096.))/35. + 
		(pow(z,12)*(216384.0756225586 + 
			    (1001*(-39153380 - 94617600*log(2) - 94617600*log(z)))/32768. + 
			    (33*(-29813210 - 60318720*log(2) - 60318720*log(z)))/4096. + 
			    (1001*(-28004480 - 49029120*log(2) - 49029120*log(z)))/16384. + 
			    (-9789956 - 14636160*log(2) - 14636160*log(z))/2048. + 
			    (21879*(-9361968 - 3440640*log(2) - 3440640*log(z)))/262144. + 
			    (2431*(18803624 + 18923520*log(2) + 18923520*log(z)))/32768. + 
			    (11*(19128967 + 33707520*log(2) + 33707520*log(z)))/4096. + 
			    (3003*(35262640 + 80424960*log(2) + 80424960*log(z)))/65536. + 
			    (143*(37511110 + 85155840*log(2) + 85155840*log(z)))/8192.))/35. + 
		(pow(z,13)*(217432.8793334961 + 
			    (3003*(-39153380 - 94617600*log(2) - 94617600*log(z)))/65536. + 
			    (143*(-29813210 - 60318720*log(2) - 60318720*log(z)))/8192. + 
			    (2431*(-28004480 - 49029120*log(2) - 49029120*log(z)))/32768. + 
			    (11*(-9789956 - 14636160*log(2) - 14636160*log(z)))/4096. + 
			    (46189*(-9361968 - 3440640*log(2) - 3440640*log(z)))/524288. + 
			    (3.9833766666666665e6 + 4804800*log(2) + 4804800*log(z))/2048. + 
			    (21879*(18803624 + 18923520*log(2) + 18923520*log(z)))/262144. + 
			    (33*(19128967 + 33707520*log(2) + 33707520*log(z)))/4096. + 
			    (1001*(35262640 + 80424960*log(2) + 80424960*log(z)))/16384. + 
			    (1001*(37511110 + 85155840*log(2) + 85155840*log(z)))/32768.))/35. + 
		(pow(z,14)*(208342.44216918945 + 
			    (1001*(-39153380 - 94617600*log(2) - 94617600*log(z)))/16384. + 
			    (1001*(-29813210 - 60318720*log(2) - 60318720*log(z)))/32768. + 
			    (21879*(-28004480 - 49029120*log(2) - 49029120*log(z)))/262144. + 
			    (33*(-9789956 - 14636160*log(2) - 14636160*log(z)))/4096. + 
			    (46189*(-9361968 - 3440640*log(2) - 3440640*log(z)))/524288. + 
			    (-1.2863026666666667e6 - 1145760*log(2) - 1145760*log(z))/2048. + 
			    (11*(3.9833766666666665e6 + 4804800*log(2) + 4804800*log(z)))/4096. + 
			    (46189*(18803624 + 18923520*log(2) + 18923520*log(z)))/524288. + 
			    (143*(19128967 + 33707520*log(2) + 33707520*log(z)))/8192. + 
			    (2431*(35262640 + 80424960*log(2) + 80424960*log(z)))/32768. + 
			    (3003*(37511110 + 85155840*log(2) + 85155840*log(z)))/65536.))/35. + 
		(pow(z,15)*(191567.53406524658 + 
			    (2431*(-39153380 - 94617600*log(2) - 94617600*log(z)))/32768. + 
			    (3003*(-29813210 - 60318720*log(2) - 60318720*log(z)))/65536. + 
			    (46189*(-28004480 - 49029120*log(2) - 49029120*log(z)))/524288. + 
			    (143*(-9789956 - 14636160*log(2) - 14636160*log(z)))/8192. + 
			    (88179*(-9361968 - 3440640*log(2) - 3440640*log(z)))/1.048576e6 + 
			    (11*(-1.2863026666666667e6 - 1145760*log(2) - 1145760*log(z)))/4096. + 
			    (326910.7272727273 + 186480*log(2) + 186480*log(z))/2048. + 
			    (33*(3.9833766666666665e6 + 4804800*log(2) + 4804800*log(z)))/4096. + 
			    (46189*(18803624 + 18923520*log(2) + 18923520*log(z)))/524288. + 
			    (1001*(19128967 + 33707520*log(2) + 33707520*log(z)))/32768. + 
			    (21879*(35262640 + 80424960*log(2) + 80424960*log(z)))/262144. + 
			    (1001*(37511110 + 85155840*log(2) + 85155840*log(z)))/16384.))/35. + 
		pow(x,2)*(-8 + 36*pow(z,2) - (33*pow(z,3))/32. + 
				(pow(z,4)*(-13428.515625 + (23551808 + 3440640*log(2) + 3440640*log(z))/2048.))/
				35. + (pow(z,5)*(-42175.546875 + 
						(-43250144 - 18923520*log(2) - 18923520*log(z))/2048. + 
						(11*(23551808 + 3440640*log(2) + 3440640*log(z)))/4096.))/35. + 
				(pow(z,6)*(-95414.1015625 + 
					   (11*(-43250144 - 18923520*log(2) - 18923520*log(z)))/4096. + 
					   (33*(23551808 + 3440640*log(2) + 3440640*log(z)))/4096. + 
					   (44401560 + 47308800*log(2) + 47308800*log(z))/2048.))/35. + 
				(pow(z,7)*(-167423.115234375 + 
					   (-29674820 - 70963200*log(2) - 70963200*log(z))/2048. + 
					   (33*(-43250144 - 18923520*log(2) - 18923520*log(z)))/4096. + 
					   (143*(23551808 + 3440640*log(2) + 3440640*log(z)))/8192. + 
					   (11*(44401560 + 47308800*log(2) + 47308800*log(z)))/4096.))/35. + 
				(pow(z,8)*(-249727.0166015625 + 
					   (11*(-29674820 - 70963200*log(2) - 70963200*log(z)))/4096. + 
					   (143*(-43250144 - 18923520*log(2) - 18923520*log(z)))/8192. + 
					   (1001*(23551808 + 3440640*log(2) + 3440640*log(z)))/32768. + 
					   (33*(44401560 + 47308800*log(2) + 47308800*log(z)))/4096. + 
					   (14759360 + 70963200*log(2) + 70963200*log(z))/2048.))/35. + 
				(pow(z,9)*(-329841.23046875 + 
					   (33*(-29674820 - 70963200*log(2) - 70963200*log(z)))/4096. + 
					   (-7735840 - 49674240*log(2) - 49674240*log(z))/2048. + 
					   (1001*(-43250144 - 18923520*log(2) - 18923520*log(z)))/32768. + 
					   (3003*(23551808 + 3440640*log(2) + 3440640*log(z)))/65536. + 
					   (143*(44401560 + 47308800*log(2) + 47308800*log(z)))/8192. + 
					   (11*(14759360 + 70963200*log(2) + 70963200*log(z)))/4096.))/35. + 
				(pow(z,10)*(-396343.212890625 + 
					    (143*(-29674820 - 70963200*log(2) - 70963200*log(z)))/8192. + 
					    (11*(-7735840 - 49674240*log(2) - 49674240*log(z)))/4096. + 
					    (3003*(-43250144 - 18923520*log(2) - 18923520*log(z)))/65536. + 
					    (1001*(23551808 + 3440640*log(2) + 3440640*log(z)))/16384. + 
					    (5566400 + 24837120*log(2) + 24837120*log(z))/2048. + 
					    (1001*(44401560 + 47308800*log(2) + 47308800*log(z)))/32768. + 
					    (33*(14759360 + 70963200*log(2) + 70963200*log(z)))/4096.))/35. + 
				(pow(z,11)*(-441263.2434082031 + 
					    (1001*(-29674820 - 70963200*log(2) - 70963200*log(z)))/32768. + 
					    (33*(-7735840 - 49674240*log(2) - 49674240*log(z)))/4096. + 
					    (1001*(-43250144 - 18923520*log(2) - 18923520*log(z)))/16384. + 
					    (-3920480 - 8870400*log(2) - 8870400*log(z))/2048. + 
					    (2431*(23551808 + 3440640*log(2) + 3440640*log(z)))/32768. + 
					    (11*(5566400 + 24837120*log(2) + 24837120*log(z)))/4096. + 
					    (3003*(44401560 + 47308800*log(2) + 47308800*log(z)))/65536. + 
					    (143*(14759360 + 70963200*log(2) + 70963200*log(z)))/8192.))/35. + 
				(pow(z,12)*(-461185.2117919922 + 
					    (3003*(-29674820 - 70963200*log(2) - 70963200*log(z)))/65536. + 
					    (143*(-7735840 - 49674240*log(2) - 49674240*log(z)))/8192. + 
					    (2431*(-43250144 - 18923520*log(2) - 18923520*log(z)))/32768. + 
					    (11*(-3920480 - 8870400*log(2) - 8870400*log(z)))/4096. + 
					    (2064840 + 2217600*log(2) + 2217600*log(z))/2048. + 
					    (21879*(23551808 + 3440640*log(2) + 3440640*log(z)))/262144. + 
					    (33*(5566400 + 24837120*log(2) + 24837120*log(z)))/4096. + 
					    (1001*(44401560 + 47308800*log(2) + 47308800*log(z)))/16384. + 
					    (1001*(14759360 + 70963200*log(2) + 70963200*log(z)))/32768.))/35. + 
				(pow(z,13)*(-456945.8367919922 + 
					    (1001*(-29674820 - 70963200*log(2) - 70963200*log(z)))/16384. + 
					    (1001*(-7735840 - 49674240*log(2) - 49674240*log(z)))/32768. + 
					    (21879*(-43250144 - 18923520*log(2) - 18923520*log(z)))/262144. + 
					    (33*(-3920480 - 8870400*log(2) - 8870400*log(z)))/4096. + 
					    (-762646.6666666666 - 369600*log(2) - 369600*log(z))/2048. + 
					    (11*(2064840 + 2217600*log(2) + 2217600*log(z)))/4096. + 
					    (46189*(23551808 + 3440640*log(2) + 3440640*log(z)))/524288. + 
					    (143*(5566400 + 24837120*log(2) + 24837120*log(z)))/8192. + 
					    (2431*(44401560 + 47308800*log(2) + 47308800*log(z)))/32768. + 
					    (3003*(14759360 + 70963200*log(2) + 70963200*log(z)))/65536.))/35. + 
				(pow(z,14)*(-432472.1401977539 + 
					    (2431*(-29674820 - 70963200*log(2) - 70963200*log(z)))/32768. + 
					    (3003*(-7735840 - 49674240*log(2) - 49674240*log(z)))/65536. + 
					    (46189*(-43250144 - 18923520*log(2) - 18923520*log(z)))/524288. + 
					    (143*(-3920480 - 8870400*log(2) - 8870400*log(z)))/8192. + 
					    (11*(-762646.6666666666 - 369600*log(2) - 369600*log(z)))/4096. + 
					    (192982.66666666666 + 36960*log(2) + 36960*log(z))/2048. + 
					    (33*(2064840 + 2217600*log(2) + 2217600*log(z)))/4096. + 
					    (46189*(23551808 + 3440640*log(2) + 3440640*log(z)))/524288. + 
					    (1001*(5566400 + 24837120*log(2) + 24837120*log(z)))/32768. + 
					    (21879*(44401560 + 47308800*log(2) + 47308800*log(z)))/262144. + 
					    (1001*(14759360 + 70963200*log(2) + 70963200*log(z)))/16384.))/35. + 
				(pow(z,15)*(-393336.06422424316 + 
					    (21879*(-29674820 - 70963200*log(2) - 70963200*log(z)))/262144. + 
					    (1001*(-7735840 - 49674240*log(2) - 49674240*log(z)))/16384. + 
					    (46189*(-43250144 - 18923520*log(2) - 18923520*log(z)))/524288. + 
					    (1001*(-3920480 - 8870400*log(2) - 8870400*log(z)))/32768. + 
					    (33*(-762646.6666666666 - 369600*log(2) - 369600*log(z)))/4096. + 
					    (-32037.39393939394 - 1680*log(2) - 1680*log(z))/2048. + 
					    (11*(192982.66666666666 + 36960*log(2) + 36960*log(z)))/4096. + 
					    (143*(2064840 + 2217600*log(2) + 2217600*log(z)))/8192. + 
					    (88179*(23551808 + 3440640*log(2) + 3440640*log(z)))/1.048576e6 + 
					    (3003*(5566400 + 24837120*log(2) + 24837120*log(z)))/65536. + 
					    (46189*(44401560 + 47308800*log(2) + 47308800*log(z)))/524288. + 
					    (2431*(14759360 + 70963200*log(2) + 70963200*log(z)))/32768.))/35.);

	return X;
}

double I1_C6(double y)
{

	double X = 1 - 18*pow(y,2) + (7*pow(y,3))/8. + 
		pow(y,4)*(137.375 + (-11504 - 6144*log(2) - 6144*log(y))/128. - 24*log(y)) + 
		pow(y,5)*(251.125 + (7*(-11504 - 6144*log(2) - 6144*log(y)))/256. + 
				(8600 + 21504*log(2) + 21504*log(y))/128.) + 
		pow(y,6)*(415.3125 + (-7334 - 29184*log(2) - 29184*log(y))/128. + 
				(7*(-11504 - 6144*log(2) - 6144*log(y)))/128. + 24*log(y) + 
				(7*(8600 + 21504*log(2) + 21504*log(y)))/256.) + 
		pow(y,7)*(559.734375 + (7*(-7334 - 29184*log(2) - 29184*log(y)))/256. + 
				(21*(-11504 - 6144*log(2) - 6144*log(y)))/256. + 
				(14645 + 16128*log(2) + 16128*log(y))/128. + 
				(7*(8600 + 21504*log(2) + 21504*log(y)))/128.) + 
		pow(y,8)*(660.76171875 + (7*(-7334 - 29184*log(2) - 29184*log(y)))/128. + 
				(105*(-11504 - 6144*log(2) - 6144*log(y)))/1024. + 
				(-21476 + 2688*log(2) + 2688*log(y))/128. + 
				(7*(14645 + 16128*log(2) + 16128*log(y)))/256. + 
				(21*(8600 + 21504*log(2) + 21504*log(y)))/256.) + 
		pow(y,9)*(656.8203125 + (21*(-7334 - 29184*log(2) - 29184*log(y)))/256. + 
				(18830 - 9408*log(2) - 9408*log(y))/128. + 
				(231*(-11504 - 6144*log(2) - 6144*log(y)))/2048. + 
				(7*(-21476 + 2688*log(2) + 2688*log(y)))/256. + 
				(7*(14645 + 16128*log(2) + 16128*log(y)))/128. + 
				(105*(8600 + 21504*log(2) + 21504*log(y)))/1024.) + 
		pow(y,10)*(633.71484375 + (105*(-7334 - 29184*log(2) - 29184*log(y)))/1024. + 
				(7*(18830 - 9408*log(2) - 9408*log(y)))/256. + 
				(231*(-11504 - 6144*log(2) - 6144*log(y)))/2048. + 
				(7*(-21476 + 2688*log(2) + 2688*log(y)))/128. + 
				(-10708 + 6048*log(2) + 6048*log(y))/128. + 
				(21*(14645 + 16128*log(2) + 16128*log(y)))/256. + 
				(231*(8600 + 21504*log(2) + 21504*log(y)))/2048.);

	return X;
}

double I2_C7(double y)
{

	double X = 8*pow(y,2) + 80*pow(y,10) + pow(y,6)*(63 - 96*log(4) - 192*log(y)) + 
		pow(y,4)*(-40.25 + 48*log(4) + 96*log(y)) + 
		pow(y,8)*(-140 + 96*log(4) + 192*log(y));

	return X;
}

//Below are a couple of Taylor expanded functions of the kinetic mixing courtesy of mathematica.

double sbeta_mathematica(double mu, double chi)
{

	if(mu>0.9){ std::cout<<"ERROR: sin(beta) function assumes Zprime mass below EW scale. Take more care."<<std::endl; return 1; }
	if(chi>0.9){ std::cout<<"ERROR: sin(beta) function works perturbatively with chi (chi<<1). Take more care."<<std::endl; return 1; }

	double X = (-0.48062459362791665 - 7.4496961006249086e-6*pow(mu,2) - 
			1.1547052050072706e-10*pow(mu,4) - 1.789796647354564e-15*pow(mu,6) - 
			2.7741903517802774e-20*pow(mu,8) - 4.30000364526672e-25*pow(mu,10))*chi + 
		(-0.10469605731194782 - 5.630728636055658e-6*pow(mu,2) - 
		 1.0938907308768873e-10*pow(mu,4) - 1.4181155435872652e-15*pow(mu,6) - 
		 8.06827027309427e-21*pow(mu,8) + 2.3958186976877815e-25*pow(mu,10))*pow(chi,3) + 
		(-0.01818861693908181 - 1.5777540649311225e-6*pow(mu,2) + 
		 4.777640898433737e-13*pow(mu,4) + 1.901589210319837e-15*pow(mu,6) + 
		 8.902204607878581e-20*pow(mu,8) + 2.787498612228834e-24*pow(mu,10))*pow(chi,5) + 
		(0.0005926669430206078 + 4.249705485402503e-7*pow(mu,2) + 
		 7.534085182368923e-11*pow(mu,4) + 4.383911746670846e-15*pow(mu,6) + 
		 1.594296315536986e-19*pow(mu,8) + 4.331358350173816e-24*pow(mu,10))*pow(chi,7) + 
		(0.0021838128835075226 + 6.873627996750816e-7*pow(mu,2) + 
		 7.450230375918738e-11*pow(mu,4) + 3.81925568196591e-15*pow(mu,6) + 
		 1.1715526294257591e-19*pow(mu,8) + 2.014105090559639e-24*pow(mu,10))*pow(chi,9);


	return X;
}

double cbeta_mathematica(double mu, double chi)
{

	if(mu>0.9){ std::cout<<"ERROR: sin(beta) function assumes Zprime mass below EW scale. Take more care."<<std::endl; return 1; }
	if(chi>0.9){ std::cout<<"ERROR: sin(beta) function works perturbatively with chi (chi<<1). Take more care."<<std::endl; return 1; }

	double X = 1 + (-0.1155 - 3.580507161014322e-6*pow(mu,2) - 8.324695798749895e-11*pow(mu,4) - 
			1.72044057262279e-15*pow(mu,6) - 3.333360276177208e-20*pow(mu,8) - 
			6.200062513814633e-25*pow(mu,10))*pow(chi,2) + 
		(-0.056989625 - 3.899769049538098e-6*pow(mu,2) - 1.2263664360858385e-10*pow(mu,4) - 
		 2.830841592203082e-15*pow(mu,6) - 5.3530988235175776e-20*pow(mu,8) - 
		 8.399018018714214e-25*pow(mu,10))*pow(chi,4) + 
		(-0.02080483052083337 - 2.1377972486153324e-6*pow(mu,2) - 
		 7.380146561367667e-11*pow(mu,4) - 1.2503943203549604e-15*pow(mu,6) + 
		 3.729648108917351e-21*pow(mu,8) + 1.0956438232376082e-24*pow(mu,10))*pow(chi,6)\
		+ (-0.005646292776809762 - 6.025871130433287e-7*pow(mu,2) - 
				3.882039773693532e-12*pow(mu,4) + 1.459802013048614e-15*pow(mu,6) + 
				1.0028529146242935e-19*pow(mu,8) + 4.136868620168733e-24*pow(mu,10))*pow(chi,8)\
		+ (-0.0008915751264421084 + 7.298323096629548e-8*pow(mu,2) + 
				3.2119852891006114e-11*pow(mu,4) + 2.8544914693644827e-15*pow(mu,6) + 
				1.4387988532822444e-19*pow(mu,8) + 5.077218790410028e-24*pow(mu,10))*pow(chi,10); 

	return X;
}


// ################################################################
// ################################################################
// ####                      Old rates                        #####
// ################################################################
// ################################################################
//
//Here are some archived functions. Please only use them for comparison.

double Gamma_NUPI0_old(double mS, double chi, double mZ, double Ue4, double Um4, double Ut4)
{
	double G_nu_pi0 = 0.0; 

	G_nu_pi0 =  (1.0-pi0_mass*pi0_mass/(mS*mS))*(1.0-pi0_mass*pi0_mass/(mS*mS))*(Ue4*Ue4+Um4*Um4+Ut4*Ut4)*pow(GF,2)*pow(fPi,2)*pow(mS,3)/(64.0*M_PI);

	if(mS<pi0_mass) { G_nu_pi0 = 0.0; }

	return G_nu_pi0;
}

double Gamma_NUMUMU_old(double mS, double chi, double mZ, double Ue4, double Um4, double Ut4)
{
	double G_nu_mu_mu = 0.0;

	if(mS>2.0*mu_mass)
	{
		double gRe = xW;
		double gRm = xW;
		double gLe = -0.5 + xW;
		double gLm = -0.5 + xW;

		G_nu_mu_mu = (pow(GF,2)*pow(mS,5)/(96.0*pow(M_PI,3)))*(Um4*Um4*((gLm*gRm +gRm)*I2_3arg_mu(mS) + (gLm*gLm + gRm*gRm + 1.0 + 2.0*gLm)*I1_3arg_mu(mS)) + (Ue4*Ue4 + Ut4*Ut4)*( gLm*gRm*I2_3arg_mu(mS) + (gLm*gLm + gRm*gRm)*I1_3arg_mu(mS)  ));
	}

	return G_nu_mu_mu;
}

double Gamma_NUMUE_old(double mS, double Ue4, double Um4, double Ut4)
{

	double G = 2.0*Um4*Um4*pow(GF,2)*pow(mS,5)/(192.0*pow(M_PI,3))*I1_3arg_mue(mS);

	return G;

}

double Gamma_EE_old(double mS, double Ue4, double Um4, double Ut4)
{
	//The standard model Z + W mediators only decay from Pascoli/Han etc.
	double G_nu_e_e = 0.0; 

	if(mS>2.0*e_mass+1e-8)
	{ 
		double gRe = xW;
		double gRm = xW;
		double gLe = -0.5 + xW;
		double gLm = -0.5 + xW;

		G_nu_e_e = (pow(GF,2)*pow(mS,5)/(96.0*pow(M_PI,3)))*(Ue4*Ue4*((gLe*gRe +gRe)*I2_C7(e_mass/mS) + (gLe*gLe + gRe*gRe + 1.0 + 2.0*gLe)*I1_C6(e_mass/mS)) + (Um4*Um4 + Ut4*Ut4)*( gLe*gRe*I2_C7(e_mass/mS) + (gLe*gLe + gRe*gRe)*I1_C6(e_mass/mS)  ));
	}

	return G_nu_e_e;
}

double Gamma_EE_oldest(double mS, double Ue4, double Um4, double Ut4)
{
	double G_nu_e_e = 0.0; 

	double gRe = xW;
	double gRm = xW;
	double gLe = -0.5 + xW;
	double gLm = -0.5 + xW;

	G_nu_e_e = pow(GF,2)*pow(mS,5)*(pow(Ue4,2)+pow(Um4,2))/(192.0*pow(M_PI,3));

	if(mS<2.0*e_mass){ G_nu_e_e = 0.0; }

	return G_nu_e_e;
}

double I1_3arg_mue(double ms)
{

	return -0.237501*((ms-0.107)/(0.5-0.107)) + 4.59382*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107)) - 6.21622*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107)) +2.58557*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107));
}

double I1_3arg_mu(double ms)
{
	return 11.0*(ms-0.212)*(ms-0.212) - 80.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212) + 180.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212);
}

double I1_3arg_e(double ms)
{
	if(ms<0.001022)
	{
		return 1e-100;
	}
	else
	{
		double temp = pow(1.0-(0.001022/ms)*(0.001022/ms),1.0)*(1.0-(0.00222/ms)*(0.00222/ms));
		if(temp<0.0) {return 0.0;}
		else {return temp;} 
	}
}

double I2_3arg_mu(double ms)
{
	return 7.0*(ms-0.212)*(ms-0.212) - 130.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212) + 30.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212) + 670.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212);
}

double I2_3arg_e(double ms)
{
	if(ms<0.001022)
	{
		return 0.0;
	}
	else
	{
		return 5200.0*(ms/0.04)*((ms - 0.001022)/(0.5-0.001022))*exp( -230.0*pow((ms-0.001022)/(0.5-0.001022),0.8)) + 0.045*(1.0-0.001022/ms)*exp(-ms/0.01); 
	}
}


