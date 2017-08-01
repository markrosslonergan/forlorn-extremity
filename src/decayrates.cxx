#include "decayrates.h"

double Gamma_EE(decay_params * params)
{
  	int comp, nregions, neval, fail;
	cubareal integral, error, prob;

  	double function = 1;

//  	double input_params[] = {p.mS, p.mZprime, p.chi, Ue4, Um4, Ut4};

  	Vegas(NDIM, NCOMP, Integrand_ee, params, NVEC,
    		EPSREL, EPSABS, VERBOSE, SEED,
    		MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    		GRIDNO, STATEFILE, SPIN,
    		&neval, &fail, &integral, &error, &prob);

return (1e-21)*(double)integral;
}

static int Integrand_ee(const int *ndim, const cubareal xx[], 
	const int *ncomp, cubareal ff[], void *userdata) {

double x=xx[0];
double y=xx[1];
double f = 0.0;

decay_params * params = (decay_params *)userdata;

double mS = (*params).mS;

if(mS<2*e_mass){ f=0.0;}
else
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
	double v_vev = 254.0;

	double tanchi = chi + pow(chi,3.0)/3.0 + (2.0/15.0)*pow(chi,5.0) + (17.0/315.0)*pow(chi,7.0) + + (62.0/2835.0)*pow(chi,9.0); //Just the tan(chi) but chi is small... works fine up to \chi=1.	
	double coschi = 1 - 0.5*pow(chi,2.0) + (1.0/24.0)*pow(chi,4.0) - (1.0/720.0)*pow(chi,6.0) + (1.0/40320.0)*pow(chi,8.0) - (1.0/3628800.0)*pow(chi,10.0); 	
	double sinbeta = sbeta_mathematica(mZprime/v_vev,chi)   ;	
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

//	double gVpgA = A/(t - mZprime*mZprime) + B/(t - Z_mass*Z_mass); 
//	double gVmgA = D/(t - mZprime*mZprime) + E/(t - Z_mass*Z_mass) + F/(u - W_mass*W_mass);
//	double gV = (gVpgA+gVmgA)/2.0;
//	double gA = (gVpgA-gVmgA)/2.0;
//
//	f  = 4.0*gVpgA*gVpgA*(t - e_mass*e_mass)*(mS*mS + e_mass*e_mass - t); 
//	f += 4.0*gVmgA*gVmgA*(u - e_mass*e_mass)*(mS*mS + e_mass*e_mass - u); 
//	f += 4.0*(gV*gV-gA*gA)*e_mass*e_mass*(t+u-2*e_mass*e_mass);
//	f *= 1e21*jacobian/(32.0*pow(2*M_PI*mS,3.0));
	
	double func_f = 4.0*(t - e_mass*e_mass)*(mS*mS + e_mass*e_mass - t); 
	double func_g = 4.0*(u - e_mass*e_mass)*(mS*mS + e_mass*e_mass - u); 
	double func_h = 4.0*e_mass*e_mass*(t+u-2*e_mass*e_mass);

	double f1 = ( A/(t - mZprime*mZprime) + B/(t - Z_mass*Z_mass)  )*((A*func_f + D*func_h)/(t - mZprime*mZprime) + (B*func_f + E*func_h)/(t - Z_mass*Z_mass));
	double f2 = pow(F/(u - W_mass*W_mass),2.0)*func_h;
	double f3 = ((A*func_h + 2.0*D*func_g)/(t - mZprime*mZprime) + (B*func_h + 2.0*E*func_g)/(t - Z_mass*Z_mass))*F/(u-W_mass*W_mass);

	f = (1-Ue4*Ue4-Um4*Um4-Ut4*Ut4)*(Ue4*Ue4*(f1+f3)+Um4*Um4*f1+Ut4*Ut4*f1) + (1-Ue4*Ue4)*Ue4*Ue4*f2;
	f *= 1e21*jacobian/(32.0*pow(2*M_PI*mS,3.0));

	ff[0] = f;
}

return 0;
}

double Gamma_NUPI0(decay_params * params)
{
	double mS = (*params).mS;
	double mZprime = (*params).mZprime;
	double chi = (*params).chi;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

	double G_nu_pi0 = 0.0; 

	if(mS>pi0_mass+1e-8)
	{ 
		double v_vev = 254.0;
		double Us4_sq = 1.0 - Ue4*Ue4 - Um4*Um4 - Ut4*Ut4;
		double tanchi = chi + pow(chi,3.0)/3.0 + (2.0/15.0)*pow(chi,5.0) + (17.0/315.0)*pow(chi,7.0) + + (62.0/2835.0)*pow(chi,9.0); //Just the tan(chi) but chi is small... works fine up to \chi=1.	
		double coschi = 1 - 0.5*pow(chi,2.0) + (1.0/24.0)*pow(chi,4.0) - (1.0/720.0)*pow(chi,6.0) + (1.0/40320.0)*pow(chi,8.0) - (1.0/3628800.0)*pow(chi,10.0); 	
		double sinbeta = sbeta_mathematica(mZprime/v_vev,chi)   ;	
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
	double mZprime = (*params).mZprime;
	double chi = (*params).chi;
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
	double mZprime = (*params).mZprime;
	double chi = (*params).chi;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

	double G = 0.0;

	if(mS>=mu_mass+e_mass+1e-6)
	{
		G = (GF*GF*pow(mS,5.0)/(192.0*pow(M_PI,3.0)))*(Ue4*Ue4*I1_C5(e_mass/mS,mu_mass/mS) + Um4*Um4*I1_C5(mu_mass/mS,e_mass/mS));
	}

return G;
}

double Gamma_EPI(decay_params * params)
{
	double mS = (*params).mS;
	double mZprime = (*params).mZprime;
	double chi = (*params).chi;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

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
	double mZprime = (*params).mZprime;
	double chi = (*params).chi;
	double Ue4 = (*params).Ue4;
	double Um4 = (*params).Um4;
	double Ut4 = (*params).Ut4;

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
	double mZprime = (*params).mZprime;
	double chi = (*params).chi;
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
	double mZprime = (*params).mZprime;
	double chi = (*params).chi;
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
double branch_plotter(double chi, double mZ, double Ue4, double Um4, double Ut4)
{
	decay_params params;
	params.mZprime=mZ;		
	params.mS=1.0;
	params.chi=chi;
	params.Ue4=1.0/sqrt(3);
	params.Um4=1.0/sqrt(3);
	params.Ut4=1.0/sqrt(3);

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

	for(log_mS=-5.0; log_mS<=log(0.490)/log(10)+1e-5; log_mS+=0.001)
	{
		params.mS=pow(10.0,log_mS);

		GT = Gamma_total(&params);
		G1 = Gamma_NUNUNU(&params); 
		G2 = Gamma_NUGAMMA(&params);
		G3 = Gamma_EE(&params);
		G4 = Gamma_NUMUE(&params);
		G5 = Gamma_NUPI0(&params);
		G6 = Gamma_EPI(&params);
		G7 = Gamma_NUMUMU(&params);
		G8 = Gamma_MUPI(&params);

		std::cout<<params.mS<<" "<<GT<<" "<<G1/GT<<" "<<G2/GT<<" "<<G3/GT<<" "<<G4/GT<<" "<<G5/GT<<" "<<G6/GT<<" "<<G7/GT<<" "<<G8/GT<<" "<<std::endl;
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

double X = pow(chi,9)*(-0.07705977726373278 - 0.5555429290052434*pow(mu,2) - 
      1.919982323225093*pow(mu,4) - 4.533147018464647*pow(mu,6) - 
      8.384273721015495*pow(mu,8) - 12.892135107956593*pow(mu,10)) + 
   pow(chi,7)*(-0.10092592904333923 - 0.6014363002828582*pow(mu,2) - 
      1.744315960766308*pow(mu,4) - 3.614361867827966*pow(mu,6) - 
      6.151350233131595*pow(mu,8) - 9.165526790764725*pow(mu,10)) + 
   pow(chi,5)*(-0.12516182698466563 - 0.5662431833043466*pow(mu,2) - 
      1.283554937941228*pow(mu,4) - 2.209156735509329*pow(mu,6) - 
      3.277656996254506*pow(mu,8) - 4.426212916054446*pow(mu,10)) + 
   pow(chi,3)*(-0.1529060225 - 0.4413224585*pow(mu,2) - 
      0.7112493080000001*pow(mu,4) - 0.9626865709999999*pow(mu,6) - 
      1.1956342474999997*pow(mu,8) - 1.4100923375*pow(mu,10)) + 
   chi*(-0.231 - 0.231*pow(mu,2) - 0.231*pow(mu,4) - 0.231*pow(mu,6) - 
      0.231*pow(mu,8) - 0.231*pow(mu,10));

return X;
}

double cbeta_mathematica(double mu, double chi)
{

if(mu>0.9){ std::cout<<"ERROR: sin(beta) function assumes Zprime mass below EW scale. Take more care."<<std::endl; return 1; }
if(chi>0.9){ std::cout<<"ERROR: sin(beta) function works perturbatively with chi (chi<<1). Take more care."<<std::endl; return 1; }

double X = 1 + pow(chi,10)*(-0.043727475807767556 - 0.37840408886083987*pow(mu,2) - 
      1.6350823810090083*pow(mu,4) - 4.899026019354672*pow(mu,6) - 
      11.65884220818534*pow(mu,8) - 23.65388866897464*pow(mu,10)) + 
   pow(chi,8)*(-0.04419701060848131 - 0.3174410748404108*pow(mu,2) - 
      1.1547778014994647*pow(mu,4) - 2.9875454720476946*pow(mu,6) - 
      6.289005646058018*pow(mu,8) - 11.524811327530113*pow(mu,10)) + 
   pow(chi,6)*(-0.04155439384633073 - 0.23279952233373588*pow(mu,2) - 
      0.6807500046241757*pow(mu,4) - 1.4729811791306302*pow(mu,6) - 
      2.678284344400904*pow(mu,8) - 4.3473220000014265*pow(mu,10)) + 
   pow(chi,4)*(-0.03567721573762501 - 0.13869047727150002*pow(mu,2) - 
      0.30512461466025004*pow(mu,4) - 0.5310644579625001*pow(mu,6) - 
      0.8125948372368751*pow(mu,8) - 1.1458005825420001*pow(mu,10)) + 
   pow(chi,2)*(-0.026680500000000003 - 0.053361000000000006*pow(mu,2) - 
      0.08004150000000002*pow(mu,4) - 0.10672200000000001*pow(mu,6) - 
      0.1334025*pow(mu,8) - 0.16008300000000003*pow(mu,10));

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


//##########################################################
//######    Crazy functions for threebody rate        ######
//##########################################################

double threebody_I(double alpha, double beta, double gamma)
{

double res = c0(alpha,beta) + gamma*c1(alpha,beta) + gamma*gamma*c2(alpha,beta) + gamma*gamma*gamma*c3(alpha,beta) + gamma*gamma*gamma*gamma*c4(alpha,beta) + gamma*gamma*gamma*gamma*gamma*c5(alpha,beta);

return res;
}

double c0(double alpha, double beta)
{

double X=0.08333333333333333 - (7*alpha)/120. + (2*pow(alpha,2))/45. - pow(alpha,3)/28. + 
   (5*pow(alpha,4))/168. - (11*pow(alpha,5))/432. + pow(alpha,6)/45. - 
   (13*pow(alpha,7))/660. + (7*pow(alpha,8))/396. - (5*pow(alpha,9))/312. + 
   (4*pow(alpha,10))/273. + (-0.058333333333333334 + (2*alpha)/45. - 
      pow(alpha,2)/28. + (5*pow(alpha,3))/168. - (11*pow(alpha,4))/432. + 
      pow(alpha,5)/45. - (13*pow(alpha,6))/660. + (7*pow(alpha,7))/396. - 
      (5*pow(alpha,8))/312. + (4*pow(alpha,9))/273. - (17*pow(alpha,10))/1260.)*
    beta + (0.044444444444444446 - alpha/28. + (5*pow(alpha,2))/168. - 
      (11*pow(alpha,3))/432. + pow(alpha,4)/45. - (13*pow(alpha,5))/660. + 
      (7*pow(alpha,6))/396. - (5*pow(alpha,7))/312. + (4*pow(alpha,8))/273. - 
      (17*pow(alpha,9))/1260. + pow(alpha,10)/80.)*pow(beta,2) + 
   (-0.03571428571428571 + (5*alpha)/168. - (11*pow(alpha,2))/432. + 
      pow(alpha,3)/45. - (13*pow(alpha,4))/660. + (7*pow(alpha,5))/396. - 
      (5*pow(alpha,6))/312. + (4*pow(alpha,7))/273. - (17*pow(alpha,8))/1260. + 
      pow(alpha,9)/80. - (19*pow(alpha,10))/1632.)*pow(beta,3) + 
   (0.02976190476190476 - (11*alpha)/432. + pow(alpha,2)/45. - 
      (13*pow(alpha,3))/660. + (7*pow(alpha,4))/396. - (5*pow(alpha,5))/312. + 
      (4*pow(alpha,6))/273. - (17*pow(alpha,7))/1260. + pow(alpha,8)/80. - 
      (19*pow(alpha,9))/1632. + (5*pow(alpha,10))/459.)*pow(beta,4) + 
   (-0.02546296296296296 + alpha/45. - (13*pow(alpha,2))/660. + 
      (7*pow(alpha,3))/396. - (5*pow(alpha,4))/312. + (4*pow(alpha,5))/273. - 
      (17*pow(alpha,6))/1260. + pow(alpha,7)/80. - (19*pow(alpha,8))/1632. + 
      (5*pow(alpha,9))/459. - (7*pow(alpha,10))/684.)*pow(beta,5) + 
   (0.022222222222222223 - (13*alpha)/660. + (7*pow(alpha,2))/396. - 
      (5*pow(alpha,3))/312. + (4*pow(alpha,4))/273. - (17*pow(alpha,5))/1260. + 
      pow(alpha,6)/80. - (19*pow(alpha,7))/1632. + (5*pow(alpha,8))/459. - 
      (7*pow(alpha,9))/684. + (11*pow(alpha,10))/1140.)*pow(beta,6) + 
   (-0.019696969696969695 + (7*alpha)/396. - (5*pow(alpha,2))/312. + 
      (4*pow(alpha,3))/273. - (17*pow(alpha,4))/1260. + pow(alpha,5)/80. - 
      (19*pow(alpha,6))/1632. + (5*pow(alpha,7))/459. - (7*pow(alpha,8))/684. + 
      (11*pow(alpha,9))/1140. - (23*pow(alpha,10))/2520.)*pow(beta,7) + 
   (0.017676767676767676 - (5*alpha)/312. + (4*pow(alpha,2))/273. - 
      (17*pow(alpha,3))/1260. + pow(alpha,4)/80. - (19*pow(alpha,5))/1632. + 
      (5*pow(alpha,6))/459. - (7*pow(alpha,7))/684. + (11*pow(alpha,8))/1140. - 
      (23*pow(alpha,9))/2520. + (2*pow(alpha,10))/231.)*pow(beta,8) + 
   (-0.016025641025641024 + (4*alpha)/273. - (17*pow(alpha,2))/1260. + 
      pow(alpha,3)/80. - (19*pow(alpha,4))/1632. + (5*pow(alpha,5))/459. - 
      (7*pow(alpha,6))/684. + (11*pow(alpha,7))/1140. - (23*pow(alpha,8))/2520. + 
      (2*pow(alpha,9))/231. - (25*pow(alpha,10))/3036.)*pow(beta,9) + 
   (0.014652014652014652 - (17*alpha)/1260. + pow(alpha,2)/80. - 
      (19*pow(alpha,3))/1632. + (5*pow(alpha,4))/459. - (7*pow(alpha,5))/684. + 
      (11*pow(alpha,6))/1140. - (23*pow(alpha,7))/2520. + (2*pow(alpha,8))/231. - 
      (25*pow(alpha,9))/3036. + (13*pow(alpha,10))/1656.)*pow(beta,10);

return X;
}

double c1(double alpha, double beta)
{
double X= -0.058333333333333334 + (7*alpha)/180. - pow(alpha,2)/35. + (5*pow(alpha,3))/224. - 
   (55*pow(alpha,4))/3024. + (11*pow(alpha,5))/720. - (13*pow(alpha,6))/990. + 
   (91*pow(alpha,7))/7920. - (35*pow(alpha,8))/3432. + (5*pow(alpha,9))/546. - 
   (34*pow(alpha,10))/4095. + (0.03888888888888889 - alpha/35. + 
      (5*pow(alpha,2))/224. - (55*pow(alpha,3))/3024. + (11*pow(alpha,4))/720. - 
      (13*pow(alpha,5))/990. + (91*pow(alpha,6))/7920. - (35*pow(alpha,7))/3432. + 
      (5*pow(alpha,8))/546. - (34*pow(alpha,9))/4095. + (17*pow(alpha,10))/2240.)*
    beta + (-0.02857142857142857 + (5*alpha)/224. - (55*pow(alpha,2))/3024. + 
      (11*pow(alpha,3))/720. - (13*pow(alpha,4))/990. + (91*pow(alpha,5))/7920. - 
      (35*pow(alpha,6))/3432. + (5*pow(alpha,7))/546. - (34*pow(alpha,8))/4095. + 
      (17*pow(alpha,9))/2240. - (19*pow(alpha,10))/2720.)*pow(beta,2) + 
   (0.022321428571428572 - (55*alpha)/3024. + (11*pow(alpha,2))/720. - 
      (13*pow(alpha,3))/990. + (91*pow(alpha,4))/7920. - (35*pow(alpha,5))/3432. + 
      (5*pow(alpha,6))/546. - (34*pow(alpha,7))/4095. + (17*pow(alpha,8))/2240. - 
      (19*pow(alpha,9))/2720. + (95*pow(alpha,10))/14688.)*pow(beta,3) + 
   (-0.018187830687830687 + (11*alpha)/720. - (13*pow(alpha,2))/990. + 
      (91*pow(alpha,3))/7920. - (35*pow(alpha,4))/3432. + (5*pow(alpha,5))/546. - 
      (34*pow(alpha,6))/4095. + (17*pow(alpha,7))/2240. - 
      (19*pow(alpha,8))/2720. + (95*pow(alpha,9))/14688. - (35*pow(alpha,10))/5814.
      )*pow(beta,4) + (0.015277777777777777 - (13*alpha)/990. + 
      (91*pow(alpha,2))/7920. - (35*pow(alpha,3))/3432. + (5*pow(alpha,4))/546. - 
      (34*pow(alpha,5))/4095. + (17*pow(alpha,6))/2240. - 
      (19*pow(alpha,7))/2720. + (95*pow(alpha,8))/14688. - 
      (35*pow(alpha,9))/5814. + (77*pow(alpha,10))/13680.)*pow(beta,5) + 
   (-0.013131313131313131 + (91*alpha)/7920. - (35*pow(alpha,2))/3432. + 
      (5*pow(alpha,3))/546. - (34*pow(alpha,4))/4095. + (17*pow(alpha,5))/2240. - 
      (19*pow(alpha,6))/2720. + (95*pow(alpha,7))/14688. - 
      (35*pow(alpha,8))/5814. + (77*pow(alpha,9))/13680. - 
      (253*pow(alpha,10))/47880.)*pow(beta,6) + 
   (0.01148989898989899 - (35*alpha)/3432. + (5*pow(alpha,2))/546. - 
      (34*pow(alpha,3))/4095. + (17*pow(alpha,4))/2240. - 
      (19*pow(alpha,5))/2720. + (95*pow(alpha,6))/14688. - 
      (35*pow(alpha,7))/5814. + (77*pow(alpha,8))/13680. - 
      (253*pow(alpha,9))/47880. + (23*pow(alpha,10))/4620.)*pow(beta,7) + 
   (-0.010198135198135198 + (5*alpha)/546. - (34*pow(alpha,2))/4095. + 
      (17*pow(alpha,3))/2240. - (19*pow(alpha,4))/2720. + 
      (95*pow(alpha,5))/14688. - (35*pow(alpha,6))/5814. + 
      (77*pow(alpha,7))/13680. - (253*pow(alpha,8))/47880. + 
      (23*pow(alpha,9))/4620. - (25*pow(alpha,10))/5313.)*pow(beta,8) + 
   (0.009157509157509158 - (34*alpha)/4095. + (17*pow(alpha,2))/2240. - 
      (19*pow(alpha,3))/2720. + (95*pow(alpha,4))/14688. - 
      (35*pow(alpha,5))/5814. + (77*pow(alpha,6))/13680. - 
      (253*pow(alpha,7))/47880. + (23*pow(alpha,8))/4620. - 
      (25*pow(alpha,9))/5313. + (325*pow(alpha,10))/72864.)*pow(beta,9) + 
   (-0.008302808302808303 + (17*alpha)/2240. - (19*pow(alpha,2))/2720. + 
      (95*pow(alpha,3))/14688. - (35*pow(alpha,4))/5814. + 
      (77*pow(alpha,5))/13680. - (253*pow(alpha,6))/47880. + 
      (23*pow(alpha,7))/4620. - (25*pow(alpha,8))/5313. + 
      (325*pow(alpha,9))/72864. - (39*pow(alpha,10))/9200.)*pow(beta,10);

return X;
}

double c2(double alpha, double beta)
{

double X= 0.044444444444444446 - alpha/35. + (103*pow(alpha,2))/5040. - 
   (473*pow(alpha,3))/30240. + (79*pow(alpha,4))/6300. - 
   (247*pow(alpha,5))/23760. + (7*pow(alpha,6))/792. - (263*pow(alpha,7))/34320. + 
   (304*pow(alpha,8))/45045. - (493*pow(alpha,9))/81900. + 
   (79*pow(alpha,10))/14560. + (-0.02857142857142857 + (103*alpha)/5040. - 
      (473*pow(alpha,2))/30240. + (79*pow(alpha,3))/6300. - 
      (247*pow(alpha,4))/23760. + (7*pow(alpha,5))/792. - 
      (263*pow(alpha,6))/34320. + (304*pow(alpha,7))/45045. - 
      (493*pow(alpha,8))/81900. + (79*pow(alpha,9))/14560. - 
      (1691*pow(alpha,10))/342720.)*beta + 
   (0.020436507936507937 - (473*alpha)/30240. + (79*pow(alpha,2))/6300. - 
      (247*pow(alpha,3))/23760. + (7*pow(alpha,4))/792. - 
      (263*pow(alpha,5))/34320. + (304*pow(alpha,6))/45045. - 
      (493*pow(alpha,7))/81900. + (79*pow(alpha,8))/14560. - 
      (1691*pow(alpha,9))/342720. + (83*pow(alpha,10))/18360.)*pow(beta,2) + 
   (-0.01564153439153439 + (79*alpha)/6300. - (247*pow(alpha,2))/23760. + 
      (7*pow(alpha,3))/792. - (263*pow(alpha,4))/34320. + 
      (304*pow(alpha,5))/45045. - (493*pow(alpha,6))/81900. + 
      (79*pow(alpha,7))/14560. - (1691*pow(alpha,8))/342720. + 
      (83*pow(alpha,9))/18360. - (1939*pow(alpha,10))/465120.)*pow(beta,3) + 
   (0.01253968253968254 - (247*alpha)/23760. + (7*pow(alpha,2))/792. - 
      (263*pow(alpha,3))/34320. + (304*pow(alpha,4))/45045. - 
      (493*pow(alpha,5))/81900. + (79*pow(alpha,6))/14560. - 
      (1691*pow(alpha,7))/342720. + (83*pow(alpha,8))/18360. - 
      (1939*pow(alpha,9))/465120. + (6743*pow(alpha,10))/1.7442e6)*pow(beta,4) + 
   (-0.010395622895622895 + (7*alpha)/792. - (263*pow(alpha,2))/34320. + 
      (304*pow(alpha,3))/45045. - (493*pow(alpha,4))/81900. + 
      (79*pow(alpha,5))/14560. - (1691*pow(alpha,6))/342720. + 
      (83*pow(alpha,7))/18360. - (1939*pow(alpha,8))/465120. + 
      (6743*pow(alpha,9))/1.7442e6 - (23*pow(alpha,10))/6384.)*pow(beta,5) + 
   (0.008838383838383838 - (263*alpha)/34320. + (304*pow(alpha,2))/45045. - 
      (493*pow(alpha,3))/81900. + (79*pow(alpha,4))/14560. - 
      (1691*pow(alpha,5))/342720. + (83*pow(alpha,6))/18360. - 
      (1939*pow(alpha,7))/465120. + (6743*pow(alpha,8))/1.7442e6 - 
      (23*pow(alpha,9))/6384. + (74*pow(alpha,10))/21945.)*pow(beta,6) + 
   (-0.0076631701631701635 + (304*alpha)/45045. - (493*pow(alpha,2))/81900. + 
      (79*pow(alpha,3))/14560. - (1691*pow(alpha,4))/342720. + 
      (83*pow(alpha,5))/18360. - (1939*pow(alpha,6))/465120. + 
      (6743*pow(alpha,7))/1.7442e6 - (23*pow(alpha,8))/6384. + 
      (74*pow(alpha,9))/21945. - (101*pow(alpha,10))/31878.)*pow(beta,7) + 
   (0.006748806748806749 - (493*alpha)/81900. + (79*pow(alpha,2))/14560. - 
      (1691*pow(alpha,3))/342720. + (83*pow(alpha,4))/18360. - 
      (1939*pow(alpha,5))/465120. + (6743*pow(alpha,6))/1.7442e6 - 
      (23*pow(alpha,7))/6384. + (74*pow(alpha,8))/21945. - 
      (101*pow(alpha,9))/31878. + (3809*pow(alpha,10))/1.27512e6)*pow(beta,8) + 
   (-0.006019536019536019 + (79*alpha)/14560. - (1691*pow(alpha,2))/342720. + 
      (83*pow(alpha,3))/18360. - (1939*pow(alpha,4))/465120. + 
      (6743*pow(alpha,5))/1.7442e6 - (23*pow(alpha,6))/6384. + 
      (74*pow(alpha,7))/21945. - (101*pow(alpha,8))/31878. + 
      (3809*pow(alpha,9))/1.27512e6 - (2859*pow(alpha,10))/1.012e6)*pow(beta,9) + 
   (0.005425824175824176 - (1691*alpha)/342720. + (83*pow(alpha,2))/18360. - 
      (1939*pow(alpha,3))/465120. + (6743*pow(alpha,4))/1.7442e6 - 
      (23*pow(alpha,5))/6384. + (74*pow(alpha,6))/21945. - 
      (101*pow(alpha,7))/31878. + (3809*pow(alpha,8))/1.27512e6 - 
      (2859*pow(alpha,9))/1.012e6 + (721*pow(alpha,10))/269100.)*pow(beta,10);

return X;
}


double c3(double alpha, double beta)
{
double X = -0.03571428571428571 + (5*alpha)/224. - (473*pow(alpha,2))/30240. + 
   (33*pow(alpha,3))/2800. - (2587*pow(alpha,4))/277200. + 
   (91*pow(alpha,5))/11880. - (133*pow(alpha,6))/20592. + 
   (167*pow(alpha,7))/30030. - (731*pow(alpha,8))/150150. + 
   (629*pow(alpha,9))/145600. - (1919*pow(alpha,10))/495040. + 
   (0.022321428571428572 - (473*alpha)/30240. + (33*pow(alpha,2))/2800. - 
      (2587*pow(alpha,3))/277200. + (91*pow(alpha,4))/11880. - 
      (133*pow(alpha,5))/20592. + (167*pow(alpha,6))/30030. - 
      (731*pow(alpha,7))/150150. + (629*pow(alpha,8))/145600. - 
      (1919*pow(alpha,9))/495040. + (361*pow(alpha,10))/102816.)*beta + 
   (-0.01564153439153439 + (33*alpha)/2800. - (2587*pow(alpha,2))/277200. + 
      (91*pow(alpha,3))/11880. - (133*pow(alpha,4))/20592. + 
      (167*pow(alpha,5))/30030. - (731*pow(alpha,6))/150150. + 
      (629*pow(alpha,7))/145600. - (1919*pow(alpha,8))/495040. + 
      (361*pow(alpha,9))/102816. - (497*pow(alpha,10))/155040.)*pow(beta,2) + 
   (0.011785714285714287 - (2587*alpha)/277200. + (91*pow(alpha,2))/11880. - 
      (133*pow(alpha,3))/20592. + (167*pow(alpha,4))/30030. - 
      (731*pow(alpha,5))/150150. + (629*pow(alpha,6))/145600. - 
      (1919*pow(alpha,7))/495040. + (361*pow(alpha,8))/102816. - 
      (497*pow(alpha,9))/155040. + (6853*pow(alpha,10))/2.3256e6)*pow(beta,3) + 
   (-0.009332611832611833 + (91*alpha)/11880. - (133*pow(alpha,2))/20592. + 
      (167*pow(alpha,3))/30030. - (731*pow(alpha,4))/150150. + 
      (629*pow(alpha,5))/145600. - (1919*pow(alpha,6))/495040. + 
      (361*pow(alpha,7))/102816. - (497*pow(alpha,8))/155040. + 
      (6853*pow(alpha,9))/2.3256e6 - (66539*pow(alpha,10))/2.44188e7)*pow(beta,4)\
    + (0.00765993265993266 - (133*alpha)/20592. + (167*pow(alpha,2))/30030. - 
      (731*pow(alpha,3))/150150. + (629*pow(alpha,4))/145600. - 
      (1919*pow(alpha,5))/495040. + (361*pow(alpha,6))/102816. - 
      (497*pow(alpha,7))/155040. + (6853*pow(alpha,8))/2.3256e6 - 
      (66539*pow(alpha,9))/2.44188e7 + (667*pow(alpha,10))/263340.)*pow(beta,5) + 
   (-0.006458818958818959 + (167*alpha)/30030. - (731*pow(alpha,2))/150150. + 
      (629*pow(alpha,3))/145600. - (1919*pow(alpha,4))/495040. + 
      (361*pow(alpha,5))/102816. - (497*pow(alpha,6))/155040. + 
      (6853*pow(alpha,7))/2.3256e6 - (66539*pow(alpha,8))/2.44188e7 + 
      (667*pow(alpha,9))/263340. - (955*pow(alpha,10))/403788.)*pow(beta,6) + 
   (0.0055611055611055615 - (731*alpha)/150150. + (629*pow(alpha,2))/145600. - 
      (1919*pow(alpha,3))/495040. + (361*pow(alpha,4))/102816. - 
      (497*pow(alpha,5))/155040. + (6853*pow(alpha,6))/2.3256e6 - 
      (66539*pow(alpha,7))/2.44188e7 + (667*pow(alpha,8))/263340. - 
      (955*pow(alpha,9))/403788. + (377*pow(alpha,10))/170016.)*pow(beta,7) + 
   (-0.004868464868464869 + (629*alpha)/145600. - (1919*pow(alpha,2))/495040. + 
      (361*pow(alpha,3))/102816. - (497*pow(alpha,4))/155040. + 
      (6853*pow(alpha,5))/2.3256e6 - (66539*pow(alpha,6))/2.44188e7 + 
      (667*pow(alpha,7))/263340. - (955*pow(alpha,8))/403788. + 
      (377*pow(alpha,9))/170016. - (14781*pow(alpha,10))/7.084e6)*pow(beta,8) + 
   (0.004320054945054945 - (1919*alpha)/495040. + (361*pow(alpha,2))/102816. - 
      (497*pow(alpha,3))/155040. + (6853*pow(alpha,4))/2.3256e6 - 
      (66539*pow(alpha,5))/2.44188e7 + (667*pow(alpha,6))/263340. - 
      (955*pow(alpha,7))/403788. + (377*pow(alpha,8))/170016. - 
      (14781*pow(alpha,9))/7.084e6 + (12957*pow(alpha,10))/6.578e6)*pow(beta,9) + 
   (-0.0038764544279250163 + (361*alpha)/102816. - (497*pow(alpha,2))/155040. + 
      (6853*pow(alpha,3))/2.3256e6 - (66539*pow(alpha,4))/2.44188e7 + 
      (667*pow(alpha,5))/263340. - (955*pow(alpha,6))/403788. + 
      (377*pow(alpha,7))/170016. - (14781*pow(alpha,8))/7.084e6 + 
      (12957*pow(alpha,9))/6.578e6 - (18067*pow(alpha,10))/9.6876e6)*pow(beta,10);

return X;
}

double c4(double alpha, double beta)
{
double X = 0.02976190476190476 - (55*alpha)/3024. + (79*pow(alpha,2))/6300. - 
   (2587*pow(alpha,3))/277200. + (3043*pow(alpha,4))/415800. - 
   (859*pow(alpha,5))/144144. + (1574*pow(alpha,6))/315315. - 
   (2567*pow(alpha,7))/600600. + (213*pow(alpha,8))/57200. - 
   (12217*pow(alpha,9))/3.7128e6 + (983*pow(alpha,10))/334152. + 
   (-0.018187830687830687 + (79*alpha)/6300. - (2587*pow(alpha,2))/277200. + 
      (3043*pow(alpha,3))/415800. - (859*pow(alpha,4))/144144. + 
      (1574*pow(alpha,5))/315315. - (2567*pow(alpha,6))/600600. + 
      (213*pow(alpha,7))/57200. - (12217*pow(alpha,8))/3.7128e6 + 
      (983*pow(alpha,9))/334152. - (8647*pow(alpha,10))/3.25584e6)*beta + 
   (0.01253968253968254 - (2587*alpha)/277200. + (3043*pow(alpha,2))/415800. - 
      (859*pow(alpha,3))/144144. + (1574*pow(alpha,4))/315315. - 
      (2567*pow(alpha,5))/600600. + (213*pow(alpha,6))/57200. - 
      (12217*pow(alpha,7))/3.7128e6 + (983*pow(alpha,8))/334152. - 
      (8647*pow(alpha,9))/3.25584e6 + (19679*pow(alpha,10))/8.1396e6)*pow(beta,2)\
    + (-0.009332611832611833 + (3043*alpha)/415800. - (859*pow(alpha,2))/144144. + 
      (1574*pow(alpha,3))/315315. - (2567*pow(alpha,4))/600600. + 
      (213*pow(alpha,5))/57200. - (12217*pow(alpha,6))/3.7128e6 + 
      (983*pow(alpha,7))/334152. - (8647*pow(alpha,8))/3.25584e6 + 
      (19679*pow(alpha,9))/8.1396e6 - (757781*pow(alpha,10))/3.418632e8)*
    pow(beta,3) + (0.007318422318422319 - (859*alpha)/144144. + 
      (1574*pow(alpha,2))/315315. - (2567*pow(alpha,3))/600600. + 
      (213*pow(alpha,4))/57200. - (12217*pow(alpha,5))/3.7128e6 + 
      (983*pow(alpha,6))/334152. - (8647*pow(alpha,7))/3.25584e6 + 
      (19679*pow(alpha,8))/8.1396e6 - (757781*pow(alpha,9))/3.418632e8 + 
      (11443*pow(alpha,10))/5.595975e6)*pow(beta,4) + 
   (-0.005959318459318459 + (1574*alpha)/315315. - (2567*pow(alpha,2))/600600. + 
      (213*pow(alpha,3))/57200. - (12217*pow(alpha,4))/3.7128e6 + 
      (983*pow(alpha,5))/334152. - (8647*pow(alpha,6))/3.25584e6 + 
      (19679*pow(alpha,7))/8.1396e6 - (757781*pow(alpha,8))/3.418632e8 + 
      (11443*pow(alpha,9))/5.595975e6 - (4595*pow(alpha,10))/2.422728e6)*
    pow(beta,5) + (0.004991833563262135 - (2567*alpha)/600600. + 
      (213*pow(alpha,2))/57200. - (12217*pow(alpha,3))/3.7128e6 + 
      (983*pow(alpha,4))/334152. - (8647*pow(alpha,5))/3.25584e6 + 
      (19679*pow(alpha,6))/8.1396e6 - (757781*pow(alpha,7))/3.418632e8 + 
      (11443*pow(alpha,8))/5.595975e6 - (4595*pow(alpha,9))/2.422728e6 + 
      (2379*pow(alpha,10))/1.34596e6)*pow(beta,6) + 
   (-0.004274059274059274 + (213*alpha)/57200. - (12217*pow(alpha,2))/3.7128e6 + 
      (983*pow(alpha,3))/334152. - (8647*pow(alpha,4))/3.25584e6 + 
      (19679*pow(alpha,5))/8.1396e6 - (757781*pow(alpha,6))/3.418632e8 + 
      (11443*pow(alpha,7))/5.595975e6 - (4595*pow(alpha,8))/2.422728e6 + 
      (2379*pow(alpha,9))/1.34596e6 - (837*pow(alpha,10))/506000.)*pow(beta,7) + 
   (0.003723776223776224 - (12217*alpha)/3.7128e6 + (983*pow(alpha,2))/334152. - 
      (8647*pow(alpha,3))/3.25584e6 + (19679*pow(alpha,4))/8.1396e6 - 
      (757781*pow(alpha,5))/3.418632e8 + (11443*pow(alpha,6))/5.595975e6 - 
      (4595*pow(alpha,7))/2.422728e6 + (2379*pow(alpha,8))/1.34596e6 - 
      (837*pow(alpha,9))/506000. + (53663*pow(alpha,10))/3.45345e7)*pow(beta,8) + 
   (-0.0032905085110967462 + (983*alpha)/334152. - (8647*pow(alpha,2))/3.25584e6 + 
      (19679*pow(alpha,3))/8.1396e6 - (757781*pow(alpha,4))/3.418632e8 + 
      (11443*pow(alpha,5))/5.595975e6 - (4595*pow(alpha,6))/2.422728e6 + 
      (2379*pow(alpha,7))/1.34596e6 - (837*pow(alpha,8))/506000. + 
      (53663*pow(alpha,9))/3.45345e7 - (5462759*pow(alpha,10))/3.729726e9)*
    pow(beta,9) + (0.00294177500059853 - (8647*alpha)/3.25584e6 + 
      (19679*pow(alpha,2))/8.1396e6 - (757781*pow(alpha,3))/3.418632e8 + 
      (11443*pow(alpha,4))/5.595975e6 - (4595*pow(alpha,5))/2.422728e6 + 
      (2379*pow(alpha,6))/1.34596e6 - (837*pow(alpha,7))/506000. + 
      (53663*pow(alpha,8))/3.45345e7 - (5462759*pow(alpha,9))/3.729726e9 + 
      (21911*pow(alpha,10))/1.582308e7)*pow(beta,10);

return X;
}

double c5(double alpha, double beta)
{
double X = -0.02546296296296296 + (11*alpha)/720. - (247*pow(alpha,2))/23760. + 
   (91*pow(alpha,3))/11880. - (859*pow(alpha,4))/144144. + 
   (7297*pow(alpha,5))/1.513512e6 - (2533*pow(alpha,6))/630630. + 
   (2193*pow(alpha,7))/640640. - (24263*pow(alpha,8))/8.16816e6 + 
   (437*pow(alpha,9))/167076. - (2819*pow(alpha,10))/1.209312e6 + 
   (0.015277777777777777 - (247*alpha)/23760. + (91*pow(alpha,2))/11880. - 
      (859*pow(alpha,3))/144144. + (7297*pow(alpha,4))/1.513512e6 - 
      (2533*pow(alpha,5))/630630. + (2193*pow(alpha,6))/640640. - 
      (24263*pow(alpha,7))/8.16816e6 + (437*pow(alpha,8))/167076. - 
      (2819*pow(alpha,9))/1.209312e6 + (40997*pow(alpha,10))/1.953504e7)*beta + 
   (-0.010395622895622895 + (91*alpha)/11880. - (859*pow(alpha,2))/144144. + 
      (7297*pow(alpha,3))/1.513512e6 - (2533*pow(alpha,4))/630630. + 
      (2193*pow(alpha,5))/640640. - (24263*pow(alpha,6))/8.16816e6 + 
      (437*pow(alpha,7))/167076. - (2819*pow(alpha,8))/1.209312e6 + 
      (40997*pow(alpha,9))/1.953504e7 - (26059*pow(alpha,10))/1.3674528e7)*
    pow(beta,2) + (0.00765993265993266 - (859*alpha)/144144. + 
      (7297*pow(alpha,2))/1.513512e6 - (2533*pow(alpha,3))/630630. + 
      (2193*pow(alpha,4))/640640. - (24263*pow(alpha,5))/8.16816e6 + 
      (437*pow(alpha,6))/167076. - (2819*pow(alpha,7))/1.209312e6 + 
      (40997*pow(alpha,8))/1.953504e7 - (26059*pow(alpha,9))/1.3674528e7 + 
      (437023*pow(alpha,10))/2.5069968e8)*pow(beta,3) + 
   (-0.005959318459318459 + (7297*alpha)/1.513512e6 - (2533*pow(alpha,2))/630630. + 
      (2193*pow(alpha,3))/640640. - (24263*pow(alpha,4))/8.16816e6 + 
      (437*pow(alpha,5))/167076. - (2819*pow(alpha,6))/1.209312e6 + 
      (40997*pow(alpha,7))/1.953504e7 - (26059*pow(alpha,8))/1.3674528e7 + 
      (437023*pow(alpha,9))/2.5069968e8 - (99145*pow(alpha,10))/6.1779564e7)*
    pow(beta,4) + (0.004821236964094107 - (2533*alpha)/630630. + 
      (2193*pow(alpha,2))/640640. - (24263*pow(alpha,3))/8.16816e6 + 
      (437*pow(alpha,4))/167076. - (2819*pow(alpha,5))/1.209312e6 + 
      (40997*pow(alpha,6))/1.953504e7 - (26059*pow(alpha,7))/1.3674528e7 + 
      (437023*pow(alpha,8))/2.5069968e8 - (99145*pow(alpha,9))/6.1779564e7 + 
      (28795*pow(alpha,10))/1.9381824e7)*pow(beta,5) + 
   (-0.004016618302332588 + (2193*alpha)/640640. - (24263*pow(alpha,2))/8.16816e6 + 
      (437*pow(alpha,3))/167076. - (2819*pow(alpha,4))/1.209312e6 + 
      (40997*pow(alpha,5))/1.953504e7 - (26059*pow(alpha,6))/1.3674528e7 + 
      (437023*pow(alpha,7))/2.5069968e8 - (99145*pow(alpha,8))/6.1779564e7 + 
      (28795*pow(alpha,9))/1.9381824e7 - (18603*pow(alpha,10))/1.34596e7)*
    pow(beta,6) + (0.0034231393606393608 - (24263*alpha)/8.16816e6 + 
      (437*pow(alpha,2))/167076. - (2819*pow(alpha,3))/1.209312e6 + 
      (40997*pow(alpha,4))/1.953504e7 - (26059*pow(alpha,5))/1.3674528e7 + 
      (437023*pow(alpha,6))/2.5069968e8 - (99145*pow(alpha,7))/6.1779564e7 + 
      (28795*pow(alpha,8))/1.9381824e7 - (18603*pow(alpha,9))/1.34596e7 + 
      (1699*pow(alpha,10))/1.3156e6)*pow(beta,7) + 
   (-0.0029704364263187792 + (437*alpha)/167076. - (2819*pow(alpha,2))/1.209312e6 + 
      (40997*pow(alpha,3))/1.953504e7 - (26059*pow(alpha,4))/1.3674528e7 + 
      (437023*pow(alpha,5))/2.5069968e8 - (99145*pow(alpha,6))/6.1779564e7 + 
      (28795*pow(alpha,7))/1.9381824e7 - (18603*pow(alpha,8))/1.34596e7 + 
      (1699*pow(alpha,9))/1.3156e6 - (903611*pow(alpha,10))/7.459452e8)*pow(beta,8)
     + (0.00261557614498791 - (2819*alpha)/1.209312e6 + 
      (40997*pow(alpha,2))/1.953504e7 - (26059*pow(alpha,3))/1.3674528e7 + 
      (437023*pow(alpha,4))/2.5069968e8 - (99145*pow(alpha,5))/6.1779564e7 + 
      (28795*pow(alpha,6))/1.9381824e7 - (18603*pow(alpha,7))/1.34596e7 + 
      (1699*pow(alpha,8))/1.3156e6 - (903611*pow(alpha,9))/7.459452e8 + 
      (396923*pow(alpha,10))/3.4810776e8)*pow(beta,9) + 
   (-0.0023310775052261122 + (40997*alpha)/1.953504e7 - 
      (26059*pow(alpha,2))/1.3674528e7 + (437023*pow(alpha,3))/2.5069968e8 - 
      (99145*pow(alpha,4))/6.1779564e7 + (28795*pow(alpha,5))/1.9381824e7 - 
      (18603*pow(alpha,6))/1.34596e7 + (1699*pow(alpha,7))/1.3156e6 - 
      (903611*pow(alpha,8))/7.459452e8 + (396923*pow(alpha,9))/3.4810776e8 - 
      (2964251*pow(alpha,10))/2.75321592e9)*pow(beta,10);

return X;
}
