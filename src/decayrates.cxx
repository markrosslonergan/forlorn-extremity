#include "decayrates.h"

double Gamma_full_EE(double chi, double mZ, double mS, double Ue4, double Um4, double Ut4)
{
//This is the ee decay with both Z and Zp mediators. It should interpolate between the two cases already coded up, but also lead to an interference enhancement.

//It's not working yet.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	double G_Zp=0.0;

	if(mS<2.0*e_mass){ G_Zp=0.0;}

return 0.0;
}

double Gamma_ZP_EE(double chi, double mZ, double mS, double Ue4, double Um4, double Ut4)
{
	double alpha = mS*mS/(mZ*mZ-mS*mS);
	double gtilde = chi*sqrt(Ue4*Ue4+Um4*Um4+Ut4*Ut4)*sqrt(4.0*M_PI/137.0)*xW/(tan(2.0*asin(xW)));
	double N = 1.0;

	if(alpha<0.01)
	{  
		N = -1.0 - alpha/(1+alpha);
		int n=0;	
		for(n=0;n<11;n++)
		{
			N += 6.0*(1+alpha)*pow(-alpha, (double)n)/((double)n+4.0);
		}
	}
	else
	{
		N = (6.0*alpha + 9*alpha*alpha + 2.0*pow(alpha,3.0) - 6.0*pow((1+alpha),2.0)*log(1+alpha)/log(exp(1.0)))/(pow(alpha,4.0)*(1+alpha)) ;
	} 

	double G_Zp = (pow(gtilde,2.0)/pow(4.0*M_PI,3.0))*(alpha*alpha*mS/3.0)*N; 

	if(mS<2.0*e_mass){ G_Zp=0.0;}

return G_Zp;
}


double Gamma_EE(double mS, double Ue4, double Um4, double Ut4)
{
	double G_nu_e_e = 0.0; 

	double gRe = xW;
	double gRm = xW;
	double gLe = -0.5 + xW;
	double gLm = -0.5 + xW;

//	G_nu_e_e = (pow(GF,2)*pow(mS,5)/(96.0*pow(M_PI,3)))*(Ue4*Ue4*((gLe*gRe +gRe)*I2_3arg_e(mS) + (gLe*gLe + gRe*gRe + 1.0 + 2.0*gLe)*I1_3arg_e(mS)) + (Um4*Um4 + Ut4*Ut4)*( gLe*gRe*I2_3arg_e(mS) + (gLe*gLe + gRe*gRe)*I1_3arg_e(mS)  ));

//	std::cout<<"old "<<G_nu_e_e<<std::endl;

	G_nu_e_e = pow(GF,2)*pow(mS,5)*(pow(Ue4,2)+pow(Um4,2))/(192.0*pow(M_PI,3));

//	std::cout<<"theirs "<<G_nu_e_e<<std::endl;

//	std::cout<<"other "<< pow(GF,2)*pow(mS,5)*(pow(Ue4,2)+pow(Um4,2))/(192.0*pow(M_PI,3))*(1/4.0-xW+2*xW*xW)<<" "<<(1/4.0-xW+2*xW*xW)<<std::endl;

	if(mS<2.0*e_mass){ G_nu_e_e = 0.0; }

return G_nu_e_e;
}

double Gamma_NUMUMU(double mS, double Ue4, double Um4, double Ut4)
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
	
double Gamma_EPI(double mS, double Ue4, double Um4, double Ut4)
{
	double G_e_pi = 0.0; 

	G_e_pi = 2.0*I1(e_mass*e_mass/(mS*mS), pi_mass*pi_mass/(mS*mS))*Ue4*Ue4*pow(Vud,2)*pow(GF,2)*pow(fPi,2)*pow(mS,3)/(16.0*M_PI);

	if(mS<pi_mass+e_mass)
	{
		G_e_pi=0.0;
	}

return G_e_pi;
}

double Gamma_MUPI(double mS, double Ue4, double Um4, double Ut4)
{
	double G_mu_pi = 0.0; 


	G_mu_pi = 2.0*I1(mu_mass*mu_mass/(mS*mS), pi_mass*pi_mass/(mS*mS))*Um4*Um4*pow(Vud,2)*pow(GF,2)*pow(fPi,2)*pow(mS,3)/(16.0*M_PI);

	if(mS<mu_mass+pi_mass){ G_mu_pi = 0.0; }

return G_mu_pi;
}


double Gamma_NUPI0(double mS, double Ue4, double Um4, double Ut4)
{
	double G_nu_pi0 = 0.0; 

	G_nu_pi0 =  (1.0-pi0_mass*pi0_mass/(mS*mS))*(1.0-pi0_mass*pi0_mass/(mS*mS))*(Ue4*Ue4+Um4*Um4+Ut4*Ut4)*pow(GF,2)*pow(fPi,2)*pow(mS,3)/(64.0*M_PI);

	if(mS<pi0_mass) { G_nu_pi0 = 0.0; }

return G_nu_pi0;
}

double Gamma_NUGAMMA(double ms,double Ue4, double Um4 ,double Ut4){
	
	double G_nu_gamma = 0.0;
		
	G_nu_gamma = GF*GF/(96*pow(M_PI,3))*(Ue4*Ue4+Um4*Um4+Ut4*Ut4)*pow(ms,5)*(27/(32*M_PI*137));

return G_nu_gamma;
}

double Gamma_NUNUNU(double mS, double Ue4, double Um4, double Ut4)
{
	/* taken from 0901.3589 */
	double G_3nu = (Ue4*Ue4+Um4*Um4+Ut4*Ut4)*pow(GF,2)*pow(mS,5)/(96.0*pow(M_PI,3)); 

return G_3nu;
}


//This should report the total decay width in GeV.
double Gamma_total(double chi, double mZ, double mS, double Ue4, double Um4, double Ut4)
{

	/* Total decay width over |U|^2 assuming only coupling to muons */
	double Gamma_without_Zp=0.0;
	double Gamma_with_Zp=0.0;
	double G_3nu = 0.0;
	double G_zp_e_e = 0.0; /* Zprime driven e^+ e^- */
	double G_nu_e_e = 0.0; 
	double G_mu_e_nue = 0.0; /* mu^- e^+ nu_e *AND* mu^+ e^- bar(nu_e) */ 
	double G_e_mu_numu = 0.0; /* e^- mu^+ nu_mu *AND* e^+ mu^- bar(nu_mu) */ 
	double G_nu_pi0 = 0.0; /* nu_e and nu_mu and nu_tau */
	double G_e_pi = 0.0; /* e^+ pi^- *AND* e^- pi^+ */
	double G_nu_mu_mu = 0.0;
	double G_mu_pi = 0.0; /* mu^+ pi^- *AND* mu^- pi^+ */

	double gRe = xW;
	double gRm = xW;
	double gLe = -0.5 + xW;
	double gLm = -0.5 + xW;

	G_3nu = Gamma_NUNUNU(mS, Ue4, Um4, Ut4); 

	/* Now we add the relevant contributions according to the mass of nu_s */
	Gamma_without_Zp = G_3nu;
	Gamma_with_Zp = G_3nu;
	
	//Add also radiative channel, why not, 
	Gamma_without_Zp += Gamma_NUGAMMA(mS,Ue4,Um4,Ut4);
	Gamma_with_Zp += Gamma_NUGAMMA(mS,Ue4,Um4,Ut4);

if(mS>=2.0*e_mass+1e-3)
{
	G_nu_e_e = Gamma_EE(mS,Ue4,Um4,Ut4);

	if(mZ>0.0)// If mZ is negative, the rate is compute for pure nuSM.
	{ 
		//WARNING, well not necessarily warning but mark changed this on 15th jul
		// OLD// G_zp_e_e = Gamma_ZP_EE(chi,mZ,mS,Ue4,Um4,Ut4);
		G_zp_e_e = Gamma_ZP_EE(chi,mZ,mS,0.0,0.0,1.0);
	}
	else 
	{
 		G_zp_e_e = 0.0; 
	}

	Gamma_with_Zp += G_zp_e_e + G_nu_e_e; /* all relevant for ms >1.02 MeV */
	Gamma_without_Zp += G_nu_e_e; /* all relevant for ms >1.02 MeV */

	if(mS>=mu_mass+e_mass)
	{
		G_mu_e_nue = 2.0*Um4*Um4*pow(GF,2)*pow(mS,5)/(192.0*pow(M_PI,3))*I1_3arg_mue(mS);//notexact.
		G_e_mu_numu = 2.0*Ue4*Ue4*pow(GF,2)*pow(mS,5)/(192.0*pow(M_PI,3))*I1_3arg_mue(mS);//not quite exact

		Gamma_with_Zp += G_mu_e_nue + G_e_mu_numu;
		Gamma_without_Zp += G_mu_e_nue + G_e_mu_numu;

		if(mS>=pi0_mass)
		{
			G_nu_pi0 = Gamma_NUPI0(mS,Ue4,Um4,Ut4);
			Gamma_with_Zp += Gamma_NUPI0(mS,Ue4,Um4,Ut4);
			Gamma_without_Zp += Gamma_NUPI0(mS,Ue4,Um4,Ut4);

			if(mS>=pi_mass+e_mass)
			{
				Gamma_without_Zp += Gamma_EPI(mS,Ue4,Um4,Ut4);
				Gamma_with_Zp += Gamma_EPI(mS,Ue4,Um4,Ut4);
				G_e_pi += Gamma_EPI(mS,Ue4,Um4,Ut4);

				if(mS>=2.0*mu_mass)
				{ 
					Gamma_without_Zp += Gamma_NUMUMU(mS,Ue4,Um4,Ut4);
					Gamma_with_Zp += Gamma_NUMUMU(mS,Ue4,Um4,Ut4);
					G_nu_mu_mu += Gamma_NUMUMU(mS,Ue4,Um4,Ut4);

					if(mS>mu_mass + pi_mass)
					{
					      Gamma_without_Zp += Gamma_MUPI(mS,Ue4,Um4,Ut4);
					      Gamma_with_Zp += Gamma_MUPI(mS,Ue4,Um4,Ut4);
					      G_mu_pi += Gamma_MUPI(mS,Ue4,Um4,Ut4);
					}
				}
			}
		}
	
	}
}
//	std::cout<<log(mS)/log(10)<<" "<<Gamma_with_Zp<<" "<<Gamma_without_Zp<<" "<<G_3nu<<" "<<G_nu_e_e<<" "<<G_zp_e_e<<" "<<G_mu_e_nue+G_e_mu_numu<<" "<<G_nu_pi0<<" "<<G_e_pi<<" "<<G_nu_mu_mu<<" "<<G_mu_pi<<std::endl;

return Gamma_with_Zp;
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


double I1_3arg_mue(double ms)
{

return -0.237501*((ms-0.107)/(0.5-0.107)) + 4.59382*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107)) - 6.21622*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107)) +2.58557*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107))*((ms-0.107)/(0.5-0.107));
}

double I3(double x, double y)
{
return (1+2.0*y)*(1-y)*sqrt(1.0 + x*x + y*y - 2.0*x - 2.0*y -2.0*x*y);
}

double I2_3arg_mu(double ms)
{
return 7.0*(ms-0.212)*(ms-0.212) - 130.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212) + 30.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212) + 670.0*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212)*(ms-0.212);
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

//This can plot the branching ratios... if you uncomment some line from the Gamma_total function.
double branch_plotter(double chi, double mZ)
{

double Ue4 = 1/sqrt(3);
double Um4 = 1/sqrt(3);
double Ut4 = 1.0/sqrt(3);
double log_mS = 0.0;

for(log_mS=-5.0; log_mS<=log(0.490)/log(10)+1e-5; log_mS+=0.001)
{
	Gamma_total(chi, mZ, pow(10,log_mS), Ue4, Um4, Ut4);
} 
return 0;
}


