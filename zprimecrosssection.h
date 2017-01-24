#include <gsl/gsl_integration.h>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <sstream>
#include <iomanip>

#define FLUX_CUTOFF 10.0
#define POSNU  0
#define POSNUBAR 1
#define NEGNU  2
#define NEGNUBAR 3

#define MP 0.13957018
#define MMu  0.1056583715
#define Me  0.0005109989
#define MKao  0.493677
#define PI 3.14159
#define Mw  80.385
#define MZ  91.1876
#define Gf  1.16637e-5
#define Mnuc 1


/**********************************************************************
 *                                                                    *
 * Flux Temp						              *
 *                                                                    *
 **********************************************************************/
double fluxPosNuMu (double En){
	std::vector<double > list = {4.544e-11,1.7132e-10,2.224e-10,2.67e-10,3.316e-10,3.64e-10,3.892e-10,4.09e-10,4.322e-10,4.482e-10,4.558e-10,4.584e-10,4.55e-10,4.506e-10,4.428e-10,4.312e-10,4.156e-10,3.984e-10,3.788e-10,3.578e-10,3.354e-10,3.116e-10,2.878e-10,2.636e-10,2.386e-10,2.138e-10,1.9006e-10,1.6712e-10,1.4556e-10,1.2584e-10,1.0792e-10,9.202e-11,7.804e-11,6.57e-11,5.52e-11,4.624e-11,3.864e-11,3.232e-11,2.71e-11,2.276e-11,1.9178e-11,1.63e-11,1.3856e-11,1.1874e-11,1.0294e-11,8.956e-12,7.87e-12,7.e-12,6.3e-12,5.734e-12,5.23e-12,4.818e-12,4.546e-12,4.22e-12,3.99e-12,3.84e-12,3.63e-12,3.452e-12,3.33e-12,3.202e-12,3.108e-12,2.986e-12,2.884e-12,2.824e-12,2.726e-12,2.646e-12,2.53e-12,2.434e-12,2.366e-12,2.28e-12,2.204e-12,2.12e-12,2.028e-12,1.94e-12,1.868e-12,1.8002e-12,1.7282e-12,1.638e-12,1.5734e-12,1.4928e-12,1.4292e-12,1.3624e-12,1.2998e-12,1.237e-12,1.1716e-12,1.1228e-12,1.064e-12,1.0032e-12,9.53e-13,9.122e-13,8.562e-13,8.174e-13,7.682e-13,7.264e-13,6.864e-13,6.526e-13,6.032e-13,5.714e-13,5.378e-13,5.058e-13,4.744e-13,4.454e-13,4.206e-13,3.914e-13,3.668e-13,3.46e-13,3.23e-13,3.026e-13,2.812e-13,2.606e-13,2.428e-13,2.258e-13,2.094e-13,1.9138e-13,1.774e-13,1.6296e-13,1.4858e-13,1.353e-13,1.2194e-13,1.0984e-13,9.954e-14,8.89e-14,7.934e-14,6.984e-14,6.074e-14,5.19e-14,4.45e-14,3.708e-14,3.074e-14,2.44e-14,1.956e-14,1.5684e-14,1.2396e-14,9.572e-15,6.668e-15,3.942e-15,1.8782e-15,5.476e-16,1.213e-16,8.27e-17,3.866e-17,1.9776e-17,8.988e-19,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double ans = 0;
	if(En <= 10.0 && En >=0){
		ans = list[floor(En/0.05)];
	}
	return ans;
}
double fluxPosNuMuBar (double En){
	std::vector<double > list ={5.12e-11,1.1342e-10,6.6e-11,4.056e-11,3.246e-11,2.79e-11,2.602e-11,2.498e-11,2.342e-11,2.108e-11,1.916e-11,1.739e-11,1.6196e-11,1.4868e-11,1.382e-11,1.2628e-11,1.181e-11,1.1008e-11,1.0158e-11,9.416e-12,8.694e-12,8.042e-12,7.406e-12,6.886e-12,6.346e-12,5.744e-12,5.194e-12,4.674e-12,4.202e-12,3.806e-12,3.436e-12,3.014e-12,2.682e-12,2.346e-12,2.106e-12,1.8482e-12,1.6376e-12,1.423e-12,1.2698e-12,1.1094e-12,9.598e-13,8.142e-13,7.184e-13,6.164e-13,5.276e-13,4.496e-13,3.756e-13,3.246e-13,2.782e-13,2.324e-13,2.02e-13,1.7382e-13,1.4764e-13,1.1998e-13,1.0008e-13,8.408e-14,7.142e-14,6.094e-14,5.194e-14,4.276e-14,3.912e-14,3.168e-14,2.454e-14,2.042e-14,1.6712e-14,1.5554e-14,1.3624e-14,1.4772e-14,1.2256e-14,1.2502e-14,1.1038e-14,7.872e-15,8.282e-15,6.79e-15,6.004e-15,5.004e-15,4.546e-15,4.598e-15,2.858e-15,3.148e-15,2.436e-15,2.56e-15,3.224e-15,1.7208e-15,1.854e-15,1.0742e-15,1.099e-15,8.552e-16,7.386e-16,1.3184e-15,1.2522e-15,4.532e-16,7.848e-16,1.0072e-15,6.102e-16,1.597e-15,3.26e-15,3.574e-15,1.1458e-15,1.2766e-16,1.0514e-16,1.0444e-16,8.738e-17,6.372e-17,7.83e-17,4.394e-17,3.38e-17,2.354e-17,1.9926e-17,1.8394e-17,1.358e-17,1.139e-17,1.0468e-17,6.418e-18,5.618e-18,5.4e-18,3.248e-18,2.766e-18,2.384e-18,1.8048e-18,1.8884e-18,1.0152e-18,1.278e-18,9.39e-19,5.468e-19,7.88e-19,4.134e-19,4.654e-19,4.588e-19,2.77e-19,3.864e-20,1.6598e-19,1.1708e-19,3.686e-20,0.,3.566e-20,8.98e-20,8.41e-20,0.,0.,0.,0.,0.,2.422e-20,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double ans = 0;

	if(En <= 10.0 && En >=0){
		ans = list[floor(En/0.05)];
	}
	return ans;
}
double fluxNegNuMu (double En){
	std::vector<double > list = {4.596e-11,1.1806e-10,7.452e-11,4.676e-11,5.14e-11,3.594e-11,3.552e-11,3.71e-11,3.668e-11,3.54e-11,3.402e-11,3.236e-11,3.11e-11,2.986e-11,2.85e-11,2.714e-11,2.604e-11,2.472e-11,2.384e-11,2.282e-11,2.174e-11,2.088e-11,1.9934e-11,1.887e-11,1.7652e-11,1.664e-11,1.5472e-11,1.436e-11,1.3218e-11,1.2106e-11,1.1066e-11,1.0116e-11,9.154e-12,8.268e-12,7.45e-12,6.672e-12,6.006e-12,5.326e-12,4.75e-12,4.252e-12,3.772e-12,3.338e-12,2.972e-12,2.62e-12,2.342e-12,2.06e-12,1.8558e-12,1.6398e-12,1.4706e-12,1.3154e-12,1.166e-12,1.0636e-12,9.644e-13,8.634e-13,7.994e-13,7.238e-13,6.75e-13,6.1e-13,5.852e-13,5.41e-13,5.318e-13,5.004e-13,4.782e-13,4.594e-13,4.382e-13,4.236e-13,4.104e-13,3.89e-13,3.772e-13,3.656e-13,3.518e-13,3.276e-13,3.138e-13,3.02e-13,2.942e-13,2.77e-13,2.698e-13,2.602e-13,2.526e-13,2.38e-13,2.344e-13,2.176e-13,2.126e-13,2.024e-13,1.909e-13,1.8716e-13,1.7328e-13,1.6778e-13,1.59e-13,1.5006e-13,1.4584e-13,1.4388e-13,1.3602e-13,1.2604e-13,1.1878e-13,1.1336e-13,1.0544e-13,1.0198e-13,9.62e-14,9.242e-14,8.728e-14,8.078e-14,7.682e-14,7.092e-14,6.776e-14,6.402e-14,5.928e-14,5.606e-14,5.302e-14,5.042e-14,4.678e-14,4.322e-14,4.024e-14,3.756e-14,3.558e-14,3.254e-14,3.04e-14,2.838e-14,2.57e-14,2.414e-14,2.18e-14,1.971e-14,1.827e-14,1.617e-14,1.4652e-14,1.313e-14,1.0906e-14,9.296e-15,7.814e-15,6.368e-15,5.026e-15,4.09e-15,3.34e-15,2.564e-15,1.882e-15,1.0478e-15,4.794e-16,1.3788e-16,1.3246e-17,4.494e-18,2.564e-18,1.5318e-18,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double ans = 0;

	if(En <= 10.0 && En >=0){
		ans = list[floor(En/0.05)];
	}
	return ans;
}
double fluxNegNuMuBar (double En){
	std::vector<double > list = {4.314e-11,1.568e-10,1.9462e-10,2.282e-10,2.638e-10,2.876e-10,2.954e-10,2.958e-10,3.e-10,2.97e-10,2.894e-10,2.812e-10,2.69e-10,2.574e-10,2.442e-10,2.304e-10,2.15e-10,1.996e-10,1.8354e-10,1.6822e-10,1.5316e-10,1.3814e-10,1.236e-10,1.101e-10,9.754e-11,8.538e-11,7.372e-11,6.302e-11,5.356e-11,4.524e-11,3.796e-11,3.16e-11,2.622e-11,2.166e-11,1.7834e-11,1.457e-11,1.1882e-11,9.668e-12,7.874e-12,6.36e-12,5.154e-12,4.132e-12,3.33e-12,2.692e-12,2.162e-12,1.7674e-12,1.4272e-12,1.1414e-12,9.24e-13,7.556e-13,6.056e-13,4.824e-13,3.954e-13,3.276e-13,2.646e-13,2.076e-13,1.7414e-13,1.3962e-13,1.2156e-13,1.0222e-13,7.838e-14,6.656e-14,5.722e-14,4.764e-14,4.59e-14,4.538e-14,3.656e-14,3.226e-14,3.074e-14,2.75e-14,2.494e-14,2.088e-14,1.9064e-14,1.6296e-14,1.497e-14,1.5674e-14,1.2404e-14,1.3228e-14,1.099e-14,1.0264e-14,9.8e-15,9.968e-15,8.958e-15,4.812e-15,4.392e-15,3.946e-15,3.524e-15,2.678e-15,2.512e-15,1.9338e-15,1.5092e-15,1.196e-15,1.2582e-15,9.48e-16,7.162e-16,6.578e-16,5.616e-16,5.412e-16,4.878e-16,2.954e-16,2.48e-16,2.07e-16,1.8002e-16,1.374e-16,1.715e-16,8.784e-17,9.472e-17,6.076e-17,4.826e-17,4.144e-17,4.122e-17,2.226e-17,1.8512e-17,1.3066e-17,1.6158e-17,1.3324e-17,6.95e-18,7.616e-18,6.99e-18,5.042e-18,3.934e-18,2.528e-18,1.3784e-18,1.4774e-18,6.934e-19,2.052e-18,3.206e-19,3.24e-19,2.26e-19,2.79e-19,4.228e-20,7.886e-20,3.776e-20,1.3424e-19,1.1888e-19,0.,2.872e-20,0.,0.,0.,0.,0.,0.,2.872e-20,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double ans = 0;

	if(En <= 10.0 && En >=0){
		ans = list[floor(En/0.05)];
	}
	return ans;
}

double flux(double En, int selector){
	double ans = 0;
	switch (selector){
		case POSNU:
			ans = fluxPosNuMu(En);
			break;
		case POSNUBAR:
			ans = fluxPosNuMuBar(En);
			break;
		case NEGNU:
			ans = fluxNegNuMu(En);
			break;
		case NEGNUBAR:
			ans = fluxNegNuMuBar(En);
			break;	
	}
	return ans;
}




double fluxPosNuMuNormed (double En){
	return fluxPosNuMu(En)/5.19358e-10;
}
double fluxPosNuMuBarNormed (double En){
	return fluxPosNuMuBar(En)/3.26002e-11;
}
double fluxNegNuMuNormed (double En){
	return fluxNegNuMu(En)/5.42276e-11;
}
double fluxNegNuMuBarNormed (double En){
	return fluxNegNuMuBar(En)/2.9282e-10;
}

double fluxNormed(double En, int selector){
	double ans = 0;
	switch (selector){
		case POSNU:
			ans = fluxPosNuMuNormed(En);
			break;
		case POSNUBAR:
			ans = fluxPosNuMuBarNormed(En);
			break;
		case NEGNU:
			ans = fluxNegNuMuNormed(En);
			break;
		case NEGNUBAR:
			ans = fluxNegNuMuBarNormed(En);
			break;	
	}
	return ans;
}



double fluxTempOldOBSOLETE(double En) {
	double result =0;
	if (0<=En && En <= 12.55){
		double kflux = 1000/50.0;
		result = kflux*(4.427824623217513e-14*pow(En,7)-9.971799853312788e-13*pow(En,6) + 9.231207886825193e-12*pow(En,5) - 4.493562336824403e-11*pow(En,4) + 1.2126620243628683e-10*pow(En,3) - 1.7206489225632898e-10*pow(En,2)+ 1.0055899381298395e-10*En - 5.096969325003961e-13);
	}

return result;
}

double fluxTemp(double En){
	double result = 0 ;
	double kflux = 1000/50.0;
	double heavi = 1.0;
	if(En > 4){
		heavi = 0.0;
	}
	double ktop=(1.51364530023684e-9*pow(En,2.5608335856419675));
	double kbot=   (pow(0.04173905329386966 + pow(En,1.942453984231271),3.7059504555184652)*(3986.923260185046/pow(En,4.646872292396358) + 0.14401684718576835*pow(En,2.5608335856419675))); 
	double ptop=(5.479159036783009e-7*pow(En,0.29615835741178265));
	double pbot=(pow(1.7304776909979156 + pow(En,1.215969360903414),11.5953242599689)*
     (0.05068826655476809*pow(En,0.29615835741178265) + 33176.287610674655/pow(1.7304776909979156 + En,11.5953242599689)));

	result = kflux*(ktop/kbot+heavi*ptop/pbot);
return result;
}



double fluxTempNormed (double En) {
	return fluxTemp(En)/4.698e-10;
}


/**********************************************************************
 *                                                                    *
 * KINEMATICS						              *
 *                                                                    *
 **********************************************************************/

double mandelT (double Es, double Ev, double mN) {
	return 2*mN*(Es-Ev);
}

double mandelS  (double Ev, double mV, double mN) {
	return mN*mN+mV*mV+2*mN*Ev;
}

double LABcos (double Es, double Ms, double Ev, double Mv, double Mn){

double top = (2*Es*Ev+2*Es*Mn-2*Ev*Mn-Ms*Ms-Mv*Mv);
double bot = 2*sqrt(Es*Es-Ms*Ms)*sqrt(Ev*Ev-Mv*Mv);

return top/bot;
}

double LABesMin (double Ev, double ms, double mv, double mn){
	double b1 = (Ev+mn)*(2*Ev*mn+ms*ms+mv*mv);
	double b2= (Ev-mv)*(Ev+mv)*(2*Ev*mn + 2*mn*ms - ms*ms + mv*mv)*(2*Ev*mn -ms*(2*mn + ms) + mv*mv);
	double b3= 2*(2*Ev*mn+mn*mn+mv*mv);

	return (b1-sqrt(b2))/b3; 
}

double LABesMax (double Ev, double ms, double mv, double mn){
	double b1 = (Ev+mn)*(2*Ev*mn+ms*ms+mv*mv);
	double b2= (Ev-mv)*(Ev+mv)*(2*Ev*mn + 2*mn*ms - ms*ms + mv*mv)*(2*Ev*mn -ms*(2*mn + ms) + mv*mv);
	double b3= 2*(2*Ev*mn+mn*mn+mv*mv);

	return (b1+sqrt(b2))/b3; 
}

double mandelTmin (double ms, double Ev, double mv, double mn) {
return mandelT(LABesMax(Ev,ms,mv,mn), Ev, mn);
}

double mandelTmax (double ms, double Ev, double mv, double mn) {
return mandelT(LABesMin(Ev,ms,mv,mn), Ev, mn);
}

double LABevMin(double ms, double mv, double mn){

return -(ms*ms+mv*mv-2*mn*ms)/(2*mn-2*ms);
}

/**********************************************************************
 *                                                                    *
 * 		Form Factor Cross Section		              *
 *                                                                    *
 **********************************************************************/



double CSffactor (double Qs, double Ev, double mz, double ms, double X, double U) {	
	//My Differentia CSection.
	double s2th= 0.23149;
	
	double ans = 0;
	double gp = 1;
	double MA = 1.39;
	double gA = 1.2671;
	double MV = 0.84;
	double Ds = 0.08;
	double muN = -1.91304273;
	double muP = 2.792847351;
	double t3  = -1;
	double tau = Qs/(4*Mnuc*Mnuc);
	double W = 4*Mnuc*Ev-Qs-ms*ms;
	//Sachs form factors

	double GD = 1/((1+Qs/(MV*MV))*(1+Qs/(MV*MV)));
	double GEn = 0;
	double GEp = GD;
	double GMn = muN*GD;
	double GMp = muP*GD;

	//Dirac and Pauly FormFac
	double	F1EMN = (GEn+tau*GMn)/(1+tau);
	double	F2EMN = (GMn-GEn)/(1+tau);
	double	F1EMP = (GEp+tau*GMp)/(1+tau);
	double	F2EMP = (GMp-GEp)/(1+tau);

	// according to 1007:4730 appendix isoscalar constsitant with 0
	double F1S =0;
	double F2S = 0;


	//write weak through dirac and pauli form fac

	double F1Z = (0.5-s2th)*(F1EMP-F1EMN)*t3-s2th*(F1EMP+F1EMN)-0.5*F1S;
	double F2Z = (0.5-s2th)*(F2EMP-F2EMN)*t3-s2th*(F2EMP+F2EMN)-0.5*F2S;
	
	double FA = gA/((1+Qs/(MA*MA))*(1+Qs/(MA*MA)));
	double FAS = Ds/((1+Qs/(MA*MA))*(1+Qs/(MA*MA)));

	double FAZ = t3*0.5*FA-0.5*FAS;
	double FPZ = 2*Mnuc*Mnuc*FAZ/(MP*MP+Qs);

	//Write A,B,C in terms of all nucelon, dirac, pauli and axial weak

	double A = (ms*ms+Qs)/(Mnuc*Mnuc)*((FAZ*FAZ*(1+tau)+(F1Z*F1Z-tau*F2Z*F2Z)*(1-tau)+4*tau*F1Z*F2Z));
	double B = 4*tau*FAZ*(F1Z+F2Z);
	double C = 0.25*(FAZ*FAZ+F1Z*F1Z+tau*F2Z*F2Z);

	
	ans = Gf*Gf*MZ*MZ*MZ*MZ*U*U*X*X*Mnuc*Mnuc*Mnuc*Mnuc*(A+B*W/(Mnuc*Mnuc)+C*W*W/(Mnuc*Mnuc*Mnuc*Mnuc))/(8*PI*Ev*Ev*(mz*mz+Qs)*(mz*mz+Qs));


	return ans;
}
double CSProton(double t, double Ev, double mz, double ms, double X, double U) {	
	//My Differentia CSection.
	double Qs = -t;
	
	double mN=1;
	
	double ee = sqrt(4*3.14148/137.0);
	double gX=1.0;
	double MV = 0.84;
	double cw = 80.385/91.1876;
	//Sachs form factors

	double GD = 1/((1+Qs/(MV*MV))*(1+Qs/(MV*MV)));
	double GE = GD;

	//0.938272;
	double muN = 2.792847351; // -1.91304273
		//0.939565;

	double GM = muN*GD;

	double tau = Qs/(4*mN*mN);




	//Dirac and Pauly FormFac
	double	F1 = (GE+tau*GM)/(1+tau);
	double	F2 = (GM-GE)/(1+tau);
	double A = 8*F1*F1*(8*Ev*Ev/(t-ms*ms)+4*Ev/mN+t/(mN*mN)+2);
	double B = 4*F2*F2*(-4*Ev*Ev*t/(mN*mN*(t-ms*ms))-2*Ev*t/(mN*mN*mN)+(ms*ms+t)/(mN*mN)+ ms*ms*t/(4*mN*mN*mN));
	double C = 8*F1*F2*(ms*ms +2*t)/(mN*mN);
	
	double Anew = 8*F1*F1*(8*Ev*Ev/(ms*ms-t)-4*Ev/mN-t/(mN*mN)+2);
	double Bnew = 4*t*F2*F2/(mN*mN)*(-4*Ev*Ev/(ms*ms-t)+2*Ev/mN-ms*ms/(4*mN*mN)+1);
	double Cnew = 8*F1*F2*(4*Ev*ms*ms/(mN*(ms*ms-t))-ms*ms/(mN*mN));


	double ans = pow(ee*gX*cw*U*X,2)*(t-ms*ms)/(32*3.14159*Ev*Ev*pow(mz*mz-t,2))*(A+B+C);
	double ansnew = pow(ee*gX*cw*U*X,2)*(ms*ms-t)/(256*3.14159*Ev*Ev*pow(mz*mz-t,2))*(Anew+Bnew+Cnew);
	return ansnew;
}
double CSNeutron(double t, double Ev, double mz, double ms, double X, double U) {	
	//My Differentia CSection.
	double Qs = -t;
	
	double mN=1;
	
	double ee = sqrt(4*3.14148/137.0);
	double gX=1.0;
	double MV = 0.84;
	double cw = 80.385/91.1876;
	//Sachs form factors

	double GD = 1/((1+Qs/(MV*MV))*(1+Qs/(MV*MV)));
	double GE = 0;

	//0.938272;
	double muN = -1.91304273;
		//0.939565;

	double GM = muN*GD;

	double tau = Qs/(4*mN*mN);

	//Dirac and Pauly FormFac
	double	F1 = (GE+tau*GM)/(1+tau);
	double	F2 = (GM-GE)/(1+tau);
	double A = 8*F1*F1*(8*Ev*Ev/(t-ms*ms)+4*Ev/mN+t/(mN*mN)+2);
	double B = 4*F2*F2*(-4*Ev*Ev*t/(mN*mN*(t-ms*ms))-2*Ev*t/(mN*mN*mN)+(ms*ms+t)/(mN*mN)+ ms*ms*t/(4*mN*mN*mN));
	double C = 8*F1*F2*(ms*ms +2*t)/(mN*mN);

	double Anew = 8*F1*F1*(8*Ev*Ev/(ms*ms-t)-4*Ev/mN-t/(mN*mN)+2);
	double Bnew = 4*t*F2*F2/(mN*mN)*(-4*Ev*Ev/(ms*ms-t)+2*Ev/mN-ms*ms/(4*mN*mN)+1);
	double Cnew = 8*F1*F2*(4*Ev*ms*ms/(mN*(ms*ms-t))-ms*ms/(mN*mN));


	double ansnew = pow(ee*gX*cw*U*X,2)*(ms*ms-t)/(256*3.14159*Ev*Ev*pow(mz*mz-t,2))*(Anew+Bnew+Cnew);
	
	double ans = pow(ee*gX*cw*U*X,2)*(t-ms*ms)/(32*3.14159*Ev*Ev*pow(mz*mz-t,2))*(A+B+C);

	return ansnew;
}

double CSCarbon(double t, double Ev, double mz, double ms, double X, double U) {	
	return 6*CSNeutron(t,Ev,mz,ms,X,U)+6*CSProton(t,Ev,mz,ms,X,U);
}


struct TcsParams { double Ev; double mz; double ms; double X; double U;};

double csToBeInt  (double tvar,  void * p) {
	struct TcsParams * params  = (struct TcsParams *)p;
	return CSffactor(-tvar, params->Ev, params->mz, params->ms,params->X, params->U); 
}

double totalCS (double Evin, double Mz, double Ms, double X, double U) {
	double ans = 0;
	// Function to integrate the diff cs from its Max and Min Mandel T bounds for any given Neutrino Energy
	
	if(Evin >= LABevMin(Ms,0.0,Mnuc)){
		double tint1 = mandelTmax(Ms,Evin,0.0,Mnuc);
		double tint2 = mandelTmin(Ms,Evin,0.0,Mnuc);
	
		
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		struct TcsParams params = {Evin, Mz, Ms,X,U};
		double result,error;
		gsl_function F;
		F.function = &csToBeInt;
		F.params = &params;
	
		//std::cout<<std::setprecision(12)<<" Evin: "<<Evin<<" ms: "<<Ms<<" mz: "<<Mz<<" Evin2: "<<params.Ev<<" ms: "<<params.ms<<" mz: "<<params.mz<<std::endl;
		gsl_integration_qags (&F,tint1, tint2 ,0, 1e-4, 1000, w, &result, &error); 
		gsl_integration_workspace_free (w);

		ans = result;
	}


return ans;

}

double csToBeIntProton  (double tvar,  void * p) {
	struct TcsParams * params  = (struct TcsParams *)p;
	return CSProton(tvar, params->Ev, params->mz, params->ms,params->X, params->U); 
}
double csToBeIntNeutron  (double tvar,  void * p) {
	struct TcsParams * params  = (struct TcsParams *)p;
	return CSNeutron(tvar, params->Ev, params->mz, params->ms,params->X, params->U); 
}

double totalCSProton (double Evin, double Mz, double Ms, double X, double U) {
	double ans = 0;
	// Function to integrate the diff cs from its Max and Min Mandel T bounds for any given Neutrino Energy
	
	if(Evin >= LABevMin(Ms,0.0,Mnuc)){
		double tint1 = mandelTmax(Ms,Evin,0.0,Mnuc);
		double tint2 = mandelTmin(Ms,Evin,0.0,Mnuc);
	
		
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		struct TcsParams params = {Evin, Mz, Ms,X,U};
		double result,error;
		gsl_function F;
		F.function = &csToBeIntProton;
		F.params = &params;
	
		//std::cout<<std::setprecision(12)<<" Evin: "<<Evin<<" ms: "<<Ms<<" mz: "<<Mz<<" Evin2: "<<params.Ev<<" ms: "<<params.ms<<" mz: "<<params.mz<<std::endl;
		gsl_integration_qags (&F,tint1, tint2 ,0, 1e-4, 1000, w, &result, &error); 
		gsl_integration_workspace_free (w);

		ans = result;
	}


return ans;

}
double totalCSNeutron (double Evin, double Mz, double Ms, double X, double U) {
	double ans = 0;
	// Function to integrate the diff cs from its Max and Min Mandel T bounds for any given Neutrino Energy
	
	if(Evin >= LABevMin(Ms,0.0,Mnuc)){
		double tint1 = mandelTmax(Ms,Evin,0.0,Mnuc);
		double tint2 = mandelTmin(Ms,Evin,0.0,Mnuc);
	
		
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		struct TcsParams params = {Evin, Mz, Ms,X,U};
		double result,error;
		gsl_function F;
		F.function = &csToBeIntNeutron;
		F.params = &params;
	
		//std::cout<<std::setprecision(12)<<" Evin: "<<Evin<<" ms: "<<Ms<<" mz: "<<Mz<<" Evin2: "<<params.Ev<<" ms: "<<params.ms<<" mz: "<<params.mz<<std::endl;
		gsl_integration_qags (&F,tint1, tint2 ,0, 1e-4, 1000, w, &result, &error); 
		gsl_integration_workspace_free (w);

		ans = result;
	}


return ans;

}

	
std::vector<double > MaxTotalCs(double Mz, double Ms, double X, double U){
// Very Simple 1-D raster scan to find maximum Total CS and its corresponding Ev 
	double emax = 0;
	double csmax = 0;
	double tempCS = 0;

	for( double ei = LABevMin(Ms,0.0,Mnuc); ei<=FLUX_CUTOFF ;ei=ei+0.0001){
		tempCS = 6*totalCSProton(ei, Mz, Ms, X, U)+6*totalCSNeutron(ei,Mz,Ms,X,U);
		if(tempCS > csmax){
			csmax = tempCS;
			emax = ei;
		}
	}
std::vector<double > ans = {emax,csmax};
return ans;

}

std::vector<double > MaxDiffCs(double Mz, double Ms, double X, double U){
// Very Simple 1-D raster scan to find maximum Diff CS and its corresponding Ev. Making assumption that for each Ev, the max will be at the tmin or tmax (to within 20%)
	double emax = 0;
	double csmax = 0;
	double tempCS = 0;
	double tmint = 0;
	double tmaxt = 0;
	
	for( double ei = LABevMin(Ms,0.0,Mnuc); ei<=FLUX_CUTOFF;ei=ei+0.01){
	         tmaxt = mandelTmax(Ms,ei,0.0,Mnuc);
		 tmint = mandelTmin(Ms,ei,0.0,Mnuc);
		 tempCS = std::max(CSCarbon(tmaxt, ei, Mz,Ms, 1,1) ,CSCarbon(tmint, ei, Mz,Ms, 1,1)  );
				
		if(tempCS >= csmax){
			csmax = tempCS;
			emax = ei;
		}
	}
std::vector<double > ans = {emax,1.2*csmax};
return ans;
}
/**********************************************************************
 *                                                                    *
 * 		Events!					              *
 *                                                                    *
 **********************************************************************/

/*double totalevents (double X, double U, double mz, double ms) {
double result=0;

	double binw = 0.01;
	double POT = 6.46e20;
	double ATOMS = 3.5e31*pow(5/6.1,3);
	double GEV2cm = pow(1e9*5.068e4,-2);
	double eE = 0.9;
	double eC = 0.75;
	
	double lowBound = LABevMin(ms,0.0,Mnuc);
	
	for(double en = lowBound; en <=FLUX_CUTOFF; en=en+binw){
	//	std::cout<<result<<std::endl;
		result = result +12*binw*GEV2cm*eC*eE*totalCS(en,mz,ms,X,U)*fluxTemp(en)*POT*ATOMS;
	}

return result;
}

*/

double POT(int selector){
	double ans = 1;
	switch (selector){
		case POSNU:
			ans = 6.46e20;
			break;
		case POSNUBAR:
		       	ans= 6.46e20;
			break;
		case NEGNU:
			ans = 11.27e20;
			break;
		case NEGNUBAR:
			ans= 11.27e20;
			break;
	}
	return ans;
}

double CSmodifier(int selector){
	double ans = 1;
	switch (selector){
		case POSNU:
			ans = 1.0;
			break;
		case POSNUBAR:
		       	ans= 1.0;//was 0.33
			break;
		case NEGNU:
			ans = 1.0;
			break;
		case NEGNUBAR:
			ans= 1.0;//was 0.33
			break;
	}
	return ans;
}


double totalevents2 (double X, double U, double mz, double ms, int selector) {
double result=0;

	double binw = 0.01;
	double ATOMS = 3.5e31*pow(5/6.1,3);
	double GEV2cm = pow(1e9*5.068e4,-2);
	double eE = 0.9;
	double eC = 0.75;
	
	double lowBound = LABevMin(ms,0.0,Mnuc);
	
	for(double en = lowBound; en <=FLUX_CUTOFF; en=en+binw){
	//	std::cout<<result<<std::endl;
		result = result +CSmodifier(selector)*binw*GEV2cm*eC*eE*(6*totalCSNeutron(en,mz,ms,X,U)+6*totalCSProton(en,mz,ms,X,U))*flux(en,selector)*POT(selector)*ATOMS;
	}

return result;
}
/**********************************************************************
 *                                                                    *
 * 		Subsequent Decay t body			              *
 *                                                                    *
 **********************************************************************/
double RFegamma(double ms, double mgam, double m5){

return (ms*ms+mgam*mgam-m5*m5)/(2*ms);
}

double egamma(double es, double ms, double mgam, double m5, gsl_rng * r){
	double gamma = es/ms;
	double beta = sqrt(1-pow(gamma,-2));
	double RFcosth = 2*gsl_rng_uniform(r)-1;
	double result = gamma*RFegamma(ms, mgam, m5)+gamma*beta*sqrt(pow(RFegamma(ms, mgam,m5),2)-mgam*mgam)*RFcosth;

return result;

}

/**********************************************************************
 *                                                                    *
 * 		 MC calc	           			      *
 *                                                                    *
 **********************************************************************/
double SampleNuFlux(double ms, gsl_rng * r, int selector){
		double Evin=0;
		double mc2 = 0;

		Evin = (FLUX_CUTOFF-LABevMin(ms,0.0,Mnuc))*gsl_rng_uniform(r)+LABevMin(ms,0.0,Mnuc);
		mc2  = gsl_rng_uniform(r); 
		//howmanyrngs++;

		while(mc2 > fluxNormed(Evin,selector)){
			Evin = (FLUX_CUTOFF-LABevMin(ms,0.0,Mnuc))*gsl_rng_uniform(r)+LABevMin(ms,0.0,Mnuc);
			mc2  = gsl_rng_uniform(r);
			//howmanyrngs++;
		}
return Evin;
}



std::vector<double > MCkin (double mz, double ms, gsl_rng * r, double CSmax, int selector){
	
	double howmanyrngs=2;  //Currently Broken

	double Evin = SampleNuFlux(ms,r,selector);
	double tmax= mandelTmax(ms,Evin,0.0,Mnuc);
	double tmin = mandelTmin(ms,Evin,0.0,Mnuc);

	double tout = (tmax-tmin)*gsl_rng_uniform(r)+tmin;
	double func = CSCarbon(tout, Evin,  mz, ms,  1.0,  1.0);
	double mc =   CSmax*gsl_rng_uniform(r);

	// Loop until an event is sampled from the distributions
	while(mc > func){
		Evin = SampleNuFlux(ms, r, selector);
		tmax= mandelTmax(ms,Evin,0.0,Mnuc);
		tmin = mandelTmin(ms,Evin,0.0,Mnuc);
		tout = (tmax-tmin)*gsl_rng_uniform(r)+tmin;
		func = CSCarbon(tout, Evin,  mz, ms,  1.0,  1.0);///totcs;
		mc = CSmax*gsl_rng_uniform(r);	
		//howmanyrngs++;
	}	

	// Calculate end observables for this event
	double Esout = (2*Evin*Mnuc + tout)/(2*Mnuc);
	double costh = LABcos(Esout, ms, Evin, 0.0 , Mnuc);
	double Edecay = egamma(Esout, ms, 0.0, 0, r);


std::vector<double > event = {Evin,tout,Esout,costh,Edecay,Esout/ms, howmanyrngs};
return event;


}



int MCmassrun( double mz, double ms, int numEvents, gsl_rng * r,int selector) {
	// Takes an input pair of Ms and Mz, along with a initilised random number generator, and a number of events to generate.

	std::vector<double > event ;
	std::stringstream str; str.str( std::string() ); str.clear();
	std::stringstream str2; str2.str( std::string() ); str2.clear();
	std::ofstream myfile;
	std::ofstream myfileHist;
	
	
 	//What I want to name the files,
	str << std::fixed<<std::setprecision(3) << "MC_"<<ms<<"_"<<mz<< ".dat";
	str2 << std::fixed<<std::setprecision(3) << "HIST_"<<ms<<"_"<<mz<< ".dat";
	myfile.open (str.str().c_str());		//Send as a standard C-style string
	myfileHist.open (str2.str().c_str());		//Send as a standard C-style string
	
	myfileHist<<"########## Ms | Mz | Total Num Sterile | MC points ##########"<<std::endl;
	myfileHist<<ms<<"\t"<<mz<<"\t"<<totalevents2(1,1,mz,ms, selector)<<"\t"<<numEvents<<std::endl;


	//Maximum Total CS for normalisation sampling. 
	std::vector<double > CSmax = MaxDiffCs(mz,ms,1,1);

	for(int j=0; j<numEvents; j++){
		event = MCkin(mz,ms,r,CSmax[1],selector); //Get Event

		//for(int i=0; i<event.size();i++){ //Print events to fil
		//                              Ev tout Esout Costh Edecay Esout/ms how many rngs
		
	
		//	myfile<<std::fixed<<std::setprecision(6)<<event[0]<<"\t"<<event[1]<<"\t"<<event[2]<<"\t"<<event[3]<<"\t"<<event[4]<<"\t"<<event[5]<<std::endl;
		
			myfile<<std::fixed<<std::setprecision(6)<<event[2]<<"\t"<<event[3]<<std::endl;
	
		//}
		
	}

	


	myfile.close();
	return 0;
}










