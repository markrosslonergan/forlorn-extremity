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

#include "channel.h"

twoIP_channel::twoIP_channel(gsl_rng * g, std::vector<double> input_params){

	std::vector<double> p  = {0.0,0.0,0.0};
//	p.push_back(0.0);	
//	p.push_back(0.0);	
//	p.push_back(0.0);	

	IP1.populate(1,p);
	IP2.populate(1,p);

	model_params = input_params;
	chan_identifier = model_params[2];	
	r = g;
}

int twoIP_channel::decayfunction(initial_sterile nuS, decay_params *)
{
	std::cout<<"You've somehow managed to call the decayfunction of the parent class (twoIP_channel). Don't do that."<<std::endl;
	return 0;
}
int twoIP_channel::observables(OBSERVABLES * output, gsl_rng *g)
{
	//OBSERVABLES { double E_sum; double Th_sum; double AngSep; double E_sterile; double Th_sterile; double E_high; double Th_high; double E_low; double Th_low; double FS_AngSep; } OBSERVABLES;

	int WHICHDIR =0;//was 0 for main

	fourmomentum sum;
	std::vector<double> p1 = {IP1.p.at(0),IP1.p.at(1),IP1.p.at(2)};
	std::vector<double> p2 = {IP2.p.at(0),IP2.p.at(1),IP2.p.at(2)};

	bool IS_ORDER = false;

	output->P_low=p2;
	output->P_high=p1;

	//		std::cout<<"test "<<(180.0/M_PI)*IP2.direction().at(0)<<" "<<acos(p2[2]/( pow(p2[0],2)+pow(p2[1],2)+pow(p2[2],2) ) )*180/3.14145<<

	//			std::cout<<Obs.Th_low<<" "<<atan2(Obs.P_low[0],Obs.P_low[2])*180/3.14159<<" "<<atan2(Obs.P_low[1],Obs.P_low[2])*180/3.14159<<std::endl;
	//			std::cout<<Obs.Th_low<<" "<<atan2(Obs.P_low[1],Obs.P_low[0])*180/3.14159<<" "<<atan2(Obs.P_low[2],Obs.P_low[0])*180/3.14159<<std::endl;
	//			std::cout<<Obs.Th_low<<" "<<atan2(Obs.P_low[0],Obs.P_low[1])*180/3.14159<<" "<<atan2(Obs.P_low[1],Obs.P_low[2])*180/3.14159<<std::endl;
	//			std::cout<<Obs.Th_low<<" "<<acos(Obs.P_low[2]/( pow(Obs.P_low[0],2)+pow(Obs.P_low[1],2)+pow(Obs.P_low[2],2) ) )*180/3.14145<<std::endl;


	std::vector<double> sum_p;
	sum_p.push_back(IP1.p.at(0) + IP2.p.at(0));
	sum_p.push_back(IP1.p.at(1) + IP2.p.at(1));
	sum_p.push_back(IP1.p.at(2) + IP2.p.at(2));

	sum.populate(IP1.E + IP2.E, sum_p);

	output->E_sum = sum.E;	
	output->Th_sum =(180.0/M_PI)*sum.direction().at(0);	

	output->AngSep =(180.0/M_PI)*acos((IP1.p.at(0)*IP2.p.at(0) + IP1.p.at(1)*IP2.p.at(1) + IP1.p.at(2)*IP2.p.at(2))/(IP1.modp*IP2.modp)); 

	//	IP1.print("IP1");
	//	IP2.print("IP2");

	//	if(IP1.p.at(2) > 0 && IP2.p.at(2) > 0)
	//	{
	//		output->FS_AngSep = (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
	//	}
	//	else 
	//	{
	//		output->FS_AngSep = 180 - (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
	//	}

	if(IP1.p.at(2)*IP2.p.at(2) >= 0 && IP2.p.at(0)*IP2.p.at(0) >= 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(fabs(atan(IP1.p.at(0)/IP1.p.at(2))) - fabs(atan(IP2.p.at(0)/IP2.p.at(2))));
	}
	else if(IP1.p.at(2)*IP2.p.at(2) >= 0 && IP2.p.at(0)*IP2.p.at(0) < 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) + atan(IP2.p.at(0)/IP2.p.at(2)));
	}
	else if(IP1.p.at(2)*IP2.p.at(2) < 0 && IP2.p.at(0)*IP2.p.at(0) >= 0)
	{
		output->FS_AngSep = 180.0- (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) + atan(IP2.p.at(0)/IP2.p.at(2)));
	}
	else if(IP1.p.at(2)*IP2.p.at(2) < 0 && IP2.p.at(0)*IP2.p.at(0) < 0)
	{
		output->FS_AngSep = 180.0- (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
	}


	output->E_high = IP1.E;	
	//	output->Th_high = (180.0/M_PI)*IP1.direction().at(0);	
	output->E_low = IP2.E;	
	//	output->Th_low = (180.0/M_PI)*IP2.direction().at(0);	


	output->Th_low = acos( p2[2]/sqrt( pow(p2[0],2)+pow(p2[1],2)+pow(p2[2],2)) )*180/M_PI; //my quick recalculation
	output->Th_high = acos( p1[2]/sqrt( pow(p1[0],2)+pow(p1[1],2)+pow(p1[2],2)) )*180/M_PI; //my quick recalculation

	output->Th_low_2 = acos( p2[1]/sqrt( pow(p2[0],2)+pow(p2[1],2)+pow(p2[2],2)) )*180/M_PI; //my quick recalculation
	output->Th_high_2 = acos( p1[1]/sqrt( pow(p1[0],2)+pow(p1[1],2)+pow(p1[2],2)) )*180/M_PI; //my quick recalculation

	output->Th_low_3 = acos( p2[0]/sqrt( pow(p2[0],2)+pow(p2[1],2)+pow(p2[2],2)) )*180/M_PI; //my quick recalculation
	output->Th_high_3 = acos( p1[0]/sqrt( pow(p1[0],2)+pow(p1[1],2)+pow(p1[2],2)) )*180/M_PI; //my quick recalculation



	double temp = 0.0;
	/*
	//Only want to do this in e+e- scenario, NOT mu pi scenario! hmm
	if(IS_ORDER && output->E_high < output->E_low)
	{ 	
	temp = output->E_low; 
	output->E_low = output->E_high; 
	output->E_high = temp; 

	temp = output->Th_low;
	output->Th_low = output->Th_high;
	output->Th_high = temp;
	}
	*/
	double mlow = sqrt(pow(IP2.E,2)-p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]);
	double mhigh = sqrt(pow(IP1.E,2)-p1[0]*p1[0]-p1[1]*p1[1]-p1[2]*p1[2]);

	output->Minvar = sqrt(mlow*mlow+mhigh*mhigh+2.0*(IP2.E*IP1.E-p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]));
	// Lets smear all of the above defined variables by preapproved gaussians

	double res_low = 0.0;
	double res_high = 0.0;
	double res_ang = 1.0;

	switch(output->chan_identifier)
	{
		case CHAN_ELECPOSI:
			res_low = 0.15;
			res_high = 0.15;
			break;
		case CHAN_ELECPI:
			res_low = 0.06;
			res_high = 0.15;
			break;
		case CHAN_MUONPI:
			res_low = 0.06;
			res_high = 0.06;
			break;
		case CHAN_NUPI0:
			res_low = 0.15;
			res_high = 0.15;
			break;
		case CHAN_GAMMA:
			res_low = 0.15;
			res_high = 0.15;
			break;
	}

	output->E_high_smear = smear_energy(output->E_high, mhigh, res_high, g);
	output->E_low_smear = smear_energy(output->E_low,  mlow, res_low, g);
	output->E_sum_smear = output->E_high_smear+output->E_low_smear;
	output->Th_high_smear = smear_angle(output->Th_high,res_ang,g);
	output->Th_low_smear = smear_angle(output->Th_low,res_ang,g);
	output->Th_sterile_smear = smear_angle(output->Th_sterile,res_ang,g);
	output->AngSep_smear = smear_angle(output->AngSep,res_ang,g);	
	output->FS_AngSep_smear = smear_angle(output->FS_AngSep,res_ang,g);	
	output->E_sterile_smear = output->E_sterile;
	output->Th_sum_smear = output->Th_sum_smear;

	output->Minvar_smear = sqrt(mlow*mlow+mhigh*mhigh+2.0*(output->E_low_smear*output->E_high-sqrt(pow(output->E_high_smear,2)-mhigh*mhigh)* sqrt(pow(output->E_low_smear,2)-mlow*mlow)*cos(output->AngSep_smear* 3.14159/180.0)  ) );

//	std::cout<<output->Minvar<<" "<<sqrt(mlow*mlow+mhigh*mhigh+2.0*(output->E_low_smear*output->E_high-sqrt(pow(output->E_high_smear,2)-mhigh*mhigh)* sqrt(pow(output->E_low_smear,2)-mlow*mlow)*cos(output->AngSep_smear* 3.14159/180.0)  ) )<<std::endl;

	//	IP1.print("IP1");
	//	IP2.print("IP2");

	return 0;
}


/* ###############################

   Below here we have a section for each channel.

############################### */

//This is the nu_s \to \nu e+ e- channel (off-shell Zprime).

threebody::threebody(gsl_rng *g, std::vector<double> input) : twoIP_channel(g, input)
{
	//Nothing to do.

}

int threebody::decayfunction(initial_sterile nuS, decay_params * params)
{
	double p0[] = {0.0,0.0,0.0,0.0};
	double p1[] = {0.0,0.0,0.0,0.0};
 	//Populate the momentum arrays in the rest frame.
	drawRestFrameDist(r,params,p0,p1);	

	//Then boost them in to the frame of nuS.
	computeLabFrameVariables(nuS,p0,p1);

return 0;
}

int threebody::computeLabFrameVariables(initial_sterile nuS, double p0[4], double p1[4])
{
	double mS = nuS.mass;
	double Es = nuS.energy;

	double beta = sqrt(1-mS*mS/(Es*Es));
	double gamma = 1.0/sqrt(1.0-beta*beta);

	double Px =nuS.labframeP.p.at(0);
	double Py =nuS.labframeP.p.at(1);
	double Pz =nuS.labframeP.p.at(2);

	std::vector< double > Vec_P0_lab = generic_boost(Es,-Px,-Py,-Pz, p0[0], p0[1],p0[2],p0[3] );
	std::vector< double > Vec_P1_lab = generic_boost(Es,-Px,-Py,-Pz, p1[0], p1[1],p1[2],p1[3] );
	std::vector< double > silly0 = {Vec_P0_lab[1],Vec_P0_lab[2],Vec_P0_lab[3]};
	std::vector< double > silly1 = {Vec_P1_lab[1],Vec_P1_lab[2],Vec_P1_lab[3]};

	IP1.populate(Vec_P0_lab[0], silly0);
	IP2.populate(Vec_P1_lab[0], silly1);

	return 0;
}


int threebody::drawRestFrameDist(gsl_rng * r, decay_params * params, double p2[4], double p3[4])
{
//This samples the decay including the ME. Another option is to pass the ME as
//a weight to the event. But that would take more thought.

	double p1[] = {1.0,0.0,0.0,0.0};
	double mS = params->mS;
	double m1 = 0.0;
	double m2 = e_mass;
	double m3 = e_mass;
	double Gamma = Gamma_EE(params);

	double PDF_MAX = 4.0;
	double x = gsl_rng_uniform(r);
	double y = gsl_rng_uniform(r);
	double z = PDF_MAX*gsl_rng_uniform(r);

	while(matrix_element_ee(x,y,params)/Gamma<=z)
	{
		x = gsl_rng_uniform(r);
		y = gsl_rng_uniform(r);
		z = PDF_MAX*gsl_rng_uniform(r);
	}	

	//We map our scaled variables back into normal parameters.
	double u_min = (m1+m2)*(m1+m2);
	double u_max = (mS-m3)*(mS-m3);
	double u = u_min + (u_max - u_min)*x;

	//The following energies and momenta are in a weird frame, but they give the bounds on t. See PDG kinematics. 
	double Estar_2 = (u + m2*m2-m1*m1)/(2.0*sqrt(u));
	double Estar_3 = (mS*mS - u - m3*m3)/(2.0*sqrt(u));
	double pstar_2 = sqrt(Estar_2*Estar_2 - e_mass*e_mass);
	double pstar_3 = sqrt(Estar_3*Estar_3 - e_mass*e_mass);

	double t_min = pow(Estar_2 + Estar_3,2.0) - (pstar_2+pstar_3)*(pstar_2+pstar_3);  
	double t_max = pow(Estar_2 + Estar_3,2.0) - (pstar_2-pstar_3)*(pstar_2-pstar_3);  
	double t = t_min + (t_max - t_min)*y;

	double s = mS*mS + m1*m1 + m2*m2 + m3*m3 - u - t;
	double phi_1 = 2*M_PI*gsl_rng_uniform(r);
	double theta_1 = M_PI*gsl_rng_uniform(r);
	double phi_2 = 2*M_PI*gsl_rng_uniform(r);

	//Then they get turned into useful variables.
	double E1 = ((mS*mS + m1*m1 - t)/(2.0*mS));	
	double E2 = ((mS*mS + m2*m2 - s)/(2.0*mS));	
	double E3 = ((mS*mS + m3*m3 - u)/(2.0*mS));	
	double mod_p1 = sqrt( E1*E1 - m1*m1);
	double mod_p2 = sqrt( E2*E2 - m2*m2);
	double mod_p3 = sqrt( E3*E3 - m3*m3);

	//Then update the final observables.
	//First we pretend p1 is aligned with the z-axis.
	p1[0]=E1;
	p1[1]=0.0; 
	p1[2]=0.0;
	p1[3]=mod_p1;

	//Which allows us to find the azimuthal angle of p2.
	double cos_theta_2 = (m1*m1 + m2*m2 + 2.0*E1*E2 - u)/(2.0*mod_p2*mod_p1);	
	p2[0]=E2;
	p2[1]=mod_p2*sqrt(1- cos_theta_2*cos_theta_2);
	p2[2]=0.0;
	p2[3]=mod_p2*cos_theta_2;

	//Then we rotate everything by phi_2 around p1, and then rotate
	//everything so that p1 is pointing along phi_1,theta_1.
	double R11 = cos(theta_1)*cos(phi_1)*cos(phi_2) + sin(phi_1)*sin(phi_2);
	double R12 = cos(theta_1)*cos(phi_1)*sin(phi_2) - sin(phi_1)*cos(phi_2);
	double R13 = sin(theta_1)*cos(phi_2); 	
	double R21 = cos(theta_1)*sin(phi_1)*cos(phi_2) - cos(phi_1)*sin(phi_2);
	double R22 = cos(theta_1)*sin(phi_1)*sin(phi_2) + cos(phi_1)*cos(phi_2);
	double R23 = sin(theta_1)*sin(phi_2); 	
	double R31 = -sin(theta_1)*cos(phi_2);
	double R32 = -sin(theta_1)*sin(phi_2);
	double R33 = cos(theta_1);

	p1[1] = R11*p1[1] +  R12*p1[2] +  R13*p1[3]; 
	p1[2] = R21*p1[1] +  R22*p1[2] +  R23*p1[3]; 
	p1[3] = R31*p1[1] +  R32*p1[2] +  R33*p1[3]; 
	p2[1] = R11*p2[1] +  R12*p2[2] +  R13*p2[3]; 
	p2[2] = R21*p2[1] +  R22*p2[2] +  R23*p2[3]; 
	p2[3] = R31*p2[1] +  R32*p2[2] +  R33*p2[3]; 

	//We construct the other momentum.
	p3[0]= mS - p1[0] - p2[0];
	p3[1]= -p1[1]-p2[1];
	p3[2]= -p1[2]-p2[2];
	p3[3]= -p1[3]-p2[3];
	
return 0;
}

std::vector<double > threebody::generic_boost(double Ep, double px, double py, double pz, double Er, double rx, double ry, double rz)
{		
	double pnorm = sqrt(px*px+py*py+pz*pz);
	double nx = px/pnorm;
	double ny = py/pnorm;
	double nz = pz/pnorm;

	double B = pnorm/Ep;
        double g = 1.0/sqrt(1-B*B);
	
	std::vector<double > boosted;

	boosted.push_back( Er*g-B*g*nx*rx-B*g*ny*ry-B*g*nz*rz);
	boosted.push_back( -B*Er*g*nx+(1+(-1+g)*nx*nx)*rx+(-1+g)*nx*ny*ry+(-1+g)*nx*nz*rz);
	boosted.push_back( -B*Er*g*ny+(-1+g)*nx*ny*rx+(1+(-1+g)*ny*ny)*ry+(-1+g)*ny*nz*rz);
	boosted.push_back( -B*Er*g*nz+(-1+g)*nx*nz*rx+(-1+g)*ny*nz*ry+(1+(-1+g)*nz*nz)*rz);
		
return boosted;
}


