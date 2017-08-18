#include "appellced.h"

double AppellF1(double alpha, double beta, double betap, double gamma, double x, double y, double* perr)
{
  double Snow = 0.0; 
  double Sprev = 0.0; 
  int n = 0;
  double test = 100.0;
  double prefactor = 1.0;

//std::cout<<"a = "<<alpha<<", b1 = "<<beta<<", b2 = "<<betap<<", c = "<<gamma<<", x = "<<x<<", y = "<<y<<std::endl;

if(beta == 0.0 || x == 0.0)
{
  Snow = G2F1(alpha,betap,gamma,y);
}
else if(betap==0.0 || y == 0.0)
{
  Snow = G2F1(alpha,beta,gamma,x);
}
else
{
  //I think these checks improve convergence for special cases.
  if(fabs((x-y)/(x-1))<y){ beta = gamma-beta-betap; y = (x-y)/(x-1); x = x/(x-1); prefactor = pow(1-x,-alpha); }
  else if (fabs((y-x)/(y-1))<x) { betap = gamma-beta - betap; x = (y-x)/(y-1); y = y/(y-1); prefactor = pow(1-y,-alpha); }

  //Then we hit the main loops.
  if(y<10.94)//I only care about y<1, and it never seems beneficial to use the other series. Need to actually implement the analytic cont. 
  {
      while(test > 1e-7)
    {
  	Sprev = Snow;
	Snow += (poch(betap, (int)n)*poch(alpha,(int)n)/(tgamma(n+1.0)*poch(gamma,(int)n)))*G2F1(alpha+n,beta,gamma+n,x)*pow(y,n);
	test = fabs((Snow-Sprev)/Snow);
	n++;
    }
  }
  else
  {
     while(test > 1e-7)
    {
	Sprev = Snow;
	Snow += (tgamma(gamma)*tgamma(alpha-gamma-betap)/(tgamma(alpha)*tgamma(betap)))*pow(1-y,gamma-alpha-betap)*(poch(beta,n)/tgamma(n+1))*G2F1(gamma-alpha,gamma+n-betap,gamma-alpha-betap+1,1-y)*pow(x,n)
 	+ (tgamma(gamma)*tgamma(gamma-alpha-betap)/(tgamma(gamma-alpha)))*(poch(beta,n)*poch(alpha,n)/(tgamma(n+1)*tgamma(gamma+n-betap)))*G2F1(alpha+n,betap,alpha-gamma+betap+1,1-y)*pow(x,n);
	test = fabs((Snow-Sprev)/Snow);
	n++;
    }
  }
}

if(perr!=NULL){
  *perr = fabs((Snow-Sprev)/Snow);
}

//std::cout<<"Iterations: "<<n<<std::endl;

return prefactor*Snow;
}

//double AppellF1_sum(double alpha, double beta, double betap, double gamma, double x, double y, double* perr)
//{
//  double low_y = 0.0;
//  double Snow = 0.0; 
//  double Sprev = 0.0; 
//  double test = 100.0; 
//  int n = 0;
//
//  //Mathematically what follows is a meaningless way to ensure convergence, but I
//  //don't know how else to do it.
//  while(test > 1e-6)
//  {
//	Sprev = Snow;
//	Snow += (poch(betap, (int)n)*poch(alpha,(int)n)/(tgamma(n+1.0)*poch(gamma,(int)n)))*G2F1(alpha+n,beta,gamma+n,x)*pow(y,n);
//	test = fabs((Snow-Sprev)/Snow);
//	n++;
//  }
//
//if(perr!=NULL){
//  *perr = fabs((Snow-Sprev)/Snow);
//}
//
//return Snow;
//}


//The pochhammer symbol.
double poch(double x, int k)
{
return tgamma(x+(double)k)/tgamma(x);
}

double G2F1(double a, double b, double c, double z)
{
  std::complex<double> ac(a,0.0);
  std::complex<double> bc(b,0.0);
  std::complex<double> cc(c,0.0);
  std::complex<double> zc(z,0.0);

return std::real(hyp_2F1(a,b,c,z));
}

double G2F1_Forrey(double a, double b, double c, double z)
{
//This was my own implementation based on https://cran.r-project.org/web/packages/hypergeo/vignettes/hypergeometric.pdf
//Which I believe is based on Forrey (1997)
//It was replaced with the Michel/Stoitsov algorithm.
//
//If the answer should be complex, it nans. I consider this a feature.

//std::cout<<"2F1: z = "<<z<<" a = "<<a<<" b = "<<b<<" c = "<<c<<std::endl;

  std::vector<double> args { fabs(z), fabs(z/(z-1.0)), fabs(1.0-z), fabs(1.0/z), fabs(1.0/(1.0-z)),fabs(1.0- 1.0/z)};

  int argMin = std::distance(args.begin(), std::min_element(args.begin(),args.end()));

//  std::cout<<"Options: "<<args.at(0)<<" "<<args.at(1)<<" "<<args.at(2)<<" "<<args.at(3)<<" "<<args.at(4)<<" "<<args.at(5)<<std::endl;
//  std::cout<<"Running with: "<<argMin<<std::endl;

  switch(argMin)
  {
    case 0: //normal
	return G2F1_sum(a,b,c,z);
    case 1: //15.3.4
	return pow(1-z,-a)*G2F1_sum(a,c-b,c,z/(z-1.0));
    case 2: //15.3.6
	return (tgamma(c)*tgamma(c-b-a)/(tgamma(c-a)*tgamma(c-b)))*G2F1_sum(a,b,a+b-c+1,1.0-z) + pow(1-z,c-a-b)*(tgamma(c)*tgamma(a+b-c)/(tgamma(a)*tgamma(b)))*G2F1_sum(c-a,c-b,c-a-b+1,1.0-z);
    case 3: //15.3.7
	return (tgamma(c)*tgamma(b-a)/(tgamma(b)*tgamma(c-a)))*pow(-z,-a)*G2F1_sum(a,1-c+a,1-b+a,1.0/z) + (tgamma(c)*tgamma(a-b)/(tgamma(a)*tgamma(c-b)))*pow(-z,-b)*G2F1_sum(b,1-c+b,1-a+b,1.0/z);
    case 4: //15.3.8
	return (tgamma(c)*tgamma(b-a)/(tgamma(b)*tgamma(c-a)))*pow(1.0-z,-a)*G2F1_sum(a,c-b,a-b+1.0,1.0/(1.0-z)) + (tgamma(c)*tgamma(a-b)/(tgamma(a)*tgamma(c-b)))*pow(1-z,-b)*G2F1_sum(b,c-a,b-a+1.0,1.0/(1.0-z));
    case 5: //15.3.9
	return (tgamma(c)*tgamma(c-b-a)/(tgamma(c-a)*tgamma(c-b)))*pow(z,-a)*G2F1_sum(a,a-c+1,a+b-c+1,1.0 - 1.0/z) + (tgamma(c)*tgamma(a+b-c)/(tgamma(a)*tgamma(b)))*pow(1-z,c-a-b)*pow(z,a-c)*G2F1_sum(c-a,1-a,c-a-b+1,1.0 - 1.0/z);
  }

//shouldn't get here.
}

double G2F1_sum(double a, double b, double c, double z)
{
  double Snow = 0.0;
  double Sprev = 0.0;
  double test = 100.0;
  int n=0;

  while(test > 1e-8)
  {
	Sprev = Snow;
	Snow += (poch(a,n)*poch(b,n)/(poch(c,n)*tgamma(n+1)))*pow(z,(double)n);
	test = fabs((Snow-Sprev)/Snow);
	n++;
  }

return Snow;
}

// What follows are actually redundant in the current implementaton. But should be used in the later 
// version.

double AppellG2(double beta1, double beta2, double gamma1, double gamma2, double x, double y) 
{
  return pow(1+x,-beta1)*pow(1+y,-beta2)*AppellF2(1-gamma1-gamma2,beta1, beta2, 1-gamma1, 1-gamma2,x/(x+1), y/(y+1)); 
}

double AppellF2(double alpha, double beta1, double beta2, double gamma1, double gamma2, double x, double y)
{

  double sum = 0.0;
  int r = 0;
  double Snow =0.0;
  double Sprev =0.0;
  double test = 100.0;

std::cout<<"F2: x = "<<x<<", y = "<<y<<std::endl;

  while(test >1e-6)
  {
    Sprev=Snow;
    Snow += (poch(alpha,r)*poch(beta1,r)*poch(beta2,r)/(poch(gamma1,r)*poch(gamma2,r)*tgamma(r+1)))*pow(x*y,r)*G2F1(alpha+r,beta1+r,gamma1+2.0*r,x)*G2F1(alpha+r,beta2+r,gamma2+2.0*r,y);
    test = fabs((Snow-Sprev)/Snow);
    r++;
  }

return Snow;
}
