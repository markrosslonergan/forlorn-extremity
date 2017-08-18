#ifndef COMPLEX_FUNCTIONS_H_
#define COMPLEX_FUNCTIONS_H_

#include <complex>
#include <cmath>
#include <iostream>

#define SIGN(a) (((a) < 0) ? (-1) : (1))

using namespace std;

// Infinite norm of a complex number.
// ----------------------------------
// It is max(|Re[z]|,|Im[z]|)

inline double inf_norm (const complex<double> &z)
{
  return max (abs (real (z)),abs (imag (z)));
}

// Test of finiteness of a complex number
// --------------------------------------
// If real or imaginary parts are finite, true is returned.
// Otherwise, false is returned

inline bool isfinite (const complex<double> &z)
{
  const double x = real (z), y = imag (z);

  return (finite (x) && finite (y));
}

// Usual operator overloads of complex numbers with integers
// ---------------------------------------------------------
// Recent complex libraries do not accept for example z+n or z==n with n integer, signed or unsigned.
// The operator overload is done here, by simply putting a cast on double to the integer.

inline complex<double> operator + (const complex<double> &z,const int n)
{
  return (z+static_cast<double> (n));
}

inline complex<double> operator - (const complex<double> &z,const int n)
{
  return (z-static_cast<double> (n));
}

inline complex<double> operator * (const complex<double> &z,const int n)
{
  return (z*static_cast<double> (n));
}

inline complex<double> operator / (const complex<double> &z,const int n)
{
  return (z/static_cast<double> (n));
}

inline complex<double> operator + (const int n,const complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline complex<double> operator - (const int n,const complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline complex<double> operator * (const int n,const complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline complex<double> operator / (const int n,const complex<double> &z)
{
  return (static_cast<double> (n)/z);
}








inline complex<double> operator + (const complex<double> &z,const unsigned int n)
{
  return (z+static_cast<double> (n));
}

inline complex<double> operator - (const complex<double> &z,const unsigned int n)
{
  return (z-static_cast<double> (n));
}

inline complex<double> operator * (const complex<double> &z,const unsigned int n)
{
  return (z*static_cast<double> (n));
}

inline complex<double> operator / (const complex<double> &z,const unsigned int n)
{
  return (z/static_cast<double> (n));
}

inline complex<double> operator + (const unsigned int n,const complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline complex<double> operator - (const unsigned int n,const complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline complex<double> operator * (const unsigned int n,const complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline complex<double> operator / (const unsigned int n,const complex<double> &z)
{
  return (static_cast<double> (n)/z);
}




inline bool operator == (const complex<double> &z,const int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const complex<double> &z,const int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const int n,const complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const int n,const complex<double> &z)
{
  return (static_cast<double> (n) != z);
}







inline bool operator == (const complex<double> &z,const unsigned int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const complex<double> &z,const unsigned int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const unsigned int n,const complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const unsigned int n,const complex<double> &z)
{
  return (static_cast<double> (n) != z);
}






complex<double> expm1 (const complex<double> &z);
complex<double> log1p (const complex<double> &z);
complex<double> log_Gamma (const complex<double> &z);
complex<double> Gamma_inv (const complex<double> &z);
complex<double> sigma_l_calc (const complex<double> &l,const complex<double> &eta);
complex<double> log_Cl_eta_calc (const complex<double> &l,const complex<double> &eta);
complex<double> log_cut_constant_AS_calc (const int omega,const complex<double> &l,const complex<double> &eta);
complex<double> log_cut_constant_CFa_calc (const bool is_it_normalized,const int omega,const complex<double> &l,const complex<double> &eta);
complex<double> log_cut_constant_CFb_calc (const bool is_it_normalized,const int omega,const complex<double> &l,const complex<double> &eta);
complex<double> sin_chi_calc (const complex<double> &l,const complex<double> &eta);
complex<double> exp_I_omega_chi_calc (const int omega,const complex<double> &l,const complex<double> &eta);

#endif
