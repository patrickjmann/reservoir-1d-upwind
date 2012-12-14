//#############################################################################
/**\file basic.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/basic.h,v $ 
$Revision: 1.4 $ 
$Date: 2012/10/21 18:32:43 $ 

\author P.J.Mann

\brief various basic quantities

*/
//#############################################################################
#if ! defined CONS1_BASIC_H
#define CONS1_BASIC_H

using namespace std;

extern int verbose;

//===================================================================
/** signum: returns the sign of x

\retval -1  x<0
\retval 0   x=0
\retval 1   x>0
*/

inline double signum( double x )
{
  if( x<0.0 ){
    return -1.0;
  } else if( x==0.0 ){
    return 0.0;
  } else if( x>0.0 ){
    return 1.0;
  } else {
    return x;    /* must be that x is Not-a-Number */
  }
}
//=============================================================================
inline double min( double x1, double x2, double x3 )
{
  return MIN( MIN(x1, x2), x3 );
}
inline double max( double x1, double x2, double x3 )
{
  return MAX( MAX(x1, x2), x3 );
}
//===================================================================
/// Not sure if this will be used.  Probably not!
struct ArtVisControl
{
  double k1,k2;
};

class Element;
//=============================================================================
/** A structure for the basic evolution variables.
 * 
 * Generally these are kept on two time slices for the evolution.
 * Some however are calculated (Darcy's law for instance) rather than evolved.
 * */

class BaseVariables
{
 public:
  double W, sigma_W;     /** Porous medium density and slope */
  double rho, sigma_rho; /** Fluid density and slope */
  double u, sigma_u;     /** Fluid velocity and slope */
  double p, sigma_p;     /** pressure and slope */
  double phi, sigma_phi; /** porosity and slope.  This is cell-dependent! */

  void MakeHighResConstant( Element*, double x  );
  void MakeHighResLinear( Element*, double x  );
  void MakeHighResCubicHermite( Element*, Element*, double x );
  
  inline void CopyFrom( BaseVariables& b ){
    W = b.W;  sigma_W = b.sigma_W;
    rho = b.rho; sigma_rho = b.sigma_rho;
    u = b.u; sigma_u = b.sigma_u;
    p = b.p; sigma_p = b.sigma_p;
    phi = b.phi; sigma_phi = b.sigma_phi;
  }
   inline void ZeroSigma(){
    sigma_W = 0.0;
    sigma_rho = 0.0;
    sigma_u = 0.0;
    sigma_p = 0.0;
    sigma_phi = 0.0;
  };
  inline void ZeroBase(){
    W = 0.0; rho = 0.0; u = 0.0; p = 0.0; phi = 0.0;
  };
  inline void ZeroAll(){
    ZeroBase(); ZeroSigma();
  }
};

inline ostream& operator << (ostream& s, const BaseVariables& b )
{
  s << "BaseVariables:\n"
    << "  W   = " << b.W   << " sigma_W   = " << b.sigma_W << '\n'
    << "  u   = " << b.u   << " sigma_u   = " << b.sigma_u << '\n'
    << "  rho = " << b.rho << " sigma_rho = " << b.sigma_rho << '\n'
    << "  p   = " << b.p   << " sigma_p   = " << b.sigma_p << '\n'
    << "  phi = " << b.phi << " sigma_phi = " << b.sigma_phi;
  return s;
};
//===================================================================
/// Flux class for nodes

class Flux
{
 public:
  double W;   /// The W component of the flux.

  inline void Make( BaseVariables& b, double q ){
    W = b.rho*b.u;
  }
  inline void Make( double aRho, double au ){
    W = aRho*au;
  }
};
inline ostream& operator << ( ostream& s, const Flux& f )
{
  s << "Flux: .W = " << f.W;
  return s;
};
#endif
