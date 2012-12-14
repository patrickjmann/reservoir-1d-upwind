//#############################################################################
/**\file Element.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Element.C,v $ 
$Revision: 1.8 $ 
$Date: 2012/11/25 00:39:23 $ 

\author P.J.Mann

Element class
*/
//#############################################################################
#include <iostream>
#include <fstream>
#include <iomanip>

#include <cstdlib>
#include <cmath>

#include "grid.h"
#include "proto.h"

using namespace std;

extern Fluid fluid;

//=============================================================================
Element::Element()
{
  cur.ZeroAll();
  old.ZeroAll();

  Exchange();  // load the "old" variables also.

  id = -1;
}
//=============================================================================
/// Very handy for testing the high-resolution methods.

void Element::ZeroCurrentSigma()
{
  cur.ZeroSigma();
}
//=============================================================================
/** load the primitive variables

Generally used for the initial data.

Note: this simply loads the primitive variables.  It does not make the
      evolution variables or the artificial viscosity.
 */
void Element::Make( const double u_in, const double rho_in )
{
  cur.u   = u_in;
  cur.rho = rho_in;
}
//======================================================================
void Element::MakePFromRho()
{
  cur.p = fluid.p( cur.rho );
}
void Element::MakeRhoFromP()
{
  cur.rho = fluid.rho( cur.p );
}
//======================================================================
/// Darcy's equations used to compute U from P

void Element::MakeUFromDarcy()
{
  const double ctmp = -(rock.kappa/fluid.mu);
  Element* eleft = left->left;
  Element* eright = right->right;
  if( eleft == 0 ){
      cur.u = ctmp * (eright->cur.p - cur.p)/(eright->x - x);
  } else if ( eright == 0 ){
      cur.u = ctmp * (cur.p - eleft->cur.p)/(x - eleft->x);
  } else {
    cur.u = ctmp * (eright->cur.p - eleft->cur.p)/(eright->x - eleft->x);
  }
}
//======================================================================
/** Make the fluid quantities rho, phi and p from W
 * See the documentation for the derivation of the equations.
*/
void Element::MakeFluidFromW()
{
  const double c = cur.W * fluid.c_f/(rock.c_r * rock.phi_r * fluid.rho_r);
  const double tmp1 = fluid.c_f/rock.c_r - 1.0;
  const double D = tmp1*tmp1 + 4.0*c;
  if( D < 0.0 ){
      cerr << "Element::MakeFluidFromW: ERROR: discriminant < 0\n"
           << "  D = " << D << '\n';
      exit(1);
  }
  cur.rho = fluid.rho_r * 0.5 * ( -tmp1 + sqrt(D) );
  
  // Now calculate the other quantities
  
  cur.phi = cur.W / cur.rho;
  cur.p = fluid.p( cur.rho ); 
}
//======================================================================
// Loops for various Make routines

void ElementList::MakeW()
{
  for( int ie=0; ie<n; ++ie ){
    e[ie].MakeW();
  }
}
//======================================================================
void ElementList::MakeUFromDarcy()
{
  for( int ie=0; ie<n; ++ie ){
    e[ie].MakeUFromDarcy();
  }
}
//======================================================================
void ElementList::MakeRhoFromP()
{
  for( int ie=0; ie<n; ++ie ){
    e[ie].MakeRhoFromP();
  }
}
//======================================================================
void ElementList::SetPhi()
{
  for( int ie=0; ie<n; ++ie ){
    e[ie].SetPhi();
  }
}
//======================================================================
void ElementList::MakeFluidFromW(){
  for( int i=0; i<n; ++i ){
    e[i].MakeFluidFromW();
  }
}
//======================================================================
void ElementList::MakeSigma( CommandLineOptions& cloption ){
  for( int i=0; i<n; ++i ){
    e[i].MakeSigma( cloption );
  }
}
void ElementList::Make( const int n_elements )
{
  if( n_elements < 4 ){
    cerr << "ElementList::Make: ERROR: n_elements = " << n_elements << '\n';
    exit(1);
  }

  n = n_elements;
  e = new Element[n];
}
//=============================================================================
/** Make the slopes of the basic variables.
These are used for the high-resolution approach to the advection terms.
*/
void Element::MakeSigma( CommandLineOptions& cloption )
{
  switch( cloption.limiter_type )
    {
    case MINMOD: MakeSigmaMinMod(); break;
    case MC: MakeSigmaMC(); break;
    case SUPERBEE: 
      cerr << "Element::MakeSigma: ERROR: SUPERBEE not supported yet.\n";
      exit(1);
      break;
    default:
      ZeroCurrentSigma();
      break;
    }
  //## ZeroCurrentSigma();
}
//-----------------------------------------------------------------
void Element::MakeSigmaMC()
{
  Element* eleft = left->left;
  Element* eright = right->right;
  
  if( eleft == 0 ){
    cur.ZeroSigma();
    return;
  }
  if( eright == 0 ){
    cur.ZeroSigma();
    return;
  }

  const double dx = eright->x - eleft->x;
  const double dxleft = x - eleft->x;
  const double dxright = eright->x - x;

  // Compute a minimum gradient.

  double dW_mid = (eright->cur.W - eleft->cur.W)/dx;
  double dW_left = (cur.W - eleft->cur.W)/dxleft;
  double dW_right = (eright->cur.W - cur.W)/dxright;
  double dW_mag = min( fabs(2.0*dW_left), fabs(dW_mid), fabs(2.0*dW_right) );
  cur.sigma_W = dW_mag * signum(dW_mid);
  if( dW_left*dW_right <= 0.0 ) cur.sigma_W = 0.0;

  double drho_mid = (eright->cur.rho - eleft->cur.rho)/dx;
  double drho_left = (cur.rho - eleft->cur.rho)/dxleft;
  double drho_right = (eright->cur.rho - cur.rho)/dxright;
  double drho_mag = min( fabs(2.0*drho_left), fabs(drho_mid), fabs(2.0*drho_right) );
  cur.sigma_rho = drho_mag * signum(drho_mid);
  if( drho_left*drho_right <= 0.0 ) cur.sigma_rho = 0.0;

  double du_mid = (eright->cur.u - eleft->cur.u)/dx;
  double du_left = (cur.u - eleft->cur.u)/dxleft;
  double du_right = (eright->cur.u - cur.u)/dxright;
  double du_mag = min( fabs(2.0*du_left), fabs(du_mid), fabs(2.0*du_right) );
  cur.sigma_u = du_mag * signum(du_mid);
  if( du_left*du_right <= 0.0 ) cur.sigma_u = 0.0;

  double dp_mid = (eright->cur.p - eleft->cur.p)/dx;
  double dp_left = (cur.p - eleft->cur.p)/dxleft;
  double dp_right = (eright->cur.p - cur.p)/dxright;
  double dp_mag = min( fabs(2.0*dp_left), fabs(dp_mid), fabs(2.0*dp_right) );
  cur.sigma_p = dp_mag * signum(dp_mid);
  if( dp_left*dp_right <= 0.0 ) cur.sigma_p = 0.0;

  double dphi_mid = (eright->cur.phi - eleft->cur.phi)/dx;
  double dphi_left = (cur.phi - eleft->cur.phi)/dxleft;
  double dphi_right = (eright->cur.phi - cur.phi)/dxright;
  double dphi_mag = min( fabs(2.0*dphi_left), fabs(dphi_mid), fabs(2.0*dphi_right) );
  cur.sigma_phi = dphi_mag * signum(dphi_mid);
  if( dphi_left*dphi_right <= 0.0 ) cur.sigma_phi = 0.0;
}
//-----------------------------------------------------------------
void Element::MakeSigmaMinMod()
{
  Element* eleft = left->left;
  Element* eright = right->right;

  if( eleft == 0 ){
    cur.ZeroSigma();
    return;
  }
  if( eright == 0 ){
    cur.ZeroSigma();
    return;
  }

  const double dx = eright->x - eleft->x;
  const double dxleft = x - eleft->x;
  const double dxright = eright->x - x;

  // Compute a minimum gradient.

  double dW_mid = (eright->cur.W - eleft->cur.W)/dx;
  double dW_left = (cur.W - eleft->cur.W)/dxleft;
  double dW_right = (eright->cur.W - cur.W)/dxright;
  double dW_mag = min( fabs(dW_left), fabs(dW_mid), fabs(dW_right) );
  cur.sigma_W = dW_mag * signum(dW_mid);
  if( dW_left*dW_right <= 0.0 ) cur.sigma_W = 0.0;

  double drho_mid = (eright->cur.rho - eleft->cur.rho)/dx;
  double drho_left = (cur.rho - eleft->cur.rho)/dxleft;
  double drho_right = (eright->cur.rho - cur.rho)/dxright;
  double drho_mag = min( fabs(drho_left), fabs(drho_mid), fabs(drho_right) );
  cur.sigma_rho = drho_mag * signum(drho_mid);
  if( drho_left*drho_right <= 0.0 ) cur.sigma_rho = 0.0;

  double du_mid = (eright->cur.u - eleft->cur.u)/dx;
  double du_left = (cur.u - eleft->cur.u)/dxleft;
  double du_right = (eright->cur.u - cur.u)/dxright;
  double du_mag = min( fabs(du_left), fabs(du_mid), fabs(du_right) );
  cur.sigma_u = du_mag * signum(du_mid);
  if( du_left*du_right <= 0.0 ) cur.sigma_u = 0.0;
 
  double dp_mid = (eright->cur.p - eleft->cur.p)/dx;
  double dp_left = (cur.p - eleft->cur.p)/dxleft;
  double dp_right = (eright->cur.p - cur.p)/dxright;
  double dp_mag = min( fabs(dp_left), fabs(dp_mid), fabs(dp_right) );
  cur.sigma_p = dp_mag * signum(dp_mid);
  if( dp_left*dp_right <= 0.0 ) cur.sigma_p = 0.0;

  double dphi_mid = (eright->cur.phi - eleft->cur.phi)/dx;
  double dphi_left = (cur.phi - eleft->cur.phi)/dxleft;
  double dphi_right = (eright->cur.phi - cur.phi)/dxright;
  double dphi_mag = min( fabs(dphi_left), fabs(dphi_mid), fabs(dphi_right) );
  cur.sigma_phi = dphi_mag * signum(dphi_mid);
  if( dphi_left*dphi_right <= 0.0 ) cur.sigma_phi = 0.0;
}
