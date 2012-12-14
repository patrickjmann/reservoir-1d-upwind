//##############################################################################
/**\file Boundary.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Boundary.C,v $
$Revision: 1.2 $
$Date: 2012/11/23 22:54:08 $

\author P.J.Mann: November, 2012

Compute fluxes on boundary nodes
*/
//##############################################################################
#include <iostream>
#include <iomanip>
#include <cmath>

#include "proto.h"
#include "cloption.h"
#include "grid.h"
#include "element.h"
#include "ptime.h"
#include "boundary.h"

using namespace std;

//==============================================================================
/** A utility function where I can put all the actual boundary
flux calculations, but still call it up for the old (predictor)
time slice and the new (corrector) time slice.

\param[in] cloption Options
\param[in] bdy_cond Boundary Conditions class
\param[in] ptime Time slice data usually for input into boundary conditions
\param[in] sign  Sign of v*deltat/2 term.  This is required
             for the forward or backward differencing term.
             forward:  sign=-1
             backward: sign=+1
\param[out] flux  output values of the calculated flux
*/
void Node::MakeFluxBoundary( CommandLineOptions& cloption,
          BoundaryCondition& bdy_cond, TimeData& ptime, int sign,
          Flux& flux )
{
  // Identify the inside element and sign for a generic calculation of
  // of the boundary conditions (left or right).
  
  Element* inside;
  double usign;
  double rho_bdy;
  if( left == 0 ){  // left-side boundary
    inside = right; // inside element is to the right
    usign = 1.0;    // test if u >= 0.0 for incoming
    rho_bdy = bdy_cond.RhoLeft( ptime.Coord() );
  } else if( right == 0 ){ // right-side boundary
    inside = left;
    usign = -1.0;
    rho_bdy = bdy_cond.RhoRight( ptime.Coord() );
  } else {  // not a boundary node!
    cerr << "Node::MakeFluxBoundary: ERROR: this is not a boundary node!\n";
    exit(1);
  }
  
  // Calculate the highres values from the inside node
  // BaseVariables is used in Element but is handy to use here although
  // various element-specific members are unused.

  BaseVariables highres_inside;
  
  if( cloption.reconstruction_type == CONSTANT ){
    highres_inside.MakeHighResConstant( inside, x );
  } else if( cloption.reconstruction_type == LINEAR ){
    highres_inside.MakeHighResLinear( inside, x );
  } else {
    cerr << "Node::MakeLeftFluxBoundary: ERROR: unknown reconstruction type.\n"
         << "      cloption.reconstruction_type=" << cloption.reconstruction_type << '\n';
    exit( 1 );
  }

  //============================================================================
  /* Now calculate the boundary conditions.
   * Velocity will be just interpolated from within as usual.  This
   * just means using the highres_inside value.  It is NOT being advected
   * so I think this is reasonable.  Just using Darcy's law which is
   * a constraint equation.
   * 
   * Density should be given at the boundary if the velocity is inward.
   * May be more "physical" to give pressure and then calculate density.
   * For this test give density.
   * 
   * If outgoing then just use the upwinded high-res values from the
   * inside element as usual.
   * */

  double u, rho;
  if( highres_inside.u*usign >= 0.0 ){
    u = highres_inside.u;   // High-res should be good interpolation
    rho = rho_bdy;          // The actual boundary condition
  } else {                  // Incoming so nothing special needed
    u = highres_inside.u;
    rho = highres_inside.rho;
  }

  flux.Make( rho, u );
  return;
}
