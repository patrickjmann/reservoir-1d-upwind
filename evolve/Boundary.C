//######################################################################
/**\file Boundary.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Boundary.C,v $
$Revision: 1.2 $
$Date: 2012/11/23 22:54:08 $

\author P.J.Mann: November, 2012

Compute fluxes on boundary nodes
*/
//######################################################################
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

//======================================================================
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
  // Use the high-resolution fitting approach
  // -----------------------
  // BaseVariables is actually used for the Element class
  // (current and old) but it's handy to use it here even
  // though other element-specific members are unused.
  
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
    cerr << "Node::MakeFluxBoundary: ERROR: this is not a boundary element!\n";
    exit(1);
  }

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

  //====================================================================
  // Simple u-based Upwind.
  /*
   * Velocity will be just interpolated from within as usual.  This
   * just means using the highres_right value.  It is NOT being advected
   * so I think this is reasonable.  Just using Darcy's law which is
   * a constraint equation.
   * 
   * Density should be given at the boundary.  May be more "physical" to
   * give pressure and then calculate density.  For this test give
   * density
   * 
   * If outgoing (v<0) then just use the upwinded values from the
   * adjacent node as usual.
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
