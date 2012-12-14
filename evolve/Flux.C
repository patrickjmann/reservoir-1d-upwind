//######################################################################
/**\file Flux.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Flux.C,v $
$Revision: 1.8 $
$Date: 2012/11/25 00:39:23 $

\author P.J.Mann: July, 2006

Compute fluxes on nodes (interface between elements) and
add to adjacent elements.
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

using namespace std;

//======================================================================
void NodeList::MakeFluxPredictor( CommandLineOptions& cloption, TimeData& ptime,
				  ArtVisControl& artvis, BoundaryCondition& bdy_cond )
{
  for( int inode=1; inode<n-1; ++inode ){
    node[inode].MakeFluxPredictor( cloption, ptime, artvis );
  }
  int sign = -1;
  node[0].MakeFluxBoundary( cloption, bdy_cond, ptime, sign, node[0].oldflux );
  node[n-1].MakeFluxBoundary( cloption, bdy_cond, ptime, sign, node[n-1].oldflux );
}
//======================================================================
/** Make the flux on each node for the corrector steps
 * \param[in] cloption Various options
 * \param[in] ptime time step data
 * \param[in] artvis Artificial viscosity control
 * */
void NodeList::MakeFluxCorrector( CommandLineOptions& cloption, TimeData& ptime,
				  ArtVisControl& artvis, BoundaryCondition& bdy_cond )
{
  for( int inode=1; inode<n-1; ++inode ){
    node[inode].MakeFluxCorrector( cloption, ptime, artvis );
  }
  int sign = +1;
  node[0].MakeFluxBoundary( cloption, bdy_cond, ptime, sign, node[0].flux );
  node[n-1].MakeFluxBoundary( cloption, bdy_cond, ptime, sign, node[n-1].flux );
}
//======================================================================
/** A utility function where I can put all the actual
flux calculations, but still call it up for the old (predictor)
time slice and the new (corrector) time slice.

\param[in] left  BaseVariables on left side of interface
\param[in] right BaseVariables on right side of interface
\param[in] leftx left element position
\param[in] rightx right element position
\param[in] artvis Artificial Viscosity control parameters
\param[in] sign  Sign of v*deltat/2 term.  This is required
             for the forward or backward differencing term.
             forward:  sign=-1
             backward: sign=+1
\param[out] flux  output values of the calculated flux
*/
void Node::MakeFlux( CommandLineOptions& cloption,
          BaseVariables& left, BaseVariables& right,
          double leftx, double rightx,
          TimeData& ptime, ArtVisControl& artvis, int sign,
          Flux& flux )
{
  //const double deltat = ptime.Delta();

  // Averages are used for various tests and upwind fallbacks,
  // so compute them all here. 

  const double u = 0.5 * ( left.u + right.u );
  //const double rho = 0.5 * ( left.rho + right.rho );

  const double du = right.u - left.u;

  // HIGH RESOLUTION FITTING
  // -----------------------
  // BaseVariables is actually used for the Element class
  // (current and old) but it's handy to use it here even
  // though other element-specific members are unused.

  BaseVariables highres_left, highres_right;
  MakeHighRes( cloption.reconstruction_type, highres_left, highres_right );

  //====================================================================
  // Simple u-based Upwind.

  double u_up, rho_up;
  if( u >= 0 ){
    u_up = highres_left.u;
    rho_up = highres_left.rho;
  } else {
    u_up = highres_right.u;
    rho_up = highres_right.rho;
  }

  // Artificial viscosity.  Necessary if straight upwind is used.

  if( du < 0.0 ){
    const double absdu = fabs(du);
    //q =  rho_up * absdu * ( artvis.k1*cs + artvis.k2*absdu );
    q =  artvis.k2 * rho_up * absdu * absdu;
      } else {
    q = 0.0;
  }

  flux.Make( rho_up, u_up );
  return;
}
//======================================================================
/** Compute flux for predictor.  Just a wrapper which calls the
 * generic MakeFlux function.
 * Note: should probably be moved into the "Flux" or "BasicVariable" 
 * classes!
 * The flux is computed and stored in the "oldflux" member.  It will
 * be used for the predictor and then for the corrector averages.
*/
void Node::MakeFluxPredictor( CommandLineOptions& cloption, TimeData& ptime, ArtVisControl& artvis )
{
  if( left == 0 ){
    // Do nothing.  Handled explicitly.
    return;
  }
  if( right == 0 ){
    // Do nothing.  Handled explicitly.
    return;
  }

  const int sign = -1;  // forward differencing (predictor)
  MakeFlux( cloption, left->old, right->old, left->x, right->x,
	    ptime, artvis, sign, oldflux );
}
//======================================================================
/** Compute flux corrector.  Just a wrapper which calls the generic
 * MakeFlux function.
 * 
 * Note: should probably be moved into the "Flux" or "BasicVariable"
 * classes!
 * 
 * The flux is put in the "flux" member, and is an average of the old
 * and new fluxes.
*/
void Node::MakeFluxCorrector( CommandLineOptions& cloption, TimeData& ptime, ArtVisControl& artvis )
{
  if( left == 0 ) return;
  if( right == 0 ) return;

  const int sign = +1;   // backward differencing (corrector)
  Flux curflux;

  // Compute fluxes

  MakeFlux( cloption, left->cur, right->cur, left->x, right->x,
	    ptime, artvis, sign, curflux );  

  flux.W = 0.5 * ( curflux.W + oldflux.W );
}
