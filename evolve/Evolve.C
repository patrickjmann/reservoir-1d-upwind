//#############################################################################
/**\file Evolve.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Evolve.C,v $ 
$Revision: 1.4 $ 
$Date: 2012/11/23 22:54:08 $ 

\author P.J.Mann: June, 2006

Evolve
*/
//#############################################################################
#include <cmath>
#include <cstdlib>

#include "putil.h"

#include "grid.h"
#include "ptime.h"

extern logstream plog;

//=============================================================================
void Grid::Evolve( CommandLineOptions& cloption, TimeData& ptime )
{
  //double dt = ptime.delta;

  corrector_level = 0;

  // exchange old/current data (only the elements need this)

  elist.Exchange();

  // Predict
  // Note:  MakeFlux uses "cur" values.  Therefore the "cur" values need to have
  // the old time slice values at this stage.  This has been forced by
  // the previous "Exchange" so both cur and old have the same values.
  // Means that "MakeFlux" only needs to address the "cur" values, even on the predictor step

  nlist.MakeFluxPredictor( cloption, ptime, artvis, bdy_cond );
  elist.Predict( cloption, ptime );
  elist.MakeFluidFromW();
  elist.MakeUFromDarcy();
  elist.MakeSigma( cloption );

  // Correct

  for( corrector_level=1; corrector_level<=n_corrector; ++corrector_level ){
    nlist.MakeFluxCorrector( cloption, ptime, artvis, bdy_cond );
    elist.Correct( cloption, ptime );
    elist.MakeFluidFromW();
    elist.MakeUFromDarcy();
    elist.MakeSigma( cloption );
  }

  // There can be oscillations and difficulties
  // so fixup as necessary here.
}
//==================================================================
/** Predictor: Evolve the element quantities using the fluxes at the
interfaces.
*/
void ElementList::Predict( CommandLineOptions& cloption, TimeData& ptime )
{
  for( int i=0; i<n; ++i ){
    e[i].EvolveSimplePredictor( cloption, ptime );
  }
}
//==================================================================
/** Corrector: Evolve the element quantities using the fluxes at the
interfaces.
*/
void ElementList::Correct( CommandLineOptions& cloption, TimeData& ptime )
{
  for( int i=0; i<n; ++i ){
    e[i].EvolveSimpleCorrector( cloption, ptime );
  }
}
//==================================================================
/** Predictor

Generally a simple 1-st order forward differencing

\param[in] cloption options and parameters
\param[in] ptime time step data
*/
void Element::EvolveSimplePredictor( CommandLineOptions& cloption, TimeData& ptime )
{
  const double dtdx  = ptime.Delta()/H();
  cur.W = old.W - dtdx * ( right->oldflux.W - left->oldflux.W );
}
//==================================================================
/** Corrector

Note that all the time-averaging for the corrector is built-in to
the Node::flux members.

\param[in] cloption options and parameters
\param[in] ptime time step data
*/
void Element::EvolveSimpleCorrector( CommandLineOptions& cloption, TimeData& ptime )
{
  const double dtdx  = ptime.Delta()/H();
  cur.W = old.W - dtdx * ( right->flux.W - left->flux.W );
}

