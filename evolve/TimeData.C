//#############################################################################
/**\file TimeData.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/TimeData.C,v $ 
$Revision: 1.6 $ 
$Date: 2012/11/25 00:39:23 $ 

\author P.J.Mann

Routines for the physical times
*/
//#############################################################################
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "pi.h"
#include "putil.h"

#include "grid.h"
#include "element.h"
#include "ptime.h"
#include "param.h"
#include "eqs.h"

#include "pexit.h"

using namespace std;

extern logstream plog;

//=============================================================================
/** Make the class

Parameters are read and data initialized.
*/
void TimeData::Make( string& param_file )
{
  // Get parameters from files

  ParameterSet param;
  param.Append( step_control.x, "time_step_control.x" );
  param.Append( step_control.c, "time_step_control.c" );
  param.Append( step_control.p1, "time_step_control.p1" );
  param.Append( step_control.p2, "time_step_control.p2" );
  param.Append( n_max, "evolve.n" );
  param.Append( coord_max, "tmax" );
  
  ifstream s;
  s.open( param_file.c_str() );
  if( !s ){
    cerr << "TimeData::Make: ERROR: unable to open \"" << param_file << "\"\n";
    exit(1);
  }
  param.Get( s );
  s.close();
  
  // Test the input
  
  if( n_max <= 0 ) {
    cerr << "TimeData: ERROR: n = " << n << endl;
    exit(1);
  }
  if( coord_max <= 0.0 ) {
    cerr << "TimeData: ERROR: coord_max = " << coord_max << endl;
    exit(1);
  }
  if( step_control.c <= 0.0 ){
    cerr << "TimeData: ERROR: step_control.c = "
         << step_control.c << endl;
    exit(1);
  }
  if( step_control.p1 <= 0.0 ){
    cerr << "TimeData: ERROR: step_control.p1 = "
         << step_control.p1 << endl;
    exit(1);
  }
  
  // Initialize basic times
  
  n = 0;
  
  coord = 0.0;
  
  // initial time step set to negative value initially for testing
  
  delta = -1.0;
  delta_min = -1.0;
  delta_max = -1.0;
}
//==============================================================================
/** Set min/max test values for the time step.
 * 
 * The time step must have been previously calculated somehow.
 * Usually the Courant condition solver will have been called during
 * initialization and then these min/max values can be set.
 * */
void TimeData::SetMinMax()
{
  if( delta <= 0.0 ){
    plog.s << "TimeData::SetMinMax: ERROR: delta still has initial value.\n"
           << "  delta = " << delta 
           << "  delta_min = " << delta_min
           << "  delta_max = " << delta_max << endl;
    ErrorExit( MAJOR_ERROR, "TimeData::SetMinMax: ERROR: delta still has initial value.\n" );
    exit(1);
  }
  delta_min = 0.01 * delta;
  delta_max = 100.0 * delta;
}
//==============================================================================
/// Test for any maximum allowed times

int TimeData::TestForMax()
{
  if( coord >= coord_max ) return OUT_OF_TIME;
  if( n     >= n_max     ) return OUT_OF_TIME_STEPS;
  return LOTS_OF_TIME;
}
//=============================================================================
/** Time step choice: (new time step put in dtnew)

Time step is controlled by:
<UL>
 <LI> sound speed courant limit
 <LI> velocity limit
 <LI> viscosity limit
</UL>

*/
void TimeData::MakeDelta( Grid& grid )
{
  // Maximum value given as multiple of previous dt and as multiple
  // of coordinate "freefall time".
  // Note: tunit is usually some characteristic crossing time.

  if( (step_control.p1 <= 0.0) ){
    plog << "TimeDate::MakeDelta: ERROR: incorrect parameters\n"
         << "  step_control.p1 = " << step_control.p1 << '\n';
    exit(1);
  }
  
  //delta_new = step_control.p1 * coord_max;
  delta_new = step_control.c * grid.dt_characteristic;
    
  // Courant limit and tests
  // go through each ELEMENT
  
  cs_max = 0.0;   // These are stored in the class.  Useful for error dumps.
  dx_max = 0.0;

  for( int ie=0; ie<grid.elist.n; ++ie ){
    
    Element* element = &(grid.elist.e[ie]);  // convenient, but this really should be in the element class
    
    const double umag = fabs(element->cur.u);
    const double dx = element->dx;
    
    // Courant limits
    
    const double cs = element->fluid.csound( element->cur.rho );
    if( cs > cs_max ) cs_max = cs;
    if( dx > dx_max ) dx_max = dx;
    
    // "Velocity" for Courant limit includes sound speed, fluid velocity
    // PJM: Nov 15, 2012:  Took out sound for testing
    
    const double courant_v = cs + umag;
    //const double courant_v = umag;
    
    if( delta_new * courant_v > step_control.c*dx ) {
      delta_new = step_control.c * dx / courant_v;
    }

    // Artificial Viscosity Limits

    const double q = 0.5 * ( element->left->q + element->right->q );
    const double artvis_v = sqrt( q / element->cur.rho );
    const double artvis_c = 0.25 * step_control.c;
    
    if( delta_new*artvis_v > artvis_c * dx ){
      delta_new = artvis_c * dx / artvis_v;
    }
  }
  
  // update the time step
  
  delta = delta_new;
}
//==============================================================================
/** Make tests as necessary on the time step.  Sometimes worthwhile
 * */
void TimeData::CheckDelta() throw(TimeException)
{
  if( delta <= delta_min ) {
    plog.s << "TimeData::CheckDelta: ERROR: small dt: delta <= delta_min\n"
	         << "   delta = " << delta_new << " s. "
	         << " delta_min = " << delta_min << " s.\n"
           << "   cs_max = " << cs_max << " m/s, dx_max = " << dx_max << " m\n";
    ErrorExit( MINOR_ERROR, "TimeData::CheckDelta: ERROR: small dt\n" );
    //exit(1);
    TimeException e;
    throw e;
  }
}
//==============================================================================
void TimeData::Update()
{
  ++n;
  coord += delta;
}
//=============================================================================
/// Put a readable summary of the time controls

void TimeData::PutParam( ofstream& s )
{
  s.precision(4);
  s.setf( ios::left, ios::adjustfield );
  
  s << GNU_COMMENT << " TimeData:: Control Parameters:\n"
    << GNU_COMMENT << "-------------------------------\n";
  
  s << GNU_COMMENT
    << " .step_control.x  = " << setw(7) << step_control.x << '\n'
    << GNU_COMMENT
    << " .step_control.c  = " << setw(7) << step_control.c << '\n'
    << GNU_COMMENT
    << " .step_control.p1 = " << setw(7) << step_control.p1 << '\n'
    << GNU_COMMENT
    << " .step_control.p2 = " << setw(7) << step_control.p2 << '\n';
  
  s << GNU_COMMENT << " coord_max = " << setw(8) << coord_max << '\n'
    << GNU_COMMENT << " n_max = " << setw(8) << n_max << endl;
}
//=============================================================================
ostream& operator << ( ostream& s, const TimeData& ptime )
{
  s.setf( ios::left, ios::adjustfield );
  s.precision( 6 );
  s << "t#" << ptime.n << " = " << ptime.coord << endl;
  return s;
}
