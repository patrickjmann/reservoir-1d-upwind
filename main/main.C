//##############################################################################
/**\file main.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/main/RCS/main.C,v $
$Revision: 1.6 $
$Date: 2012/11/25 00:39:16 $
$Author: mann $

\author P.J.Mann

Original: Sep.19, 2012

1d Single Phase flow through porous media
*/
//##############################################################################
// System includes

#include <iostream>
#include <fstream>
#include <iomanip>

#include <cstdlib>
#include <cmath>

#include <cstdio>   // For writing to a string

// My utilities with such things as binary i/o classes,
// 1d shape functions, sorting routines, ...

#include "putil.h"        // basic utilities
#include "nutil.h"        // numerical utilities
#include "param.h"        // Parameter input from files

// Local includes.

#include "grid.h"
#include "param.h"	
#include "cpu.h"
#include "proto.h"
#include "pexit.h"
#include "ptime.h"
#include "cloption.h"

using namespace std;

//==============================================================================
// Global classes/structures needed everywhere.

//Fluid fluid;               /// Moved to element. Equation of state
//LooseSand rock;            /// Moved to element. Rock properties
TimeData ptime;            /// various (physical) times.

CpuClass cpu;              /// Cpu times

logstream plog;            /// Run log

//==============================================================================
/** Local class to handle dump times during the evolution.

  Tests various quantities (interval since last dump,
  proper time elapsed, ...) and sets a flag if a dump is indicated.
*/
class LastDump
{
public:
  enum dump_returns { DUMP = 0, NODUMP };
  
private:
  int interval;        // number of time steps since last dump
  double time;         // time of last dump
  int interval_max;    // max number of times steps between dumps
  double time_max;     // max time at centre between dumps
  
public:
  LastDump( TimeData& );
  inline void Reset( TimeData& ptime ){
    interval = 0; time = ptime.Coord();
  };
  inline void Update(){ interval++; };
  int Test( TimeData& ptime );
};

LastDump::LastDump( TimeData& ptime )
{
  ParameterSet param;
  double epsdmp;
  param.Append( epsdmp, "epsdmp" );
  param.Append( interval_max, "dump.n" );
  ifstream s( "param.dat" );
  param.Get( s );
  s.close();
  
  time_max = epsdmp * ptime.CoordMax();   // time between dumps
  
  interval = 0; time = 0.0;
};

int LastDump::Test( TimeData& ptime )
{
  if( interval >= interval_max ) return DUMP;
  if( fabs(ptime.Coord() - time) > time_max ) return DUMP;
  return NODUMP;
};

//#############################################################################
//#############################################################################
//                         The Evolution   
//#############################################################################
//#############################################################################
/** Top level evolution */

void Evolve( TimeData& ptime, CommandLineOptions& cloption,
	     Grid& grid )
{
  // "LastDump" class decides on frequency of complete time slice dumps
  
  LastDump last_dump( ptime );

  grid.PutSlice( ptime );            // Dump the initial slice
  
  // Use the first time step to set various max/min values
  
  ptime.MakeDelta( grid );
  ptime.SetMinMax();  // Set test values from the time step estimated above

  //===========================================================================
  // Evolve
  
  plog.SetLineStyle( '=' );
  plog.PutLine();
  plog << "Evolve: INFO: beginning the evolution..\n";
  
  cpu.StartTimer();
  for( int n_time=1; n_time<=ptime.NMax(); n_time++ ){
    //--------------------------------------------------------------------------
    ptime.MakeDelta( grid );
    try {
      ptime.CheckDelta();     // not necessary, but useful sometimes
    } catch ( TimeException& e ) {
      cerr << e.what() << '\n';
      exit(1);
    }
    grid.Evolve( cloption, ptime );
    ptime.Update();   // update all times
    
    if( verbose ){
      plog << "main: INFO: " << n_time
          << " tcoord = " << ptime.Coord()
          << " dt = " << ptime.Delta()
          << '\n';
    } else {
      if( n_time % 50 == 0 ){
        plog << "main: INFO: " << n_time
            << " tcoord = " << ptime.Coord()
            << " dt = " << ptime.Delta()
            << '\n';
      }
    }     

    // Check for end of run (physical or cpu time) reached.
    
    switch( ptime.TestForMax() ){
    case TimeData::OUT_OF_TIME:
      NormalExit( "main: Normal Exit: time\n", grid, ptime, cpu );
      exit(0);
      break;
    case TimeData::OUT_OF_TIME_STEPS:
      NormalExit( "main: Normal Exit: time step limit\n", grid, ptime, cpu );
      exit(0);
      break;
    }
    
    // time step too small?
    
    if( ptime.Delta() < ptime.DeltaMin() ) {
      ErrorExit( MAJOR_ERROR,"main: ERROR: time step too small\n", grid );
      exit(1);
    }
    
    // Check to see if a time slice dump is due.
    
    last_dump.Update();
    if( last_dump.Test( ptime ) == LastDump::DUMP ){
      grid.PutSlice( ptime );
      last_dump.Reset( ptime );
    }
  }   // end of evolution loop
}

//#############################################################################
//#############################################################################
//                         The Main Program
//#############################################################################
//#############################################################################
/** Does all the setup, gets initial data, and call the evolver */

int main ( const int argc, char *argv[] )
{
  // command line

  CommandLineOptions cloption;   // (sets defaults)
  
  cloption.Get( argc, argv );
  
  // id numbers and strings
  
  int run_id = GetId();
  char* host_name = GetHostName();
  
  char run_id_string[200];
  sprintf( run_id_string, "%s: #%d, %s {%s}",
	   PROGRAM_NAME, run_id, DateTime(), host_name );
  
  if( verbose ) cerr << "main: INFO: Make the output files\n";
  MakeOutFiles();
  
  // Parameter file
  
  string param_file( "param.dat" );
  
  // time control

  if( verbose) cerr << "main: INFO: ptime.Make()..\n";
  ptime.Make( param_file);
  
  // Initialize the basic classes
  
  plog << "main: INFO: New Run..\n";
  
  // Dump run info
  // (usually pasted into the system LOG file after a run)
  
  plog.SetLineStyle( '#' );
  plog.PutLine();
  plog << run_id_string << '\n';
 
  Grid grid( cloption, run_id_string, run_id );

  // Initialize
  // NOTE: initial data is read by the Grid class
  
  plog.SetLineStyle( '=' );
  plog.PutLine();
  plog << "main: INFO: MakeInitial..\n";
  grid.MakeInitial( cloption );
  
  // Dump the parameters for reference
  
  plog.SetLineStyle( '=' );

  plog.PutLine();
  grid.PutParam( plog.s );
  
  plog.PutLine();
  cloption.Put( plog.s );

  plog.PutLine();
  ptime.PutParam( plog.s );
  
  // Evolve
  
  plog << "main: INFO: Evolve..\n";
  Evolve( ptime, cloption, grid );
  
  // dump the last slice to a standard final configuration file
  
  grid.PutSlice( ptime );
}
