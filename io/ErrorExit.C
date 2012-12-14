//#############################################################################
/**\file ErrorExit.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/ErrorExit.C,v $ 
$Revision: 1.1 $ 
$Date: 2012/10/02 21:52:35 $ 

\author P.J.Mann

Error finish and clean up.

Input
-----
msg  ! message from calling program
isev ! severity level

Contents
--------
perror   ! error exit: ascii dump of state of run (BIG!!!)
puterr   ! error message

*/
//#############################################################################
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "pexit.h"
#include "ptime.h"
#include "cpu.h"
#include "proto.h"

using namespace std;

extern logstream plog;
extern CpuClass cpu;
extern TimeData ptime;

//=============================================================================
void dump_global_stuff( string msg )
{
  plog.SetLineStyle( '=' );
  plog.PutLine();
  plog << msg << '\n';
  plog.PutLine();
  plog.s << ptime;
  plog.PutLine();
  plog.s << cpu;
}
//=============================================================================
void dump_readable( Grid& grid )
{
  grid.PutTable( "grid_error_dump" );
}
//=============================================================================
/// Standard exit on error with detailed info about grid

void ErrorExit( int severity, string msg, Grid& grid )
{
  cpu.StopTimer();                     // closing cpu time
  grid.PutSlice( ptime );              // last slice
  
  dump_global_stuff( msg );
  
  plog.SetLineStyle( '=' );
  plog.PutLine();
  plog.s << "ErrorExit: " << msg << '\n';
  
  plog.close();
  
  exit(1);
}
//=============================================================================
void ErrorExit( int severity, string msg, Grid& grid, Element& element_error )
{
  cpu.StopTimer();                     // closing cpu time
  grid.PutSlice( ptime );              // last slice
  
  plog.SetLineStyle( '-' );
  plog.PutLine();
  plog.s << "ErrorExit: " << msg << '\n';
  
  plog.PutLine();
  plog.s << "The error occurred on the following element:\n"
	 << element_error;
  
  plog.close();
  
  exit(1);
}
//=============================================================================
/** Standard exit on error with just a message.
    No other info available here so severity is ignored.
*/
void ErrorExit( int severity, string msg )
{
  cpu.StopTimer();
  dump_global_stuff( msg );
  plog.close();
  exit(1);
}
//=============================================================================
void ErrorExit( int severity, string msg, Element& element_error )
{
  cpu.StopTimer();
  dump_global_stuff( msg );
  
  plog.SetLineStyle( '=' );
  plog.PutLine();
  plog.s << "The error occurred on the following element:\n"
	 << element_error;
  
  plog.close();
  exit(1);
}
