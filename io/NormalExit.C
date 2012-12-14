//#############################################################################
/**\file NormalExit.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/NormalExit.C,v $ 
$Revision: 1.1 $ 
$Date: 2012/10/02 21:52:35 $ 

\author P.J.Mann

Graceful finish and clean up.

Revisions
---------

*/
//#############################################################################
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "putil.h"

#include "cpu.h"

#include "grid.h"
#include "ptime.h"
#include "proto.h"

using namespace std;

extern logstream plog;

//=============================================================================
void NormalExit( string msg, Grid& grid, TimeData& ptime, CpuClass& cpu )
{
  cpu.StopTimer();                   // "cpu" is global

  grid.PutSlice( ptime );
  
  // messages
  
  plog.SetLineStyle( '=' );
  plog.PutLine();
  plog << "NormalExit: " << msg;

  plog.s << cpu << '\n';
  cout << cpu << '\n';
  
  plog.close();
  exit(0);
}

