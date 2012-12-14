//#############################################################################
// $Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/MakeOutFiles.C,v $ 
// $Revision: 1.2 $ 
// $Date: 2012/11/12 18:33:49 $ 

// P.J.Mann

// Construct names of output files and open them.
//#############################################################################
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "putil.h"

using namespace std;

extern logstream plog;

//=============================================================================
void MakeOutFiles()
{
  plog.open( "log.out" );
  plog.s.precision( 4 );
  plog.s.setf( ios::left, ios::adjustfield );
}
