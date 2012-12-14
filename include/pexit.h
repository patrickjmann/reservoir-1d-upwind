//#############################################################################
/**\file pexit.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/pexit.h,v $ 
$Revision: 1.1 $ 
$Date: 2012/10/02 21:52:01 $ 

\author P.J.Mann

Definitions for exit routines
*/
//#############################################################################

#if ! defined CONS1_PEXIT_H
#define CONS1_PEXIT_H

#include "grid.h"

using namespace std;

#define MINOR_ERROR 1    // Basic summary dump only
#define MAJOR_ERROR 3    // Complete dump of Grid, ElementList, ...

void NormalExit( string msg, Grid& grid, TimeData& tdata, CpuClass& );
/*inline void NormalExit( string msg, Grid& grid, TimeData& tdata ){
  NormalExit( msg.c_str(), grid, tdata );
}
*/

void ErrorExit( int severity, string msg );
/*inline void ErrorExit( int severity, string msg ){ ErrorExit(severity, msg.c_str()); };

void ErrorExit( int severity, char* msg, Node& error_node );
inline void ErrorExit( int severity, string msg, Node& error_node ){
  ErrorExit(severity, msg.c_str(), error_node );
}
*/

void ErrorExit( int severity, string msg, Element& error_element );
/*inline void ErrorExit( int severity, const string msg, Element& error_element ){
  ErrorExit(severity, msg.c_str(), error_element );
}
*/

void ErrorExit( int severity, string msg, Grid& grid );
/*
inline void ErrorExit( int severity, const string msg, Grid& grid ){
  ErrorExit(severity, msg.c_str(), grid );
}
*/
#endif
