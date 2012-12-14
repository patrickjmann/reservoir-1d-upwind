//#############################################################################
/** GetId.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/GetId.C,v $ 
$Revision: 1.1 $ 
$Date: 2012/10/02 21:52:35 $ 

Get run id and restart numbers from files.

Input
-----
newrun     ! RESTART or NEWRUN

Returns    ! id number of the run
-------

Output
------
irnum      ! restart number of run

Pretty straightforward here.

If this is a new run then get the id number
from a system data file in the home directory.  Then update the id number
and put it back in the file.  Also output id number to a restart file
kept locally.

If this is a restart then just read the local restart file.
*/
//#############################################################################
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>

#include "putil.h"

using namespace std;

//=============================================================================
int GetId( )
{
  int run_id;
  char run_id_file[100];
  
  char* home = getenv( "HOME" );
  strcpy( run_id_file, home );
  strcat( run_id_file, "/.rshock_" );
  strcat( run_id_file, PACKAGE );
  
  ifstream s_in( run_id_file );
  if( !s_in ){
    run_id = 0;
  } else {
    s_in >> run_id; s_in.ignore(100,'\n');
  }
  s_in.close();
  
  run_id++;
  
  ofstream s_out( run_id_file );
  s_out << run_id << "   # run id number\n";
  s_out.close();
  
  run_id += HostNameId();      // Adds a host name number to the id
  return run_id;
}
