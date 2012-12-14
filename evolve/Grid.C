//##############################################################################
/*\file Grid.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Grid.C,v $ 
$Revision: 1.3 $ 
$Date: 2012/11/25 00:39:23 $ 

\author P.J.Mann

constructors/destructors/... for grid classes
*/
//##############################################################################
#include <iostream>
#include <fstream>
#include <iomanip>

#include <cstdlib>
#include <cstring>

#include <cmath>

#include "ptime.h"
#include "grid.h"
#include "param.h"

using namespace std;

//==============================================================================
/** Sets defaults, then reads changed parameters from files.

Note:  This constructor simply creates the class and
       puts in the basic parameters.

       The initial data is NOT set up.  See initial data calls.
*/
Grid::Grid(  CommandLineOptions& cloption,
	     char* run_id_string_in, int run_id_in )
{
  const int MAX_N_STRING = 200;
  run_id_string = new char[MAX_N_STRING];
  strncpy( run_id_string, run_id_string_in, (MAX_N_STRING-1) );
  run_id = run_id_in;  

  // The defaults: lots of them: see "grid.h" for explanations
  
  artvis.k1 = 0.0;
  artvis.k2 = 0.0;
  
  corrector_level = -1;
  n_corrector = 1;

  // Make a generic parameter list
    
  ParameterSet param;

  param.Append( artvis.k1, "artvis.k1" );
  param.Append( artvis.k2, "artvis.k2" );

  param.Append( n_corrector, "corrector.n" );

  // Get parameters from standard data file

  ifstream s_in( "param.dat" );
  param.Get( s_in );
  s_in.close();
}
//==============================================================================
/// test some critical parameters

int Grid::test_param()
{
  if( nlist.n <= 4 ){
    cerr << "Grid::test_param: ERROR: nlist.n = " << nlist.n << endl;
    return 1;
  }
  if( elist.n <= 4 ){
    cerr << "Grid::test_param: ERROR: elist.n = " << elist.n << endl;
    return 2;
  }
  if( n_corrector <= 0 ){
    cerr << "Grid::test_param: ERROR: n_corrector = "
         << n_corrector << endl;
    return 3;
  }
  
  return 0;
}
//==============================================================================
/** Make a characteristic velocity in the grid
 * 
 * Just use Darcy's law with a pressure differential given by 
 * P_max-P_min where I just use P_min=0 over the range of the
 * reservoir.
 * */
void Grid::MakeCharacteristicQuantities()
{
  // Calculate values required by Darcy's law.
  // Just use means
  
  double p_sum = 0.0;
  double kappa_sum = 0.0;
  double mu_sum = 0.0;
  
  for( int ie=0; ie<elist.n; ++ie ){
    p_sum += elist.e[ie].cur.p;
    kappa_sum += elist.e[ie].rock.Kappa();
    mu_sum += elist.e[ie].fluid.Mu();
  }
  double p_mean = p_sum / double( elist.n );
  double kappa_mean = kappa_sum / double(elist.n);
  double mu_mean = mu_sum / double(elist.n);
  
  const double dx = elist.DX();
  
  // Apply Darcy's law
  
  const double ctmp = kappa_mean/mu_mean;  // Using +ve values
  u_characteristic = ctmp * (p_mean - 0.0)/dx;
  dt_characteristic = dx / u_characteristic;
}
