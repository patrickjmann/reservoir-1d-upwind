//#############################################################################
/*\file Initialize.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Initialize.C,v $ 
$Revision: 1.9 $ 
$Date: 2012/11/23 22:54:08 $ 

\author P.J.Mann

Initialization routines for various initial data sets.
*/
//#############################################################################
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

extern logstream plog;

//=============================================================================
/** Make initial data
 */
void Grid::MakeInitial( CommandLineOptions& cloption )
{
  switch( cloption.init_type ){
  case TUBE:
    plog << "Grid::MakeInitial: INFO: grid.MakeTube()..\n";
    MakeTube( cloption );
    break;
  default:
    plog << "Grid::MakeInitial: ERROR: unknown init_type "
         << cloption.init_type << '\n';
    exit(1);
  }

  MakeCharacteristicQuantities();
                
  elist.Exchange();
}
//=============================================================================
/** Make the grid geometry and positions.
    
A simple constant spaced grid.  Elements contain pointers to the nodes
and nodes contain pointers to left/right elements.
*/

int Grid::MakeGeometry( int n_elements, double x_left, double x_right )
{
  if( (fabs(x_left - x_right) <= 1.0e-4) ){
    cerr << "Grid::MakeGeometry: ERROR: initial left/right data\n"
	 << "  x_left = " << x_left << " x_right = " << x_right << '\n';
    exit(1);
  }
  if( n_elements < 4 ){
    cerr << "Grid::MakeGeometry: ERROR: n_elements is too small.\n"
	 << "  n_elements = " << n_elements << '\n';
    exit(1);
  }

  int n_nodes = n_elements + 1;

  dx_mean = (x_right - x_left )/double(n_elements);
  
  nlist.Make( n_nodes );
  elist.Make( n_elements );

  // Make the geometry

  for( int ie=0; ie<n_elements; ++ie ){
    elist.e[ie].left  = &(nlist.node[ie]);
    elist.e[ie].right = &(nlist.node[ie+1]);
  }

  nlist.node[0].left  = 0;  // no element outside the boundary
  nlist.node[0].right = &(elist.e[0]);

  for( int in=1; in<(n_nodes-1); ++in ){
    nlist.node[in].left  = &(elist.e[in-1]);
    nlist.node[in].right = &(elist.e[in]);
  }

  nlist.node[n_nodes-1].left  = &(elist.e[n_nodes-2]);
  nlist.node[n_nodes-1].right = 0;

  // make useful but unnecessary id's

  for( int ie=0; ie<n_elements; ++ie ){
    elist.e[ie].id = ie;
  }

  for( int in=0; in<n_nodes; ++in ){
    nlist.node[in].id = in;
  }

  // Positions
  
  // Make the nodes first and then the elements can
  // be calculated from the nodal positions.

  for( int i_node=0; i_node<n_nodes; ++i_node ){
    nlist.node[i_node].x = x_left + dx_mean * double(i_node);
  }

  for( int ie=0; ie<n_elements; ++ie ){
    elist.e[ie].x = 0.5 * ( elist.e[ie].left->x + elist.e[ie].right->x );
    elist.e[ie].dx = elist.e[ie].right->x - elist.e[ie].left->x;
  }
  
  return 0;
}
//=============================================================================
/** Make initial data for the tube.

This is a single call which does everything.
* 
* Default is just Loose Sand and Oil
*/
int Grid::MakeTube( CommandLineOptions& cloption )
{
  // Set default tube parameters and then get the rest.
  // These are in meters (stick with MKS everywhere for now).

  double x_left  = -500.0;   // meters 
  double x_right =  500.0;

  /* Pressure will also depend on the depth.  There is a formula for
   * that somewhere.  However I'll just try the reference (surface)
   * pressure for now.
   * */

  Element e_tmp;                       // Kludge to get the built-in fluid reference values.  These should be input from datafile!
  double p_left  = e_tmp.fluid.P_R(); // reference pressure on the left
  double p_right = 0.5*p_left;        // will zero work?
    
  // Make a generic parameter list
  // For now some sort of linear fit will be sufficient.
    
  ParameterSet param;

  param.Append( x_left , "x_left" );
  param.Append( x_right, "x_right" );
  param.Append( p_left , "p_left" );
  param.Append( p_right, "p_right" );

  // Get parameters from standard data file

  string param_file = "param.dat";
  if( verbose ){
    cerr << "Grid::MakeTube: INFO:  opening parameter file \"" << param_file << "\"\n";
  }
  ifstream s_in( param_file.c_str() );
  if ( ! s_in ){
    cerr << "Grid::MakeTube: WARNING: there is no data file \"" << param_file << "\" so default values are used\n";
  } else {
    param.Get( s_in );
    s_in.close();
  }

  // Make the grid and geometry

  int ierr = MakeGeometry( cloption.n_elements, x_left, x_right );
  if( ierr != 0 ){
    cerr << "Grid::MakeTube: ERROR: MakeGeometry returns " << ierr << '\n';
    exit(1);
  }

  // Make the tube

  MakeTube( cloption, x_left, x_right, p_left, p_right );
  
  //  Boundary
  
  bdy_cond.Make( elist.e[0].fluid.rho(p_left),
                elist.e[elist.n-1].fluid.rho(p_right) );

  // Dump initial data

  plog.s << "Initial Data:  Horizontal Tube\n"
         << "  x_left  = " << setw(10) << x_left
         << "  x_right = " << setw(10) << x_right << '\n'
         << "  p_left  = " << setw(10) << p_left 
         << "  p_right = " << setw(10) << p_right << '\n';
  
  return 0;
}
//=============================================================================
/** Make Reservoir Tube initial data.

Note that only the elements have data.  The nodes compute data as
required from the adjoining elements.
* 
* Note that generally the pressure p is given.
*/
int Grid::MakeTube( CommandLineOptions& cloption, double x_left, double x_right,
		    double p_left, double p_right )
{
  // should check parameters

  if( (p_left <= 0.0) ||
      (p_right <= 0.0) ){
    cerr << "Grid::MakeTube: ERROR: incorrect initial data\n"
	       << "  p_left = " << p_left << " p_right = " << p_right << '\n';
    exit(1);
  }
  
  /* Porosity
  * Just set to the reference porosity.
  * */
  
  for( int ie=0; ie<elist.n; ++ie ){
    elist.e[ie].SetPhi();
  }
  
  /* Other Fluid and Rock properties
   * In this first version the fluid and rock properties are constant
   * and created by the Constructors.
   * So nothing to do here.
   * */
  
  // Make p on each node.  Just a linear fit for now.
  
  const double slope = (p_right - p_left)/(x_right - x_left);
  for( int ie=0; ie<elist.n; ++ie ){
    const double x = elist.e[ie].x;
    elist.e[ie].cur.p = p_left + slope * ( x - x_left ); 
  }
  
  /* Make the rest of the element variables.
   * Order is important in the below.
   * For instance need to make Pressure before making U from Darcy's law.
   * */
  elist.MakeRhoFromP();   // from equation of state
  elist.MakeUFromDarcy();
  elist.MakeW();          // From formula

  return 0;
}    
