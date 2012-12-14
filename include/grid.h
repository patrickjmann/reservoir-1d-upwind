//######################################################################
/**\file grid.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/grid.h,v $ 
$Revision: 1.5 $ 
$Date: 2012/11/23 18:25:06 $ 

\author P.J.Mann

\brief the basic 3d grid structure

*/
//######################################################################
#if ! defined CONS1_GRID_H
#define CONS1_GRID_H

#include <cmath>

#include "putil.h"
#include "cpu.h"

#include "cloption.h"

#include "element.h"
#include "node.h"
#include "boundary.h"

// some prototypes

using namespace std;

class TimeData;
class InputParameterSet;

//======================================================================
/** Simply a container for the nodelist and elementlist.
*/

class Grid
{
 public:
  
 private:
  ElementList elist;
  NodeList nlist;

  // Various control parameters

  ArtVisControl artvis;         /// artificial viscosity
  double dx_mean;               /// mean step size  
  
  // boundaries
  
  BoundaryCondition bdy_cond;   /// boundary conditions

  // There are a lot of characteristic quantities and counters here

  int corrector_level;
  int n_corrector;

  char* run_id_string;          /// string describing the run
  int run_id;                   /// run id number
  
  double u_characteristic;      /// used in time step choice.
  double dt_characteristic;     /// from u_characteristic

  // private routines

  int test_param();

  // Make the geometry (containing elist and nlist);

  int MakeGeometry( int n_elements, double x_left, double x_right );

  // Make initial data (MakeGeometry has to be called first)

  int MakeTube( CommandLineOptions& );  /// make initial data (reads various parameters)
  int MakeTube( CommandLineOptions&, double x_left, double x_right, double p_left, double p_right ); /// Make initial data from the given input values.
  
  void MakeCharacteristicQuantities();  // Various useful quantities

  // Construction

 public:
  
  Grid( CommandLineOptions&, char* run_id_string_in, int run_id_in );
  virtual ~Grid(){};

  // Initial data
  /* these are single calls that do everything.  They call
     MakeGeometry and then a suitable data function.
  */

  void MakeInitial( CommandLineOptions& cloption );

  // Accessors
  
  inline int CorrectorLevel() const { return corrector_level; };
  
  inline char* RunIdString() { return run_id_string; };
  inline int RunId() const { return run_id; };
  
  inline ArtVisControl& ArtVis(){ return artvis; };
  
  // i/o
  
  friend ostream& operator << ( ostream&, const Grid& );
  void PutParam( ofstream& );
  void PutTable( const char* filename_base );
  void PutSlice( TimeData& );

  // Evolution routines
 
  void Evolve( CommandLineOptions& cloption, TimeData& ptime );

  friend class TimeData;
};

#endif
