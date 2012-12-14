//#############################################################################
/**\file node.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/node.h,v $ 
$Revision: 1.3 $ 
$Date: 2012/11/23 22:54:11 $ 

\author P.J.Mann

\brief the basic 3d grid structure

*/
//#############################################################################
#if ! defined CONS1_NODE_H
#define CONS1_NODE_H

#include <cmath>

#include "putil.h"
#include "cpu.h"

#include "cloption.h"
#include "basic.h"
#include "boundary.h"

// some prototypes

using namespace std;

class Element;
class TimeData;
class InputParameterSet;
class ElementList;

//=============================================================================
/** A Single Node in the grid
    
Note:  these are the "i+1/2" points at the intersection between nodes

Fundamentally these nodes only need to know which element is on
the left and right.  They derive the fluxes from the elements on the left
and right.

*/

class Node {
 private:
  int id;        /// id is useful

  Element* left;
  Element* right;

  double x;

  Flux flux;    /// mean Flux used to actually do the evolution (avg?, upwind?, ..)
  Flux oldflux; /// Flux on old slice used for predictor-correctors

  double q;     /// need to store this for time-step Courant limits
   
 public:

  int flag;    /// Can be used by any routine to flag or identify particular nodes

  Node();
  virtual ~Node(){};

  /* Evolution:  Note that this only needs to know
     how to make the flux, and how to add resulting integrated
     quantities to the adjoining elements.
     However this may require a Riemann solver
  */
  void MakeHighRes( int, BaseVariables& highres_left, BaseVariables& highres_right );

  void MakeFlux( CommandLineOptions&, BaseVariables& base_left, BaseVariables& base_right,
      double leftx, double rightx,
      TimeData& ptime, ArtVisControl& artvis, int sign, Flux& flux );
      
  void MakeFluxBoundary( CommandLineOptions&, BoundaryCondition&, TimeData&, int sign, Flux& );

  void MakeFluxPredictor( CommandLineOptions&, TimeData&, ArtVisControl& );
  void MakeFluxCorrector( CommandLineOptions&, TimeData&, ArtVisControl& );

  // Various useful accessors

  inline double X() const { return x; };
  inline int Id() const { return id; };
  
  // i/o

  void PutTable( ostream& );
  friend ostream& operator << ( ostream&, Node& );
  friend ostream& operator << ( ostream&, const Element& );

  void PutTable( ofstream& );
  
  // Friends: More-or-less every class eventually accesses data on the nodes.
  
  friend class Grid;
  friend class NodeList;
  friend class Element;
  friend class ElementList;
  friend class TimeData;
};

//=============================================================================
/** Basically a vector of nodes.
*/

class NodeList
{
 public:
  
 private:
  Node* node;           /// vector of nodes: (not pointers to nodes)
  int n;                /// number of nodes
      
  // private routines

  int make_initial_nodes();

  // Construction

  public:

  NodeList(){ node=0; n=-1; };  // constructor used by Grid
  virtual ~NodeList(){
    if( node != 0 ){ delete[] node; };
  };

  void Make( const int n_nodes );

  // Accessors
    
  inline int N() const { return n; };
    
  inline Node* operator []( int i ){ return &(node[i]); };

  // i/o

  friend ostream& operator << ( ostream&, const NodeList& );
  void PutTable( ofstream& );

  // Evolution

  void MakeFluxPredictor( CommandLineOptions&, TimeData&, ArtVisControl&, BoundaryCondition& );
  void MakeFluxCorrector( CommandLineOptions&, TimeData&, ArtVisControl&, BoundaryCondition& );

  // Miscellaneous

  friend class Grid;  
  friend class Element;
  friend class TimeData;
};

#endif
