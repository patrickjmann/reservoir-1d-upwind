//#############################################################################
/**\file element.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/element.h,v $ 
$Revision: 1.8 $ 
$Date: 2012/11/25 00:39:17 $ 

\author P.J.Mann

Element classes
*/
//#############################################################################
#if ! defined CONS1_ELEMENT_H
#define CONS1_ELEMENT_H

#include "putil.h"
#include "basic.h"
#include "node.h"
#include "physical.h"

using namespace std;

//=============================================================================
/** Contains data (averaged physical quantities) and links to
surrounding nodes.  In 1d only left and right nodes.
*/
class Element
{
 public:
  enum sizes { N_NODES = 2 };
  
 private:

  int id;        /// useful for debugging

  Node* left;    /// left side node
  Node* right;   /// right side node

  double x;     /// centroid of element
  double dx;    /// width of element
  
  /* Element Rock properties
   * The following are rock properties in this cell.
   * The properties are generally cell-dependent.
   *  ie) various types of rock in this cell.
   * Also might want to use high-resolution methods.
   * Note that various rock types, etc. have their own data.
   * But for parallelization each element should have it's own
   * rock properties.
  */
  
  Rock rock;            /// Element rock properties
  
  /* Basic fluid properties.  These are fluid-dependent but
   * it is very handy to have them available in an element.
   * Especially for parallelization it localizes data in the
   * elements
   * */
   
  Fluid fluid;          /// Element fluid properties

  // Basic quantities computed during evolution
  // Note: these are integrated quantities normalized by volume
    
  BaseVariables cur;   /// current
  BaseVariables old;   /// used for evolution and predictor-corrector

 public:
  Element();
  virtual ~Element(){ };

  // Makes
  // Note: p must be made by a call to the equation of state class

  void Make( const double u_in, const double rho_in );

  void MakePFromRho();    /// not inline because it uses equation of state
  void MakeRhoFromP();
  inline void MakeW(){ cur.W = cur.phi*cur.rho; };
  void MakeUFromDarcy();  /// From Darcy's equation.  Needs p to be made first.
  void MakeFluidFromW();

  void MakeSigma( CommandLineOptions& );
  void MakeSigmaMinMod();
  void MakeSigmaMC();
  void ZeroCurrentSigma();

  // Evolve

  void EvolveSimplePredictor( CommandLineOptions&, TimeData& );
  void EvolveSimpleCorrector( CommandLineOptions&, TimeData& );

  inline void Exchange(){
    old.CopyFrom( cur );
  }

  // Accessors

  inline double xmid()   const { return (0.5*(left->x + right->x)); };
  inline double X()      const { return x; };
  inline double Dx()     const { return dx; };
  inline double Volume() const { return dx; };
  inline double H()      const { return dx; };

  inline int Id(){ return id; };

  inline Node* Left() { return left; };
  inline Node* Right() { return right; };
  
  inline void SetPhi(){
    cur.phi = rock.phi_r;   // default
  };
  inline void SetPhi( double in_phi ){
    cur.phi = in_phi;
  };
  
  // i/o
  
  friend ostream& operator << ( ostream&, const Element& );
  void PutTable( ostream& );
  
  // Friendly classes
  
  friend class ElementList;
  friend class Node;
  friend class Grid;
  friend class TimeData;
  friend class BaseVariables;
};
//=============================================================================
/** The list of elements in the grid.
 */
class ElementList
{
 private:
  Element* e;            /// vector of elements (not pointers)
  int n;
  
 public:
  ElementList(){ e=0; n=-1; };
  virtual ~ElementList(){ if( e != 0 ){ delete[] e;} };

  // Makes

  void Make( const int n_elements );

  void MakePFromRho();
  void MakeRhoFromP();
  void MakeUFromDarcy();
  void MakeW();
  void MakeFluidFromW();
  void MakeConserved();
  
  void MakeSigma( CommandLineOptions& );
  
  void SetPhi();

  // Evolve

  void Predict( CommandLineOptions&, TimeData& );
  void Correct( CommandLineOptions&, TimeData& );

  inline void Exchange(){
    for( int i=0; i<n; ++i ){ e[i].Exchange(); }
  }
  
  // Accessors    
  
  inline int N() const { return n; };
  inline double RestMass();
  inline double DX(){ return (e[n-1].x - e[0].x); };
  
  // Public access to an element in the list: No tests here for efficiency.
  
  inline Element* operator []( int i ){ return &(e[i]); };
  
  // i/o
  
  friend ostream& operator << ( ostream&, const ElementList& );
  void PutTable( ostream& );
  
  // Friendly classes
  
  friend class Grid;
  friend class TimeData;
};

#endif
