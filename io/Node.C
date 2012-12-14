//#############################################################################
/**\file Node.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/Node.C,v $ 
$Revision: 1.1 $ 
$Date: 2012/10/02 21:52:35 $ 

\author P.J.Mann

April 29, 2006

i/o routines for nodes

Input
-----
s     ! input stream (my binary system)

Notes
-----
d = relativistic density D = u(t) rho0
e = specific internal energy epsilon

Always assumes that a binary data set is 64 bit (double)
*/
//#############################################################################
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "putil.h"
#include "eqs.h"

#include "grid.h"

using namespace std;

extern Fluid fluid;

//=============================================================================
/** Put a table of node values.  This table is as compact
    as possible while still maintaining readability
*/
void Node::PutTable( ofstream& s )
{
  s.precision(4);
  s.setf( ios::right, ios::adjustfield );
  
  s << setw( 6) << id << ' '
    << setw(12) << x  << ' '
    << setw(12) << flux.W << ' ';
  if( left != 0 ){
    s << setw(12) << left->id  << ' ';
  } else {
    s << setw(12) << "0 (bdy)" << ' ';
  }
  if( right != 0 ){
    s << setw(12) << right->id  << ' ';
  } else {
    s << setw(12) << "0 (bdy)" << ' ';
  }

  s << '\n';
}
//=============================================================================
/** Put a table of node values.  This table is as compact
    as possible while still maintaining readability
*/
void NodeList::PutTable( ofstream& s )
{
  s.precision(4);
  s.setf( ios::right, ios::adjustfield );
  
  s << GNU_COMMENT << " NodeList\n"
    << GNU_COMMENT << "=========\n";
  
  s << GNU_COMMENT << ' '
    << setw( 6) << "Id(1) "
    << setw(13) << "x(2) "
    << setw(13) << "flux.W(3) "
    << setw(13) << "left id(6) "
    << setw(13) << "right id(6) "
    << '\n';
  
  SetLineStyle( '-' );
  PutLine( s, 106 );
  
  for( int i=0; i<n; i++ ){
    node[i].PutTable( s );
  }
  s.flush();
}
//=============================================================================
ostream& operator << ( ostream& s, Node& node )
{
  s.precision( 4 );
  s.setf( ios::left, ios::adjustfield );
  
  s << node.id
    << " : Position = " << node.x << '\n';
  
  s << " : left  element = " << node.left->Id() << '\n'
    << " : right element = " << node.right->Id() << '\n';

  return s;
}

ostream& operator << ( ostream& s, Node* node )
{
  s << *(node);
  return s;
}
//=============================================================================
ostream& operator << ( ostream& s, const NodeList& nlist )
{
  s.setf( ios::left, ios::adjustfield );
  
  SetLineStyle( '=' );
  PutLine( s );
  s << GNU_COMMENT << " NodeList\n";
  SetLineStyle( '-' );
  PutLine( s );
  s << GNU_COMMENT << " .n = " << nlist.n << '\n';
  for( int i=0; i<nlist.n; ++i ){
    PutLine( s );
    s << "Node #" << i << '\n' << nlist.node[i];
  }
  return s;
}
