//#############################################################################
/**\file Element.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/Element.C,v $ 
$Revision: 1.4 $ 
$Date: 2012/11/25 00:39:19 $ 

\file P.J.Mann

get/put element data

Input
-----
s    ! input stream  (my binary/ascii stream)

Notes
-----
This routine assumes that the ElementList has been constructed with
some given sizes.  It does, however, check to see that the size is
consistent.
*/
//#############################################################################
#include <cstdlib>
#include <iostream>

#include "putil.h"
#include "element.h"

using namespace std;

//=============================================================================
void Element::PutTable( ostream& s )
{
  s.setf( ios::left, ios::adjustfield );
  s.precision( 4 );
  const int width = 10;
  s << "  "
    << setw(6)  << id << ' '
    << setw(width) << xmid() << ' '
    << setw(width) << Volume()  << ' ' 
    << setw(width) << cur.W << ' '
    << setw(width) << cur.rho << ' '
    << setw(width) << cur.phi << ' '
    << setw(width) << cur.u << ' '
    << setw(width) << cur.p << ' '
    << setw(width) << fluid.c_f << ' '
    << setw(width) << fluid.p_r << ' '
    << setw(width) << fluid.rho_r << ' '
    << setw(width) << fluid.mu << ' '
    << setw(width) << left->id << ' '
    << setw(width) << right->id << ' '
    << '\n';
}
//=============================================================================
void ElementList::PutTable( ostream& s )
{
  s.setf( ios::left, ios::adjustfield );
  const int width = 11;
  
  s << GNU_COMMENT << " ElementList: element quantities\n"
    << GNU_COMMENT << "--------------------------------\n";
  
  s << GNU_COMMENT << ' '
    << setw( 6) << "id(1) "
    << setw(width) << "x(2) "
    << setw(width) << "Volume(3) "
    << setw(width) << "W(4) "
    << setw(width) << "rho(5) " 
    << setw(width) << "phi(6) "
    << setw(width) << "u(7) "
    << setw(width) << "p(8) "
    << setw(width) << "fluid.c_f(9) "
    << setw(width) << "fluid.p_r(10) "
    << setw(width) << "fluid.rho_r(11) "
    << setw(width) << "fluid.mu(12) "
    << setw(width) << "left(13) "
    << setw(width) << "right(14) "
    << '\n';
  
  s.precision( 4 );
  for( int ie=0; ie<n; ie++ ){
    e[ie].PutTable( s );
  }
  s.flush();
}
//=============================================================================
ostream& operator << ( ostream& s, const Element& element )
{
  s.precision( 4 );
  s.setf( ios::right, ios::adjustfield );
  
  s << "Element:  id=" << element.id
    << "  xmid=" << setw(10) << element.xmid()      << ' '
    << "  volume=" << setw(10) << element.Volume()  << '\n';

  s << "  old.W = "       << setw(10) << element.old.W           << ' '
    << "  old.sigma_W = " << setw(10) << element.old.sigma_W << '\n'
    << "  old.rho = "     << setw(10) << element.old.rho       << ' '
    << "  old.sigma_rho = " << setw(10) << element.old.sigma_rho << '\n'
    << "  old.phi = "     << setw(10) << element.old.phi       << ' '
    << "  old.sigma_phi = " << setw(10) << element.old.sigma_phi << '\n';

   s << "  cur.W = " << setw(10) << element.cur.W          << ' '
     << "  cur.sigma_W = " << setw(10) << element.cur.sigma_W << '\n'
     << "  cur.rho = " << setw(10) << element.cur.rho      << ' '
     << "  cur.sigma_rho = " << setw(10) << element.cur.sigma_rho << '\n'
     << "  cur.phi = " << setw(10) << element.cur.phi      << ' '
     << "  cur.sigma_phi = " << setw(10) << element.cur.sigma_phi << '\n';

  s << "  left=" << setw(2) << element.left->id << ' '
    << "  right=" << setw(2) << element.right->id << '\n';
  
  return s;
}
//=============================================================================
ostream& operator << ( ostream& s, const ElementList& elist )
{
  s.setf( ios::left, ios::adjustfield );
  
  SetLineStyle( '=' );
  PutLine( s );
  s << GNU_COMMENT << " THE ELEMENT LIST\n";
  SetLineStyle( '-' );
  PutLine( s );
  
  s << GNU_COMMENT << " .n = " << elist.n << '\n';
  
  for( int i=0; i<elist.n; ++i ){
    PutLine( s );
    s << elist.e[i];
  }
  s.flush();
  return s;
}
