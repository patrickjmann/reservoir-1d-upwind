//#############################################################################
/**\file Grid.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/Grid.C,v $ 
$Revision: 1.4 $ 
$Date: 2012/11/25 00:39:19 $ 

\author P.J.Mann

April 29, 2006

i/o routines for Grid

*/
//#############################################################################
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "putil.h"

#include "grid.h"
#include "ptime.h"

using namespace std;

//=============================================================================
/** Put the elementlist and nodelist to separate files with time information
 */
void Grid::PutSlice( TimeData& ptime )
{
  // make filename

  const int N_FILENAME = 64;
  char* filename_base = new char[N_FILENAME];

  sprintf( filename_base, "t%05i", ptime.n );

  int n_base = strlen( filename_base );
  int n_filename = n_base + 20;
  char* elist_filename = new char[n_filename];
  char* nlist_filename = new char[n_filename];

  sprintf( elist_filename, "%s.%s", filename_base, "elements.gnu" );
  sprintf( nlist_filename, "%s.%s", filename_base, "nodes.gnu" );

  // Print the data in tabular form

  ofstream se( elist_filename );
  se << "# Time step number: " << ptime.n << '\n'
     << "# Coordinate Time:  " << ptime.coord << " s. = "
     << ptime.coord/(60.0*60.0) << " days.\n";

  elist.PutTable( se );
  se.close();

  ofstream sn( nlist_filename );
  sn << "# Time step number: " << ptime.n << '\n'
     << "# Coordinate Time:  " << ptime.coord << " s.\n";
  nlist.PutTable( sn );
  sn.close();

  delete[] elist_filename;
  delete[] nlist_filename;
}
//=============================================================================
/** Put the elementlist and nodelist to separate files in tabular, readable form
 */
void Grid::PutTable( const char* filename_base )
{
  int n_base = strlen( filename_base );
  int n_filename = n_base + 20;
  char* elist_filename = new char[n_filename];
  char* nlist_filename = new char[n_filename];

  sprintf( elist_filename, "%s.%s", filename_base, "nodes.gnu" );
  sprintf( nlist_filename, "%s.%s", filename_base, "elements.gnu" );

  ofstream se( elist_filename );
  elist.PutTable( se );
  se.close();

  ofstream sn( nlist_filename );
  nlist.PutTable( sn );
  sn.close();

  delete[] elist_filename;
  delete[] nlist_filename;
}
//=============================================================================
void Grid::PutParam( ofstream& s )
{
  s.setf( ios::left, ios::adjustfield );

  s << GNU_COMMENT << " Grid:: Control Parameters:\n"
    << GNU_COMMENT << "---------------------------\n"
    << GNU_COMMENT << " ElementList .elist.n = " << setw(7) << elist.n << '\n'
    << GNU_COMMENT << " NodeList    .nlist.n = " << setw(7) << nlist.n << '\n'
    << GNU_COMMENT << " .n_corrector = " << setw(2) << n_corrector << '\n';

  s.precision(3);
  s << GNU_COMMENT << " .artvis: .k1 = " << setw(7) << artvis.k1
    << " .k2 = " << setw(7) << artvis.k2 << '\n';
  
  s.flush();
}
//=============================================================================
ostream& operator << ( ostream& s, const Grid& grid )
{
  s.setf( ios::left, ios::adjustfield );
  
  SetLineStyle( '=' );
  PutLine( s );
  s << GNU_COMMENT << " Elements\n";
  SetLineStyle( '-' );
  PutLine( s );
  s << grid.elist << '\n';

  SetLineStyle( '=' );
  PutLine( s );
  s << GNU_COMMENT << " Nodes\n";
  SetLineStyle( '-' );
  PutLine( s );
  s << grid.nlist << '\n';

  s.flush();
  return s;
}
