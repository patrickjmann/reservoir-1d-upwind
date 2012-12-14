//#############################################################################
/*\file Node.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/Node.C,v $ 
$Revision: 1.2 $ 
$Date: 2012/10/30 22:10:19 $ 

\author P.J.Mann

Node and NodeList member functions
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

//=============================================================================
Node::Node()
{
  id = -1;
  left = 0;  right = 0;
  x = 0.0;
  flux.W = 0.0;
  oldflux.W = 0.0;
  q = 0.0;
  flag = -1;
}
//=============================================================================
void NodeList::Make( const int n_nodes )
{
  if( n_nodes < 4 ){
    cerr << "Node::Make: ERROR: n_nodes = " << n_nodes << '\n';
    exit(1);
  }
  n = n_nodes;
  node = new Node[n];
}

