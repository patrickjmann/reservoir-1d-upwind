//#############################################################################
/**\file MakeHighRes.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/evolve/RCS/MakeHighRes.C,v $
$Revision: 1.2 $
$Date: 2012/10/30 22:10:19 $

\author P.J.Mann: June, 2012

Compute high-resolution approximations on left and right sides of the node.
*/
//#############################################################################
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "cloption.h"
#include "grid.h"
#include "element.h"
#include "ptime.h"

//=============================================================================
/** Cubic Hermite interpolation basis functions

\param[in] nim Linear basis function N_{i-1} on the left
\param[in] ni  Linear basis function N_{i} on the right
\param[in] h length left to right edge
\param[out] cubic  Cubic Hermite basis functions
*/
inline void CubicHermite( double nim, double ni,
                          double h,
                          double *cubic )
{
  cubic[0] = nim * nim * (3.0 - nim * 2.0);
  cubic[2] = 1.0 - cubic[0];
  
  const double a  = h * nim * ni;
  cubic[3] = -a * ni;
  cubic[1] = a*nim;
}
//============================================================================
/** Make constant extrapolations from base position to highres position.

\param[in] elem Element from which to start
\param[in] x Position at which to calculate the high-resolution extrapolation
\param[out] highres High resolution base variables
*/
void BaseVariables::MakeHighResConstant( Element* elem, double x )
{
  BaseVariables* base = &(elem->cur);   // just a handy short form
  
  W = base->W;
  rho = base->rho;
  u = base->u;
  p = base->p;
  phi = base->phi;
}

//============================================================================
/** Make linear extrapolations from base position to highres position.

\param[in] elem Element from which to start
\param[in] x Position at which to calculate the high-resolution extrapolation
\param[out] highres High resolution base variables
*/
void BaseVariables::MakeHighResLinear( Element* elem, double x )
{
  const double x_base = elem->x;
  BaseVariables* base = &(elem->cur);   // just a handy short form
  
  const double dx = x - x_base;
  W = base->W + base->sigma_W * dx;
  rho = base->rho + base->sigma_rho * dx;
  u = base->u + base->sigma_u * dx;
  p = base->p + base->sigma_p * dx;
  phi = base->phi + base->sigma_phi * dx;
}
//============================================================================
/** Make Cubic Hermite extrapolations from the elements elem1 and elem2 to x

Note: elem1 must be to the left of elem2

\param[in] elem1 Element 1
\param[in] elem2 Element 2
\param[in] x  Position at which to calculate the high-resolution (cubic hermite) extrapolation
\param[out] highres_left High resolution base variables extrapolated to x
*/
void BaseVariables::MakeHighResCubicHermite( Element* elem1, Element* elem2, double x )
{
  if( elem1->x >= elem2->x ){
    cerr << "BaseVariables::MakeHighResCubicHermite: ERROR: elem1 and elem2 are not in order\n";
    exit(1);
  }
  const double h = elem2->x - elem1->x;
  const double ni = (x - elem1->x) / h;
  //const double ni = (elem2->x - x) / h;
  const double nim = 1.0 - ni;

  double cubic[4];
  CubicHermite( nim, ni, h, cubic );

  W = elem1->cur.W * cubic[0] + elem1->cur.sigma_W * cubic[1] +
    elem2->cur.W * cubic[2] + elem2->cur.sigma_W * cubic[3];

  rho = elem1->cur.rho * cubic[0] + elem1->cur.sigma_rho * cubic[1] +
    elem2->cur.rho * cubic[2] + elem2->cur.sigma_rho * cubic[3];

  u = elem1->cur.u * cubic[0] + elem1->cur.sigma_u * cubic[1] +
    elem2->cur.u * cubic[2] + elem2->cur.sigma_u * cubic[3];

  p = elem1->cur.p * cubic[0] + elem1->cur.sigma_p * cubic[1] +
    elem2->cur.p * cubic[2] + elem2->cur.sigma_p * cubic[3];

  phi = elem1->cur.phi * cubic[0] + elem1->cur.sigma_phi * cubic[1] +
    elem2->cur.phi * cubic[2] + elem2->cur.sigma_phi * cubic[3];
}
//============================================================================
/** Make the high resolution approximations on the right and left of the current Node

\param[in] reconstruction_type   Switch for linear or cubic reconstruction
\param[out] highres_left    High resolution base variables on left
\param[out] highres_right   High resolution base variables on right 
*/
void Node::MakeHighRes( int reconstruction_type,
			BaseVariables& highres_left, BaseVariables& highres_right )
{
  if( reconstruction_type == CONSTANT ){
    highres_left.MakeHighResConstant( left, x );
    highres_right.MakeHighResConstant( right, x );
    return;
  }

  if( reconstruction_type == LINEAR ){
    highres_left.MakeHighResLinear( left,  x );
    highres_right.MakeHighResLinear( right, x );
    return;
  }

  if( reconstruction_type != CUBICHERMITE ){
    cerr << "Node::MakeHighRes: ERROR: unknown reconstruction type.  reconstruction_type="
         << reconstruction_type << '\n';
    exit( 1 );
  }

  // Elements

  Element* elem1;
  Element* elem2;  // will be ordered, elem1 to the left of elem2
    
  // Left

  elem2 = left; 
  elem1 = elem2->left->left;  // Left element pointing to lefter node pointing to lefter element
  if( elem1 == 0 ){
    highres_left.MakeHighResLinear( elem2, x );
  } else {
    highres_left.MakeHighResCubicHermite( elem1, elem2, x );
  }

  // Right

  elem1 = right;
  elem2 = elem1->right->right;
  if( elem2 == 0 ){
    highres_right.MakeHighResLinear( elem1, x );
  } else {
    highres_right.MakeHighResCubicHermite( elem1, elem2, x );
  }
}
