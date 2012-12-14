//######################################################################
/**\file boundary.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/boundary.h,v $ 
$Revision: 1.1 $ 
$Date: 2012/11/23 18:25:06 $ 

\author P.J.Mann

\brief boundary conditions

*/
//######################################################################
#if ! defined RS_BOUNDARY_H
#define RS_BOUNDARY_H

#include <cstdlib>

using namespace std;
//======================================================================
class BoundaryCondition
{
  private:
  double rho_left;
  double rho_right;
  public:
  BoundaryCondition(){
    rho_left = -1;
    rho_right = -1;
  }
  inline void Make( double rho_left_in, double rho_right_in ){
    if( rho_left_in<= 0.0 ){
      cerr << "BoundaryCondition: ERROR: rho_left_in < 0\n";
      exit(1);
    }
    if( rho_right_in<= 0.0 ){
      cerr << "BoundaryCondition: ERROR: rho_right_in < 0\n";
      exit(1);
    }
    rho_left = rho_left_in;
    rho_right = rho_right_in;
  }
  inline double RhoLeft( double t ){
    if( rho_left <= 0.0 ){
      cerr << "BoundaryCondition: ERROR: rho_left<0.0 "
           << ": rho_left = " << rho_left << '\n'
           << "      Probably not initialized correctly.\n";
      exit(1);
    }
    return rho_left;    // just constant for now
  }
 inline double RhoRight( double t ){
    if( rho_left <= 0.0 ){
      cerr << "BoundaryCondition: ERROR: rho_right<0.0 "
           << ": rho_left = " << rho_right << '\n'
           << "      Probably not initialized correctly.\n";
      exit(1);
    }
    return rho_right;    // just constant for now
  }
};

#endif
