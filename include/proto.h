//#############################################################################
/**\file proto.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/proto.h,v $ 
$Revision: 1.1 $ 
$Date: 2012/10/02 21:52:01 $ 

\author P.J.Mann

*/
//-----------------------------------------------------------------------------
#if ! defined CONS1_PROTO_H
#define CONS1_PROTO_H

using namespace std;

// Various useful utilities

int GetId();
void MakeOutFiles();
void Set_OS_Stuff( );

// Primitive variables solver

int fvalue( const double pcurrent,
	    const double S, const double D, const double E,
	    const double gamma,
	    double& f, double& fprime );

int PrimitiveSolve( const double S, const double D, const double E,
		    const double gamma,
		    const double pguess,
		    double& v, double& rho0, double& eps, double& p );


#endif
