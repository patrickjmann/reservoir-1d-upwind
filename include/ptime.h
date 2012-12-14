//#############################################################################
/**\file ptime.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/ptime.h,v $ 
$Revision: 1.3 $ 
$Date: 2012/11/25 00:39:17 $ 

\author P.J.Mann

Physical times and time step control
*/
//#############################################################################

#if ! defined CONS1_P_TIME_H
#define CONS1_P_TIME_H

#include "putil.h"
#include "grid.h"
#include <exception>

using namespace std;

//==============================================================================
/// Exception that can be thrown by any of the time step functions

class TimeException: public exception
{
  public:
  virtual const char* what() const throw()
  {
    return "TimeException: a problem occurred in the Time Step classes";
  }
};
//=============================================================================
struct TimeStepControl {
  double x;
  double c;    /// Courant parameter
  double p1;   /// initial delta = p1*freefall
  double p2;   /// unused?
};

class TimeData {
 public:
  enum TimeLimits {
    LOTS_OF_TIME = 0,
		OUT_OF_TIME,
		OUT_OF_TIME_STEPS
  };
  
 private:
  int n;                                 /// time step number
  int n_max;                             /// number of time steps required
  
  double coord;                          /// coordinate time
  double coord_max;
  
  double delta,delta_new;                /// time step size
  double delta_min,delta_max;            /// min/max allowed time steps
  double cs_max, dx_max;                 /// Courant values useful for error dumps
  
  TimeStepControl step_control;
  
 public:
  TimeData(){ n = -1; n_max=-1; coord=-1.0; coord_max=-1.0; delta=-1.0; delta_new=-1.0; delta_min=-1.0; delta_max=-1.0; cs_max=-1.0; dx_max=-1.0; };
  void Make( string& param_file );
  virtual ~TimeData(){};
  
  void MakeDelta( Grid& );
  void Update();
  void SetMinMax();
  void CheckDelta() throw(TimeException);
  int TestForMax();
  
  void TestRun(){ n_max = 3; };
  
  // accessors
  
  inline double Delta() const { return delta; };
  inline double DeltaMin() const { return delta_min; };
  
  inline int N() const { return n; };
  inline int NMax() const { return n_max; };
  
  inline double Coord() const { return coord; };
  inline double CoordMax() const { return coord_max; };
  
  // i/o
  
  void PutParam( ofstream& );
  friend ostream& operator << ( ostream&, const TimeData& );
  
  // friends
  
  friend class Grid;
  friend class LastDump;  // Used in main() to check on time for another dump.
};

#endif
