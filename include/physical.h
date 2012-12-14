//######################################################################
/**\file physical.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/physical.h,v $ 
$Revision: 1.5 $ 
$Date: 2012/11/22 23:23:30 $ 

\author P.J.Mann

\brief Properties of the physical mediums (rock and fluid).  MKS units.

*/
//######################################################################
#if ! defined RS_PHYSICAL_H
#define RS_PHYSICAL_H

using namespace std;

extern int verbose;

//======================================================================
/**Fluid properties and equation of state
 * This should be in a hierarchy or templated.
 * */

class Fluid
{
  private:
    string name;
    string source;
    double c_f;   /// fluid compressibility Pascals^{-1}
    double p_r;   /// reference pressure Pascals
    double rho_r; /// reference density kg/m^3
    double mu;    /// Dynamic viscosity pa-s

  public:
    Fluid(){
      name = "Crude Oil";
      source = "Wikipedia";
      c_f = 1.0e-6;            /// Fluid Compressibility Pa^{-1} default value using AS typical values.
      p_r = 101.325 * 1000.0;  /// Pressure Pa Atmospheric pressure is the standard reference
      rho_r = 790.0;           /// Density kg/m^3 Crude Oil 48 Degrees API
      mu = 0.319;              /// Viscosity Pa-s motor oil
    }
  
    inline const string Name(){ return name; };
    inline string Source(){ return source; };
    inline double C_F(){ return c_f; };
    inline double P_R(){ return p_r; };
    inline double Rho_R(){ return rho_r; };
    inline double Mu(){ return mu; };

    inline double p(double rho){    /// Equation of state for low compressibility
      return p_r + (1.0/c_f) * (rho/rho_r - 1.0);
    };
    inline double rho( double p ){  /// Equation of state inverse
      return rho_r * ( c_f * ( p - p_r ) + 1.0 );
    }; 
    inline double csound( double rho ){
      return ( sqrt(1.0/(c_f * rho_r)) );
    };
    
    inline void CopyFrom( Fluid& r ){
      name = r.name;
      source = r.source;
      c_f = r.c_f;
      p_r = r.p_r;
      rho_r = r.rho_r;
      mu = r.mu;
    };
    
    friend class Element;
    friend ostream& operator << ( ostream&, const Fluid& );
};
inline ostream& operator << ( ostream& s, const Fluid& e )
{
    s << "Fluid: .name = " << e.name
      << ", .c_f = " << e.c_f
      << ", .p_r = " << e.p_r
      << ", .rho_r = " << e.rho_r
      << ", .mu = " << e.mu;
    return s;
}
//===================================================================
/// Rock Properties for various rocks
/**
 * Note that porosity is calculated as it is dependent on the
 * compressibility of the rock.  So the value is kept in the
 * the cell data
 * 
 * This should be inherited from some base class
 */

class Rock
{
  private:
    string name;    /// Useful to have the name available
    string source;  /// Source of data.  Maybe a list?
    double phi_r;   /// reference porosity
    double c_r;     /// rock compressibility
    double kappa;   /// rock permeability

  public:
    Rock(){
      name = "Loose Sand";
      source = "Wikipedia";
      phi_r = 0.4;             /// dimensionless
      c_r = 1.0e-8;            /// Pa^{-1}
      kappa = 1.0e-5 * 1.0e-4; /// m^2
    };
  
    inline string Name(){ return name; };
    inline string Source(){ return source; };
    inline double C_R(){ return c_r; };
    inline double Kappa(){ return kappa; };
    inline double Phi_R(){ return phi_r; };
    
    inline void CopyFrom( Rock& r ){
      name = r.name;
      source = r.source;
      phi_r = r.phi_r;
      c_r = r.c_r;
      kappa = r.kappa;
    };
    
    friend class Element;
    friend ostream& operator << ( ostream&, const Rock& );
};
inline ostream& operator << ( ostream& s, const Rock& r )
{
    s << r.name << ": .c_r = " << r.c_r
      << ", .kappa = " << r.kappa
      << ", .phi_r = " << r.phi_r;
    return s;
}
#endif
