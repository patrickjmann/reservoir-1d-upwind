//#############################################################################
/**\file cloption.C

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/io/RCS/cloption.C,v $ 
$Revision: 1.4 $ 
$Date: 2012/11/25 00:39:19 $ 

\author P.J.Mann

Command-line options
*/
//#############################################################################
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

#include <getopt.h>

#include "putil.h"

#include "grid.h"
#include "cloption.h"

using namespace std;
int verbose;

//=============================================================================
CommandLineOptions::CommandLineOptions()     // defaults
{
  n_elements          = 100;
  init_type           = TUBE;
  evolution_type      = UPWIND;
  limiter_type        = NONE;
  reconstruction_type = CONSTANT;
  verbose             = FALSE;
}
//=============================================================================
void CommandLineOptions::give_help()
{
  cerr << '\n'
       << PACKAGE << ": "
       << "Relativistic Shock: High-resolution version\n";
  
  cerr << '\n'
       << "USAGE: " << PACKAGE << " [-v] [-g #] [-n #] [-i]\n";
  
  cerr << '\n'
       << "   -v #  : version information \n"
       << "   -n #  : number of nodes (default = " << n_elements << ") \n"
       << "   -d    : create default data file \n"
       << "   --verbose: Lots and lots of output\n";

  cerr << '\n'
       << "Initial Data\n"
       << "  --tube  : Simple Tube (DEFAULT)\n";

  cerr << '\n'
       << "Evolution Method\n"
       << "  --upwind  : Use the Upwind method (DEFAULT)\n";

  cerr << '\n'
       << "Limiter Used to limit the Slopes (default is NONE)\n"
       << "  --minmod   : Use MINMOD limiter\n"
       << "  --mc       : Use MC limiter\n"
       << "  --superbee : Use SuperBee limiter (not currently supported)\n";

  cerr << '\n'
       << "Reconstruction method (default is CONSTANT)\n"
       << "  --constant     : Use Constant reconstruction (default)\n"
       << "  --linear       : Use linear reconstruction\n"
       << "  --cubichermite : Use Cubic Hermite reconstruction\n";
  
  cerr << endl;
}
//=============================================================================
void MakeDataFile( char* datafile )
{
  time_t now = time(0);
  char* now_string = ctime(&now);
  
  ofstream s( datafile );
  if( !s ){
    cerr << "MakeDataFile: ERROR: unable to open \"" << datafile << "\"\n";
    exit(1);
  }
  s << "# Test parameter file for " << PACKAGE << '\n'
    << "# Generated " << now_string << '\n'
    << "# Units:  MKS" << '\n'
    << "#\n";
  
  const double default_time = 24.0*60.0*60.0;  // one day
  
  s << "tmax     " << default_time << " # total time (one day in seconds)\n"
    
    << "evolve.n     5000  # number of evolutions\n"
    << "epsdmp        1.0  # dump control: allowed change in evolution variables\n"
    << "dump.n       1000  # interval between dumps\n"
    << "corrector.n     2  # number of correctors\n"
    << "#\n";
  
  s << "artvis.k1     0.0  # artificial viscosity (dimensionless)\n"
    << "artvis.k2     0.0\n"
    << "#\n";
  
  s << "time_step_control.c  0.4   # courant\n"
    << "time_step_control.p1 0.4   # initial delta\n"
    << "#\n";

  s << "x_left       -500.0  # (meters)\n"
    << "x_right      +500.0  # (meters)\n"
    << "#\n";

  double p_atmos = 101325.0 ; // Atmospheric pressure
  s << "p_left     " << p_atmos << " # atmospheric pressure (pascals)\n"
    << "p_right    " << 2.0*p_atmos << " # (pascals)\n"
    << "#\n";
}
//=============================================================================
void CommandLineOptions::Get( const int argc, char* argv[] )
{
  extern int opterr;
  //extern int optind;
  extern char* optarg;
  int opt;

  // Note: turn off getopt() error checking
  
  opterr = FALSE;

  // Long options

  static struct option long_options[] =
    {
      {"verbose", no_argument, &verbose, TRUE },
      
      {"tube",  no_argument, &init_type, TUBE },
      
      {"upwind", no_argument, &evolution_type, UPWIND },

      {"minmod",   no_argument, &limiter_type, MINMOD },
      {"mc",       no_argument, &limiter_type, MC },
      {"superbee", no_argument, &limiter_type, SUPERBEE },

      {"constant", no_argument, &reconstruction_type, CONSTANT },
      {"linear",   no_argument, &reconstruction_type, LINEAR },
      {"cubic",    no_argument, &reconstruction_type, CUBICHERMITE },

      {"elements", required_argument, 0, 'n' },
      {"datafile", no_argument, 0, 'd' },
      {"version", no_argument, 0, 'v' },
      {"help", no_argument, 0, 'h' },
      {0,0,0,0}
    };
  
  while( 1 ){

    int option_index = 0;  // getopt_long stores the option index here

    opt = getopt_long( argc, argv, "h?vn:de:", long_options, &option_index );
    
    if( opt == -1 ) break;

    switch( opt ){
      case 0:
        if( long_options[option_index].flag != 0 ) break;
        break;
      case '?':         // getopt() does not recognize the argument
        cerr << "CommandLineOptions: ERROR: invalid command line option\n";
        give_help();
        exit( 0 );
      case 'h':
        give_help();
        exit( 0 );
      case 'n':
        n_elements = atoi( optarg );
        break;
      case 'v':
        cerr << PACKAGE << " version " << VERSION << '\n';
        exit( 0 );
      case 'd':{
        char* datafile = (char*)"param.dat.test";
        cerr << "cloption: INFO: "
             << "Test parameter file generated in \"" << datafile << "\"\n"
             << "       (should be moved to \"param.dat\")\n";
        MakeDataFile( datafile );
        exit( 0 );
      }
      default:
        cerr << "CommandLineOptions: ERROR: invalid command line option\n";
        give_help();
        exit( 0 );
    }
  }
}
//=============================================================================
void CommandLineOptions::Put( ostream& s )
{
  s << "# CommandLineOptions:\n"
    << "#--------------------\n";

  s << "#   .n_elements  = " << n_elements << '\n';

  s << "#   .init_type   = " << init_type;
  switch( init_type ){
    case TUBE: s << " Standard TUBE"; break;
  }
  s << '\n';

  s << "#   .evolution_type = " << evolution_type;
  switch( evolution_type ){
    case UPWIND: s << " (UPWIND)"; break;
    default: s << " (UNKNOWN)"; break;
  }
  s << '\n';

  s <<"#   .limiter_type = " << limiter_type;
  switch( limiter_type ){
    case NONE:
      s << " (NONE)  This is dangerous.  No method is stable without limiting of some sort. ******";
      cerr << "**** WARNING: Using NO limiter, which is ALWAYS UNSTABLE ***********\n";
      break;
    case MINMOD:   s << " (MINMOD)"; break;
    case MC:       s << " (MC)"; break;
    case SUPERBEE:
      s << " (SUPERBEE).  WARNING: Not currently supported";
      cerr << "CommandLineOptions::Put: ERROR: SUPERBEE is not currently supported\n";
      exit(1);
      break;
  }
  s << '\n';

  s << "#   .reconstruction_type = " << reconstruction_type;
  switch( reconstruction_type ){
    case CONSTANT: s << " (CONSTANT)"; break;
    case LINEAR: s << "  (LINEAR)"; break;
    case CUBICHERMITE: s << " (CUBICHERMITE)"; break;
  }
  s << '\n';
}
