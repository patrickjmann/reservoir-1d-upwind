//#############################################################################
/**\file cloption.h

$Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/include/RCS/cloption.h,v $ 
$Revision: 1.2 $ 
$Date: 2012/10/30 22:10:13 $ 

\author P.J.Mann

\brief Command-line options
*/
//#############################################################################
#if ! defined CONS_CLOPTIONS_H
#define CONS_CLOPTIONS_H

using namespace std;

#define PROGRAM_NAME PACKAGE

enum INIT_TYPE { TUBE };
enum EVOLUTION_TYPE { UPWIND };
enum LIMITER_TYPE { NONE, MINMOD, MC, SUPERBEE };
enum RECONSTRUCTION_TYPE { CONSTANT, LINEAR, CUBICHERMITE };

//=============================================================================
class CommandLineOptions {
  public:
    int n_elements;

    int evolution_type;
    int limiter_type;
    int reconstruction_type;
    int init_type;

    void give_help();

    CommandLineOptions();
    virtual ~CommandLineOptions(){};

    void Get( const int argc, char* argv[] );
    void Put( ostream& );
};
    
#endif
