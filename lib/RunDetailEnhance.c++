#include "DetailEnhance.h"
#include <tclap/CmdLine.h>

using namespace TCLAP;

int main( int argc, char **argv )
{
  try
  {
    CmdLine cmd("Computes edge aware smoothing.", ' ', "2.0");

    ValueArg<string> outputArg( "o", "output", "Output name", true, "out", "string", cmd);

    ValueArg<string> pictureArg( "p", "picture", "Picture name", true, "in", "string", cmd);

    cmd.parse( argc, argv );

    string inputName                  = pictureArg.getValue();
    string outputName                 = outputArg.getValue();

    if( VIPS_INIT( argv[0] ) ) return( -1 );

    RunDetailEnhance( inputName, outputName );
  }
  catch (ArgException &e)  // catch any exceptions
  {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  return 0;
}
