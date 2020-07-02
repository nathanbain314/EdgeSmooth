#include "EdgeSmooth.h"
#include <tclap/CmdLine.h>

using namespace TCLAP;

int main( int argc, char **argv )
{
  try
  {
    CmdLine cmd("Computes edge aware smoothing.", ' ', "2.0");

    ValueArg<float> decArg( "d", "dec", "Dec value", false, 2.0, "float", cmd);

    ValueArg<int> KArg( "k", "k", "K value", false, 3, "int", cmd);

    ValueArg<float> sigma2Arg( "t", "sigma2", "Sigma2 value", false, 3.0, "float", cmd);

    ValueArg<float> sigma1Arg( "s", "sigma1", "Sigma1 value", false, 1.0, "float", cmd);

    ValueArg<float> lambdaArg( "l", "lambda", "Lambda value", false, 0.01f, "float", cmd);

    ValueArg<string> outputArg( "o", "output", "Output name", true, "out", "string", cmd);

    ValueArg<string> pictureArg( "p", "picture", "Picture name", true, "in", "string", cmd);

    cmd.parse( argc, argv );

    string inputName                  = pictureArg.getValue();
    string outputName                 = outputArg.getValue();
    float lambda                      = lambdaArg.getValue();
    float sigma1                      = sigma1Arg.getValue();
    float sigma2                      = sigma2Arg.getValue();
    int K                             = KArg.getValue();
    float dec                         = decArg.getValue();

    if( VIPS_INIT( argv[0] ) ) return( -1 );

    RunEdgeSmooth( inputName, outputName, lambda, sigma1, sigma2, K, dec );
  }
  catch (ArgException &e)  // catch any exceptions
  {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  return 0;
}
