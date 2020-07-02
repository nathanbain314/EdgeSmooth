#include "EdgeSmooth.h"

int numberOfCPUS()
{
  int numCPU;

#ifdef WIN32
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  numCPU = sysinfo.dwNumberOfProcessors;
#else
  numCPU = sysconf(_SC_NPROCESSORS_ONLN);
#endif

  return numCPU;
}

void normalize( double * g, int n )
{
  double sum = 0;
  for( int i = 0; i < n; ++i )
  {
    sum += g[i];
  }
  for( int i = 0; i < n; ++i )
  {
    g[i] /= sum;
  }
}

void fastBlur( double * img, double sigma, int width, int height, int dim )
{
  int ksize = round( 5.0f*sigma );
  ksize = ksize | 1;

  double * g = new double[ksize];

  int hk = ksize/2;

  for( int k = 0; k < ksize; ++k )
  {
    double x = abs(hk-k);
    g[k] = exp(-(x*x)/(2.0f*sigma*sigma));
  }

  normalize(g,ksize);

  double * ret = new double[3*width*height];

  for( int i = 0, p = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      for( int k = 0; k < dim; ++k, ++p )
      {
        double gsum = 0.0;
        for( int l = 0; l < ksize; ++l )
        {
          int l1 = l - hk;
          int j2 = j + l1;

          if( j2 >= 0 && j2 < width )
          {
            gsum += g[l]*img[dim*(i*width+j2)+k];
          }
        }
        ret[p] = gsum;
      }
    }
  }

  for( int i = 0, p = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      for( int k = 0; k < dim; ++k, ++p )
      {
        double gsum = 0.0;
        for( int l = 0; l < ksize; ++l )
        {
          int l1 = l - hk;
          int i2 = i + l1;

          if( i2 >= 0 && i2 < height )
          {
            gsum += g[l]*ret[dim*(i2*width+j)+k];
          }
        }
        img[p] = gsum;
      }
    }
  }

  delete [] g;
  delete [] ret;
}

void computeReWeights( double * s, double sigma1, double sigma2, double * wx, double * wy, int width, int height )
{
  double eps = 0.00001;

  double * dx = new double[3*width*height];
  double * dy = new double[3*width*height];

  for( int i = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      int p = 3*(i*width+j);
      for( int k = 0; k < 3; ++k )
      {
        if( j == width-1 )
        {
          dx[p+k]=0;
        }
        else
        {
          dx[p+k] = s[p+k+3] - s[p+k];
        }

        if( i == height-1 )
        {
          dy[p+k]=0;
        }
        else
        {
          dy[p+k] = s[p+k+3*width] - s[p+k];
        }
      }
    }
  }

  double * gd1 = new double[3*width*height];

  for( int p = 0; p < 3*width*height; ++p )
  {
    gd1[p] = sqrt(dx[p]*dx[p]+dy[p]*dy[p]);
  }

  fastBlur( gd1, sigma1, width, height, 3 );
  fastBlur( dx, sigma2, width, height, 3 );
  fastBlur( dy, sigma2, width, height, 3 );

  for( int i = 0, p = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j, ++p )
    {
      double mg = 0.0;
      double mx = 0.0;
      double my = 0.0;

      for( int k = 0; k < 3; ++k )
      {
        mg += fabs(gd1[3*p+k]);
        mx += fabs(dx[3*p+k]);
        my += fabs(dy[3*p+k]);
      }

      mg /= 3.0;
      mx /= 3.0;
      my /= 3.0;

      wx[p] = 1.0/max(mg*mx,eps);
      wy[p] = 1.0/max(mg*my,eps);
    }
  }

  fastBlur( wx, sigma1/2.0, width, height, 1 );
  fastBlur( wy, sigma1/2.0, width, height, 1 );

  for( int i = 0; i < height; ++i )
  {
    wx[i*width+width-1] = 0.0;
  }

  for( int j = 0; j < width; ++j )
  {
    wy[(height-1)*width+j] = 0.0;
  }

  delete [] gd1;
  delete [] dx;
  delete [] dy;
}

void columnize( double * a, int width, int height )
{
  double * t = new double[width*height];

  for( int p = 0; p < width*height; ++p )
  {
    t[p] = a[p];
  }


  for( int j = 0, p = 0; j < width; ++j )
  {
    for( int i = 0; i < height; ++i, ++p )
    {
      a[p] = t[i*width+j];
    }
  }

  delete [] t;
}

void decolumnize( double * a, int width, int height )
{
  double * t = new double[width*height];

  for( int p = 0; p < width*height; ++p )
  {
    t[p] = a[p];
  }

  for( int j = 0, p = 0; j < width; ++j )
  {
    for( int i = 0; i < height; ++i, ++p )
    {
      a[i*width+j] = t[p];
    }
  }

  delete [] t;
}

void solveLinearEquation( double * s, double * in, double * wx, double * wy, double lambda, int width, int height )
{
  int n = width * height;

  double * D = new double[n];

  columnize( wx, width, height );
  columnize( wy, width, height );

  for( int i = 0; i < n; ++i )
  {
    wx[i] *= lambda;
    wy[i] *= lambda;
    D[i] = wx[i] + wy[i];
    if( i >= height ) D[i] += wx[i-height];
    if( i > 0 ) D[i] += wy[i-1];
  }

  std::vector<T> B;

  for( int i = 0; i < n; ++i )
  {
    if( i < n - height )
    {
      if( abs(wx[i]) > 0.00001 )
      {
        B.push_back( T(i+height,i,-wx[i]));
//        B.push_back( T(i,i+height,-wx[i]));
      }
    }
    if( i < n - 1 )
    {
      if( abs(wy[i]) > 0.00001 )
      {
        B.push_back( T(i+1,i,-wy[i]));
///        B.push_back( T(i,i+1,-wy[i]));
      }
    }

    B.push_back( T(i,i,D[i]+1.0f));
  }

  SpMat A1(n,n);
  A1.setFromTriplets(B.begin(), B.end());
 
//  Eigen::SimplicialLLT<SpMat> chol(A1);
  Eigen::SimplicialCholesky<SpMat> chol(A1);
//  Eigen::SimplicialLDLT<SpMat> chol;//(A1);
//  Eigen::SparseLU<SpMat> chol(A1);

//  chol.analyzePattern(A1);
//  chol.factorize(A1);

  Eigen::VectorXd b(n);

  double * tin = new double[n];


  for( int k = 0; k < 3; ++k )
  {
    for( int i = 0; i < n; ++i )
    {
      tin[i] = in[3*i+k];
    }

    columnize(tin,width,height);

    for( int i = 0; i < n; ++i )
    {
      b[i] = tin[i];
    }

    Eigen::VectorXd x = chol.solve(b); 

    for( int i = 0; i < n; ++i )
    {
      tin[i] = x(i);
    }

    decolumnize( tin, width, height );

    for( int i = 0; i < n; ++i )
    {
      s[3*i+k] = tin[i];
    }
  }

  delete [] D;
  delete [] tin;
}

void generateEdgeWeights( string inputName )
{
  VImage edgeImage = VImage::vipsload( (char *)inputName.c_str() ).canny();//VImage::option()->set( "sigma", 1 ));//.colourspace(VIPS_INTERPRETATION_B_W);

  float * edgeData1 = ( float * )edgeImage.data();

  int width = edgeImage.width();
  int height = edgeImage.height();

  unsigned char * edgeData2 = new unsigned char[width*height];

  float m = 0;

  for( int p = 0; p < width*height; ++p )
  {
    edgeData2[p] = 0;
    if( edgeData1[3*p+0] > 1.0 ) edgeData2[p] = 255;
    if( edgeData1[3*p+1] > 1.0 ) edgeData2[p] = 255;
    if( edgeData1[3*p+2] > 1.0 ) edgeData2[p] = 255;
  }

  VImage::new_from_memory( edgeData2, width*height, width, height, 1, VIPS_FORMAT_UCHAR ).vipssave("edge.png");

}

void RunEdgeSmooth( string inputName, string outputName, double lambda = 0.01f, double sigma1 = 1, double sigma2 = 3, int K = 1, double dec = 2.0f )
{
  VImage image = VImage::vipsload( (char *)inputName.c_str() ).autorot();//.colourspace(VIPS_INTERPRETATION_B_W);

  // Convert to a three band image
  if( image.bands() == 1 )
  {
    image = image.bandjoin(image).bandjoin(image);
  }
  if( image.bands() == 4 )
  {
    image = image.flatten();
  }

  int width = image.width();
  int height = image.height();

  unsigned char * c = ( unsigned char * )image.data();
  double * c2 = new double[3*width*height];
  double * s = new double[3*width*height];

  for( int p = 0; p < 3*width * height; ++p )
  {
    c2[p] = c[p];
    c2[p] /= 255.0f;
    s[p] = c2[p];
  }

  double * wx = new double[width*height];
  double * wy = new double[width*height];

  ProgressBar *ComputeEdge = new ProgressBar(K, "Computing Edge Smoothing");

  ComputeEdge->Progressed(0);

  for( int k1 = 0; k1 < K; ++k1 )
  {
    computeReWeights( s, sigma1, sigma2, wx, wy, width, height );
    solveLinearEquation( s, c2, wx, wy, lambda, width, height );
    sigma1 /= dec;
    sigma2 /= dec;
    ComputeEdge->Increment();
  }

  ComputeEdge->Finish();

  for( int p = 0; p < 3*width * height; ++p )
  {
    int i = round(s[p] * 255.0);
    i = min( i, 255 );
    i = max( i, 0 );
    c[p] = i;
  }

  VImage::new_from_memory( c, 3*width*height, width, height, 3, VIPS_FORMAT_UCHAR ).vipssave((char *)outputName.c_str());
}