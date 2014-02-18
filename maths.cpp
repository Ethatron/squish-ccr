/* -----------------------------------------------------------------------------

	Copyright (c) 2006 Simon Brown                          si@sjbrown.co.uk
        Copyright (c) 2006 Ignacio Castano                   icastano@nvidia.com
	Copyright (c) 2012 Niels Fröhling              niels@paradice-insight.us

	Permission is hereby granted, free of charge, to any person obtaining
	a copy of this software and associated documentation files (the
	"Software"), to	deal in the Software without restriction, including
	without limitation the rights to use, copy, modify, merge, publish,
	distribute, sublicense, and/or sell copies of the Software, and to
	permit persons to whom the Software is furnished to do so, subject to
	the following conditions:

	The above copyright notice and this permission notice shall be included
	in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
	OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
	CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
	TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

   -------------------------------------------------------------------------- */

/*! @file

	The symmetric eigensystem solver algorithm is from
	http://www.geometrictools.com/Documentation/EigenSymmetric3x3.pdf

	The power-method is from nvtt
	http://code.google.com/p/nvidia-texture-tools/
*/

#include "maths.h"
#include "simd.h"

#undef max
#undef min

namespace squish {

/* *****************************************************************************
 */

#ifdef FEATURE_METRIC_COVARIANCE
#define CoVar(a,b)  (a - b) * metric
#else
#define CoVar(a,b)  (a - b)
#endif

#if	!defined(SQUISH_USE_PRE)
void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec3 &centroid, int n, Vec3 const* points, Vec3 const &metric) {
  // compute the centroid
  Vec3 center = Vec3(0.0f);

  for (int i = 0; i < n; ++i)
    center += points[i];

  center /= n;

  // accumulate the covariance smatrix
  Vec3 covariance_035 = Vec3(0.0f);
  Vec3 covariance_14 = Vec3(0.0f);
  Vec3 covariance_2 = Vec3(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec3 a = CoVar(points[i], center);
    Vec3 b = a;
    Vec3 c = a * b;
    Vec3 d = a * RotateLeft<1>(b);
    Vec3 e = a * RotateLeft<2>(b);

    covariance_035 += c;
    covariance_14 += d;
    covariance_2 += e;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_035.StoreX(&covariance[0]);
  covariance_035.StoreY(&covariance[3]);
  covariance_035.StoreZ(&covariance[5]);
  covariance_14.StoreX(&covariance[1]);
  covariance_14.StoreY(&covariance[4]);
  covariance_2.StoreX(&covariance[2]);
}

void ComputeWeightedCovariance2(Sym2x2 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric) {
  // compute the centroid
  Vec4 center = Vec4(0.0f);

  for (int i = 0; i < n; ++i)
    center += points[i];

  center /= n;

  // accumulate the covariance smatrix
  Vec4 covariance_02 = Vec4(0.0f);
  Vec4 covariance_1 = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec4 a = CoVar(points[i], center);
    Vec4 b = a;
    Vec4 c = a * b;
    Vec4 d = a * RotateLeft<1>(b);

    covariance_02 += c;
    covariance_1 += d;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_02.StoreX(&covariance[0]);
  covariance_02.StoreY(&covariance[2]);
  covariance_1.StoreX(&covariance[1]);
}

void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric) {
  // compute the centroid
  Vec4 center = Vec4(0.0f);

  for (int i = 0; i < n; ++i)
    center += points[i];

  center /= n;

  // accumulate the covariance smatrix
  Vec4 covariance_035 = Vec4(0.0f);
  Vec4 covariance_14 = Vec4(0.0f);
  Vec4 covariance_2 = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec4 a = CoVar(points[i], center);
    Vec4 b = a;
    Vec4 c = a * b;
    Vec4 d = a * RotateLeft<1>(b);
    Vec4 e = a * RotateLeft<2>(b);

    covariance_035 += c;
    covariance_14 += d;
    covariance_2 += e;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_035.StoreX(&covariance[0]);
  covariance_035.StoreY(&covariance[3]);
  covariance_035.StoreZ(&covariance[5]);
  covariance_14.StoreX(&covariance[1]);
  covariance_14.StoreY(&covariance[4]);
  covariance_2.StoreX(&covariance[2]);
}

void ComputeWeightedCovariance4(Sym4x4 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric) {
  // compute the centroid
  Vec4 center = Vec4(0.0f);

  for (int i = 0; i < n; ++i)
    center += points[i];

  center /= n;

  // accumulate the covariance smatrix
  Vec4 covariance_0479 = Vec4(0.0f);
  Vec4 covariance_158 = Vec4(0.0f);
  Vec4 covariance_26 = Vec4(0.0f);
  Vec4 covariance_3 = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec4 a = CoVar(points[i], center);
    Vec4 b = a;
    Vec4 c = a * b;
    Vec4 d = a * RotateLeft<1>(b);
    Vec4 e = a * RotateLeft<2>(b);
    Vec4 f = a * RotateLeft<3>(b);

    covariance_0479 += c;
    covariance_158 += d;
    covariance_26 += e;
    covariance_3 += f;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_0479.StoreX(&covariance[0]);
  covariance_0479.StoreY(&covariance[4]);
  covariance_0479.StoreZ(&covariance[7]);
  covariance_0479.StoreW(&covariance[9]);
  covariance_158.StoreX(&covariance[1]);
  covariance_158.StoreY(&covariance[5]);
  covariance_158.StoreZ(&covariance[8]);
  covariance_26.StoreX(&covariance[2]);
  covariance_26.StoreY(&covariance[6]);
  covariance_3.StoreX(&covariance[3]);
}

void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec3 &centroid, int n, Vec3 const* points, Vec3 const &metric, Scr3 const* weights) {
  // compute the centroid
  Scr3 total = Scr3(0.0f);
  Vec3 center = Vec3(0.0f);

  for (int i = 0; i < n; ++i) {
    total  += weights[i];
    center += weights[i] * points[i];
  }

  center /= total;

  // accumulate the covariance smatrix
  Vec3 covariance_035 = Vec3(0.0f);
  Vec3 covariance_14 = Vec3(0.0f);
  Vec3 covariance_2 = Vec3(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec3 a = CoVar(points[i], center);
    Vec3 b = weights[i] * a;
    Vec3 c = a * b;
    Vec3 d = a * RotateLeft<1>(b);
    Vec3 e = a * RotateLeft<2>(b);

    covariance_035 += c;
    covariance_14 += d;
    covariance_2 += e;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_035.StoreX(&covariance[0]);
  covariance_035.StoreY(&covariance[3]);
  covariance_035.StoreZ(&covariance[5]);
  covariance_14.StoreX(&covariance[1]);
  covariance_14.StoreY(&covariance[4]);
  covariance_2.StoreX(&covariance[2]);
}

void ComputeWeightedCovariance2(Sym2x2 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, Vec4 const* weights) {
  // compute the centroid
  Vec4 total = Vec4(0.0f);
  Vec4 center = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    total  += weights[i];
    center += weights[i] * points[i];
  }

  center /= total;

  // accumulate the covariance smatrix
  Vec4 covariance_02 = Vec4(0.0f);
  Vec4 covariance_1 = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec4 a = CoVar(points[i], center);
    Vec4 b = weights[i] * a;
    Vec4 c = a * b;
    Vec4 d = a * RotateLeft<1>(b);
    Vec4 e = a * RotateLeft<2>(b);

    covariance_02 += c;
    covariance_1 += d;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_02.StoreX(&covariance[0]);
  covariance_02.StoreY(&covariance[2]);
  covariance_1.StoreX(&covariance[1]);
}

void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, Vec4 const* weights) {
  // compute the centroid
  Vec4 total  = Vec4(0.0f);
  Vec4 center = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    total  += weights[i];
    center += weights[i] * points[i];
  }

  center /= total;

  // accumulate the covariance smatrix
  Vec4 covariance_035 = Vec4(0.0f);
  Vec4 covariance_14 = Vec4(0.0f);
  Vec4 covariance_2 = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec4 a = CoVar(points[i], center);
    Vec4 b = weights[i] * a;
    Vec4 c = a * b;
    Vec4 d = a * RotateLeft<1>(b);
    Vec4 e = a * RotateLeft<2>(b);

    covariance_035 += c;
    covariance_14 += d;
    covariance_2 += e;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_035.StoreX(&covariance[0]);
  covariance_035.StoreY(&covariance[3]);
  covariance_035.StoreZ(&covariance[5]);
  covariance_14.StoreX(&covariance[1]);
  covariance_14.StoreY(&covariance[4]);
  covariance_2.StoreX(&covariance[2]);
}

void ComputeWeightedCovariance4(Sym4x4 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, Vec4 const* weights) {
  // compute the centroid
  Vec4 total  = Vec4(0.0f);
  Vec4 center = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    total    += weights[i];
    center += weights[i] * points[i];
  }

  center /= total;

  // accumulate the covariance smatrix
  Vec4 covariance_0479 = Vec4(0.0f);
  Vec4 covariance_158 = Vec4(0.0f);
  Vec4 covariance_26 = Vec4(0.0f);
  Vec4 covariance_3 = Vec4(0.0f);

  for (int i = 0; i < n; ++i) {
    Vec4 a = CoVar(points[i], center);
    Vec4 b = weights[i] * a;
    Vec4 c = a * b;
    Vec4 d = a * RotateLeft<1>(b);
    Vec4 e = a * RotateLeft<2>(b);
    Vec4 f = a * RotateLeft<3>(b);

    covariance_0479 += c;
    covariance_158 += d;
    covariance_26 += e;
    covariance_3 += f;
  }

  // save the centroid
  centroid = center;

  // save the covariance smatrix (TODO: swizzled store)
  covariance_0479.StoreX(&covariance[0]);
  covariance_0479.StoreY(&covariance[4]);
  covariance_0479.StoreZ(&covariance[7]);
  covariance_0479.StoreW(&covariance[9]);
  covariance_158.StoreX(&covariance[1]);
  covariance_158.StoreY(&covariance[5]);
  covariance_158.StoreZ(&covariance[8]);
  covariance_26.StoreX(&covariance[2]);
  covariance_26.StoreY(&covariance[6]);
  covariance_3.StoreX(&covariance[3]);
}

/* .............................................................................
 */

template<class Vec>
static void GetMultiplicity1Evector(Sym3x3 const& smatrix, Vec &out, float evalue)
{
  // compute M
  Sym3x3 m;
  m[0] = smatrix[0] - evalue;
  m[1] = smatrix[1];
  m[2] = smatrix[2];
  m[3] = smatrix[3] - evalue;
  m[4] = smatrix[4];
  m[5] = smatrix[5] - evalue;

  // compute U
  Sym3x3 u;
  u[0] = m[3]*m[5] - m[4]*m[4];
  u[1] = m[2]*m[4] - m[1]*m[5];
  u[2] = m[1]*m[4] - m[2]*m[3];
  u[3] = m[0]*m[5] - m[2]*m[2];
  u[4] = m[1]*m[2] - m[4]*m[0];
  u[5] = m[0]*m[3] - m[1]*m[1];

  // find the largest component
  float mc = std::fabs( u[0] );
  int mi = 0;
  for( int i = 1; i < 6; ++i )
  {
    float c = std::fabs( u[i] );
    if( c > mc )
    {
      mc = c;
      mi = i;
    }
  }

  // pick the column with this component
  switch( mi )
  {
    case 0:
      out = Vec( u[0], u[1], u[2] ); break;

    case 1:
    case 3:
      out = Vec( u[1], u[3], u[4] ); break;

    default:
      out = Vec( u[2], u[4], u[5] ); break;
  }
}

template<class Vec>
static void GetMultiplicity2Evector(Sym3x3 const& smatrix, Vec &out, float evalue)
{
  // compute M
  Sym3x3 m;
  m[0] = smatrix[0] - evalue;
  m[1] = smatrix[1];
  m[2] = smatrix[2];
  m[3] = smatrix[3] - evalue;
  m[4] = smatrix[4];
  m[5] = smatrix[5] - evalue;

  // find the largest component
  float mc = std::fabs( m[0] );
  int mi = 0;
  for( int i = 1; i < 6; ++i )
  {
    float c = std::fabs( m[i] );
    if( c > mc )
    {
      mc = c;
      mi = i;
    }
  }

  // pick the first eigenvector based on this index
  switch( mi )
  {
    case 0:
    case 1:
      out = Vec( -m[1], m[0], 0.0f ); break;

    case 2:
      out = Vec( m[2], 0.0f, -m[0] ); break;

    case 3:
    case 4:
      out = Vec( 0.0f, -m[4], m[3] ); break;

    default:
      out = Vec( 0.0f, -m[5], m[4] ); break;
  }
}

void ComputePrincipleComponent(Sym3x3 const& smatrix, Vec3 &out)
{
  // compute the cubic coefficients
  float c0 =
           smatrix[0] * smatrix[3] * smatrix[5]
  + 2.0f * smatrix[1] * smatrix[2] * smatrix[4]
  -        smatrix[0] * smatrix[4] * smatrix[4]
  -        smatrix[3] * smatrix[2] * smatrix[2]
  -        smatrix[5] * smatrix[1] * smatrix[1];
  float c1 =
           smatrix[0] * smatrix[3]
  +        smatrix[0] * smatrix[5]
  +        smatrix[3] * smatrix[5]
  -        smatrix[1] * smatrix[1]
  -        smatrix[2] * smatrix[2]
  -        smatrix[4] * smatrix[4];
  float c2 =
           smatrix[0]
  +        smatrix[3]
  +        smatrix[5];

  // compute the quadratic coefficients
  float a = c1 - ( 1.0f /  3.0f) * c2 * c2;
  float b =      (-2.0f / 27.0f) * c2 * c2 * c2
  +              ( 1.0f /  3.0f) * c1 * c2 - c0;

  // compute the root count check
  float Q = 0.25f * b * b + (1.0f / 27.0f) * a * a * a;

  // test the multiplicity
  if (FLT_EPSILON < Q) {
    // only one root, which implies we have a multiple of the identity
    out = Vec3(1.0f); return; }

  if (Q < -FLT_EPSILON) {
    // three distinct roots
    float theta = std::atan2(math::sqrt(-Q), -0.5f * b);
    float rho = math::sqrt(0.25f * b * b - Q);

//  float rt = std::pow( rho, 1.0f/3.0f );
    float rt = math::cbrt( rho );
    float ct = std::cos( theta/3.0f );
    float st = std::sin( theta/3.0f );

    float l1 = ( 1.0f/3.0f )*c2 + 2.0f*rt*ct;
    float l2 = ( 1.0f/3.0f )*c2 - rt*( ct + ( float )sqrt( 3.0f )*st );
    float l3 = ( 1.0f/3.0f )*c2 - rt*( ct - ( float )sqrt( 3.0f )*st );

    // pick the larger
    if( std::fabs( l2 ) > std::fabs( l1 ) )
      l1 = l2;
    if( std::fabs( l3 ) > std::fabs( l1 ) )
      l1 = l3;

    // get the eigenvector
    GetMultiplicity1Evector( smatrix, out, l1 );
  }
  else // if( -FLT_EPSILON <= Q && Q <= FLT_EPSILON )
  {
    // two roots
    float rt;
    if( b < 0.0f )
//    rt = -std::pow( -0.5f*b, 1.0f/3.0f );
      rt = -math::cbrt( -0.5f*b );
    else
//    rt = std::pow( 0.5f*b, 1.0f/3.0f );
      rt = math::cbrt( 0.5f*b );

    float l1 = ( 1.0f/3.0f )*c2 + rt;		// repeated
    float l2 = ( 1.0f/3.0f )*c2 - 2.0f*rt;

    // get the eigenvector
    if( std::fabs( l1 ) > std::fabs( l2 ) )
      GetMultiplicity2Evector( smatrix, out, l1 );
    else
      GetMultiplicity1Evector( smatrix, out, l2 );
  }
}

void ComputePrincipleComponent(Sym2x2 const& smatrix, Vec4 &out)
{
  out = Vec4(smatrix[0]); abort();
}

void ComputePrincipleComponent(Sym3x3 const& smatrix, Vec4 &out)
{
  // compute the cubic coefficients
  float c0 = 
    smatrix[0] * smatrix[3] * smatrix[5]
  + smatrix[1] * smatrix[2] * smatrix[4] * 2.0f
  - smatrix[0] * smatrix[4] * smatrix[4]
  - smatrix[3] * smatrix[2] * smatrix[2]
  - smatrix[5] * smatrix[1] * smatrix[1];
  float c1 = 
    smatrix[0] * smatrix[3] + smatrix[0] * smatrix[5] + smatrix[3] * smatrix[5]
  - smatrix[1] * smatrix[1] - smatrix[2] * smatrix[2] - smatrix[4] * smatrix[4];
  float c2 = 
    smatrix[0] + smatrix[3] + smatrix[5];

  // compute the quadratic coefficients
  float a = c1 - ( 1.0f/3.0f )*c2*c2;
  float b = ( -2.0f/27.0f )*c2*c2*c2 + ( 1.0f/3.0f )*c1*c2 - c0;

  // compute the root count check
  float Q = 0.25f*b*b + ( 1.0f/27.0f )*a*a*a;

  // test the multiplicity
  if( FLT_EPSILON < Q )
  {
    // only one root, which implies we have a multiple of the identity
    out = Vec4( 1.0f ); return;
  }

  if( Q < -FLT_EPSILON )
  {
    // three distinct roots
    float theta = std::atan2( math::sqrt( -Q ), -0.5f*b );
    float rho = math::sqrt( 0.25f*b*b - Q );

//  float rt = std::pow( rho, 1.0f/3.0f );
    float rt = math::cbrt( rho );
    float ct = std::cos( theta/3.0f );
    float st = std::sin( theta/3.0f );

    float l1 = ( 1.0f/3.0f )*c2 + 2.0f*rt*ct;
    float l2 = ( 1.0f/3.0f )*c2 - rt*( ct + ( float )sqrt( 3.0f )*st );
    float l3 = ( 1.0f/3.0f )*c2 - rt*( ct - ( float )sqrt( 3.0f )*st );

    // pick the larger
    if( std::fabs( l2 ) > std::fabs( l1 ) )
      l1 = l2;
    if( std::fabs( l3 ) > std::fabs( l1 ) )
      l1 = l3;

    // get the eigenvector
    GetMultiplicity1Evector( smatrix, out, l1 );
  }
  else // if( -FLT_EPSILON <= Q && Q <= FLT_EPSILON )
  {
    // two roots
    float rt;
    if( b < 0.0f )
//    rt = -std::pow( -0.5f*b, 1.0f/3.0f );
      rt = -math::cbrt( -0.5f*b );
    else
//    rt = std::pow( 0.5f*b, 1.0f/3.0f );
      rt = math::cbrt( 0.5f*b );

    float l1 = ( 1.0f/3.0f )*c2 + rt;		// repeated
    float l2 = ( 1.0f/3.0f )*c2 - 2.0f*rt;

    // get the eigenvector
    if( std::fabs( l1 ) > std::fabs( l2 ) )
      GetMultiplicity2Evector( smatrix, out, l1 );
    else
      GetMultiplicity1Evector( smatrix, out, l2 );
  }
}

/* .............................................................................
 *         16bit      8bit      8bit
 *
 * [ 0]        0         0         0
 * [ 1]    10855     10855         0
 * [ 2]    13404    792634  37707174
 * [ 3]   147487    571689  30087713
 * [ 4]   301039    341389  18500805
 *   =    472785   1716567  86295692
 * [ 5]   335015    194660  10656029
 * [ 6]   297186    112792   6337250
 * [ 7]   238783     68232   3904835
 * [ 8]   184111     42843   2561885
 *   =   1527880   2135094 109755691
 *
 * [ 9]   140932     28428   1718595
 * [10]   107016     19683   1207264
 * [11]    82490     14083    861012
 * [12]    64155     10192    633528
 * [13]    50592      7341    484082
 * [14]    40557      5815    368445
 * [15]    32226      4433    288080
 * [16]    26161      3554    239111
 * [17]    21705      2867    182662
 * [18]    18087      2360    152417
 * [19]    15273      1889    125116
 * [20]    12975      1575    104299
 * [21]    10988      1343     87505
 * [22]     9442      1135     75360
 * [23]     8345       976     63749
 * [24]     7246       803     56561
 * [25]     6122       734     52284
 * [26]     5474       605     42680
 * [27]     4805       556     36671
 * [28]     4327       464     33079
 * [29]     3875       411     32825
 * [30]     3491       369     25787
 * [31]     3116       334     23838
 * [32]     2816       286     21195
 * [33]     2530       239     19282
 * [34]     2258       261     18245
 * [35]     2064       226     16229
 * [36]     1843       182     14444
 * [37]     1774       189     13264
 * [38]     1599       172     12364
 * [39]     1453       150     10995
 * [40]     1369       129     11121
 * [41]     1197       128      9621
 * [42]     1129       106      8961
 * [43]     1131       113      7806
 * [44]     1036       105      7406
 * [45]      930        93      6796
 * [46]      895        72      6156
 * [47]      814        90      6051
 * [48]      718        86      5969
 * [49]      693        71      5102
 * [50]      710        57      4841
 * [51]      634        67      4426
 * [52]      594        58      4358
 * [53]      576        61      4100
 * [54]      533        54      3510
 * [55]      458        53      3579
 * [56]      451        35      3164
 * [57]      434        46      2948
 * [58]      441        42      2851
 * [59]      355        44      2767
 * [60]      373        52      2713
 * [61]      396        37      2415
 * [62]      327        48      2204
 * [63]     9401       816     47776
 *   =    721332    114118   7187609
 */
#define POWER_ITERATION_COUNT	-1
#define POWER_ITERATION_PREC	1.0f / 256

void EstimatePrincipleComponent(Sym2x2 const& matrix, Vec4 &out)
{
//Vec4 const row0(matrix[0], matrix[1], 0.0f, 0.0f);
//Vec4 const row1(matrix[1], matrix[2], 0.0f, 0.0f);
  Vec4 row0; 
  Vec4 row1; 
  Vec4 v;
  
  LoadUnaligned(row0, matrix(0));
  
  row1 = Merge<1,2,1,2>(row0, row0) & Vec4(true, true, false, false);
  row0 =                row0        & Vec4(true, true, false, false);
  
  assert(row0.GetO(0) == matrix[0]); assert(row0.GetO(1) == matrix[1]); assert(row0.GetO(2) == 0.0f); assert(row0.GetO(3) == 0.0f);
  assert(row1.GetO(0) == matrix[1]); assert(row1.GetO(1) == matrix[2]); assert(row1.GetO(2) == 0.0f); assert(row1.GetO(3) == 0.0f);

  Scr4 r0 = LengthSquared(row0);
  Scr4 r1 = LengthSquared(row1);

  if (r0 > r1) v = row0;
  else v = row1;

#if POWER_ITERATION_COUNT > 0
  for (int i = 0; i < POWER_ITERATION_COUNT; i++) {
#else
  int i = 0; Vec4 d; do { d = v; i++;
#endif
    Scr4 x = Dot(v, row0);
    Scr4 y = Dot(v, row1);

    v  = Vec4(x, y);
    v *= Reciprocal(HorizontalMax(Abs(v)));
  }
#if POWER_ITERATION_COUNT <= 0
  while (CompareAnyGreaterThan(AbsoluteDifference(v, d), Vec4(POWER_ITERATION_PREC)) && (i < 64));
#endif

#if defined(TRACK_STATISTICS)
  gstat.num_poweritrs[std::min(i, 63)]++;
#endif
  
  assert(v.GetZ() == 0.0f);
  assert(v.GetW() == 0.0f);
  out = v;
}

void EstimatePrincipleComponent(Sym3x3 const& matrix, Vec3 &out)
{
  Vec4 tmp; EstimatePrincipleComponent(matrix, tmp); out = tmp.GetVec3();
}

void EstimatePrincipleComponent(Sym3x3 const& matrix, Vec4 &out)
{
//Vec4 const row0(matrix[0], matrix[1], matrix[2], 0.0f);
//Vec4 const row1(matrix[1], matrix[3], matrix[4], 0.0f);
//Vec4 const row2(matrix[2], matrix[4], matrix[5], 0.0f);
  Vec4 row0;
  Vec4 row1;
  Vec4 row2;
  Vec4 v;

  LoadUnaligned(row0, matrix(0));
  LoadUnaligned(row1, matrix(2));

  row2 = Merge<0,2,3,3>(row1, row1) & Vec4(true, true, true, false);
  row1 = Merge<1,3,2,2>(row0, row1) & Vec4(true, true, true, false);
  row0 =                row0        & Vec4(true, true, true, false);

  assert(row0.GetO(0) == matrix[0]); assert(row0.GetO(1) == matrix[1]); assert(row0.GetO(2) == matrix[2]); assert(row0.GetO(3) == 0.0f);
  assert(row1.GetO(0) == matrix[1]); assert(row1.GetO(1) == matrix[3]); assert(row1.GetO(2) == matrix[4]); assert(row1.GetO(3) == 0.0f);
  assert(row2.GetO(0) == matrix[2]); assert(row2.GetO(1) == matrix[4]); assert(row2.GetO(2) == matrix[5]); assert(row2.GetO(3) == 0.0f);

  Scr4 r0 = LengthSquared(row0);
  Scr4 r1 = LengthSquared(row1);
  Scr4 r2 = LengthSquared(row2);

  if (r0 > r1 && r0 > r2) v = row0;
  else if (r1 > r2) v = row1;
  else v = row2;

#if POWER_ITERATION_COUNT > 0
  for (int i = 0; i < POWER_ITERATION_COUNT; i++) {
#else
  int i = 0; Vec4 d; do { d = v; i++;
#endif
    Scr4 x = Dot(v, row0);
    Scr4 y = Dot(v, row1);
    Scr4 z = Dot(v, row2);

    v  = Vec4(x, y, z);
    v *= Reciprocal(HorizontalMax(Abs(v)));
  }
#if POWER_ITERATION_COUNT <= 0
  while (CompareAnyGreaterThan(AbsoluteDifference(v, d), Vec4(POWER_ITERATION_PREC)) && (i < 64));
#endif

#if defined(TRACK_STATISTICS)
  gstat.num_poweritrs[std::min(i, 63)]++;
#endif

  assert(v.GetW() == 0.0f);
  out = v;
}

void EstimatePrincipleComponent(Sym4x4 const& matrix, Vec4 &out)
{
  Vec4 const row0(matrix[0], matrix[1], matrix[2], matrix[3]);
  Vec4 const row1(matrix[1], matrix[4], matrix[5], matrix[6]);
  Vec4 const row2(matrix[2], matrix[5], matrix[7], matrix[8]);
  Vec4 const row3(matrix[3], matrix[6], matrix[8], matrix[9]);
  Vec4 v;

  Scr4 r0 = LengthSquared(row0);
  Scr4 r1 = LengthSquared(row1);
  Scr4 r2 = LengthSquared(row2);
  Scr4 r3 = LengthSquared(row3);

  if (r0 > r1 && r0 > r2 && r0 > r3) v = row0;
  else if (r1 > r2 && r1 > r3) v = row1;
  else if (r2 > r3) v = row2;
  else v = row3;

#if POWER_ITERATION_COUNT > 0
  for (int i = 0; i < POWER_ITERATION_COUNT; i++) {
#else
  int i = 0; Vec4 d; do { d = v; i++;
#endif
    Scr4 x = Dot(v, row0);
    Scr4 y = Dot(v, row1);
    Scr4 z = Dot(v, row2);
    Scr4 w = Dot(v, row3);

    v  = Vec4(x, y, z, w);
    v *= Reciprocal(HorizontalMax(Abs(v)));
  }
#if POWER_ITERATION_COUNT <= 0
  while (CompareAnyGreaterThan(AbsoluteDifference(v, d), Vec4(POWER_ITERATION_PREC)) && (i < 64));
#endif

#if defined(TRACK_STATISTICS)
  gstat.num_poweritrs[std::min(i, 63)]++;
#endif

  out = v;
}

/* -----------------------------------------------------------------------------
 */
template<class VecX, class ScrX>
void GetPrincipleProjection(VecX &enter, VecX &leave, VecX const &principle, VecX const &centroid, int n, VecX const* points)
{
  ScrX div = Reciprocal(Dot(principle, principle));
  VecX rec = Reciprocal(    principle            );
  ScrX len, min, max;
  VecX chk;

  // compute the projection
  min = max = Dot(points[0] - centroid, principle);

  for (int i = 1; i < n; ++i) {
    len = Dot(points[i] - centroid, principle);
    min = Min(min, len);
    max = Max(max, len);
  }

  VecX start = centroid + principle * min * div;
  VecX end   = centroid + principle * max * div;
    
  // take care that low magnitude overshoots can have very high
  // derivatives on the other axis (-0.0039 in R may be +1 in G,
  // thus we at least go to -0.0039/255 to get +1/255 -> 1/255²)

  // intersect negative undershoot with axis-plane(s), clamp to 0.0
  chk = start;
  while (CompareAnyLessThan(chk, VecX(-1.0f / (255.0f * 255.0f)))) {
    VecX fct = chk * rec;
    VecX hin = Select(fct, chk, HorizontalMin(chk));

    start -= principle * hin;
    chk = start;
  }

  // intersect negative undershoot with axis-plane(s), clamp to 0.0
  chk = end;
  while (CompareAnyLessThan(chk, VecX(-1.0f / (255.0f * 255.0f)))) {
    VecX fct = chk * rec;
    VecX hin = Select(fct, chk, HorizontalMin(chk));

    end -= principle * hin;
    chk = end;
  }

  // intersect positive overshoot with axis-plane(s), clamp to 1.0
  chk = start - VecX(1.0f);
  while (CompareAnyGreaterThan(chk, VecX(1.0f / (255.0f * 255.0f)))) {
    VecX fct = chk * rec;
    VecX hax = Select(fct, chk, HorizontalMax(chk));

    start -= principle * hax;
    chk = start - VecX(1.0f);
  }

  // intersect positive overshoot with axis-plane(s), clamp to 1.0
  chk = end - VecX(1.0f);
  while (CompareAnyGreaterThan(chk, VecX(1.0f / (255.0f * 255.0f)))) {
    VecX fct = chk * rec;
    VecX hax = Select(fct, chk, HorizontalMax(chk));

    end -= principle * hax;
    chk = end - VecX(1.0f);
  }

/*assert(HorizontalMin(start).X() > -0.0001);
  assert(HorizontalMin(end  ).X() > -0.0001);
  assert(HorizontalMax(start).X() <  1.0001);
  assert(HorizontalMax(end  ).X() <  1.0001);	 */

  enter = start;
  leave = end;
}

void GetPrincipleProjection(Vec3 &enter, Vec3 &leave, Vec3 const &principle, Vec3 const &centroid, int n, Vec3 const* points)
{
  GetPrincipleProjection<Vec3, Scr3>(enter, leave, principle, centroid, n, points);
}

void GetPrincipleProjection(Vec4 &enter, Vec4 &leave, Vec4 const &principle, Vec4 const &centroid, int n, Vec4 const* points)
{
  GetPrincipleProjection<Vec4, Scr4>(enter, leave, principle, centroid, n, points);
}

/* -----------------------------------------------------------------------------
 * float q = LUTindex / ((1 << b) - 1);
 * int t = (int)floor(q * 255.0f) >> (8 - b);
 * int d = (t << (8 - b)); d |= (d >> (b)); d |= (d >> (2 * b)); d |= (d >> (3 * b));
 * LUTvalue = d / 255.0f;
 */
const a16 float qLUT_1all[2] = {
    0.0f / 255.0f, 255.0f / 255.0f};
const a16 float qLUT_1clr[2] = {
    0.0f / 255.0f,   0.0f / 255.0f};
const a16 float qLUT_1set[2] = {
  255.0f / 255.0f, 255.0f / 255.0f};

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
const a16 float qLUT_2all[4] = {
    0.0f / 255.0f,  85.0f / 255.0f, 170.0f / 255.0f, 255.0f / 255.0f};
const a16 float qLUT_2clr[4] = {
    0.0f / 255.0f,   0.0f / 255.0f, 170.0f / 255.0f, 170.0f / 255.0f};
const a16 float qLUT_2set[4] = {
   85.0f / 255.0f,  85.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f};

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
const a16 float qLUT_3all[8] = {
    0.0f / 255.0f,  36.0f / 255.0f,  73.0f / 255.0f, 109.0f / 255.0f, 146.0f / 255.0f, 182.0f / 255.0f, 219.0f / 255.0f, 255.0f / 255.0f,
};

const a16 float qLUT_3clr[8] = {
    0.0f / 255.0f,  73.0f / 255.0f,  73.0f / 255.0f,  73.0f / 255.0f, 146.0f / 255.0f, 146.0f / 255.0f, 219.0f / 255.0f, 219.0f / 255.0f,
};

const a16 float qLUT_3set[8] = {
   36.0f / 255.0f,  36.0f / 255.0f, 109.0f / 255.0f, 109.0f / 255.0f, 182.0f / 255.0f, 182.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f,
};

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
const a16 float qLUT_4all[16] = {
    0.0f / 255.0f,  17.0f / 255.0f,  34.0f / 255.0f,  51.0f / 255.0f,  68.0f / 255.0f,  85.0f / 255.0f, 102.0f / 255.0f, 119.0f / 255.0f,
  136.0f / 255.0f, 153.0f / 255.0f, 170.0f / 255.0f, 187.0f / 255.0f, 204.0f / 255.0f, 221.0f / 255.0f, 238.0f / 255.0f, 255.0f / 255.0f,
};

const a16 float qLUT_4clr[16] = {
    0.0f / 255.0f,  34.0f / 255.0f,  34.0f / 255.0f,  68.0f / 255.0f,  68.0f / 255.0f,  68.0f / 255.0f, 102.0f / 255.0f, 102.0f / 255.0f,
  136.0f / 255.0f, 136.0f / 255.0f, 170.0f / 255.0f, 170.0f / 255.0f, 204.0f / 255.0f, 204.0f / 255.0f, 238.0f / 255.0f, 238.0f / 255.0f,
};

const a16 float qLUT_4set[16] = {
   17.0f / 255.0f,  17.0f / 255.0f,  51.0f / 255.0f,  51.0f / 255.0f,  85.0f / 255.0f,  85.0f / 255.0f, 119.0f / 255.0f, 119.0f / 255.0f,
  153.0f / 255.0f, 153.0f / 255.0f, 187.0f / 255.0f, 187.0f / 255.0f, 221.0f / 255.0f, 221.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f,
};

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
const a16 float qLUT_5all[32] = {
    0.0f / 255.0f,   8.0f / 255.0f,  16.0f / 255.0f,  24.0f / 255.0f,  33.0f / 255.0f,  41.0f / 255.0f,  49.0f / 255.0f,  57.0f / 255.0f,
   66.0f / 255.0f,  74.0f / 255.0f,  82.0f / 255.0f,  90.0f / 255.0f,  99.0f / 255.0f, 107.0f / 255.0f, 115.0f / 255.0f, 123.0f / 255.0f,
  132.0f / 255.0f, 140.0f / 255.0f, 148.0f / 255.0f, 156.0f / 255.0f, 165.0f / 255.0f, 173.0f / 255.0f, 181.0f / 255.0f, 189.0f / 255.0f,
  198.0f / 255.0f, 206.0f / 255.0f, 214.0f / 255.0f, 222.0f / 255.0f, 231.0f / 255.0f, 239.0f / 255.0f, 247.0f / 255.0f, 255.0f / 255.0f,
};

const a16 float qLUT_5clr[32] = {
    0.0f / 255.0f,  16.0f / 255.0f,  16.0f / 255.0f,  33.0f / 255.0f,  33.0f / 255.0f,  49.0f / 255.0f,  49.0f / 255.0f,  66.0f / 255.0f,
   66.0f / 255.0f,  82.0f / 255.0f,  82.0f / 255.0f,  99.0f / 255.0f,  99.0f / 255.0f, 115.0f / 255.0f, 115.0f / 255.0f, 115.0f / 255.0f,
  132.0f / 255.0f, 132.0f / 255.0f, 148.0f / 255.0f, 148.0f / 255.0f, 165.0f / 255.0f, 165.0f / 255.0f, 181.0f / 255.0f, 181.0f / 255.0f,
  198.0f / 255.0f, 198.0f / 255.0f, 214.0f / 255.0f, 214.0f / 255.0f, 231.0f / 255.0f, 231.0f / 255.0f, 247.0f / 255.0f, 247.0f / 255.0f,
};

const a16 float qLUT_5set[32] = {
    8.0f / 255.0f,   8.0f / 255.0f,  24.0f / 255.0f,  24.0f / 255.0f,  41.0f / 255.0f,  41.0f / 255.0f,  57.0f / 255.0f,  57.0f / 255.0f,
   74.0f / 255.0f,  74.0f / 255.0f,  90.0f / 255.0f,  90.0f / 255.0f, 107.0f / 255.0f, 107.0f / 255.0f, 123.0f / 255.0f, 123.0f / 255.0f,
  140.0f / 255.0f, 140.0f / 255.0f, 156.0f / 255.0f, 156.0f / 255.0f, 173.0f / 255.0f, 173.0f / 255.0f, 189.0f / 255.0f, 189.0f / 255.0f,
  206.0f / 255.0f, 206.0f / 255.0f, 222.0f / 255.0f, 222.0f / 255.0f, 239.0f / 255.0f, 239.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f,
};

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
const a16 float qLUT_6all[64] = {
    0.0f / 255.0f,   4.0f / 255.0f,   8.0f / 255.0f,  12.0f / 255.0f,  16.0f / 255.0f,  20.0f / 255.0f,  24.0f / 255.0f,  28.0f / 255.0f,
   32.0f / 255.0f,  36.0f / 255.0f,  40.0f / 255.0f,  44.0f / 255.0f,  48.0f / 255.0f,  52.0f / 255.0f,  56.0f / 255.0f,  60.0f / 255.0f,
   65.0f / 255.0f,  69.0f / 255.0f,  73.0f / 255.0f,  77.0f / 255.0f,  81.0f / 255.0f,  85.0f / 255.0f,  89.0f / 255.0f,  93.0f / 255.0f,
   97.0f / 255.0f, 101.0f / 255.0f, 105.0f / 255.0f, 109.0f / 255.0f, 113.0f / 255.0f, 117.0f / 255.0f, 121.0f / 255.0f, 125.0f / 255.0f,
  130.0f / 255.0f, 134.0f / 255.0f, 138.0f / 255.0f, 142.0f / 255.0f, 146.0f / 255.0f, 150.0f / 255.0f, 154.0f / 255.0f, 158.0f / 255.0f,
  162.0f / 255.0f, 166.0f / 255.0f, 170.0f / 255.0f, 174.0f / 255.0f, 178.0f / 255.0f, 182.0f / 255.0f, 186.0f / 255.0f, 190.0f / 255.0f,
  195.0f / 255.0f, 199.0f / 255.0f, 203.0f / 255.0f, 207.0f / 255.0f, 211.0f / 255.0f, 215.0f / 255.0f, 219.0f / 255.0f, 223.0f / 255.0f,
  227.0f / 255.0f, 231.0f / 255.0f, 235.0f / 255.0f, 239.0f / 255.0f, 243.0f / 255.0f, 247.0f / 255.0f, 251.0f / 255.0f, 255.0f / 255.0f,
};

const a16 float qLUT_6clr[64] = {
    0.0f / 255.0f,   8.0f / 255.0f,   8.0f / 255.0f,  16.0f / 255.0f,  16.0f / 255.0f,  24.0f / 255.0f,  24.0f / 255.0f,  32.0f / 255.0f,
   32.0f / 255.0f,  40.0f / 255.0f,  40.0f / 255.0f,  48.0f / 255.0f,  48.0f / 255.0f,  56.0f / 255.0f,  56.0f / 255.0f,  65.0f / 255.0f,
   65.0f / 255.0f,  65.0f / 255.0f,  73.0f / 255.0f,  73.0f / 255.0f,  81.0f / 255.0f,  89.0f / 255.0f,  89.0f / 255.0f,  97.0f / 255.0f,
   97.0f / 255.0f, 105.0f / 255.0f, 105.0f / 255.0f, 105.0f / 255.0f, 113.0f / 255.0f, 113.0f / 255.0f, 121.0f / 255.0f, 121.0f / 255.0f,
  130.0f / 255.0f, 130.0f / 255.0f, 138.0f / 255.0f, 138.0f / 255.0f, 146.0f / 255.0f, 146.0f / 255.0f, 154.0f / 255.0f, 154.0f / 255.0f,
  162.0f / 255.0f, 162.0f / 255.0f, 170.0f / 255.0f, 170.0f / 255.0f, 178.0f / 255.0f, 178.0f / 255.0f, 186.0f / 255.0f, 186.0f / 255.0f,
  195.0f / 255.0f, 195.0f / 255.0f, 203.0f / 255.0f, 203.0f / 255.0f, 211.0f / 255.0f, 211.0f / 255.0f, 219.0f / 255.0f, 219.0f / 255.0f,
  227.0f / 255.0f, 227.0f / 255.0f, 235.0f / 255.0f, 235.0f / 255.0f, 243.0f / 255.0f, 243.0f / 255.0f, 251.0f / 255.0f, 251.0f / 255.0f,
};

const a16 float qLUT_6set[64] = {
    4.0f / 255.0f,   4.0f / 255.0f,  12.0f / 255.0f,  12.0f / 255.0f,  20.0f / 255.0f,  20.0f / 255.0f,  28.0f / 255.0f,  28.0f / 255.0f,
   36.0f / 255.0f,  36.0f / 255.0f,  44.0f / 255.0f,  44.0f / 255.0f,  52.0f / 255.0f,  52.0f / 255.0f,  60.0f / 255.0f,  60.0f / 255.0f,
   69.0f / 255.0f,  69.0f / 255.0f,  77.0f / 255.0f,  77.0f / 255.0f,  85.0f / 255.0f,  85.0f / 255.0f,  93.0f / 255.0f,  93.0f / 255.0f,
  101.0f / 255.0f, 101.0f / 255.0f, 109.0f / 255.0f, 109.0f / 255.0f, 117.0f / 255.0f, 117.0f / 255.0f, 125.0f / 255.0f, 125.0f / 255.0f,
  134.0f / 255.0f, 134.0f / 255.0f, 142.0f / 255.0f, 142.0f / 255.0f, 150.0f / 255.0f, 150.0f / 255.0f, 158.0f / 255.0f, 158.0f / 255.0f,
  166.0f / 255.0f, 166.0f / 255.0f, 174.0f / 255.0f, 174.0f / 255.0f, 182.0f / 255.0f, 182.0f / 255.0f, 190.0f / 255.0f, 190.0f / 255.0f,
  199.0f / 255.0f, 199.0f / 255.0f, 207.0f / 255.0f, 207.0f / 255.0f, 215.0f / 255.0f, 215.0f / 255.0f, 223.0f / 255.0f, 223.0f / 255.0f,
  231.0f / 255.0f, 231.0f / 255.0f, 239.0f / 255.0f, 239.0f / 255.0f, 247.0f / 255.0f, 247.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f,
};

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
const a16 float qLUT_7all[128] = {
    0.0f / 255.0f,   2.0f / 255.0f,   4.0f / 255.0f,   6.0f / 255.0f,   8.0f / 255.0f,  10.0f / 255.0f,  12.0f / 255.0f,  14.0f / 255.0f,
   16.0f / 255.0f,  18.0f / 255.0f,  20.0f / 255.0f,  22.0f / 255.0f,  24.0f / 255.0f,  26.0f / 255.0f,  28.0f / 255.0f,  30.0f / 255.0f,
   32.0f / 255.0f,  34.0f / 255.0f,  36.0f / 255.0f,  38.0f / 255.0f,  40.0f / 255.0f,  42.0f / 255.0f,  44.0f / 255.0f,  46.0f / 255.0f,
   48.0f / 255.0f,  50.0f / 255.0f,  52.0f / 255.0f,  54.0f / 255.0f,  56.0f / 255.0f,  58.0f / 255.0f,  60.0f / 255.0f,  62.0f / 255.0f,
   64.0f / 255.0f,  66.0f / 255.0f,  68.0f / 255.0f,  70.0f / 255.0f,  72.0f / 255.0f,  74.0f / 255.0f,  76.0f / 255.0f,  78.0f / 255.0f,
   80.0f / 255.0f,  82.0f / 255.0f,  84.0f / 255.0f,  86.0f / 255.0f,  88.0f / 255.0f,  90.0f / 255.0f,  92.0f / 255.0f,  94.0f / 255.0f,
   96.0f / 255.0f,  98.0f / 255.0f, 100.0f / 255.0f, 102.0f / 255.0f, 104.0f / 255.0f, 106.0f / 255.0f, 108.0f / 255.0f, 110.0f / 255.0f,
  112.0f / 255.0f, 114.0f / 255.0f, 116.0f / 255.0f, 118.0f / 255.0f, 120.0f / 255.0f, 122.0f / 255.0f, 124.0f / 255.0f, 126.0f / 255.0f,
  129.0f / 255.0f, 131.0f / 255.0f, 133.0f / 255.0f, 135.0f / 255.0f, 137.0f / 255.0f, 139.0f / 255.0f, 141.0f / 255.0f, 143.0f / 255.0f,
  145.0f / 255.0f, 147.0f / 255.0f, 149.0f / 255.0f, 151.0f / 255.0f, 153.0f / 255.0f, 155.0f / 255.0f, 157.0f / 255.0f, 159.0f / 255.0f,
  161.0f / 255.0f, 163.0f / 255.0f, 165.0f / 255.0f, 167.0f / 255.0f, 169.0f / 255.0f, 171.0f / 255.0f, 173.0f / 255.0f, 175.0f / 255.0f,
  177.0f / 255.0f, 179.0f / 255.0f, 181.0f / 255.0f, 183.0f / 255.0f, 185.0f / 255.0f, 187.0f / 255.0f, 189.0f / 255.0f, 191.0f / 255.0f,
  193.0f / 255.0f, 195.0f / 255.0f, 197.0f / 255.0f, 199.0f / 255.0f, 201.0f / 255.0f, 203.0f / 255.0f, 205.0f / 255.0f, 207.0f / 255.0f,
  209.0f / 255.0f, 211.0f / 255.0f, 213.0f / 255.0f, 215.0f / 255.0f, 217.0f / 255.0f, 219.0f / 255.0f, 221.0f / 255.0f, 223.0f / 255.0f,
  225.0f / 255.0f, 227.0f / 255.0f, 229.0f / 255.0f, 231.0f / 255.0f, 233.0f / 255.0f, 235.0f / 255.0f, 237.0f / 255.0f, 239.0f / 255.0f,
  241.0f / 255.0f, 243.0f / 255.0f, 245.0f / 255.0f, 247.0f / 255.0f, 249.0f / 255.0f, 251.0f / 255.0f, 253.0f / 255.0f, 255.0f / 255.0f,
};

const a16 float qLUT_7clr[128] = {
    0.0f / 255.0f,   4.0f / 255.0f,   4.0f / 255.0f,   8.0f / 255.0f,   8.0f / 255.0f,  12.0f / 255.0f,  12.0f / 255.0f,  16.0f / 255.0f,
   16.0f / 255.0f,  20.0f / 255.0f,  20.0f / 255.0f,  24.0f / 255.0f,  24.0f / 255.0f,  28.0f / 255.0f,  28.0f / 255.0f,  32.0f / 255.0f,
   32.0f / 255.0f,  36.0f / 255.0f,  36.0f / 255.0f,  40.0f / 255.0f,  40.0f / 255.0f,  44.0f / 255.0f,  44.0f / 255.0f,  48.0f / 255.0f,
   48.0f / 255.0f,  52.0f / 255.0f,  52.0f / 255.0f,  56.0f / 255.0f,  56.0f / 255.0f,  60.0f / 255.0f,  60.0f / 255.0f,  64.0f / 255.0f,
   64.0f / 255.0f,  68.0f / 255.0f,  68.0f / 255.0f,  72.0f / 255.0f,  72.0f / 255.0f,  76.0f / 255.0f,  76.0f / 255.0f,  80.0f / 255.0f,
   80.0f / 255.0f,  84.0f / 255.0f,  84.0f / 255.0f,  84.0f / 255.0f,  88.0f / 255.0f,  88.0f / 255.0f,  92.0f / 255.0f,  92.0f / 255.0f,
   96.0f / 255.0f,  96.0f / 255.0f, 100.0f / 255.0f, 100.0f / 255.0f, 104.0f / 255.0f, 104.0f / 255.0f, 108.0f / 255.0f, 108.0f / 255.0f,
  112.0f / 255.0f, 112.0f / 255.0f, 116.0f / 255.0f, 116.0f / 255.0f, 120.0f / 255.0f, 120.0f / 255.0f, 124.0f / 255.0f, 124.0f / 255.0f,
  129.0f / 255.0f, 129.0f / 255.0f, 133.0f / 255.0f, 133.0f / 255.0f, 137.0f / 255.0f, 137.0f / 255.0f, 141.0f / 255.0f, 141.0f / 255.0f,
  145.0f / 255.0f, 145.0f / 255.0f, 149.0f / 255.0f, 149.0f / 255.0f, 153.0f / 255.0f, 153.0f / 255.0f, 157.0f / 255.0f, 157.0f / 255.0f,
  161.0f / 255.0f, 161.0f / 255.0f, 165.0f / 255.0f, 165.0f / 255.0f, 169.0f / 255.0f, 169.0f / 255.0f, 173.0f / 255.0f, 173.0f / 255.0f,
  177.0f / 255.0f, 177.0f / 255.0f, 181.0f / 255.0f, 181.0f / 255.0f, 185.0f / 255.0f, 185.0f / 255.0f, 189.0f / 255.0f, 189.0f / 255.0f,
  193.0f / 255.0f, 193.0f / 255.0f, 197.0f / 255.0f, 197.0f / 255.0f, 201.0f / 255.0f, 201.0f / 255.0f, 205.0f / 255.0f, 205.0f / 255.0f,
  209.0f / 255.0f, 209.0f / 255.0f, 213.0f / 255.0f, 213.0f / 255.0f, 217.0f / 255.0f, 217.0f / 255.0f, 221.0f / 255.0f, 221.0f / 255.0f,
  225.0f / 255.0f, 225.0f / 255.0f, 229.0f / 255.0f, 229.0f / 255.0f, 233.0f / 255.0f, 233.0f / 255.0f, 237.0f / 255.0f, 237.0f / 255.0f,
  241.0f / 255.0f, 241.0f / 255.0f, 245.0f / 255.0f, 245.0f / 255.0f, 249.0f / 255.0f, 249.0f / 255.0f, 253.0f / 255.0f, 253.0f / 255.0f,
};

const a16 float qLUT_7set[128] = {
    2.0f / 255.0f,   2.0f / 255.0f,   6.0f / 255.0f,   6.0f / 255.0f,  10.0f / 255.0f,  10.0f / 255.0f,  14.0f / 255.0f,  14.0f / 255.0f,
   18.0f / 255.0f,  18.0f / 255.0f,  22.0f / 255.0f,  22.0f / 255.0f,  26.0f / 255.0f,  26.0f / 255.0f,  30.0f / 255.0f,  30.0f / 255.0f,
   34.0f / 255.0f,  34.0f / 255.0f,  38.0f / 255.0f,  38.0f / 255.0f,  42.0f / 255.0f,  42.0f / 255.0f,  46.0f / 255.0f,  46.0f / 255.0f,
   50.0f / 255.0f,  50.0f / 255.0f,  54.0f / 255.0f,  54.0f / 255.0f,  58.0f / 255.0f,  58.0f / 255.0f,  62.0f / 255.0f,  62.0f / 255.0f,
   66.0f / 255.0f,  66.0f / 255.0f,  70.0f / 255.0f,  70.0f / 255.0f,  74.0f / 255.0f,  74.0f / 255.0f,  78.0f / 255.0f,  78.0f / 255.0f,
   82.0f / 255.0f,  82.0f / 255.0f,  86.0f / 255.0f,  86.0f / 255.0f,  90.0f / 255.0f,  90.0f / 255.0f,  94.0f / 255.0f,  94.0f / 255.0f,
   98.0f / 255.0f,  98.0f / 255.0f, 102.0f / 255.0f, 102.0f / 255.0f, 106.0f / 255.0f, 106.0f / 255.0f, 110.0f / 255.0f, 110.0f / 255.0f,
  114.0f / 255.0f, 114.0f / 255.0f, 118.0f / 255.0f, 118.0f / 255.0f, 122.0f / 255.0f, 122.0f / 255.0f, 126.0f / 255.0f, 126.0f / 255.0f,
  131.0f / 255.0f, 131.0f / 255.0f, 135.0f / 255.0f, 135.0f / 255.0f, 139.0f / 255.0f, 139.0f / 255.0f, 143.0f / 255.0f, 143.0f / 255.0f,
  147.0f / 255.0f, 147.0f / 255.0f, 151.0f / 255.0f, 151.0f / 255.0f, 155.0f / 255.0f, 155.0f / 255.0f, 159.0f / 255.0f, 159.0f / 255.0f,
  163.0f / 255.0f, 163.0f / 255.0f, 167.0f / 255.0f, 167.0f / 255.0f, 171.0f / 255.0f, 171.0f / 255.0f, 175.0f / 255.0f, 175.0f / 255.0f,
  179.0f / 255.0f, 179.0f / 255.0f, 183.0f / 255.0f, 183.0f / 255.0f, 187.0f / 255.0f, 187.0f / 255.0f, 191.0f / 255.0f, 191.0f / 255.0f,
  195.0f / 255.0f, 195.0f / 255.0f, 199.0f / 255.0f, 199.0f / 255.0f, 203.0f / 255.0f, 203.0f / 255.0f, 207.0f / 255.0f, 207.0f / 255.0f,
  211.0f / 255.0f, 211.0f / 255.0f, 215.0f / 255.0f, 215.0f / 255.0f, 219.0f / 255.0f, 219.0f / 255.0f, 223.0f / 255.0f, 223.0f / 255.0f,
  227.0f / 255.0f, 227.0f / 255.0f, 231.0f / 255.0f, 231.0f / 255.0f, 235.0f / 255.0f, 235.0f / 255.0f, 239.0f / 255.0f, 239.0f / 255.0f,
  243.0f / 255.0f, 243.0f / 255.0f, 247.0f / 255.0f, 247.0f / 255.0f, 251.0f / 255.0f, 251.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f,
};

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
const a16 float qLUT_8all[256] = {
    0.0f / 255.0f,   1.0f / 255.0f,   2.0f / 255.0f,   3.0f / 255.0f,   4.0f / 255.0f,   5.0f / 255.0f,   6.0f / 255.0f,   7.0f / 255.0f,
    8.0f / 255.0f,   9.0f / 255.0f,  10.0f / 255.0f,  11.0f / 255.0f,  12.0f / 255.0f,  13.0f / 255.0f,  14.0f / 255.0f,  15.0f / 255.0f,
   16.0f / 255.0f,  17.0f / 255.0f,  18.0f / 255.0f,  19.0f / 255.0f,  20.0f / 255.0f,  21.0f / 255.0f,  22.0f / 255.0f,  23.0f / 255.0f,
   24.0f / 255.0f,  25.0f / 255.0f,  26.0f / 255.0f,  27.0f / 255.0f,  28.0f / 255.0f,  29.0f / 255.0f,  30.0f / 255.0f,  31.0f / 255.0f,
   32.0f / 255.0f,  33.0f / 255.0f,  34.0f / 255.0f,  35.0f / 255.0f,  36.0f / 255.0f,  37.0f / 255.0f,  38.0f / 255.0f,  39.0f / 255.0f,
   40.0f / 255.0f,  41.0f / 255.0f,  42.0f / 255.0f,  43.0f / 255.0f,  44.0f / 255.0f,  45.0f / 255.0f,  46.0f / 255.0f,  47.0f / 255.0f,
   48.0f / 255.0f,  49.0f / 255.0f,  50.0f / 255.0f,  51.0f / 255.0f,  52.0f / 255.0f,  53.0f / 255.0f,  54.0f / 255.0f,  55.0f / 255.0f,
   56.0f / 255.0f,  57.0f / 255.0f,  58.0f / 255.0f,  59.0f / 255.0f,  60.0f / 255.0f,  61.0f / 255.0f,  62.0f / 255.0f,  63.0f / 255.0f,
   64.0f / 255.0f,  65.0f / 255.0f,  66.0f / 255.0f,  67.0f / 255.0f,  68.0f / 255.0f,  69.0f / 255.0f,  70.0f / 255.0f,  71.0f / 255.0f,
   72.0f / 255.0f,  73.0f / 255.0f,  74.0f / 255.0f,  75.0f / 255.0f,  76.0f / 255.0f,  77.0f / 255.0f,  78.0f / 255.0f,  79.0f / 255.0f,
   80.0f / 255.0f,  81.0f / 255.0f,  82.0f / 255.0f,  83.0f / 255.0f,  84.0f / 255.0f,  85.0f / 255.0f,  86.0f / 255.0f,  87.0f / 255.0f,
   88.0f / 255.0f,  89.0f / 255.0f,  90.0f / 255.0f,  91.0f / 255.0f,  92.0f / 255.0f,  93.0f / 255.0f,  94.0f / 255.0f,  95.0f / 255.0f,
   96.0f / 255.0f,  97.0f / 255.0f,  98.0f / 255.0f,  99.0f / 255.0f, 100.0f / 255.0f, 101.0f / 255.0f, 102.0f / 255.0f, 103.0f / 255.0f,
  104.0f / 255.0f, 105.0f / 255.0f, 106.0f / 255.0f, 107.0f / 255.0f, 108.0f / 255.0f, 109.0f / 255.0f, 110.0f / 255.0f, 111.0f / 255.0f,
  112.0f / 255.0f, 113.0f / 255.0f, 114.0f / 255.0f, 115.0f / 255.0f, 116.0f / 255.0f, 117.0f / 255.0f, 118.0f / 255.0f, 119.0f / 255.0f,
  120.0f / 255.0f, 121.0f / 255.0f, 122.0f / 255.0f, 123.0f / 255.0f, 124.0f / 255.0f, 125.0f / 255.0f, 126.0f / 255.0f, 127.0f / 255.0f,
  128.0f / 255.0f, 129.0f / 255.0f, 130.0f / 255.0f, 131.0f / 255.0f, 132.0f / 255.0f, 133.0f / 255.0f, 134.0f / 255.0f, 135.0f / 255.0f,
  136.0f / 255.0f, 137.0f / 255.0f, 138.0f / 255.0f, 139.0f / 255.0f, 140.0f / 255.0f, 141.0f / 255.0f, 142.0f / 255.0f, 143.0f / 255.0f,
  144.0f / 255.0f, 145.0f / 255.0f, 146.0f / 255.0f, 147.0f / 255.0f, 148.0f / 255.0f, 149.0f / 255.0f, 150.0f / 255.0f, 151.0f / 255.0f,
  152.0f / 255.0f, 153.0f / 255.0f, 154.0f / 255.0f, 155.0f / 255.0f, 156.0f / 255.0f, 157.0f / 255.0f, 158.0f / 255.0f, 159.0f / 255.0f,
  160.0f / 255.0f, 161.0f / 255.0f, 162.0f / 255.0f, 163.0f / 255.0f, 164.0f / 255.0f, 165.0f / 255.0f, 166.0f / 255.0f, 167.0f / 255.0f,
  168.0f / 255.0f, 169.0f / 255.0f, 170.0f / 255.0f, 171.0f / 255.0f, 172.0f / 255.0f, 173.0f / 255.0f, 174.0f / 255.0f, 175.0f / 255.0f,
  176.0f / 255.0f, 177.0f / 255.0f, 178.0f / 255.0f, 179.0f / 255.0f, 180.0f / 255.0f, 181.0f / 255.0f, 182.0f / 255.0f, 183.0f / 255.0f,
  184.0f / 255.0f, 185.0f / 255.0f, 186.0f / 255.0f, 187.0f / 255.0f, 188.0f / 255.0f, 189.0f / 255.0f, 190.0f / 255.0f, 191.0f / 255.0f,
  192.0f / 255.0f, 193.0f / 255.0f, 194.0f / 255.0f, 195.0f / 255.0f, 196.0f / 255.0f, 197.0f / 255.0f, 198.0f / 255.0f, 199.0f / 255.0f,
  200.0f / 255.0f, 201.0f / 255.0f, 202.0f / 255.0f, 203.0f / 255.0f, 204.0f / 255.0f, 205.0f / 255.0f, 206.0f / 255.0f, 207.0f / 255.0f,
  208.0f / 255.0f, 209.0f / 255.0f, 210.0f / 255.0f, 211.0f / 255.0f, 212.0f / 255.0f, 213.0f / 255.0f, 214.0f / 255.0f, 215.0f / 255.0f,
  216.0f / 255.0f, 217.0f / 255.0f, 218.0f / 255.0f, 219.0f / 255.0f, 220.0f / 255.0f, 221.0f / 255.0f, 222.0f / 255.0f, 223.0f / 255.0f,
  224.0f / 255.0f, 225.0f / 255.0f, 226.0f / 255.0f, 227.0f / 255.0f, 228.0f / 255.0f, 229.0f / 255.0f, 230.0f / 255.0f, 231.0f / 255.0f,
  232.0f / 255.0f, 233.0f / 255.0f, 234.0f / 255.0f, 235.0f / 255.0f, 236.0f / 255.0f, 237.0f / 255.0f, 238.0f / 255.0f, 239.0f / 255.0f,
  240.0f / 255.0f, 241.0f / 255.0f, 242.0f / 255.0f, 243.0f / 255.0f, 244.0f / 255.0f, 245.0f / 255.0f, 246.0f / 255.0f, 247.0f / 255.0f,
  248.0f / 255.0f, 249.0f / 255.0f, 250.0f / 255.0f, 251.0f / 255.0f, 252.0f / 255.0f, 253.0f / 255.0f, 254.0f / 255.0f, 255.0f / 255.0f,
};

const a16 float qLUT_8clr[256] = {
    0.0f / 255.0f,   0.0f / 255.0f,   2.0f / 255.0f,   2.0f / 255.0f,   4.0f / 255.0f,   4.0f / 255.0f,   6.0f / 255.0f,   6.0f / 255.0f,
    8.0f / 255.0f,   8.0f / 255.0f,  10.0f / 255.0f,  10.0f / 255.0f,  12.0f / 255.0f,  12.0f / 255.0f,  14.0f / 255.0f,  14.0f / 255.0f,
   16.0f / 255.0f,  16.0f / 255.0f,  18.0f / 255.0f,  18.0f / 255.0f,  20.0f / 255.0f,  20.0f / 255.0f,  22.0f / 255.0f,  22.0f / 255.0f,
   24.0f / 255.0f,  24.0f / 255.0f,  26.0f / 255.0f,  26.0f / 255.0f,  28.0f / 255.0f,  28.0f / 255.0f,  30.0f / 255.0f,  30.0f / 255.0f,
   32.0f / 255.0f,  32.0f / 255.0f,  34.0f / 255.0f,  34.0f / 255.0f,  36.0f / 255.0f,  36.0f / 255.0f,  38.0f / 255.0f,  38.0f / 255.0f,
   40.0f / 255.0f,  40.0f / 255.0f,  42.0f / 255.0f,  42.0f / 255.0f,  44.0f / 255.0f,  44.0f / 255.0f,  46.0f / 255.0f,  46.0f / 255.0f,
   48.0f / 255.0f,  48.0f / 255.0f,  50.0f / 255.0f,  50.0f / 255.0f,  52.0f / 255.0f,  52.0f / 255.0f,  54.0f / 255.0f,  54.0f / 255.0f,
   56.0f / 255.0f,  56.0f / 255.0f,  58.0f / 255.0f,  58.0f / 255.0f,  60.0f / 255.0f,  60.0f / 255.0f,  62.0f / 255.0f,  62.0f / 255.0f,
   64.0f / 255.0f,  64.0f / 255.0f,  66.0f / 255.0f,  66.0f / 255.0f,  68.0f / 255.0f,  68.0f / 255.0f,  70.0f / 255.0f,  70.0f / 255.0f,
   72.0f / 255.0f,  72.0f / 255.0f,  74.0f / 255.0f,  74.0f / 255.0f,  76.0f / 255.0f,  76.0f / 255.0f,  78.0f / 255.0f,  78.0f / 255.0f,
   80.0f / 255.0f,  80.0f / 255.0f,  82.0f / 255.0f,  82.0f / 255.0f,  84.0f / 255.0f,  84.0f / 255.0f,  86.0f / 255.0f,  86.0f / 255.0f,
   88.0f / 255.0f,  88.0f / 255.0f,  90.0f / 255.0f,  90.0f / 255.0f,  92.0f / 255.0f,  92.0f / 255.0f,  94.0f / 255.0f,  94.0f / 255.0f,
   96.0f / 255.0f,  96.0f / 255.0f,  98.0f / 255.0f,  98.0f / 255.0f, 100.0f / 255.0f, 100.0f / 255.0f, 102.0f / 255.0f, 102.0f / 255.0f,
  104.0f / 255.0f, 104.0f / 255.0f, 106.0f / 255.0f, 106.0f / 255.0f, 108.0f / 255.0f, 108.0f / 255.0f, 110.0f / 255.0f, 110.0f / 255.0f,
  112.0f / 255.0f, 112.0f / 255.0f, 114.0f / 255.0f, 114.0f / 255.0f, 116.0f / 255.0f, 116.0f / 255.0f, 118.0f / 255.0f, 118.0f / 255.0f,
  120.0f / 255.0f, 120.0f / 255.0f, 122.0f / 255.0f, 122.0f / 255.0f, 124.0f / 255.0f, 124.0f / 255.0f, 126.0f / 255.0f, 126.0f / 255.0f,
  128.0f / 255.0f, 128.0f / 255.0f, 130.0f / 255.0f, 130.0f / 255.0f, 132.0f / 255.0f, 132.0f / 255.0f, 134.0f / 255.0f, 134.0f / 255.0f,
  136.0f / 255.0f, 136.0f / 255.0f, 138.0f / 255.0f, 138.0f / 255.0f, 140.0f / 255.0f, 140.0f / 255.0f, 142.0f / 255.0f, 142.0f / 255.0f,
  144.0f / 255.0f, 144.0f / 255.0f, 146.0f / 255.0f, 146.0f / 255.0f, 148.0f / 255.0f, 148.0f / 255.0f, 150.0f / 255.0f, 150.0f / 255.0f,
  152.0f / 255.0f, 152.0f / 255.0f, 154.0f / 255.0f, 154.0f / 255.0f, 156.0f / 255.0f, 156.0f / 255.0f, 158.0f / 255.0f, 158.0f / 255.0f,
  160.0f / 255.0f, 160.0f / 255.0f, 162.0f / 255.0f, 162.0f / 255.0f, 164.0f / 255.0f, 164.0f / 255.0f, 166.0f / 255.0f, 166.0f / 255.0f,
  168.0f / 255.0f, 168.0f / 255.0f, 170.0f / 255.0f, 170.0f / 255.0f, 172.0f / 255.0f, 172.0f / 255.0f, 174.0f / 255.0f, 174.0f / 255.0f,
  176.0f / 255.0f, 176.0f / 255.0f, 178.0f / 255.0f, 178.0f / 255.0f, 180.0f / 255.0f, 180.0f / 255.0f, 182.0f / 255.0f, 182.0f / 255.0f,
  184.0f / 255.0f, 184.0f / 255.0f, 186.0f / 255.0f, 186.0f / 255.0f, 188.0f / 255.0f, 188.0f / 255.0f, 190.0f / 255.0f, 190.0f / 255.0f,
  192.0f / 255.0f, 192.0f / 255.0f, 194.0f / 255.0f, 194.0f / 255.0f, 196.0f / 255.0f, 196.0f / 255.0f, 198.0f / 255.0f, 198.0f / 255.0f,
  200.0f / 255.0f, 200.0f / 255.0f, 202.0f / 255.0f, 202.0f / 255.0f, 204.0f / 255.0f, 204.0f / 255.0f, 206.0f / 255.0f, 206.0f / 255.0f,
  208.0f / 255.0f, 208.0f / 255.0f, 210.0f / 255.0f, 210.0f / 255.0f, 212.0f / 255.0f, 212.0f / 255.0f, 214.0f / 255.0f, 214.0f / 255.0f,
  216.0f / 255.0f, 216.0f / 255.0f, 218.0f / 255.0f, 218.0f / 255.0f, 220.0f / 255.0f, 220.0f / 255.0f, 222.0f / 255.0f, 222.0f / 255.0f,
  224.0f / 255.0f, 224.0f / 255.0f, 226.0f / 255.0f, 226.0f / 255.0f, 228.0f / 255.0f, 228.0f / 255.0f, 230.0f / 255.0f, 230.0f / 255.0f,
  232.0f / 255.0f, 232.0f / 255.0f, 234.0f / 255.0f, 234.0f / 255.0f, 236.0f / 255.0f, 236.0f / 255.0f, 238.0f / 255.0f, 238.0f / 255.0f,
  240.0f / 255.0f, 240.0f / 255.0f, 242.0f / 255.0f, 242.0f / 255.0f, 244.0f / 255.0f, 244.0f / 255.0f, 246.0f / 255.0f, 246.0f / 255.0f,
  248.0f / 255.0f, 248.0f / 255.0f, 250.0f / 255.0f, 250.0f / 255.0f, 252.0f / 255.0f, 252.0f / 255.0f, 254.0f / 255.0f, 254.0f / 255.0f,
};

const a16 float qLUT_8set[256] = {
    1.0f / 255.0f,   1.0f / 255.0f,   3.0f / 255.0f,   3.0f / 255.0f,   5.0f / 255.0f,   5.0f / 255.0f,   7.0f / 255.0f,   7.0f / 255.0f,
    9.0f / 255.0f,   9.0f / 255.0f,  11.0f / 255.0f,  11.0f / 255.0f,  13.0f / 255.0f,  13.0f / 255.0f,  15.0f / 255.0f,  15.0f / 255.0f,
   17.0f / 255.0f,  17.0f / 255.0f,  19.0f / 255.0f,  19.0f / 255.0f,  21.0f / 255.0f,  21.0f / 255.0f,  23.0f / 255.0f,  23.0f / 255.0f,
   25.0f / 255.0f,  25.0f / 255.0f,  27.0f / 255.0f,  27.0f / 255.0f,  29.0f / 255.0f,  29.0f / 255.0f,  31.0f / 255.0f,  31.0f / 255.0f,
   33.0f / 255.0f,  33.0f / 255.0f,  35.0f / 255.0f,  35.0f / 255.0f,  37.0f / 255.0f,  37.0f / 255.0f,  39.0f / 255.0f,  39.0f / 255.0f,
   41.0f / 255.0f,  41.0f / 255.0f,  43.0f / 255.0f,  43.0f / 255.0f,  45.0f / 255.0f,  45.0f / 255.0f,  47.0f / 255.0f,  47.0f / 255.0f,
   49.0f / 255.0f,  49.0f / 255.0f,  51.0f / 255.0f,  51.0f / 255.0f,  53.0f / 255.0f,  53.0f / 255.0f,  55.0f / 255.0f,  55.0f / 255.0f,
   57.0f / 255.0f,  57.0f / 255.0f,  59.0f / 255.0f,  59.0f / 255.0f,  61.0f / 255.0f,  61.0f / 255.0f,  63.0f / 255.0f,  63.0f / 255.0f,
   65.0f / 255.0f,  65.0f / 255.0f,  67.0f / 255.0f,  67.0f / 255.0f,  69.0f / 255.0f,  69.0f / 255.0f,  71.0f / 255.0f,  71.0f / 255.0f,
   73.0f / 255.0f,  73.0f / 255.0f,  75.0f / 255.0f,  75.0f / 255.0f,  77.0f / 255.0f,  77.0f / 255.0f,  79.0f / 255.0f,  79.0f / 255.0f,
   81.0f / 255.0f,  81.0f / 255.0f,  83.0f / 255.0f,  83.0f / 255.0f,  85.0f / 255.0f,  85.0f / 255.0f,  87.0f / 255.0f,  87.0f / 255.0f,
   89.0f / 255.0f,  89.0f / 255.0f,  91.0f / 255.0f,  91.0f / 255.0f,  93.0f / 255.0f,  93.0f / 255.0f,  95.0f / 255.0f,  95.0f / 255.0f,
   97.0f / 255.0f,  97.0f / 255.0f,  99.0f / 255.0f,  99.0f / 255.0f, 101.0f / 255.0f, 101.0f / 255.0f, 103.0f / 255.0f, 103.0f / 255.0f,
  105.0f / 255.0f, 105.0f / 255.0f, 107.0f / 255.0f, 107.0f / 255.0f, 109.0f / 255.0f, 109.0f / 255.0f, 111.0f / 255.0f, 111.0f / 255.0f,
  113.0f / 255.0f, 113.0f / 255.0f, 115.0f / 255.0f, 115.0f / 255.0f, 117.0f / 255.0f, 117.0f / 255.0f, 119.0f / 255.0f, 119.0f / 255.0f,
  121.0f / 255.0f, 121.0f / 255.0f, 123.0f / 255.0f, 123.0f / 255.0f, 125.0f / 255.0f, 125.0f / 255.0f, 127.0f / 255.0f, 127.0f / 255.0f,
  129.0f / 255.0f, 129.0f / 255.0f, 131.0f / 255.0f, 131.0f / 255.0f, 133.0f / 255.0f, 133.0f / 255.0f, 135.0f / 255.0f, 135.0f / 255.0f,
  137.0f / 255.0f, 137.0f / 255.0f, 139.0f / 255.0f, 139.0f / 255.0f, 141.0f / 255.0f, 141.0f / 255.0f, 143.0f / 255.0f, 143.0f / 255.0f,
  145.0f / 255.0f, 145.0f / 255.0f, 147.0f / 255.0f, 147.0f / 255.0f, 149.0f / 255.0f, 149.0f / 255.0f, 151.0f / 255.0f, 151.0f / 255.0f,
  153.0f / 255.0f, 153.0f / 255.0f, 155.0f / 255.0f, 155.0f / 255.0f, 157.0f / 255.0f, 157.0f / 255.0f, 159.0f / 255.0f, 159.0f / 255.0f,
  161.0f / 255.0f, 161.0f / 255.0f, 163.0f / 255.0f, 163.0f / 255.0f, 165.0f / 255.0f, 165.0f / 255.0f, 167.0f / 255.0f, 167.0f / 255.0f,
  169.0f / 255.0f, 169.0f / 255.0f, 171.0f / 255.0f, 171.0f / 255.0f, 173.0f / 255.0f, 173.0f / 255.0f, 175.0f / 255.0f, 175.0f / 255.0f,
  177.0f / 255.0f, 177.0f / 255.0f, 179.0f / 255.0f, 179.0f / 255.0f, 181.0f / 255.0f, 181.0f / 255.0f, 183.0f / 255.0f, 183.0f / 255.0f,
  185.0f / 255.0f, 185.0f / 255.0f, 187.0f / 255.0f, 187.0f / 255.0f, 189.0f / 255.0f, 189.0f / 255.0f, 191.0f / 255.0f, 191.0f / 255.0f,
  193.0f / 255.0f, 193.0f / 255.0f, 195.0f / 255.0f, 195.0f / 255.0f, 197.0f / 255.0f, 197.0f / 255.0f, 199.0f / 255.0f, 199.0f / 255.0f,
  201.0f / 255.0f, 201.0f / 255.0f, 203.0f / 255.0f, 203.0f / 255.0f, 205.0f / 255.0f, 205.0f / 255.0f, 207.0f / 255.0f, 207.0f / 255.0f,
  209.0f / 255.0f, 209.0f / 255.0f, 211.0f / 255.0f, 211.0f / 255.0f, 213.0f / 255.0f, 213.0f / 255.0f, 215.0f / 255.0f, 215.0f / 255.0f,
  217.0f / 255.0f, 217.0f / 255.0f, 219.0f / 255.0f, 219.0f / 255.0f, 221.0f / 255.0f, 221.0f / 255.0f, 223.0f / 255.0f, 223.0f / 255.0f,
  225.0f / 255.0f, 225.0f / 255.0f, 227.0f / 255.0f, 227.0f / 255.0f, 229.0f / 255.0f, 229.0f / 255.0f, 231.0f / 255.0f, 231.0f / 255.0f,
  233.0f / 255.0f, 233.0f / 255.0f, 235.0f / 255.0f, 235.0f / 255.0f, 237.0f / 255.0f, 237.0f / 255.0f, 239.0f / 255.0f, 239.0f / 255.0f,
  241.0f / 255.0f, 241.0f / 255.0f, 243.0f / 255.0f, 243.0f / 255.0f, 245.0f / 255.0f, 245.0f / 255.0f, 247.0f / 255.0f, 247.0f / 255.0f,
  249.0f / 255.0f, 249.0f / 255.0f, 251.0f / 255.0f, 251.0f / 255.0f, 253.0f / 255.0f, 253.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f,
};

//cQuantizer3<5,6,5  > q3_565;
//cQuantizer4<5,6,5,0> q4_565;

/* -----------------------------------------------------------------------------
 */

// sRGB spec
float basefpartition = 0.0031308f;
float baseipartition = 0.004045f;
float basefslope = 12.92f / 1.0f;
float baseislope = 1.0f / 12.92f;
float basefgamma = 2.4f / 1.0f;
float baseigamma = 1.0f / 2.4f;
float baseoffset = 0.055f;
float baseLUT_sRGB[256] = {
   0.0f,0.000303527f,0.00114819f,0.00132772f,0.00152264f,0.00173331f,0.00196007f,0.00220325f,
   0.00246318f,0.00274017f,0.00303452f,0.00334654f,0.00367651f,0.00402472f,0.00439144f,0.00477695f,
   0.00518152f,0.00560539f,0.00604883f,0.00651209f,0.00699541f,0.00749903f,0.00802319f,0.00856813f,
   0.00913406f,0.00972122f,0.0103298f,0.0109601f,0.0116122f,0.0122865f,0.012983f,0.0137021f,
   0.0144438f,0.0152085f,0.0159963f,0.0168074f,0.017642f,0.0185002f,0.0193824f,0.0202886f,
   0.021219f,0.0221739f,0.0231534f,0.0241576f,0.0251869f,0.0262412f,0.0273209f,0.028426f,
   0.0295568f,0.0307134f,0.031896f,0.0331048f,0.0343398f,0.0356013f,0.0368894f,0.0382044f,
   0.0395462f,0.0409152f,0.0423114f,0.043735f,0.0451862f,0.0466651f,0.0481718f,0.0497066f,
   0.0512695f,0.0528606f,0.0544803f,0.0561285f,0.0578054f,0.0595112f,0.0612461f,0.06301f,
   0.0648033f,0.0666259f,0.0684782f,0.0703601f,0.0722719f,0.0742136f,0.0761854f,0.0781874f,
   0.0802198f,0.0822827f,0.0843762f,0.0865005f,0.0886556f,0.0908417f,0.093059f,0.0953075f,
   0.0975873f,0.0998987f,0.102242f,0.104616f,0.107023f,0.109462f,0.111932f,0.114435f,
   0.116971f,0.119538f,0.122139f,0.124772f,0.127438f,0.130136f,0.132868f,0.135633f,
   0.138432f,0.141263f,0.144128f,0.147027f,0.14996f,0.152926f,0.155926f,0.158961f,
   0.162029f,0.165132f,0.168269f,0.171441f,0.174647f,0.177888f,0.181164f,0.184475f,
   0.187821f,0.191202f,0.194618f,0.198069f,0.201556f,0.205079f,0.208637f,0.212231f,
   0.215861f,0.219526f,0.223228f,0.226966f,0.23074f,0.234551f,0.238398f,0.242281f,
   0.246201f,0.250158f,0.254152f,0.258183f,0.262251f,0.266356f,0.270498f,0.274677f,
   0.278894f,0.283149f,0.287441f,0.291771f,0.296138f,0.300544f,0.304987f,0.309469f,
   0.313989f,0.318547f,0.323143f,0.327778f,0.332452f,0.337164f,0.341914f,0.346704f,
   0.351533f,0.3564f,0.361307f,0.366253f,0.371238f,0.376262f,0.381326f,0.386429f,
   0.391572f,0.396755f,0.401978f,0.40724f,0.412543f,0.417885f,0.423268f,0.42869f,
   0.434154f,0.439657f,0.445201f,0.450786f,0.456411f,0.462077f,0.467784f,0.473531f,
   0.47932f,0.48515f,0.491021f,0.496933f,0.502886f,0.508881f,0.514918f,0.520996f,
   0.527115f,0.533276f,0.539479f,0.545724f,0.552011f,0.55834f,0.564712f,0.571125f,
   0.57758f,0.584078f,0.590619f,0.597202f,0.603827f,0.610496f,0.617207f,0.62396f,
   0.630757f,0.637597f,0.64448f,0.651406f,0.658375f,0.665387f,0.672443f,0.679542f,
   0.686685f,0.693872f,0.701102f,0.708376f,0.715693f,0.723055f,0.730461f,0.73791f,
   0.745404f,0.752942f,0.760525f,0.768151f,0.775822f,0.783538f,0.791298f,0.799103f,
   0.806952f,0.814847f,0.822786f,0.83077f,0.838799f,0.846873f,0.854993f,0.863157f,
   0.871367f,0.879622f,0.887923f,0.896269f,0.904661f,0.913099f,0.921582f,0.930111f,
   0.938686f,0.947307f,0.955973f,0.964686f,0.973445f,0.982251f,0.991102f,1.0f};
float baseLUT_Linear[256] = {
   0.0f,0.00392157f,0.00784314f,0.0117647f,0.0156863f,0.0196078f,0.0235294f,0.027451f,
   0.0313726f,0.0352941f,0.0392157f,0.0431373f,0.0470588f,0.0509804f,0.054902f,0.0588235f,
   0.0627451f,0.0666667f,0.0705882f,0.0745098f,0.0784314f,0.0823529f,0.0862745f,0.0901961f,
   0.0941176f,0.0980392f,0.101961f,0.105882f,0.109804f,0.113725f,0.117647f,0.121569f,
   0.12549f,0.129412f,0.133333f,0.137255f,0.141176f,0.145098f,0.14902f,0.152941f,
   0.156863f,0.160784f,0.164706f,0.168627f,0.172549f,0.176471f,0.180392f,0.184314f,
   0.188235f,0.192157f,0.196078f,0.2f,0.203922f,0.207843f,0.211765f,0.215686f,
   0.219608f,0.223529f,0.227451f,0.231373f,0.235294f,0.239216f,0.243137f,0.247059f,
   0.25098f,0.254902f,0.258824f,0.262745f,0.266667f,0.270588f,0.27451f,0.278431f,
   0.282353f,0.286275f,0.290196f,0.294118f,0.298039f,0.301961f,0.305882f,0.309804f,
   0.313726f,0.317647f,0.321569f,0.32549f,0.329412f,0.333333f,0.337255f,0.341176f,
   0.345098f,0.34902f,0.352941f,0.356863f,0.360784f,0.364706f,0.368627f,0.372549f,
   0.376471f,0.380392f,0.384314f,0.388235f,0.392157f,0.396078f,0.4f,0.403922f,
   0.407843f,0.411765f,0.415686f,0.419608f,0.423529f,0.427451f,0.431373f,0.435294f,
   0.439216f,0.443137f,0.447059f,0.45098f,0.454902f,0.458824f,0.462745f,0.466667f,
   0.470588f,0.47451f,0.478431f,0.482353f,0.486275f,0.490196f,0.494118f,0.498039f,
   0.501961f,0.505882f,0.509804f,0.513726f,0.517647f,0.521569f,0.52549f,0.529412f,
   0.533333f,0.537255f,0.541176f,0.545098f,0.54902f,0.552941f,0.556863f,0.560784f,
   0.564706f,0.568627f,0.572549f,0.576471f,0.580392f,0.584314f,0.588235f,0.592157f,
   0.596078f,0.6f,0.603922f,0.607843f,0.611765f,0.615686f,0.619608f,0.623529f,
   0.627451f,0.631373f,0.635294f,0.639216f,0.643137f,0.647059f,0.65098f,0.654902f,
   0.658824f,0.662745f,0.666667f,0.670588f,0.67451f,0.678431f,0.682353f,0.686275f,
   0.690196f,0.694118f,0.698039f,0.701961f,0.705882f,0.709804f,0.713726f,0.717647f,
   0.721569f,0.72549f,0.729412f,0.733333f,0.737255f,0.741176f,0.745098f,0.74902f,
   0.752941f,0.756863f,0.760784f,0.764706f,0.768627f,0.772549f,0.776471f,0.780392f,
   0.784314f,0.788235f,0.792157f,0.796078f,0.8f,0.803922f,0.807843f,0.811765f,
   0.815686f,0.819608f,0.823529f,0.827451f,0.831373f,0.835294f,0.839216f,0.843137f,
   0.847059f,0.85098f,0.854902f,0.858824f,0.862745f,0.866667f,0.870588f,0.87451f,
   0.878431f,0.882353f,0.886275f,0.890196f,0.894118f,0.898039f,0.901961f,0.905882f,
   0.909804f,0.913725f,0.917647f,0.921569f,0.92549f,0.929412f,0.933333f,0.937255f,
   0.941176f,0.945098f,0.94902f,0.952941f,0.956863f,0.960784f,0.964706f,0.968627f,
   0.972549f,0.976471f,0.980392f,0.984314f,0.988235f,0.992157f,0.996078f,1.0f};

void SetGamma() {
  ;
}

const float *ComputeGammaLUT(bool sRGB) {
  return (sRGB ? baseLUT_sRGB : baseLUT_Linear);
}
#endif

} // namespace squish

#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#if	!defined(SQUISH_USE_COMPUTE)
// this one works on a vector<6>
#include "maths_vector.cpp"
#else
// this one works on two float4s
#include "maths_packed.cpp"
#endif
#endif
