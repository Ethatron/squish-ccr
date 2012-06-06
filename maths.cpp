/* -----------------------------------------------------------------------------

	Copyright (c) 2006 Simon Brown                          si@sjbrown.co.uk
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
*/

#include "maths.h"

namespace squish {

#if	!defined(USE_PRE)
Sym3x3 ComputeWeightedCovariance( int n, Vec3 const* points, float const* weights )
{
	// compute the centroid
	float total = 0.0f;
	Vec3 centroid( 0.0f );
	for( int i = 0; i < n; ++i )
	{
		total += weights[i];
		centroid += weights[i]*points[i];
	}
	centroid /= total;

	// accumulate the covariance smatrix
	Sym3x3 covariance( 0.0f );
	for( int i = 0; i < n; ++i )
	{
		Vec3 a = points[i] - centroid;
		Vec3 b = weights[i]*a;

		covariance[0] += a.X()*b.X();
		covariance[1] += a.X()*b.Y();
		covariance[2] += a.X()*b.Z();
		covariance[3] += a.Y()*b.Y();
		covariance[4] += a.Y()*b.Z();
		covariance[5] += a.Z()*b.Z();
	}

	// return it
	return covariance;
}

static Vec3 GetMultiplicity1Evector( Sym3x3 const& smatrix, float evalue )
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
      return Vec3( u[0], u[1], u[2] );

    case 1:
    case 3:
      return Vec3( u[1], u[3], u[4] );

    default:
      return Vec3( u[2], u[4], u[5] );
  }
}

static Vec3 GetMultiplicity2Evector( Sym3x3 const& smatrix, float evalue )
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
      return Vec3( -m[1], m[0], 0.0f );

    case 2:
      return Vec3( m[2], 0.0f, -m[0] );

    case 3:
    case 4:
      return Vec3( 0.0f, -m[4], m[3] );

    default:
      return Vec3( 0.0f, -m[5], m[4] );
  }
}

Vec3 ComputePrincipleComponent( Sym3x3 const& smatrix )
{
  // compute the cubic coefficients
  float c0 = smatrix[0]*smatrix[3]*smatrix[5]
  + 2.0f*smatrix[1]*smatrix[2]*smatrix[4]
  - smatrix[0]*smatrix[4]*smatrix[4]
  - smatrix[3]*smatrix[2]*smatrix[2]
  - smatrix[5]*smatrix[1]*smatrix[1];
  float c1 = smatrix[0]*smatrix[3] + smatrix[0]*smatrix[5] + smatrix[3]*smatrix[5]
  - smatrix[1]*smatrix[1] - smatrix[2]*smatrix[2] - smatrix[4]*smatrix[4];
  float c2 = smatrix[0] + smatrix[3] + smatrix[5];

  // compute the quadratic coefficients
  float a = c1 - ( 1.0f/3.0f )*c2*c2;
  float b = ( -2.0f/27.0f )*c2*c2*c2 + ( 1.0f/3.0f )*c1*c2 - c0;

  // compute the root count check
  float Q = 0.25f*b*b + ( 1.0f/27.0f )*a*a*a;

  // test the multiplicity
  if( FLT_EPSILON < Q )
  {
    // only one root, which implies we have a multiple of the identity
    return Vec3( 1.0f );
  }
  else if( Q < -FLT_EPSILON )
  {
    // three distinct roots
    float theta = std::atan2( std::sqrt( -Q ), -0.5f*b );
    float rho = std::sqrt( 0.25f*b*b - Q );

    float rt = std::pow( rho, 1.0f/3.0f );
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
    return GetMultiplicity1Evector( smatrix, l1 );
  }
  else // if( -FLT_EPSILON <= Q && Q <= FLT_EPSILON )
  {
    // two roots
    float rt;
    if( b < 0.0f )
      rt = -std::pow( -0.5f*b, 1.0f/3.0f );
    else
      rt = std::pow( 0.5f*b, 1.0f/3.0f );

    float l1 = ( 1.0f/3.0f )*c2 + rt;		// repeated
    float l2 = ( 1.0f/3.0f )*c2 - 2.0f*rt;

    // get the eigenvector
    if( std::fabs( l1 ) > std::fabs( l2 ) )
      return GetMultiplicity2Evector( smatrix, l1 );
    else
      return GetMultiplicity1Evector( smatrix, l2 );
  }
}
#endif

#if	defined(USE_AMP) || defined(USE_COMPUTE)
#if	defined(USE_COMPUTE)
  tile_static float  _totals[8];
  tile_static float3 _centroids[8];
  tile_static Sym3x3 covariance[8];
#endif

Sym3x3 ComputeWeightedCovariance(tile_barrier barrier, const int thread, int n, point16 points, weight16 weights) amp_restricted
{
  // compute the centroid (AMP: reduction, O(count) vs. O(ln(16)) = O(4),
  //                            unused coeffs have been set to zero and
  //                            remain neutral in the covariance)
#if	!defined(USE_COMPUTE)
  tile_static float  _totals[8];
  tile_static float3 _centroids[8];
#endif

  threaded_for(tc1, 8) {
    const int dn = (tc1 << 1) + 0;
    const int up = (tc1 << 1) + 1;

    _totals   [tc1] = weights[dn]              + weights[up];
    _centroids[tc1] = weights[dn] * points[dn] + weights[up] * points[up];
  }

  threaded_for(tc2, 4) {
    const int dn = (tc2 << 1) + 0;
    const int up = (tc2 << 1) + 1;

    _totals   [tc2] = _totals   [dn] + _totals   [up];
    _centroids[tc2] = _centroids[dn] + _centroids[up];
  }

  threaded_for(tc3, 2) {
    const int dn = (tc3 << 1) + 0;
    const int up = (tc3 << 1) + 1;

    _totals   [tc3] = _totals   [dn] + _totals   [up];
    _centroids[tc3] = _centroids[dn] + _centroids[up];
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // compute the centroid (AMP: pull into register)
  float3 _centroid;

  _centroid = (_centroids[0] + _centroids[1]) /
              (_totals   [0] + _totals   [1]);

  // accumulate the covariance smatrix
#if	!defined(USE_COMPUTE)
  tile_static Sym3x3 covariance[8];
#endif

  threaded_for(co1, 8) {
    const int dn = (co1 << 1) + 0;
    const int up = (co1 << 1) + 1;

    float3 a =  points[dn] - _centroid;
    float3 A =  points[up] - _centroid;
    float3 b = weights[dn] * a;
    float3 B = weights[up] * A;

    covariance[co1][0] = a.x * b.x + A.x * B.x;
    covariance[co1][1] = a.x * b.y + A.x * B.y;
    covariance[co1][2] = a.x * b.z + A.x * B.z;
    covariance[co1][3] = a.y * b.y + A.y * B.y;
    covariance[co1][4] = a.y * b.z + A.y * B.z;
    covariance[co1][5] = a.z * b.z + A.z * B.z;
  }

  threaded_for(co2, 2*4) {
    const int c2 = co2 & 3;
    const int dn = ((co2 & 3) << 1) + 0;
    const int up = ((co2 & 3) << 1) + 1;
    const int of = (co2 >> 2) * 3;

    // AMP: prefer 3-wide vectorized op
    covariance[c2][0 + of] = covariance[dn][0 + of] + covariance[up][0 + of];
    covariance[c2][1 + of] = covariance[dn][1 + of] + covariance[up][1 + of];
    covariance[c2][2 + of] = covariance[dn][2 + of] + covariance[up][2 + of];
  }

  threaded_for(co3, 2*2) {
    const int c3 = co2 & 1;
    const int dn = ((co3 & 1) << 1) + 0;
    const int up = ((co3 & 1) << 1) + 1;
    const int of = (co3 >> 1) * 3;

    // AMP: prefer 3-wide vectorized op
    covariance[c3][0 + of] = covariance[dn][0 + of] + covariance[up][0 + of];
    covariance[c3][1 + of] = covariance[dn][1 + of] + covariance[up][1 + of];
    covariance[c3][2 + of] = covariance[dn][2 + of] + covariance[up][2 + of];
  }

  threaded_for(co4, 2) {
    const int c4 = 0;
    const int dn = 0;
    const int up = 1;
    const int of = (co4 >> 0) * 3;

    // AMP: prefer 3-wide vectorized op
    covariance[c4][0 + of] = covariance[dn][0 + of] + covariance[up][0 + of];
    covariance[c4][1 + of] = covariance[dn][1 + of] + covariance[up][1 + of];
    covariance[c4][2 + of] = covariance[dn][2 + of] + covariance[up][2 + of];
  }

  return covariance[0];

#if 0
  // compute the centroid
  float  _total = 0.0f;
  float3 _centroid = 0.0f;

  for (int i = 0; i < n; ++i) {
    _total    += weights[i];
    _centroid += weights[i] * points[i];
  }

  _centroid /= _total;

  // accumulate the covariance smatrix
  Sym3x3 covariance;

  covariance[0] = 0;
  covariance[1] = 0;
  covariance[2] = 0;
  covariance[3] = 0;
  covariance[4] = 0;
  covariance[5] = 0;

  for (int i = 0; i < n; ++i) {
    float3 a =  points[i] - _centroid;
    float3 b = weights[i] * a;

    covariance[0] += a.x * b.x;
    covariance[1] += a.x * b.y;
    covariance[2] += a.x * b.z;
    covariance[3] += a.y * b.y;
    covariance[4] += a.y * b.z;
    covariance[5] += a.z * b.z;
  }

  // return it
  return covariance;
#endif
}

static float3 GetMultiplicity1Evector(tile_barrier barrier, const int thread, Sym3x3r smatrix, float evalue) amp_restricted
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

  u[0] = m[3] * m[5] - m[4] * m[4];
  u[1] = m[2] * m[4] - m[1] * m[5];
  u[2] = m[1] * m[4] - m[2] * m[3];
  u[3] = m[0] * m[5] - m[2] * m[2];
  u[4] = m[1] * m[2] - m[4] * m[0];
  u[5] = m[0] * m[3] - m[1] * m[1];

  // find the largest component
  float mc = fabsf(u[0]);
  int mi = 0;

  for (int i = 1; i < 6; ++i) {
    float c = fabsf(u[i]);
    if (c > mc) {
      mc = c;
      mi = i;
    }
  }

  // pick the column with this component
  switch (mi) {
    case 0:
      return float3(u[0], u[1], u[2]);

    case 1:
    case 3:
      return float3(u[1], u[3], u[4]);

    default:
      return float3(u[2], u[4], u[5]);
  }
}

static float3 GetMultiplicity2Evector(tile_barrier barrier, const int thread, Sym3x3r smatrix, float evalue) amp_restricted
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
  float mc = fabsf(m[0]);
  int mi = 0;

  for (int i = 1; i < 6; ++i) {
    float c = fabsf(m[i]);
    if (c > mc) {
      mc = c;
      mi = i;
    }
  }

  // pick the first eigenvector based on this index
  switch (mi) {
    case 0:
    case 1:
      return float3(-m[1], m[0],  0.0f);

    case 2:
      return float3( m[2], 0.0f, -m[0]);

    case 3:
    case 4:
      return float3( 0.0f, -m[4], m[3]);

    default:
      return float3( 0.0f, -m[5], m[4]);
  }
}

float3 ComputePrincipleComponent(tile_barrier barrier, const int thread, Sym3x3r smatrix) amp_restricted
{
  // compute the cubic coefficients
  float c0 = smatrix[0] * smatrix[3] * smatrix[5]
    + 2.0f * smatrix[1] * smatrix[2] * smatrix[4]
    -        smatrix[0] * smatrix[4] * smatrix[4]
    -        smatrix[3] * smatrix[2] * smatrix[2]
    -        smatrix[5] * smatrix[1] * smatrix[1];
  float c1 = smatrix[0] * smatrix[3]
    +        smatrix[0] * smatrix[5]
    +        smatrix[3] * smatrix[5]
    -        smatrix[1] * smatrix[1]
    -        smatrix[2] * smatrix[2]
    -        smatrix[4] * smatrix[4];
  float c2 = smatrix[0]
    +        smatrix[3]
    +        smatrix[5];

  // compute the quadratic coefficients
  float a = c1 - ( 1.0f /  3.0f) * c2 * c2;
  float b =      (-2.0f / 27.0f) * c2 * c2 * c2
               + ( 1.0f /  3.0f) * c1 * c2 - c0;

  // compute the root count check
  float Q = 0.25f * b * b + (1.0f / 27.0f) * a * a * a;

  // test the multiplicity
  if (FLT_EPSILON < Q) {
    // only one root, which implies we have a multiple of the identity
    return 1.0f;
  }
  else if (Q < -FLT_EPSILON) {
    // three distinct roots
    float theta = atan2f( sqrtf( -Q ), -0.5f*b );
    float rho = sqrtf( 0.25f*b*b - Q );

    float rt = powf(rho, 1.0f / 3.0f);
    float ct = cosf(theta / 3.0f);
    float st = sinf(theta / 3.0f);

    float l1 = (1.0f / 3.0f) * c2 + 2.0f * rt * ct;
    float l2 = (1.0f / 3.0f) * c2 - rt * (ct + (float)sqrtf(3.0f) * st);
    float l3 = (1.0f / 3.0f) * c2 - rt * (ct - (float)sqrtf(3.0f) * st);

    // pick the larger
    if (fabsf(l2) > fabsf(l1))
      l1 = l2;
    if (fabsf(l3) > fabsf(l1))
      l1 = l3;

    // get the eigenvector
    return GetMultiplicity1Evector(barrier, thread, smatrix, l1);
  }
  else { // if( -FLT_EPSILON <= Q && Q <= FLT_EPSILON )
    // two roots
    float rt;
    if (b < 0.0f)
      rt = -powf(-0.5f * b, 1.0f / 3.0f);
    else
      rt =  powf( 0.5f * b, 1.0f / 3.0f);

    float l1 = (1.0f / 3.0f) * c2 +        rt;		// repeated
    float l2 = (1.0f / 3.0f) * c2 - 2.0f * rt;

    // get the eigenvector
    if (fabsf(l1) > fabsf(l2))
      return GetMultiplicity2Evector(barrier, thread, smatrix, l1);
    else
      return GetMultiplicity1Evector(barrier, thread, smatrix, l2);
  }
}
#endif

} // namespace squish
