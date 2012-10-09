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

namespace squish {

#if	defined(SQUISH_USE_COMPUTE)
  tile_static float4 centroid_t[8];
  tile_static Sym3x3 covariance[8];
#endif

Sym3x3 ComputeWeightedCovariance3(tile_barrier barrier, const int thread, int n, point16 points, weight16 weights) amp_restricted
{
  // compute the centroid (AMP: reduction, O(count) vs. O(ln(16)) = O(4),
  //                            unused coeffs have been set to zero and
  //                            remain neutral in the covariance)
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static float4 centroid_t[8];
#endif

  // same for all
  const int dn = (thread << 1) + 0;
  const int up = (thread << 1) + 1;

  threaded_for(tc1, 8) {
    centroid_t[tc1].w   = weights[dn]              + weights[up];
    centroid_t[tc1].xyz = weights[dn] * points[dn] + weights[up] * points[up]; }

  // AMP: prefer 4-wide vectorized op
  threaded_for(tc2, 4) {
    centroid_t[tc2]     = centroid_t[dn] + centroid_t[up]; }
  threaded_for(tc3, 2) {
    centroid_t[tc3]     = centroid_t[dn] + centroid_t[up]; }
  threaded_for(tc4, 1) {
    centroid_t[tc4]     = centroid_t[dn] + centroid_t[up]; }

  // make it visible
  tile_static_memory_fence(barrier);

  // compute the centroid (AMP: pull into register)
  float3 _centroid = centroid_t[0].xyz / centroid_t[0].w;

  // accumulate the covariance smatrix
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static Sym3x3 covariance[8];
#endif

  // same for all
  const int tf = (thread >> 1);
  const int of = (thread & 1);
  const int dm = (tf << 1) + 0;
  const int uq = (tf << 1) + 1;

  threaded_for(co1, 8) {
    float3 a =  points[dn] - _centroid;
    float3 A =  points[up] - _centroid;
    float3 b = weights[dn] * a;
    float3 B = weights[up] * A;

    covariance[co1][0] = float4(
      a.x * b.x + A.x * B.x,
      a.x * b.y + A.x * B.y,
      a.x * b.z + A.x * B.z, 0
    );

    covariance[co1][1] = float4(
      a.y * b.y + A.y * B.y,
      a.y * b.z + A.y * B.z,
      a.z * b.z + A.z * B.z, 0
    );
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // AMP: prefer 3-wide vectorized op
  threaded_for(co2, 4*2) {
    covariance[tf][of] = covariance[dm][of] + covariance[uq][of]; }
  threaded_for(co3, 2*2) {
    covariance[tf][of] = covariance[dm][of] + covariance[uq][of]; }
  threaded_for(co4, 1*2) {
    covariance[tf][of] = covariance[dm][of] + covariance[uq][of]; }

  // make it visible
  tile_static_memory_fence(barrier);

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

  m[0].x = smatrix[0].x - evalue;
  m[0].y = smatrix[0].y;
  m[0].z = smatrix[0].z;
  m[1].x = smatrix[1].x - evalue;
  m[1].y = smatrix[1].y;
  m[1].z = smatrix[1].z - evalue;

  // compute U
  Sym3x3 u;

  u[0].x = m[1].x * m[1].z - m[1].y * m[1].y;
  u[0].y = m[0].z * m[1].y - m[0].y * m[1].z;
  u[0].z = m[0].y * m[1].y - m[0].z * m[1].x;
  u[1].x = m[0].x * m[1].z - m[0].z * m[0].z;
  u[1].y = m[0].y * m[0].z - m[1].y * m[0].x;
  u[1].z = m[0].x * m[1].x - m[0].y * m[0].y;

  // find the largest component
  float
      mc = fabsf(u[0].x); int mi = 0;
  if (mc < fabsf(u[0].y)) {
      mc = fabsf(u[0].y); mi = 1; }
  if (mc < fabsf(u[0].z)) {
      mc = fabsf(u[0].z); mi = 2; }
  if (mc < fabsf(u[1].x)) {
      mc = fabsf(u[1].x); mi = 3; }
  if (mc < fabsf(u[1].y)) {
      mc = fabsf(u[1].y); mi = 4; }
  if (mc < fabsf(u[1].z)) {
      mc = fabsf(u[1].z); mi = 5; }

  // pick the column with this component
  switch (mi) {
    case 0:
      return float3(u[0].x, u[0].y, u[0].z);

    case 1:
    case 3:
      return float3(u[0].y, u[1].x, u[1].y);

    default:
      return float3(u[0].z, u[1].y, u[1].z);
  }
}

static float3 GetMultiplicity2Evector(tile_barrier barrier, const int thread, Sym3x3r smatrix, float evalue) amp_restricted
{
  // compute M
  Sym3x3 m;

  m[0].x = smatrix[0].x - evalue;
  m[0].y = smatrix[0].y;
  m[0].z = smatrix[0].z;
  m[1].x = smatrix[1].x - evalue;
  m[1].y = smatrix[1].y;
  m[1].z = smatrix[1].z - evalue;

  // find the largest component
  float
      mc = fabsf(m[0].x); int mi = 0;
  if (mc < fabsf(m[0].y)) {
      mc = fabsf(m[0].y); mi = 1; }
  if (mc < fabsf(m[0].z)) {
      mc = fabsf(m[0].z); mi = 2; }
  if (mc < fabsf(m[1].x)) {
      mc = fabsf(m[1].x); mi = 3; }
  if (mc < fabsf(m[1].y)) {
      mc = fabsf(m[1].y); mi = 4; }
  if (mc < fabsf(m[1].z)) {
      mc = fabsf(m[1].z); mi = 5; }

  // pick the first eigenvector based on this index
  switch (mi) {
    case 0:
    case 1:
      return float3(-m[0].y, m[0].x,  0.0f);

    case 2:
      return float3( m[0].z, 0.0f, -m[0].x);

    case 3:
    case 4:
      return float3( 0.0f, -m[1].y, m[1].x);

    default:
      return float3( 0.0f, -m[1].z, m[1].y);
  }
}

float3 ComputePrincipleComponent(tile_barrier barrier, const int thread, Sym3x3r smatrix) amp_restricted
{
  // compute the cubic coefficients
  float c0 = smatrix[0].x * smatrix[1].x * smatrix[1].z
    + 2.0f * smatrix[0].y * smatrix[0].z * smatrix[1].y
    -        smatrix[0].x * smatrix[1].y * smatrix[1].y
    -        smatrix[1].x * smatrix[0].z * smatrix[0].z
    -        smatrix[1].z * smatrix[0].y * smatrix[0].y;
  float c1 = smatrix[0].x * smatrix[1].x
    +        smatrix[0].x * smatrix[1].z
    +        smatrix[1].x * smatrix[1].z
    -        smatrix[0].y * smatrix[0].y
    -        smatrix[0].z * smatrix[0].z
    -        smatrix[1].y * smatrix[1].y;
  float c2 = smatrix[0].x
    +        smatrix[1].x
    +        smatrix[1].z;

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

} // namespace squish
