/* -----------------------------------------------------------------------------

	Copyright (c) 2006 Simon Brown                          si@sjbrown.co.uk
	Copyright (c) 2007 Ignacio Castano                   icastano@nvidia.com
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

#include <assert.h>

#include "bitoneclusterfit.h"
#include "bitoneset.h"
#include "bitoneblock.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
BitoneClusterFit::BitoneClusterFit(BitoneSet const* bitones, int flags)
  : BitoneFit(bitones, flags)
{
  // set the iteration count
  m_iterationCount = (m_flags & kColourIterativeClusterFits) / kColourClusterFit;

  // initialize the best error
  m_besterror = Scr4(FLT_MAX);

  // cache some values
  bool const unweighted = m_bitones->IsUnweighted();
  int const count = m_bitones->GetCount();
  Vec3 const* values = m_bitones->GetPoints();

  Sym3x3 covariance;
  Vec3 centroid;

  // get the covariance matrix
  if (unweighted)
    ComputeWeightedCovariance3(covariance, centroid, count, values, Vec3(1.0f));
  else
    ComputeWeightedCovariance3(covariance, centroid, count, values, Vec3(1.0f), m_bitones->GetWeights());

  // compute the principle component
  GetPrincipleComponent(covariance, m_principle);

  // we have tables for this
  m_optimizable = unweighted & (count == 16);
}

bool BitoneClusterFit::ConstructOrdering(Vec3 const& axis, int iteration)
{
  // cache some values
  int const count = m_bitones->GetCount();
  Vec3 const* values = m_bitones->GetPoints();

  // build the list of dot products
  float dps[16];
  u8* order = (u8*)m_order + 16 * iteration;
  for (int i = 0; i < count; ++i) {
    Dot(values[i], axis, dps + i);
    order[i] = (u8)i;
  }

  // stable sort using them
  for (int i = 0; i < count; ++i) {
    for (int j = i; j > 0 && dps[j] < dps[j - 1]; --j) {
      std::swap(dps  [j], dps  [j - 1]);
      std::swap(order[j], order[j - 1]);
    }
  }

  // check this ordering is unique
  Col4 curr; LoadAligned(curr, m_order + 16 * iteration);
  for (int it = 0; it < iteration; ++it) {
    Col4 prev; LoadAligned(prev, m_order + 16 * it);

    if (CompareAllEqualTo(curr, prev))
      return false;
  }

  // copy the ordering and weight all the points
  Vec3 const* unweighted = m_bitones->GetPoints();
  Scr3 const* weights = m_bitones->GetWeights();
  m_xsum_wsum = VEC4_CONST(0.0f);
  for(int i = 0; i < count; ++i) {
    int j = order[i];

    Vec4 p = Vec4(unweighted[j], 1.0f);
    Vec4 w(weights[j]);
    Vec4 x = p * w;

    m_points_weights[i] = x;
    m_xsum_wsum        += x;
  }

  return true;
}

void BitoneClusterFit::ClusterFit4(void* block)
{
  // declare variables
  int const count = m_bitones->GetCount();
  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const one = VEC4_CONST(1.0f);

  Vec4 const onethird_onethird2  (1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 9.0f);
  Vec4 const twothirds_twothirds2(2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f, 4.0f / 9.0f);
  Vec4 const twonineths                                     = VEC4_CONST(2.0f / 9.0f);

  assume((count > 0) && (count <= 16));

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle, 0);

  // check all possible clusters and iterate on the total order
  Vec4 beststart = VEC4_CONST(0.0f);
  Vec4 bestend = VEC4_CONST(0.0f);
  Scr4 besterror = m_besterror;
  u8 bestindices[16];
  int bestiteration = 0;
  int besti = 0, bestj = 0, bestk = 0;

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {
      // second cluster [i,j) is one third along
      Vec4 part1 = VEC4_CONST(0.0f);
      for (int j = i;;) {
	// third cluster [j,k) is two thirds along
	Vec4 part2 = (j == 0) ? m_points_weights[0] : VEC4_CONST(0.0f);
	int kmin = (j == 0) ? 1 : j;
	for (int k = kmin;;) {
	  // last cluster [k,count) is at the end
	  Vec4 part3 = m_xsum_wsum - part2 - part1 - part0;

	  // compute least squares terms directly
	  Vec4 const alphax_sum = MultiplyAdd(part2, onethird_onethird2, MultiplyAdd(part1, twothirds_twothirds2, part0));
	  Vec4 const  betax_sum = MultiplyAdd(part1, onethird_onethird2, MultiplyAdd(part2, twothirds_twothirds2, part3));

	  Vec4 const alpha2_sum = alphax_sum.SplatW();
	  Vec4 const  beta2_sum =  betax_sum.SplatW();

	  Vec4 const alphabeta_sum = twonineths * (part1 + part2).SplatW();

	  // compute the least-squares optimal points
	  Vec4 factor = Reciprocal(NegativeMultiplySubtract(alphabeta_sum, alphabeta_sum, alpha2_sum * beta2_sum));
	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_sum, alphax_sum *  beta2_sum) * factor;
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_sum,  betax_sum * alpha2_sum) * factor;

	  // snap floating-point-values to the integer-lattice
	  a = Truncate(a * 255.0f) * (1.0f / 255.0f);
	  b = Truncate(b * 255.0f) * (1.0f / 255.0f);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = e4;

	  // keep the solution if it wins
	  if (besterror > eS) {
	    besterror = eS;
	    beststart = a;
	    bestend = b;
	    besti = i;
	    bestj = j;
	    bestk = k;
	    bestiteration = iterationIndex;
	  }

	  // advance
	  if (k == count) break;
	  part2 += m_points_weights[k]; ++k;
	}

	// advance
	if (j == count) break;
	part1 += m_points_weights[j]; ++j;
      }

      // advance
      part0 += m_points_weights[i];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;

    // advance if possible
    ++iterationIndex;
    if (iterationIndex == m_iterationCount)
      break;

    // stop if a new iteration is an ordering that has already been tried
    Vec3 axis = (bestend - beststart).GetVec3();
    if (!ConstructOrdering(axis, iterationIndex))
      break;
  }

  // save the block if necessary
  if (besterror < m_besterror) {
    // save the error
    m_besterror = besterror;

    // remap the indices
    u8 const* order = (u8*)m_order + 16 * bestiteration;

    u8 unordered[16];
    for (int m =     0; m < besti; ++m)
      unordered[order[m]] = 0;
    for (int m = besti; m < bestj; ++m)
      unordered[order[m]] = 2;
    for (int m = bestj; m < bestk; ++m)
      unordered[order[m]] = 3;
    for (int m = bestk; m < count; ++m)
      unordered[order[m]] = 1;

    m_bitones->RemapIndices(unordered, bestindices);

    // save the block
    WriteBitoneBlock4(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
}

#include "bitoneclusterfit.inl"

void BitoneClusterFit::ClusterFit4Constant(void* block)
{
  // declare variables
  int const count = m_bitones->GetCount();
  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const one = VEC4_CONST(1.0f);

  Vec4 const onethird_onethird2  (1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 9.0f);
  Vec4 const twothirds_twothirds2(2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f, 4.0f / 9.0f);
  Vec4 const twonineths                                     = VEC4_CONST(2.0f / 9.0f);

  assume((count > 0) && (count <= 16));

  // check all possible clusters and iterate on the total order
  Vec4 beststart = VEC4_CONST(0.0f);
  Vec4 bestend = VEC4_CONST(0.0f);
  Scr4 besterror = m_besterror;
  u8 bestindices[16];
  int bestiteration = 0;
  int besti = 0, bestj = 0, bestk = 0;

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle, 0);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_wsum = m_xsum_wsum;

    // constants if weights == 1
    Vec4 alphabeta_dltas  = *((Vec4 *)part2delta[0]);
    Vec4 *alphabeta_inits = (Vec4 *)part2inits[0];
    float *alphabeta_factors = (float *)part2factors;

#if 0
  Vec4 lasta = Vec4(0.0f);
  Vec4 lastb = xsum_wsum;
  Vec4 lastc = Vec4(0.0f);
#endif

    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

    // second cluster [i,j) is one third along
    Vec4 part1 = VEC4_CONST(0.0f);
    for (int j = i;;) {

    // third cluster [j,k) is two thirds along
    Vec4 part2 = (j == 0) ? m_points_weights[0] : VEC4_CONST(0.0f);
    Vec4 alphabeta_val = *alphabeta_inits++;
    int kmin = (j == 0) ? 1 : j;
    for (int k = kmin;;) {
	  // TODO: the inner alphabeta_sum seems always to be the same sequence
	  Vec4 alphabeta_factor = alphabeta_val * Vec4(*alphabeta_factors++);

	  // compute least squares terms directly
	  Vec4 const alphax_sum =   MultiplyAdd(part2, onethird_onethird2, MultiplyAdd(part1, twothirds_twothirds2, part0));
	  Vec4 const  betax_sum = /*MultiplyAdd(part1, onethird_onethird2, MultiplyAdd(part2, twothirds_twothirds2, part3))*/ xsum_wsum - alphax_sum;

	  Vec4 const    alpha2_sum = alphabeta_val.SplatX();
	  Vec4 const     beta2_sum = alphabeta_val.SplatY();
	  Vec4 const alphabeta_sum = alphabeta_val.SplatZ();

	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_factor.SplatZ(), alphax_sum * alphabeta_factor.SplatY());
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_factor.SplatZ(),  betax_sum * alphabeta_factor.SplatX());

#if 0
	  // last cluster [k,count) is at the end
	  Vec4 part3 = xsum_wsum - part2 - part1 - part0;

	  // compute least squares terms directly
	  Vec4 const _alphax_sum = MultiplyAdd(part2, onethird_onethird2, MultiplyAdd(part1, twothirds_twothirds2, part0));
	  Vec4 const  _betax_sum = MultiplyAdd(part1, onethird_onethird2, MultiplyAdd(part2, twothirds_twothirds2, part3));
//	  Vec4 const  _betac_sum = xsum_wsum - _alphax_sum;

	  Vec4 const _alpha2_sum = _alphax_sum.SplatW();
	  Vec4 const  _beta2_sum =  _betax_sum.SplatW();

	  Vec4 const _alphabeta_sum = twonineths * (part1 + part2).SplatW();

	  // compute the least-squares optimal points
	  Vec4 _factor = Reciprocal(NegativeMultiplySubtract(_alphabeta_sum, _alphabeta_sum, _alpha2_sum * _beta2_sum));
	  Vec4 _a = NegativeMultiplySubtract( _betax_sum, _alphabeta_sum, _alphax_sum *  _beta2_sum) * _factor;
	  Vec4 _b = NegativeMultiplySubtract(_alphax_sum, _alphabeta_sum,  _betax_sum * _alpha2_sum) * _factor;

#undef	limit
#define	limit 2e-5
	  assert(fabs(_alpha2_sum.W() - alpha2_sum.X()) < limit);
	  assert(fabs(_beta2_sum.W() - beta2_sum.X()) < limit);
	  assert(fabs(_alphabeta_sum.W() - alphabeta_sum.X()) < limit);

	  if (alphabeta_factors[-1] != FLT_MAX) {
	    assert(fabs(_factor.W() - alphabeta_factors[-1]) < limit);

	    assert(fabs(_alpha2_sum.W()    * _factor.W() - alphabeta_factor.X()) < limit);
	    assert(fabs(_beta2_sum.W()     * _factor.W() - alphabeta_factor.Y()) < limit);
	    assert(fabs(_alphabeta_sum.W() * _factor.W() - alphabeta_factor.Z()) < limit);

	    assert(fabs(a.X() - _a.X()) < limit);
	    assert(fabs(a.Y() - _a.Y()) < limit);
	    assert(fabs(a.Z() - _a.Z()) < limit);

	    assert(fabs(b.X() - _b.X()) < limit);
	    assert(fabs(b.Y() - _b.Y()) < limit);
	    assert(fabs(b.Z() - _b.Z()) < limit);
	  }

#if 0
	  fprintf(stderr, "{%.9ff},", _factor.W());
	  if (k == kmin)
	    fprintf(stderr, "{%.9f, %.9f, %.9f},\n", alpha2_sum.W(), beta2_sum.W(), alphabeta_sum.W());
	  fprintf(stderr, "{%.9f/*%.9f*/,%.9f/*%.9f*/,%.9f,%.9f},\n",
	    alpha2_sum.W(), lasta.W() - alpha2_sum.W(),
	    beta2_sum.W(), lastb.W() - beta2_sum.W(),
	    alphabeta_sum.W(), lastc.W() - alphabeta_sum.W(),
	    factor.W());

	  lasta = alpha2_sum;
	  lastb = beta2_sum;
	  lastc = alphabeta_sum;
#endif
#endif

	  // snap floating-point-values to the integer-lattice
	  a = Truncate(a * 255.0f) * (1.0f / 255.0f);
	  b = Truncate(b * 255.0f) * (1.0f / 255.0f);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = e4;

	  // keep the solution if it wins
	  if (besterror > eS) {
	    besterror = eS;

	    beststart = a;
	    bestend = b;
	    bestiteration = iterationIndex;

	    besti = i;
	    bestj = j;
	    bestk = k;
	  }

      alphabeta_val += alphabeta_dltas;

      // advance
      if (k == count) break;
      part2 += m_points_weights[k]; ++k; }

      // advance
      if (j == count) break;
      part1 += m_points_weights[j]; ++j; }

      // advance
      part0 += m_points_weights[i];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;

    // advance if possible
    ++iterationIndex;
    if (iterationIndex == m_iterationCount)
      break;

    // stop if a new iteration is an ordering that has already been tried
    Vec3 axis = (bestend - beststart).GetVec3();
    if (!ConstructOrdering(axis, iterationIndex))
      break;
  }

  // save the block if necessary
  if (besterror < m_besterror) {
    // save the error
    m_besterror = besterror;

    // remap the indices
    u8 const* order = (u8*)m_order + 16 * bestiteration;

    u8 unordered[16];
    for (int m =     0; m < besti; ++m)
      unordered[order[m]] = 0;
    for (int m = besti; m < bestj; ++m)
      unordered[order[m]] = 2;
    for (int m = bestj; m < bestk; ++m)
      unordered[order[m]] = 3;
    for (int m = bestk; m < count; ++m)
      unordered[order[m]] = 1;

    m_bitones->RemapIndices(unordered, bestindices);

    // save the block
    WriteBitoneBlock4(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
}

void BitoneClusterFit::Compress4(void* block)
{
#if defined(TRACK_STATISTICS)
  /* there is a clear skew towards unweighted clusterfit (all weights == 1)
   *
   * C == 4, numset ==
   *  [0]	0x00f796e0 {13856, 5184}
   */
  if (m_optimizable & (m_bitones->GetCount() == 16))
    gstat.has_noweightsets[1][0][0]++;
  else
    gstat.has_noweightsets[1][0][1]++;
#endif
  
  if (m_optimizable & (m_bitones->GetCount() == 16))
    ClusterFit4Constant(block);
  else
    ClusterFit4        (block);
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
