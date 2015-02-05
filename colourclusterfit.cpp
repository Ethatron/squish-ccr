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

#include "colourclusterfit.h"
#include "colourset.h"
#include "colourblock.h"

#include "coloursinglefit.h"
#include "coloursinglesnap.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
ColourClusterFit::ColourClusterFit(ColourSet const* colours, int flags)
  : ColourFit(colours, flags)
{
  // set the iteration count
  m_iterationCount = (m_flags & kColourIterativeClusterFits) / kColourClusterFit;
  
  // initialize endpoints
  ComputeEndPoints();
}

void ColourClusterFit::ComputeEndPoints()
{
  // cache some values
  bool const unweighted = m_colours->IsUnweighted();
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();

  Sym3x3 covariance;
  Vec3 centroid;

  // get the covariance matrix
  if (unweighted)
    ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric);
  else
    ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric, m_colours->GetWeights());

  // compute the principle component
  GetPrincipleComponent(covariance, m_principle);

  // we have tables for this
  m_optimizable = unweighted & ((count == 16) | ((m_flags & kBtcp) == kBtc1));
}

void ColourClusterFit::SumError3(u8 (&closest)[16], Vec4 &beststart, Vec4 &bestend, Scr4 &besterror) {
  cQuantizer4<5,6,5,0> q = cQuantizer4<5,6,5,0>();
  
  // snap floating-point-values to the integer-lattice
  Vec3 start = q.SnapToLattice(beststart);
  Vec3 end   = q.SnapToLattice(bestend);
  
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();

  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[3]; Codebook3(codes, m_metric * start, m_metric * end);

  // error-sum
  for (int i = 0; i < count; ++i)
    besterror += LengthSquared(m_metric * values[i] - codes[closest[i]]) * freq[i];
}

void ColourClusterFit::SumError4(u8 (&closest)[16], Vec4 &beststart, Vec4 &bestend, Scr4 &besterror) {
  cQuantizer4<5,6,5,0> q = cQuantizer4<5,6,5,0>();
  
  // snap floating-point-values to the integer-lattice
  Vec3 start = q.SnapToLattice(beststart);
  Vec3 end   = q.SnapToLattice(bestend);
  
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();

  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[4]; Codebook4(codes, m_metric * start, m_metric * end);

  // error-sum
  for (int i = 0; i < count; ++i)
    besterror += LengthSquared(m_metric * values[i] - codes[closest[i]]) * freq[i];
}

bool ColourClusterFit::ConstructOrdering(Vec3 const& axis, int iteration)
{
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();

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
  Vec3 const* unweighted = m_colours->GetPoints();
  Scr3 const* weights = m_colours->GetWeights();

  m_xsum_wsum = VEC4_CONST(0.0f);
  for (int i = 0; i < count; ++i) {
    int j = order[i];

    Vec4 p = Vec4(unweighted[j], 1.0f);
    Scr4 w(weights[j]);
    Vec4 x = p * w;

    m_points_weights[i] = x;
    m_xsum_wsum        += x;
  }

  return true;
}

#ifdef	FEATURE_METRIC_SQUARED
#define CMetric(m)  m * m
#else
#define CMetric(m)  m
#endif

void ColourClusterFit::ClusterFit3(void* block)
{
  cQuantizer4<5,6,5,0> q = cQuantizer4<5,6,5,0>();

  // declare variables
  int const count = m_colours->GetCount();
  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const one = VEC4_CONST(1.0f);
  Vec4 const half_half2(0.5f, 0.5f, 0.5f, 0.25f);

  assume((count > 0) && (count <= 16));

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle, 0);

  // check all possible clusters and iterate on the total order
  Vec4 beststart = VEC4_CONST(0.0f);
  Vec4 bestend = VEC4_CONST(0.0f);
  Scr4 besterror = m_besterror;
  a16 u8 bestindices[16];
  int bestiteration = 0;
  int besti = 0, bestj = 0;

  // metric is squared as well
  Vec4 cmetric = CMetric(m_metric);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {
      // second cluster [i,j) is half along
      Vec4 part1 = (i == 0) ? m_points_weights[0] : VEC4_CONST(0.0f);
      int jmin = (i == 0) ? 1 : i;
      for (int j = jmin;;) {
	// last cluster [j,count) is at the end
	Vec4 part2 = m_xsum_wsum - part1 - part0;

	// compute least squares terms directly
	Vec4 alphax_sum = MultiplyAdd(part1, half_half2, part0);
	Vec4 alpha2_sum = alphax_sum.SplatW();

	Vec4 betax_sum = MultiplyAdd(part1, half_half2, part2);
	Vec4 beta2_sum = betax_sum.SplatW();

	Vec4 alphabeta_sum = (part1 * half_half2).SplatW();

	// compute the least-squares optimal points
	Vec4 factor = Reciprocal( NegativeMultiplySubtract(alphabeta_sum, alphabeta_sum, alpha2_sum * beta2_sum));
	Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_sum, alphax_sum *  beta2_sum) * factor;
	Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_sum,  betax_sum * alpha2_sum) * factor;

	// snap floating-point-values to the integer-lattice
	a = q.SnapToLattice(a);
	b = q.SnapToLattice(b);

	// compute the error (we skip the constant xxsum)
	Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	Vec4 e4 = MultiplyAdd(two, e3, e1);

	// apply the metric to the error term
	Scr4 eS = Dot(e4, cmetric);

	// keep the solution if it wins
	if (besterror > eS) {
	  besterror = eS;
	  beststart = a;
	  bestend = b;
	  besti = i;
	  bestj = j;
	  bestiteration = iterationIndex;
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
  
#if 1
  u8 const* order = (u8*)m_order + 16 * bestiteration;

  a16 u8 unordered[16];
  for (int m =     0; m < besti; ++m)
    unordered[order[m]] = 0;
  for (int m = besti; m < bestj; ++m)
    unordered[order[m]] = 2;
  for (int m = bestj; m < count; ++m)
    unordered[order[m]] = 1;
    
  // save the block if necessary
  besterror = Scr4(0.0f);
  SumError3(unordered, beststart, bestend, besterror);
  if (Scr3(besterror) < m_besterror) {
    // save the error
    m_besterror = besterror;
    
    // remap the indices
    m_colours->RemapIndices(unordered, bestindices);
  
    // save the block
    WriteColourBlock3(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#else
  Vec4 m_xxsum_wwsum = VEC4_CONST(0.0f);
  for (int i = 0; i < count; ++i)
    m_xxsum_wwsum += m_points_weights[i] * m_points_weights[i];

  besterror += Dot(m_xxsum_wwsum, cmetric);

  // save the block if necessary
  if (besterror < m_besterror) {
    // save the error
    m_besterror = besterror;

    // remap the indices
    u8 const* order = (u8*)m_order + 16 * bestiteration;

    a16 u8 unordered[16];
    for (int m =     0; m < besti; ++m)
      unordered[order[m]] = 0;
    for (int m = besti; m < bestj; ++m)
      unordered[order[m]] = 2;
    for (int m = bestj; m < count; ++m)
      unordered[order[m]] = 1;

    m_colours->RemapIndices(unordered, bestindices);
    
    // save the block
    WriteColourBlock3(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#endif
}

void ColourClusterFit::ClusterFit4(void* block)
{
  cQuantizer4<5,6,5,0> q = cQuantizer4<5,6,5,0>();

  // declare variables
  int const count = m_colours->GetCount();
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

  // metric is squared as well
  Vec4 cmetric = CMetric(m_metric);

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
	  a = q.SnapToLattice(a);
	  b = q.SnapToLattice(b);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = Dot(e4, cmetric);

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
  
#if 1
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

  // save the block if necessary
  besterror = Scr4(0.0f);
  SumError4(unordered, beststart, bestend, besterror);
  if (Scr3(besterror) < m_besterror) {
    // save the error
    m_besterror = besterror;
    
    // remap the indices
    m_colours->RemapIndices(unordered, bestindices);
  
    // save the block
    WriteColourBlock4(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#else
  Vec4 m_xxsum_wwsum = VEC4_CONST(0.0f);
  for (int i = 0; i < count; ++i)
    m_xxsum_wwsum += m_points_weights[i] * m_points_weights[i];

  besterror += Dot(m_xxsum_wwsum, cmetric);

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

    m_colours->RemapIndices(unordered, bestindices);

    // save the block
    WriteColourBlock4(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#endif
}

#include "colourclusterfit.inl"

void ColourClusterFit::ClusterFit3Constant(void* block)
{
  cQuantizer4<5,6,5,0> q = cQuantizer4<5,6,5,0>();

  // declare variables
  int const count = m_colours->GetCount();
  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const one = VEC4_CONST(1.0f);
  Vec4 const half_half2(0.5f, 0.5f, 0.5f, 0.25f);

  assume((count > 0) && (count <= 16));

  // check all possible clusters and iterate on the total order
  Vec4 beststart = VEC4_CONST(0.0f);
  Vec4 bestend = VEC4_CONST(0.0f);
  Scr4 besterror = m_besterror;
  a16 u8 bestindices[16];
  int bestiteration = 0;
  int besti = 0, bestj = 0;

  // metric is squared as well
  Vec4 cmetric = CMetric(m_metric);

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle, 0);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_wsum = m_xsum_wsum;

    // constants if weights == 1
    Vec4 alphabeta_dltas  = *((Vec4 *)part1delta[count - 1][0]);
    Vec4 *alphabeta_inits = (Vec4 *)part1inits[count - 1][0];
    float *alphabeta_factors = (float *)part1factors[count - 1];

#if 0
  Vec4 lasta = Vec4(0.0f);
  Vec4 lastb = xsum_wsum;
  Vec4 lastc = Vec4(0.0f);
#endif

    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

    // second cluster [i,j) is half along
    Vec4 part1 = (i == 0) ? m_points_weights[0] : VEC4_CONST(0.0f);
    Vec4 alphabeta_val = *alphabeta_inits++;
    int jmin = (i == 0) ? 1 : i;
    for (int j = jmin;;) {
	// TODO: the inner alphabeta_sum seems always to be the same sequence
  	Vec4 alphabeta_factor = alphabeta_val * Vec4(*alphabeta_factors++);

	// compute least squares terms directly
	Vec4 const alphax_sum =   MultiplyAdd(part1, half_half2, part0);
	Vec4 const  betax_sum = /*MultiplyAdd(part1, half_half2, part2)*/ xsum_wsum - alphax_sum;

	Vec4 const    alpha2_sum = alphabeta_val.SplatX();
	Vec4 const     beta2_sum = alphabeta_val.SplatY();
	Vec4 const alphabeta_sum = alphabeta_val.SplatZ();

	Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_factor.SplatZ(), alphax_sum * alphabeta_factor.SplatY());
	Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_factor.SplatZ(),  betax_sum * alphabeta_factor.SplatX());

#if 0
	// last cluster [j,count) is at the end
	Vec4 part2 = xsum_wsum - part1 - part0;

	// compute least squares terms directly
	Vec4 const _alphax_sum = MultiplyAdd(part1, half_half2, part0);
	Vec4 const  _betax_sum = MultiplyAdd(part1, half_half2, part2);

	Vec4 const _alpha2_sum = _alphax_sum.SplatW();
	Vec4 const  _beta2_sum =  _betax_sum.SplatW();

	Vec4 const _alphabeta_sum = (part1 * half_half2).SplatW();

	// compute the least-squares optimal points
	Vec4 const _factor = Reciprocal(NegativeMultiplySubtract(_alphabeta_sum, _alphabeta_sum, _alpha2_sum * _beta2_sum));
	Vec4 const _a = NegativeMultiplySubtract( _betax_sum, _alphabeta_sum, _alphax_sum *  _beta2_sum) * _factor;
	Vec4 const _b = NegativeMultiplySubtract(_alphax_sum, _alphabeta_sum,  _betax_sum * _alpha2_sum) * _factor;

#define	limit 1e-5
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
	  fprintf(stderr, "{%.9ff},", factor.W());
	  if (j == jmin)
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
	a = q.SnapToLattice(a);
	b = q.SnapToLattice(b);

	// compute the error (we skip the constant xxsum)
	Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	Vec4 e4 = MultiplyAdd(two, e3, e1);

	// apply the metric to the error term
	Scr4 eS = Dot(e4, cmetric);

	// keep the solution if it wins
	if (besterror > eS) {
	  besterror = eS;
	  beststart = a;
	  bestend = b;
	  besti = i;
	  bestj = j;
	  bestiteration = iterationIndex;
	}

      alphabeta_val += alphabeta_dltas;

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
  
#if 1
  u8 const* order = (u8*)m_order + 16 * bestiteration;
  
  a16 u8 unordered[16];
  for (int m =     0; m < besti; ++m)
    unordered[order[m]] = 0;
  for (int m = besti; m < bestj; ++m)
    unordered[order[m]] = 2;
  for (int m = bestj; m < count; ++m)
    unordered[order[m]] = 1;

  // save the block if necessary
  besterror = Scr4(0.0f);
  SumError3(unordered, beststart, bestend, besterror);
  if (Scr3(besterror) < m_besterror) {
    // save the error
    m_besterror = besterror;
    
    // remap the indices
    m_colours->RemapIndices(unordered, bestindices);
  
    // save the block
    WriteColourBlock3(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#else
  Vec4 m_xxsum_wwsum = VEC4_CONST(0.0f);
  for (int i = 0; i < count; ++i)
    m_xxsum_wwsum += m_points_weights[i] * m_points_weights[i];

  besterror += Dot(m_xxsum_wwsum, cmetric);

  // save the block if necessary
  if (besterror < m_besterror) {
    // save the error
    m_besterror = besterror;

    // remap the indices
    u8 const* order = (u8*)m_order + 16 * bestiteration;

    a16 u8 unordered[16];
    for (int m =     0; m < besti; ++m)
      unordered[order[m]] = 0;
    for (int m = besti; m < bestj; ++m)
      unordered[order[m]] = 2;
    for (int m = bestj; m < count; ++m)
      unordered[order[m]] = 1;

    m_colours->RemapIndices(unordered, bestindices);

    // save the block
    WriteColourBlock3(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#endif
}

void ColourClusterFit::ClusterFit4Constant(void* block)
{
  cQuantizer4<5,6,5,0> q = cQuantizer4<5,6,5,0>();

  // declare variables
  int const count = m_colours->GetCount();
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

  // metric is squared as well
  Vec4 cmetric = CMetric(m_metric);

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
	  a = q.SnapToLattice(a);
	  b = q.SnapToLattice(b);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = Dot(e4, cmetric);

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
  
#if 1
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

  // save the block if necessary
  besterror = Scr4(0.0f);
  SumError4(unordered, beststart, bestend, besterror);
  if (Scr3(besterror) < m_besterror) {
    // save the error
    m_besterror = besterror;
    
    // remap the indices
    m_colours->RemapIndices(unordered, bestindices);
  
    // save the block
    WriteColourBlock4(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#else
  Vec4 m_xxsum_wwsum = VEC4_CONST(0.0f);
  for (int i = 0; i < count; ++i)
    m_xxsum_wwsum += m_points_weights[i] * m_points_weights[i];

  besterror += Dot(m_xxsum_wwsum, cmetric);

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

    m_colours->RemapIndices(unordered, bestindices);

    // save the block
    WriteColourBlock4(beststart.GetVec3(), bestend.GetVec3(), bestindices, block);
  }
#endif
}

void ColourClusterFit::Compress3b(void* block)
{
  ColourSet copy = *m_colours;
  m_colours = &copy;
  
  Scr3 m_destroyed = Scr3(0.0f);
  while (copy.RemoveBlack(m_metric, m_destroyed) && !(m_besterror < m_destroyed)) {
    m_besterror -= m_destroyed;
    
    if (copy.GetCount() == 1) {
      // always do a single colour fit
      ColourSingleMatch fit(m_colours, m_flags);

      fit.SetError(m_besterror);
      fit.Compress(block);

      m_besterror = fit.GetError();
    }
    else {
      ComputeEndPoints();
      Compress3(block);
    }

    m_besterror += m_destroyed;
  }
}

void ColourClusterFit::Compress3(void* block)
{
#if defined(TRACK_STATISTICS)
  /* there is a clear skew towards unweighted clusterfit (all weights == 1)
   *
   * C == 3, numset ==
   *  [0]	0x00f796e0 {124800, 15616}
   */
  if (m_optimizable)
    gstat.has_noweightsets[0][0][0]++;
  else
    gstat.has_noweightsets[0][0][1]++;
#endif

  if (m_optimizable)
    ClusterFit3Constant(block);
  else
    ClusterFit3        (block);
}

void ColourClusterFit::Compress4(void* block)
{
#if defined(TRACK_STATISTICS)
  /* there is a clear skew towards unweighted clusterfit (all weights == 1)
   *
   * C == 4, numset ==
   *  [0]	0x00f796e0 {13856, 5184}
   */
  if (m_optimizable & (m_colours->GetCount() == 16))
    gstat.has_noweightsets[1][0][0]++;
  else
    gstat.has_noweightsets[1][0][1]++;
#endif

  if (m_optimizable & (m_colours->GetCount() == 16))
    ClusterFit4Constant(block);
  else
    ClusterFit4        (block);
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
void ClusterFit_CCR::AssignSet(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, const int metric, const int fit ) amp_restricted
{
  ColourFit_CCR::AssignSet(barrier, thread, m_colours, metric, fit);

  threaded_cse(0) {
    // set the iteration count
    m_iterationCount = (fit == SQUISH_FIT_CLUSTER) ? 1 : kMaxIterations;

    // initialize the best error
    m_besterror = FLT_MAX;

    // initialize the metric
    if (metric == SQUISH_METRIC_UNIT)
      m_metric4 = float4( 1.0f, 1.0f, 0.0f, 0.0f );
    else if (metric == SQUISH_METRIC_PERCEPTUAL)
      m_metric4 = float4( 0.2126f, 0.7152f, 0.0722f, 0.0f );
    else
      m_metric4 = float4( 1.0f, 1.0f, 1.0f, 0.0f );
  }

  // cache some values
  // AMP: causes a full copy, for indexibility
  const int count = m_colours.GetCount();
  point16 values = m_colours.GetPoints();
  weight16 weights = m_colours.GetWeights();

  // get the covariance smatrix
  Sym3x3 covariance = ComputeWeightedCovariance3(barrier, thread, count, values, weights);

  // compute the principle component
  float3 principle = GetPrincipleComponent(barrier, thread, covariance);

  threaded_cse(0) {
    // compute the principle component
    m_principle = principle;
  }
}

bool ClusterFit_CCR::ConstructOrdering(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, float3r axis, int iteration) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // cache some values (AMP: causes a full copy, for indexibility)
  const int count = m_colours.GetCount();
  point16 values = m_colours.GetPoints();

  // build the list of dot products
  float dps[16];
  index16 order = m_order[iteration];

  for (int i = 0; i < count; ++i) {
    dps  [i] = dot(values[i], axis);
    order[i] = (ccr8)i;
  }

  // stable sort using them
  for (int i = 0; i < count; ++i) {
    for (int j = i; j > 0 && dps[j] < dps[j - 1]; --j) {
      float c =   dps[j];   dps[j] =   dps[j - 1];   dps[j - 1] = c;
      ccr8  d = order[j]; order[j] = order[j - 1]; order[j - 1] = d;
    }
  }

  // check this ordering is unique
  for (int it = 0; it < iteration; ++it) {
    index16 prev = m_order[it];
    bool same = true;

    for (int i = 0; i < count; ++i) {
      if (order[i] != prev[i]) {
	same = false;
	break;
      }
    }

    if (same)
      return false;
  }

  // copy the ordering and weight all the points
  // AMP: causes a full copy, for indexibility
  point16 unweighted = m_colours.GetPoints();
  weight16 weights = m_colours.GetWeights();

  m_xsum_wsum = 0.0f;
  for (int i = 0; i < count; ++i) {
    int j = order[i];

    float4 p = float4(unweighted[j].x, unweighted[j].y, unweighted[j].z, 1.0f       );
    float4 w = float4(   weights[j]  ,    weights[j]  ,    weights[j]  , weights[j] );
    float4 x = p * w;

    m_points_weights[i] = x;
    m_xsum_wsum += x;
  }

  return true;
}

void ClusterFit_CCR::Compress(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block, const bool trans,
			        IndexBlockLUT yArr) amp_restricted
{
  /* all or nothing branches, OK, same for all threads */
  bool isBtc1 = (trans);
  if (isBtc1) {
    Compress3(barrier, thread, m_colours, block, yArr);
    if (!m_colours.IsTransparent())
      Compress4(barrier, thread, m_colours, block, yArr);
  }
  else
    Compress4(barrier, thread, m_colours, block, yArr);
}

void ClusterFit_CCR::Compress3(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
			         IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // declare variables
  const int count = m_colours.GetCount();

  const float4 two = 2.0;
  const float4 half_half2 = float4(0.5f, 0.5f, 0.5f, 0.25f);
  const float4 half = 0.5f;
  const float4 grid = float4(31.0f, 63.0f, 31.0f, 0.0f);
  const float4 gridrcp = float4(1.0f/31.0f, 1.0f/63.0f, 1.0f/31.0f, 0.0f);

  // prepare an ordering using the principle axis
  float3 bestline[CVALS];

  bestline[CSTRT] = 0.0f;
  bestline[CSTOP] = m_principle;

  // check all possible clusters and iterate on the total order
  float besterror = m_besterror;
  int bestiteration = 0;
  int besti = 0, bestj = 0;

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0; iterationIndex < m_iterationCount; iterationIndex++) {
    float3 axis = (m_line[CSTOP] - m_line[CSTRT]);

    // stop if a new iteration is an ordering that has already been tried
    if (!ConstructOrdering(barrier, thread, m_colours, axis, iterationIndex) && iterationIndex)
      break;

    // first cluster [0,i) is at the start
    float4 part0 = 0.0f;
    for (int i = 0; i < count; ++i) {
      // second cluster [i,j) is half along
      float4 part1 = (i == 0) ? m_points_weights[0] : 0.0f;
      int jmin = (i == 0) ? 1 : i;
      for (int j = jmin;;) {
	// last cluster [j,count) is at the end
	float4 part2 = m_xsum_wsum - part1 - part0;

	// compute least squares terms directly
	float4 alphax_sum = muladd( part1, half_half2, part0 );
	float4 alpha2_sum; alpha2_sum = alphax_sum.w;

	float4 betax_sum = muladd( part1, half_half2, part2 );
	float4 beta2_sum; beta2_sum = betax_sum.w;

	float4 alphabeta_sum; alphabeta_sum = ( part1*half_half2 ).w;

	// compute the least-squares optimal points
	float4 factor = recip( submul( alphabeta_sum, alphabeta_sum, alpha2_sum*beta2_sum ) );
	float4 a = submul( betax_sum, alphabeta_sum, alphax_sum*beta2_sum )*factor;
	float4 b = submul( alphax_sum, alphabeta_sum, betax_sum*alpha2_sum )*factor;

	// snap floating-point-values to the integer-lattice
	a = saturate( a );
	b = saturate( b );
	a = truncate( muladd( grid, a, half ) )*gridrcp;
	b = truncate( muladd( grid, b, half ) )*gridrcp;

	// compute the error (we skip the constant xxsum)
	float4 e1 = muladd( a*a, alpha2_sum, b*b*beta2_sum );
	float4 e2 = submul( a, alphax_sum, a*b*alphabeta_sum );
	float4 e3 = submul( b, betax_sum, e2 );
	float4 e4 = muladd( two, e3, e1 );

	// apply the metric to the error term
	float4 e5 = e4 * m_metric4;
	float error; error = e5.x + e5.y + e5.z;

	// keep the solution if it wins
	if (less(error, besterror)) {
#ifdef USE_AMP_DEBUG
	  bestline[CSTRT] = a.GetVec3();
	  bestline[CSTOP] = b.GetVec3();
#else
	  bestline[CSTRT] = a.xyz;
	  bestline[CSTOP] = b.xyz;
#endif

	  besti = i;
	  bestj = j;
	  besterror = error;
	  bestiteration = iterationIndex;
	}

	// advance
	if (j == count)
	  break;

	part1 += m_points_weights[j];
	++j;
      }

      // advance
      part0 += m_points_weights[i];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;
  }

  // save the block if necessary
  if (less(besterror, m_besterror)) {
    // save the error
    m_besterror = besterror;

    // remap the indices
    index16 order = m_order[bestiteration];
    for (int m =     0; m < besti; ++m)
      m_matches[1][order[m]] = 0;
    for (int m = besti; m < bestj; ++m)
      m_matches[1][order[m]] = 2;
    for (int m = bestj; m < count; ++m)
      m_matches[1][order[m]] = 1;

    m_line[CSTRT] = bestline[CSTRT];
    m_line[CSTOP] = bestline[CSTOP];

    // save the block
    ColourFit_CCR::Compress3(barrier, thread, m_colours, block, yArr);
  }
}

void ClusterFit_CCR::Compress4(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
			         IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // declare variables
  const int count = m_colours.GetCount();

  const float4 two = 2.0f;
  const float4 onethird_onethird2 = float4( 1.0f/3.0f, 1.0f/3.0f, 1.0f/3.0f, 1.0f/9.0f );
  const float4 twothirds_twothirds2 = float4( 2.0f/3.0f, 2.0f/3.0f, 2.0f/3.0f, 4.0f/9.0f );
  const float4 twonineths = 2.0f/9.0f;
  const float4 half = 0.5f;
  const float4 grid = float4( 31.0f, 63.0f, 31.0f, 0.0f );
  const float4 gridrcp = float4( 1.0f/31.0f, 1.0f/63.0f, 1.0f/31.0f, 0.0f );

  // prepare an ordering using the principle axis
  float3 bestline[CVALS];

  bestline[CSTRT] = 0.0f;
  bestline[CSTOP] = m_principle;

  // check all possible clusters and iterate on the total order
  float besterror = m_besterror;
  int bestiteration = 0;
  int besti = 0, bestj = 0, bestk = 0;

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0; iterationIndex < m_iterationCount; iterationIndex++) {
    float3 axis = (m_line[CSTOP] - m_line[CSTRT]);

    // stop if a new iteration is an ordering that has already been tried
    if (!ConstructOrdering(barrier, thread, m_colours, axis, iterationIndex) && iterationIndex)
      break;

    // first cluster [0,i) is at the start
    float4 part0 = 0.0f;
    for (int i = 0; i < count; ++i) {
      // second cluster [i,j) is one third along
      float4 part1 = 0.0f;
      for (int j = i;;) {
	// third cluster [j,k) is two thirds along
	float4 part2 = (j == 0) ? m_points_weights[0] : 0.0f;
	int kmin = (j == 0) ? 1 : j;
	for (int k = kmin;;) {
	  // last cluster [k,count) is at the end
	  float4 part3 = m_xsum_wsum - part2 - part1 - part0;

	  // compute least squares terms directly
	  float4 alphax_sum = muladd( part2, onethird_onethird2, muladd( part1, twothirds_twothirds2, part0 ) );
	  float4 alpha2_sum; alpha2_sum = alphax_sum.w;

	  float4 betax_sum = muladd( part1, onethird_onethird2, muladd( part2, twothirds_twothirds2, part3 ) );
	  float4 beta2_sum; beta2_sum = betax_sum.w;

	  float4 alphabeta_sum; alphabeta_sum = twonineths * ( part1 + part2 ).w;

	  // compute the least-squares optimal points
	  float4 factor = recip( submul( alphabeta_sum, alphabeta_sum, alpha2_sum*beta2_sum ) );
	  float4 a = submul( betax_sum, alphabeta_sum, alphax_sum*beta2_sum )*factor;
	  float4 b = submul( alphax_sum, alphabeta_sum, betax_sum*alpha2_sum )*factor;

	  // snap floating-point-values to the integer-lattice
	  a = saturate( a );
	  b = saturate( b );
	  a = truncate( muladd( grid, a, half ) )*gridrcp;
	  b = truncate( muladd( grid, b, half ) )*gridrcp;

	  // compute the error (we skip the constant xxsum)
	  float4 e1 = muladd( a*a, alpha2_sum, b*b*beta2_sum );
	  float4 e2 = submul( a, alphax_sum, a*b*alphabeta_sum );
	  float4 e3 = submul( b, betax_sum, e2 );
	  float4 e4 = muladd( two, e3, e1 );

	  // apply the metric to the error term
	  float4 e5 = e4 * m_metric4;
	  float error; error = e5.x + e5.y + e5.z;

	  // keep the solution if it wins
	  if (less(error, besterror)) {
#ifdef USE_AMP_DEBUG
	    bestline[CSTRT] = a.GetVec3();
	    bestline[CSTOP] = b.GetVec3();
#else
	    bestline[CSTRT] = a.xyz;
	    bestline[CSTOP] = b.xyz;
#endif

	    besterror = error;
	    besti = i;
	    bestj = j;
	    bestk = k;
	    bestiteration = iterationIndex;
	  }

	  // advance
	  if (k == count)
	    break;

	  part2 += m_points_weights[k];
	  ++k;
	}

	// advance
	if (j == count)
	  break;

	part1 += m_points_weights[j];
	++j;
      }

      // advance
      part0 += m_points_weights[i];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;
  }

  // save the block if necessary
  if (less(besterror, m_besterror)) {
    // save the error
    m_besterror = besterror;

    // remap the indices
    index16 order = m_order[bestiteration];
    for (int m =     0; m < besti; ++m)
      m_matches[1][order[m]] = 0;
    for (int m = besti; m < bestj; ++m)
      m_matches[1][order[m]] = 2;
    for (int m = bestj; m < bestk; ++m)
      m_matches[1][order[m]] = 3;
    for (int m = bestk; m < count; ++m)
      m_matches[1][order[m]] = 1;

    m_line[CSTRT] = bestline[CSTRT];
    m_line[CSTOP] = bestline[CSTOP];

    // save the block
    ColourFit_CCR::Compress4(barrier, thread, m_colours, block, yArr);
  }
}
#endif

} // namespace squish
