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

#include "paletteclusterfit.h"
#include "paletteset.h"
#include "paletteblock.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
PaletteClusterFit::PaletteClusterFit(PaletteSet const* palette, int flags, int swap, int shared)
  : PaletteSingleMatch(palette, flags, swap, shared)
  ,  PaletteChannelFit(palette, flags, swap, shared)
  ,         PaletteFit(palette, flags, swap, shared)
{
  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;
  bool const trns = m_palette->IsMergedAlpha() && m_palette->IsTransparent();

  assume((isets >  0) && (isets <= 3));
  assume((asets >= 0) && (asets <= 3));
  assume(((isets    +    asets) <= 3));

  // set the iteration count
  m_iterationCount = (m_flags & kColourIterativeClusterFits) / kColourClusterFit;

  // loop over all sets
  for (int s = 0; s < (isets + asets); s++) {
    // cache some values
    bool const unweighted = m_palette->IsUnweighted(s);
    int const count = m_palette->GetCount(s);
    Vec4 const* values = m_palette->GetPoints(s);
    Scr4 const* weights = m_palette->GetWeights(s);

    // we don't do this for sparse sets
    if (count != 1) {
      Vec4 centroid;

      // combined alpha
      if (trns) {
        Sym4x4 covariance;

        // get the covariance matrix
        if (unweighted)
	  ComputeWeightedCovariance4(covariance, centroid, count, values, m_metric[s]);
        else
	  ComputeWeightedCovariance4(covariance, centroid, count, values, m_metric[s], weights);

	// compute the principle component
	GetPrincipleComponent(covariance, m_principle[s]);
      }
      // no or separate alpha
      else {
        Sym3x3 covariance;

        // get the covariance matrix
        if (unweighted)
	  ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric[s]);
        else
	  ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric[s], weights);

	// compute the principle component
	GetPrincipleComponent(covariance, m_principle[s]);
      }

      // we have tables for this
      m_optimizable[s] = unweighted & (count == 16);
    }
  }

  // set all orders to zero to allow sizeless compare
  memset(m_order, 0, sizeof(m_order));
}

bool PaletteClusterFit::ConstructOrdering(Vec4 const& axis, int iteration, int set)
{
  // cache some values
  int const count = m_palette->GetCount(set);
  Vec4 const* values = m_palette->GetPoints(set);

  // build the list of dot products
  float dps[16];
  u8* order = (u8*)m_order[set] + 16 * iteration;
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
  Col4 curr; LoadAligned(curr, m_order[set] + 16 * iteration);
  for (int it = 0; it < iteration; ++it) {
    Col4 prev; LoadAligned(prev, m_order[set] + 16 * it);

    if (CompareAllEqualTo(curr, prev))
      return false;
  }

  // copy the ordering and weight all the points
  Vec4 const* unweighted = m_palette->GetPoints(set);
  Scr4 const* weights = m_palette->GetWeights(set);
  bool const trns = m_palette->IsMergedAlpha();

  if (trns) {
    m_xsum_wsum  [(set << 1) + 0] = VEC4_CONST(0.0f);
    m_xsum_wsum  [(set << 1) + 1] = VEC4_CONST(0.0f);
//  m_xxsum_wwsum[(set << 1) + 0] = VEC4_CONST(0.0f);
//  m_xxsum_wwsum[(set << 1) + 1] = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {
      int j = order[i];

      Vec4 p = unweighted[j];
      Scr4 w(weights[j]);
      Vec4 x = p * w;

      m_points_weights[set][(i << 1) + 0] = x;
      m_points_weights[set][(i << 1) + 1] = w;

      m_xsum_wsum     [(set << 1) + 0]   += x;
      m_xsum_wsum     [(set << 1) + 1]   += w;
//    m_xxsum_wwsum   [(set << 1) + 0]   += x * x;
//    m_xxsum_wwsum   [(set << 1) + 1]   += w * w;
    }
  }
  else {
    m_xsum_wsum  [set] = VEC4_CONST(0.0f);
//  m_xxsum_wwsum[set] = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {
      int j = order[i];

      Vec4 p = TransferW(unweighted[j], Vec4(1.0f));
      Scr4 w(weights[j]);
      Vec4 x = p * w;

      m_points_weights[set][i] = x;
      m_xsum_wsum     [set]   += x;
//    m_xxsum_wwsum   [set]   += x * x;
    }
  }

  return true;
}

Scr4 PaletteClusterFit::ClusterSearch4(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb)
{
  /*
  Vec4 const weight1(21.0f / 64.0f, 21.0f / 64.0f, 21.0f / 64.0f,  441.0f / 4096.0f);
  Vec4 const weight2(43.0f / 64.0f, 43.0f / 64.0f, 43.0f / 64.0f, 1849.0f / 4096.0f);
  Vec4 const twonineths                               = VEC4_CONST(882.0f / 4096.0f);
  */
  Vec4 const weight1(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 9.0f);
  Vec4 const weight2(2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f, 4.0f / 9.0f);
  Vec4 const twonineths                        = VEC4_CONST(2.0f / 9.0f);
  Vec4 const fivenineths                       = VEC4_CONST(5.0f / 9.0f);
  Vec4 const convnineths                       = VEC4_CONST(2.5f / 1.0f);

  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const half = VEC4_CONST(0.5f);

  assume((count > 0) && (count <= 16));

  // match each point to the closest code
  int besti = 0, bestj = 0, bestk = 0;
  int bestiteration = 0;
  Vec4 beststart = VEC4_CONST(0.0f);
  Vec4 bestend   = VEC4_CONST(0.0f);

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle[set], 0, set);

  // check all possible clusters and iterate on the total order
  Scr4 besterror = Scr4(FLT_MAX);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_wsum = m_xsum_wsum[set];

    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

    // second cluster [i,j) is one third along
    Vec4 part1 = VEC4_CONST(0.0f);
    for (int j = i;;) {

    // third cluster [j,k) is two thirds along
    Vec4 part2 = (j == 0) ? m_points_weights[set][0] : VEC4_CONST(0.0f);
    int kmin = (j == 0) ? 1 : j;
    for (int k = kmin;;) {

	  // last cluster [k,count) is at the end
	  Vec4 part3 = xsum_wsum - part2 - part1 - part0;
//	  Vec4 partI = xsum_wsum - part2 - part1;

	  // compute least squares terms directly
	  Vec4 const alphax_sum = MultiplyAdd(part2, weight1, MultiplyAdd(part1, weight2, part0));
	  Vec4 const  betax_sum = MultiplyAdd(part1, weight1, MultiplyAdd(part2, weight2, part3));

	  Vec4 const alpha2_sum = alphax_sum.SplatW();
	  Vec4 const  beta2_sum =  betax_sum.SplatW();

	  Vec4 const alphabeta_sum = twonineths * (part1 + part2).SplatW();

	  // equivalent
//	  Vec4 const  betac_sum =                xsum_wsum                - alphax_sum;
//	  Vec4 const beta2a_sum = (fivenineths * (part1 + part2) + partI) - alpha2_sum;
//	  Vec4 const beta2b_sum = (convnineths *  alphabeta_sum  + partI) - alpha2_sum;

	  //   alpha.x = p3 * 0/3  + p2 * 1/3  + p1 * 2/3  + p0 * 3/3
	  //    beta.x = p3 * 3/3  + p2 * 2/3  + p1 * 1/3  + p0 * 0/3
	  //   gamma.x = p3 * 3/3  + p2 * 3/3  + p1 * 3/3  + p0 * 3/3 = alpha.x + beta.x = xsum_wsum;
	  //
	  //   alpha.w = p3 * 0/3  + p2 * 1/9  + p1 * 4/9  + p0 * 9/9
	  //    beta.w = p3 * 3/3  + p2 * 2/3  + p1 * 1/3  + p0 * 0/3
	  //
	  //   alpha.w =      0/3² +      1/3² +      2/3² +      ?.?²
	  //    beta.w =      ?.?² +      2/3² +      1/3² +      0/3²
	  // alphabeta = p3 * 0/-- + p2 * 2/9  + p1 * 2/9  + p0 * 0/--

	  // compute the least-squares optimal points
	  Vec4 factor = Reciprocal(NegativeMultiplySubtract(alphabeta_sum, alphabeta_sum, alpha2_sum * beta2_sum));
	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_sum, alphax_sum *  beta2_sum) * factor;
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_sum,  betax_sum * alpha2_sum) * factor;

	  // factor = 1.0 / (alpha.w * beta.w - alphabeta * alphabeta)
	  // a = (alpha.x * beta.w -  beta.x * alphabeta) * factor
	  // b = (beta.x * alpha.w - alpha.x * alphabeta) * factor
	  // a = (alpha.x * beta.w -  beta.x * alphabeta) / (alpha.w * beta.w - alphabeta * alphabeta)
	  // b = (beta.x * alpha.w - alpha.x * alphabeta) / (alpha.w * beta.w - alphabeta * alphabeta)

	  // snap floating-point-values to the integer-lattice
	  a = q.SnapToLattice(a, sb, 1 << SBSTART);
	  b = q.SnapToLattice(b, sb, 1 << SBEND);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // e1 = a * a * alpha2_sum + b * b * beta2_sum
	  // e2 = a * b * alphabeta_sum - a * alphax_sum
	  // e3 = e2 - b * betax_sum
	  // e4 = 2 * e3 + e1
	  //
	  // 9 muls, 4 adds

	  // e =
	  //     2 * a * b * alphabeta_sum -
	  //     2 * a * alphax_sum -
	  //     a * a * alpha2_sum +
	  //     2 * b * betax_sum +
	  //     b * b * beta2_sum

	  // e =
	  //     a * b * 2 * alphabeta_sum -
	  //     a * (a * alpha2_sum + 2 * alphax_sum) +
	  //     b * (b *  beta2_sum + 2 *  betax_sum)
	  //
	  // 9 muls, 4 adds

	  // e / 2 =
	  //     a * b * alphabeta_sum -
	  //     a * (a * alpha2_sum / 2 + alphax_sum) +
	  //     b * (b *  beta2_sum / 2 +  betax_sum)
	  //
	  // 8 muls, 4 adds

	  // apply the metric to the error term
	  Scr4 eS = Dot(e4, metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (besterror > eS) {
	    besterror = eS;

	    beststart = a;
	    bestend   = b;
	    bestiteration = iterationIndex;

	    besti = i,
	    bestj = j,
	    bestk = k;
	  }

      // advance
      if (k == count) break;
      part2 += m_points_weights[set][k]; ++k; }

      // advance
      if (j == count) break;
      part1 += m_points_weights[set][j]; ++j; }

      // advance
      part0 += m_points_weights[set][i];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;

    // advance if possible
    ++iterationIndex;
    if (iterationIndex == m_iterationCount)
      break;

    // stop if a new iteration is an ordering that has already been tried
    Vec4 axis = KillW(bestend - beststart);
    if (!ConstructOrdering(axis, iterationIndex, set))
      break;
  }

  // assign closest points
  u8 const* order = (u8*)m_order[set] + 16 * bestiteration;

  for (int m =     0; m < besti; ++m)
    closest[set][order[m]] = 0;
  for (int m = besti; m < bestj; ++m)
    closest[set][order[m]] = 1;
  for (int m = bestj; m < bestk; ++m)
    closest[set][order[m]] = 2;
  for (int m = bestk; m < count; ++m)
    closest[set][order[m]] = 3;

  // copy rgb into the common start/end definition
  m_start[set] = beststart;
  m_end  [set] = bestend;

  return besterror;
}

Scr4 PaletteClusterFit::ClusterSearch4Alpha(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb)
{
  /*
  Vec4 const weight1(21.0f / 64.0f, 21.0f / 64.0f, 21.0f / 64.0f,  441.0f / 4096.0f);
  Vec4 const weight2(43.0f / 64.0f, 43.0f / 64.0f, 43.0f / 64.0f, 1849.0f / 4096.0f);
  Vec4 const twonineths                               = VEC4_CONST(882.0f / 4096.0f);
  */
  Vec4 const weight1x(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f);
  Vec4 const weight2x(2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f);
  Vec4 const weight1w(1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f);
  Vec4 const weight2w(4.0f / 9.0f, 4.0f / 9.0f, 4.0f / 9.0f, 4.0f / 9.0f);
  Vec4 const twonineths                         = VEC4_CONST(2.0f / 9.0f);
  Vec4 const fivenineths                        = VEC4_CONST(5.0f / 9.0f);
  Vec4 const convnineths                        = VEC4_CONST(2.5f / 1.0f);

  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const half = VEC4_CONST(0.5f);

  assume((count > 0) && (count <= 16));

  // match each point to the closest code
  int besti = 0, bestj = 0, bestk = 0;
  int bestiteration = 0;
  Vec4 beststart = VEC4_CONST(0.0f);
  Vec4 bestend   = VEC4_CONST(0.0f);

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle[set], 0, set);

  // check all possible clusters and iterate on the total order
  Scr4 besterror = Scr4(FLT_MAX);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_xsum = m_xsum_wsum[(set << 1) + 0];
    Vec4 const wsum_wsum = m_xsum_wsum[(set << 1) + 1];

    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    Vec4 wgts0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

    // second cluster [i,j) is one third along
    Vec4 part1 = VEC4_CONST(0.0f);
    Vec4 wgts1 = VEC4_CONST(0.0f);
    for (int j = i;;) {

    // third cluster [j,k) is two thirds along
    Vec4 part2 = (j == 0) ? m_points_weights[set][0] : VEC4_CONST(0.0f);
    Vec4 wgts2 = (j == 0) ? m_points_weights[set][1] : VEC4_CONST(0.0f);
    int kmin = (j == 0) ? 1 : j;
    for (int k = kmin;;) {

	  // last cluster [k,count) is at the end
	  Vec4 part3 = xsum_xsum - part2 - part1 - part0;
	  Vec4 wgts3 = wsum_wsum - wgts2 - wgts1 - wgts0;
//	  Vec4 partI = xsum_xsum - part2 - part1;
//	  Vec4 wgtsI = wsum_wsum - wgts2 - wgts1;

	  // compute least squares terms directly
	  Vec4 const alphax_sum = MultiplyAdd(part2, weight1x, MultiplyAdd(part1, weight2x, part0));
	  Vec4 const  betax_sum = MultiplyAdd(part1, weight1x, MultiplyAdd(part2, weight2x, part3));

	  Vec4 const alpha2_sum = MultiplyAdd(wgts2, weight1w, MultiplyAdd(wgts1, weight2w, wgts0));
	  Vec4 const  beta2_sum = MultiplyAdd(wgts1, weight1w, MultiplyAdd(wgts2, weight2w, wgts3));

	  Vec4 const alphabeta_sum = twonineths * (wgts1 + wgts2);

	  // compute the least-squares optimal points
	  Vec4 factor = Reciprocal(NegativeMultiplySubtract(alphabeta_sum, alphabeta_sum, alpha2_sum * beta2_sum));
	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_sum, alphax_sum *  beta2_sum) * factor;
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_sum,  betax_sum * alpha2_sum) * factor;

	  // snap floating-point-values to the integer-lattice
	  a = q.SnapToLattice(a, sb, 1 << SBSTART);
	  b = q.SnapToLattice(b, sb, 1 << SBEND);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = Dot(e4, metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (besterror > eS) {
	    besterror = eS;

	    beststart = a;
	    bestend   = b;
	    bestiteration = iterationIndex;

	    besti = i,
	    bestj = j,
	    bestk = k;
	  }

      // advance
      if (k == count) break;
      part2 += m_points_weights[set][(k << 1) + 0];
      wgts2 += m_points_weights[set][(k << 1) + 1]; ++k; }

      // advance
      if (j == count) break;
      part1 += m_points_weights[set][(j << 1) + 0];
      wgts1 += m_points_weights[set][(j << 1) + 1]; ++j; }

      // advance
      part0 += m_points_weights[set][(i << 1) + 0];
      wgts0 += m_points_weights[set][(i << 1) + 1];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;

    // advance if possible
    ++iterationIndex;
    if (iterationIndex == m_iterationCount)
      break;

    // stop if a new iteration is an ordering that has already been tried
    Vec4 axis = KillW(bestend - beststart);
    if (!ConstructOrdering(axis, iterationIndex, set))
      break;
  }

  // assign closest points
  u8 const* order = (u8*)m_order[set] + 16 * bestiteration;

  for (int m =     0; m < besti; ++m)
    closest[set][order[m]] = 0;
  for (int m = besti; m < bestj; ++m)
    closest[set][order[m]] = 1;
  for (int m = bestj; m < bestk; ++m)
    closest[set][order[m]] = 2;
  for (int m = bestk; m < count; ++m)
    closest[set][order[m]] = 3;

  // copy rgb into the common start/end definition
  m_start[set] = beststart;
  m_end  [set] = bestend;

  return besterror;
}

#include "paletteclusterfit.inl"

Scr4 PaletteClusterFit::ClusterSearch4Constant(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb)
{
  /*
  Vec4 const weight1(21.0f / 64.0f, 21.0f / 64.0f, 21.0f / 64.0f,  441.0f / 4096.0f);
  Vec4 const weight2(43.0f / 64.0f, 43.0f / 64.0f, 43.0f / 64.0f, 1849.0f / 4096.0f);
  Vec4 const twonineths                               = VEC4_CONST(882.0f / 4096.0f);
  */
  Vec4 const weight1(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 9.0f);
  Vec4 const weight2(2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f, 4.0f / 9.0f);
  Vec4 const twonineths                        = VEC4_CONST(2.0f / 9.0f);

  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const half = VEC4_CONST(0.5f);

  assume((count > 0) && (count <= 16));

  // match each point to the closest code
  int besti = 0, bestj = 0, bestk = 0;
  int bestiteration = 0;
  Vec4 beststart = VEC4_CONST(0.0f);
  Vec4 bestend   = VEC4_CONST(0.0f);

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle[set], 0, set);

  // check all possible clusters and iterate on the total order
  Scr4 besterror = Scr4(FLT_MAX);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_wsum = m_xsum_wsum[set];

    // constants if weights == 1
    Vec4 alphabeta_dltas = *((Vec4 *)part2delta[0]);
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
    Vec4 part2 = (j == 0) ? m_points_weights[set][0] : VEC4_CONST(0.0f);
    Vec4 alphabeta_val = *alphabeta_inits++;
    int kmin = (j == 0) ? 1 : j;
    for (int k = kmin;;) {
	  // TODO: the inner alphabeta_sum seems always to be the same sequence
	  Vec4 alphabeta_factor = alphabeta_val * Vec4(*alphabeta_factors++);

	  // compute least squares terms directly
	  Vec4 const alphax_sum =   MultiplyAdd(part2, weight1, MultiplyAdd(part1, weight2, part0));
	  Vec4 const  betax_sum = /*MultiplyAdd(part1, weight1, MultiplyAdd(part2, weight2, part3))*/ xsum_wsum - alphax_sum;

	  Vec4 const    alpha2_sum = alphabeta_val.SplatX();
	  Vec4 const     beta2_sum = alphabeta_val.SplatY();
	  Vec4 const alphabeta_sum = alphabeta_val.SplatZ();

	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_factor.SplatZ(), alphax_sum * alphabeta_factor.SplatY());
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_factor.SplatZ(),  betax_sum * alphabeta_factor.SplatX());

#if 0
	  // last cluster [k,count) is at the end
	  Vec4 part3 = xsum_wsum - part2 - part1 - part0;

	  // compute least squares terms directly
	  Vec4 const _alphax_sum = MultiplyAdd(part2, weight1, MultiplyAdd(part1, weight2, part0));
	  Vec4 const  _betax_sum = MultiplyAdd(part1, weight1, MultiplyAdd(part2, weight2, part3));
//	  Vec4 const  _betac_sum = xsum_wsum - _alphax_sum;

	  Vec4 const _alpha2_sum = _alphax_sum.SplatW();
	  Vec4 const  _beta2_sum =  _betax_sum.SplatW();

	  Vec4 const _alphabeta_sum = twonineths * (part1 + part2).SplatW();

	  //   alpha.x = p3 * 0/3  + p2 * 1/3  + p1 * 2/3  + p0 * 3/3
	  //    beta.x = p3 * 3/3  + p2 * 2/3  + p1 * 1/3  + p0 * 0/3
	  //   gamma.x = p3 * 3/3  + p2 * 3/3  + p1 * 3/3  + p0 * 3/3 = alpha.x + beta.x = xsum_wsum;
	  //
	  //   alpha.w =      0/3² +      1/3² +      2/3² +      ?.?²
	  //    beta.w =      ?.?² +      2/3² +      1/3² +      0/3²
	  // alphabeta = p3 * 0/-- + p2 * 2/9  + p1 * 2/9  + p0 * 0/--

	  // compute the least-squares optimal points
	  Vec4 _factor = Reciprocal(NegativeMultiplySubtract(_alphabeta_sum, _alphabeta_sum, _alpha2_sum * _beta2_sum));
	  Vec4 _a = NegativeMultiplySubtract( _betax_sum, _alphabeta_sum, _alphax_sum *  _beta2_sum) * _factor;
	  Vec4 _b = NegativeMultiplySubtract(_alphax_sum, _alphabeta_sum,  _betax_sum * _alpha2_sum) * _factor;

#define	limit 1e-5
	  assert(fabs(_alpha2_sum.W() - alpha2_sum.X()) < limit);
	  assert(fabs(_beta2_sum.W() - beta2_sum.X()) < limit);
	  assert(fabs(_alphabeta_sum.W() - alphabeta_sum.X()) < limit);

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
	  a = q.SnapToLattice(a, sb, 1 << SBSTART);
	  b = q.SnapToLattice(b, sb, 1 << SBEND);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b, betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = Dot(e4, metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (besterror > eS) {
	    besterror = eS;

	    beststart = a;
	    bestend   = b;
	    bestiteration = iterationIndex;

	    besti = i,
	    bestj = j,
	    bestk = k;
	  }

      alphabeta_val += alphabeta_dltas;

      // advance
      if (k == count) break;
      part2 += m_points_weights[set][k]; ++k; }

      // advance
      if (j == count) break;
      part1 += m_points_weights[set][j]; ++j; }

      // advance
      part0 += m_points_weights[set][i];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;

    // advance if possible
    ++iterationIndex;
    if (iterationIndex == m_iterationCount)
      break;

    // stop if a new iteration is an ordering that has already been tried
    Vec4 axis = KillW(bestend - beststart);
    if (!ConstructOrdering(axis, iterationIndex, set))
      break;
  }

  // assign closest points
  u8 const* order = (u8*)m_order[set] + 16 * bestiteration;

  for (int m =     0; m < besti; ++m)
    closest[set][order[m]] = 0;
  for (int m = besti; m < bestj; ++m)
    closest[set][order[m]] = 1;
  for (int m = bestj; m < bestk; ++m)
    closest[set][order[m]] = 2;
  for (int m = bestk; m < count; ++m)
    closest[set][order[m]] = 3;

  // copy rgb into the common start/end definition
  m_start[set] = beststart;
  m_end  [set] = bestend;

  return besterror;
}

Scr4 PaletteClusterFit::ClusterSearch8(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb)
{
  /*
  Vec4 const weight1( 9.0f / 64.0f,  9.0f / 64.0f,  9.0f / 64.0f,   81.0f / 4096.0f);
  Vec4 const weight2(18.0f / 64.0f, 18.0f / 64.0f, 18.0f / 64.0f,  324.0f / 4096.0f);
  Vec4 const weight3(27.0f / 64.0f, 27.0f / 64.0f, 27.0f / 64.0f,  729.0f / 4096.0f);
  Vec4 const weight4(37.0f / 64.0f, 37.0f / 64.0f, 37.0f / 64.0f, 1369.0f / 4096.0f);
  Vec4 const weight5(46.0f / 64.0f, 46.0f / 64.0f, 46.0f / 64.0f, 2116.0f / 4096.0f);
  Vec4 const weight6(55.0f / 64.0f, 55.0f / 64.0f, 55.0f / 64.0f, 3025.0f / 4096.0f);
  Vec4 const twonineths                               = VEC4_CONST(486.0f / 4096.0f);
  */
  Vec4 const weight1(1.0f / 7.0f, 1.0f / 7.0f, 1.0f / 7.0f,  1.0f / 49.0f);
  Vec4 const weight2(2.0f / 7.0f, 2.0f / 7.0f, 2.0f / 7.0f,  4.0f / 49.0f);
  Vec4 const weight3(3.0f / 7.0f, 3.0f / 7.0f, 3.0f / 7.0f,  9.0f / 49.0f);
  Vec4 const weight4(4.0f / 7.0f, 4.0f / 7.0f, 4.0f / 7.0f, 16.0f / 49.0f);
  Vec4 const weight5(5.0f / 7.0f, 5.0f / 7.0f, 5.0f / 7.0f, 25.0f / 49.0f);
  Vec4 const weight6(6.0f / 7.0f, 6.0f / 7.0f, 6.0f / 7.0f, 36.0f / 49.0f);
  Vec4 const twonineths                        = VEC4_CONST( 6.0f / 49.0f);
  Vec4 const threenineths                      = VEC4_CONST(10.0f / 49.0f);
  Vec4 const fournineths                       = VEC4_CONST(12.0f / 49.0f);
  Vec4 const twonineths_                       = VEC4_CONST(37.0f / 49.0f);
  Vec4 const threenineths_                     = VEC4_CONST(29.0f / 49.0f);
  Vec4 const fournineths_                      = VEC4_CONST(25.0f / 49.0f);
  Vec4 const twoninethsc                       = VEC4_CONST((37.0f / 49.0f) / ( 6.0f / 49.0f));
  Vec4 const threeninethsc                     = VEC4_CONST((29.0f / 49.0f) / (10.0f / 49.0f));
  Vec4 const fourninethsc                      = VEC4_CONST((25.0f / 49.0f) / (12.0f / 49.0f));
  // onenineths   = sqr(1.0 / distance) = sqr(1.0 / 7.0) = (1.0 / 49.0)
  // twonineths   = weight1 * weight6 = (6.0 / 49.0)
  // threenineths = weight2 * weight5 = (10.0 / 49.0)
  // fournineths  = weight3 * weight4 = (12.0 / 49.0)

  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const half = VEC4_CONST(0.5f);
  Vec4 const zero = VEC4_CONST(0.0f);

  assume((count > 0) && (count <= 16));

  // match each point to the closest code
  int besti = 0, bestj = 0, bestk = 0;
  int bestl = 0, bestm = 0, bestn = 0;
  int besto = 0;

  int bestiteration = 0;
  Vec4 beststart    = VEC4_CONST(0.0f);
  Vec4 bestend      = VEC4_CONST(0.0f);

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle[set], 0, set);

  // check all possible clusters and iterate on the total order
  Scr4 besterror = Scr4(FLT_MAX);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_wsum = m_xsum_wsum[set];

#if 0
    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

    // second cluster [i,j) is one seventh along
    Vec4 part1 = VEC4_CONST(0.0f);
    for (int j = i;;) {
	  Vec4 sub1  =   part0 + part1;
	  Vec4 mul1a = weight6 * part1;
	  Vec4 mul1b = weight1 * part1;

    // third cluster [j,k) is one seventh along
    Vec4 part2 = VEC4_CONST(0.0f);
    for (int k = j;;) {
	  Vec4 sub2  =    sub1 + part2;
	  Vec4 mul2a = weight5 * part2;
	  Vec4 mul2b = weight2 * part2;
	  Vec4 add2a =   mul1a + mul2a;
	  Vec4 add2b =   mul1b + mul2b;

    // fourth cluster [k,l) is one seventh along
    Vec4 part3 = VEC4_CONST(0.0f);
    for (int l = k;;) {
	  Vec4 sub3  =    sub2 + part3;
	  Vec4 mul3a = weight4 * part3;
	  Vec4 mul3b = weight3 * part3;
	  Vec4 add3a =   add2a + mul3a;
	  Vec4 add3b =   add2b + mul3b;

    // fifth cluster [l,m) is one seventh along
    Vec4 part4 = VEC4_CONST(0.0f);
    for (int m = l;;) {
	  Vec4 sub4  =    sub3 + part4;
	  Vec4 mul4a = weight3 * part4;
	  Vec4 mul4b = weight4 * part4;
	  Vec4 add4a =   add3a + mul4a;
	  Vec4 add4b =   add3b + mul4b;

	  Vec4 sum34 = fournineths  * (part3 + part4).SplatW();

    // sixths cluster [m,n) is one seventh along
    Vec4 part5 = VEC4_CONST(0.0f);
    for (int n = m;;) {
	  Vec4 sub5  =    sub4 + part5;
	  Vec4 mul5a = weight2 * part5;
	  Vec4 mul5b = weight5 * part5;
	  Vec4 add5a =   add4a + mul5a;
	  Vec4 add5b =   add4b + mul5b;

	  Vec4 sum25 = threenineths * (part2 + part5).SplatW();

    // seventh cluster [n,o) is six seventh along
    Vec4 part6 = (n == 0) ? m_points_weights[set][0] : VEC4_CONST(0.0f);
    int omin = (n == 0) ? 1 : n;
    for (int o = omin;;) {
	  Vec4 sub6  =    sub5 + part6;
	  Vec4 mul6a = weight1 * part6;
	  Vec4 mul6b = weight6 * part6;
	  Vec4 add6a =   add5a + mul6a;
	  Vec4 add6b =   add5b + mul6b;

	  Vec4 sum16 = twonineths   * (part1 + part6).SplatW();

	  // last cluster [k,count) is at the end
	  Vec4 part7 = xsum_wsum - part6 - part5 - part4 - part3 - part2 - part1 - part0;
	  Vec4 partX = xsum_wsum - sub6;

	  Vec4 diff = SummedAbsoluteDifference(part7, partX);
	  if (CompareFirstGreaterThan(diff, Vec4(2.4e-6f))) {
	  	int bugout = 1; bugout = 2;
	  }

	  Vec4 mul0a =         part0;
	  Vec4 mul7b =         part7;
	  Vec4 add0a = add6a + mul0a;
	  Vec4 add7b = add6b + mul7b;

	  // compute least squares terms directly
	  Vec4 const alphax_sum =
	  	MultiplyAdd(part6, weight1,
	  	MultiplyAdd(part5, weight2,
	  	MultiplyAdd(part4, weight3,
	  	MultiplyAdd(part3, weight4,
	  	MultiplyAdd(part2, weight5,
	  	MultiplyAdd(part1, weight6,
	  	            part0))))));
	  Vec4 const alphax_chk1 =
	        part0 + mul1a + mul2a + mul3a + mul4a + mul5a + mul6a;
	  Vec4 const alphax_chk2 =
	        add0a;

	  diff = SummedAbsoluteDifference(alphax_sum, alphax_chk2);
	  if (CompareFirstGreaterThan(diff, Vec4(3e-6f))) {
	  	int bugout = 1; bugout = 2;
	  }

	  Vec4 const  betax_sum =
	  	MultiplyAdd(part1, weight1,
	  	MultiplyAdd(part2, weight2,
	  	MultiplyAdd(part3, weight3,
	  	MultiplyAdd(part4, weight4,
	  	MultiplyAdd(part5, weight5,
	  	MultiplyAdd(part6, weight6,
	  	            part7))))));
	  Vec4 const betax_chk1 =
	        part7 + mul6b + mul5b + mul4b + mul3b + mul2b + mul1b;
	  Vec4 const betax_chk2 =
	        add7b;

	  diff = SummedAbsoluteDifference(betax_sum, betax_chk2);
	  if (CompareFirstGreaterThan(diff, Vec4(3e-6f))) {
	  	int bugout = 1; bugout = 2;
	  }

	  Vec4 const alpha2_sum = alphax_sum.SplatW();
	  Vec4 const  beta2_sum =  betax_sum.SplatW();

	  Vec4 const alphabeta_sum =
	  	twonineths   * (part1 + part6).SplatW() +
	  	threenineths * (part2 + part5).SplatW() +
	  	fournineths  * (part3 + part4).SplatW();
	  Vec4 const alphabeta_chk =
	        sum16 + sum25 + sum34;

	  diff = SummedAbsoluteDifference(alphabeta_sum, alphabeta_chk);
	  if (CompareFirstGreaterThan(diff, Vec4(0.0f))) {
	  	int bugout = 1; bugout = 2;
	  }

	  //  alphax.t = p7 * 0/7  + p6 * 1/7  + p5 *  2/7  + p4 *  3/7  + p3 *  4/7  + p2 *  5/7  + p1 * 6/7  + p0 * 7/7
	  //   betax.t = p7 * 7/7  + p6 * 6/7  + p5 *  5/7  + p4 *  4/7  + p3 *  3/7  + p2 *  2/7  + p1 * 1/7  + p0 * 0/7
	  //   gamma.x = p7 * 7/7  + p6 * 7/7  + p5 *  7/7  + p4 *  7/7  + p3 *  7/7  + p2 *  7/7  + p1 * 7/7  + p0 * 7/7 = alpha.x + beta.x = xsum_wsum;
	  //
	  //  alphax.w =      0/7² +      1/7² +       2/7² +       3.7² +       4/7² +       5/7² +      6/7² +      ?/?²
	  //   betax.w =      ?.?² +      6/7² +       5/7² +       4/7² +       3/7² +       2/7² +      1/7² +      0/7²
	  // alphabeta = p7 * 0/-- + p6 * 6/49 + p5 * 10/49 + p4 * 12/49 + p3 * 12/49 + p2 * 10/49 + p1 * 6/49 + p0 * 0/--

	  // compute the least-squares optimal points
	  Vec4 factor = Reciprocal(NegativeMultiplySubtract(alphabeta_sum, alphabeta_sum, alpha2_sum * beta2_sum));
	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_sum, alphax_sum *  beta2_sum) * factor;
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_sum,  betax_sum * alpha2_sum) * factor;

	  // factor = 1.0 / (alpha.w * beta.w - alphabeta * alphabeta)
	  // a = (alpha.x * beta.w -  beta.x * alphabeta) * factor
	  // b = (beta.x * alpha.w - alpha.x * alphabeta) * factor
	  // a = (alpha.x * beta.w -  beta.x * alphabeta) / (alpha.w * beta.w - alphabeta * alphabeta)
	  // b = (beta.x * alpha.w - alpha.x * alphabeta) / (alpha.w * beta.w - alphabeta * alphabeta)

	  // snap floating-point-values to the integer-lattice
	  a = q.SnapToLattice(a, sb, 1 << SBSTART);
	  b = q.SnapToLattice(b, sb, 1 << SBEND);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b,  betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Vec4 eS = Dot(e4, metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (besterror > eS) {
	    besterror = eS;

	    beststart = a;
	    bestend   = b;
	    bestiteration = iterationIndex;

	    besti = i,
	    bestj = j,
	    bestk = k,
	    bestl = l,
	    bestm = m,
	    bestn = n,
	    besto = o;
	  }

      // advance
      if (o == count) break;
      part6 += m_points_weights[set][o]; ++o; }

      // advance
      if (n == count) break;
      part5 += m_points_weights[set][n]; ++n; }

      // advance
      if (m == count) break;
      part4 += m_points_weights[set][m]; ++m; }

      // advance
      if (l == count) break;
      part3 += m_points_weights[set][l]; ++l; }

      // advance
      if (k == count) break;
      part2 += m_points_weights[set][k]; ++k; }

      // advance
      if (j == count) break;
      part1 += m_points_weights[set][j]; ++j; }

      // advance
      part0 += m_points_weights[set][i];
    }
#else
    // complexity:
    //  1 * (c - 0) * (c - 0) * (c - 0) * (c - 0) * (c - 0) * (c - 1) +
    //  1 * (c - 1) * (c - 1) * (c - 1) * (c - 1) * (c - 1) * (c - 1) +
    //  1 * (c - 2) * (c - 2) * (c - 2) * (c - 2) * (c - 2) * (c - 2) +
    // ...
    //  1 * (c - c) * (c - c) * (c - c) * (c - c) * (c - c) * (c - c)
    //
    //  (c - 0) ^ 6 + (c - 1) ^ 6 + (c - 2) ^ 6 + ... + (c - c) ^ 6
    //
    // worst case, c == 16
    //  1 * 16 * 16 * 16 * 16 * 16 * 15 +
    //  1 * 15 * 15 * 15 * 15 * 15 * 15 +
    //  1 * 14 * 14 * 14 * 14 * 14 * 14 +
    //  ...
    //
    // 2^24 + 11390625 + 7529536 + 4826809 + 2985984 + 1771561 + 1000000
    // 531441 + 262144 + 117649 + 46656 + 15625 + 4096 + 729 + 64 + 1
    // = 47.260.136
    //
    // if the inner loop body takes 90 cycles it takes 1 second
    // on a 4GHz machine to compress 1 block
    //
    // every single cycle removed/optimized gives us a speedup of 1ms or
    // a speed gain of 1%
    //
    // if we'd use the graphics card we could simply start 47 million threads
    // say it's the same 90 cycles, we start 64 item thread-groups, down to
    // 738.439 threads-groups, say 10 can run concurrently, we have to process
    // 90 cycles about 73.843 times, equals 6.645.956 cycles or 1/161th of
    // a second on 1 GHz card, equals 6ms or 161 blocks per second, with is
    // about one compressed row of a texture
    //
    // almost acceptable, but it's only one iteration of worst-case

    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

    // second cluster [i,j) is one seventh along
    Vec4 part1 = VEC4_CONST(0.0f);
    for (int j = i;;) {
	  Vec4 sub1  =   part0 + part1;
	  Vec4 add1a = weight6 * part1;
	  Vec4 add1b = weight1 * part1;

    // third cluster [j,k) is one seventh along
    Vec4 part2 = VEC4_CONST(0.0f);
    for (int k = j;;) {
	  Vec4 sub2  =    sub1 + part2;
	  Vec4 add2a = weight5 * part2 + add1a;
	  Vec4 add2b = weight2 * part2 + add1b;

    // fourth cluster [k,l) is one seventh along
    Vec4 part3 = VEC4_CONST(0.0f);
    for (int l = k;;) {
	  Vec4 sub3  =    sub2 + part3;
	  Vec4 add3a = weight4 * part3 + add2a;
	  Vec4 add3b = weight3 * part3 + add2b;

    // fifth cluster [l,m) is one seventh along
    Vec4 part4 = VEC4_CONST(0.0f);
    for (int m = l;;) {
	  Vec4 sub4  =    sub3 + part4;
	  Vec4 add4a = weight3 * part4 + add3a;
	  Vec4 add4b = weight4 * part4 + add3b;

	  Vec4 sum34 = fournineths  * (part3 + part4).SplatW();

    // sixths cluster [m,n) is one seventh along
    Vec4 part5 = VEC4_CONST(0.0f);
    for (int n = m;;) {
	  Vec4 sub5  =    sub4 + part5;
	  Vec4 add5a = weight2 * part5 + add4a;
	  Vec4 add5b = weight5 * part5 + add4b;

	  Vec4 sum25 = threenineths * (part2 + part5).SplatW();

    // seventh cluster [n,o) is six seventh along
    Vec4 part6 = (n == 0) ? m_points_weights[set][0] : VEC4_CONST(0.0f);
    int omin = (n == 0) ? 1 : n;
    for (int o = omin;;) {
	  Vec4 sub6  =    sub5 + part6;
	  Vec4 add6a = weight1 * part6 + add5a;
	  Vec4 add6b = weight6 * part6 + add5b;

	  Vec4 sum16 = twonineths   * (part1 + part6).SplatW();

	  // last cluster [k,count) is at the end
	  Vec4 part7 = xsum_wsum - sub6;
//	  Vec4 partI = xsum_wsum - (sub6 - part0);

	  // compute least squares terms directly
	  Vec4 const alphax_sum = add6a + part0;
	  Vec4 const  betax_sum = add6b + part7;

	  Vec4 const alpha2_sum = alphax_sum.SplatW();
	  Vec4 const  beta2_sum =  betax_sum.SplatW();

	  Vec4 const alphabeta_sum = sum16 + sum25 + sum34;

	  // equivalent
//	  Vec4 const  betac_sum =                                               xsum_wsum                                                    - alphax_sum;
//	  Vec4 const beta2a_sum = (twonineths_ * (part1 + part6) + threenineths_ * (part2 + part5) + fournineths_ * (part3 + part4) + partI) - alpha2_sum;
//	  Vec4 const beta2b_sum = (twoninethsc *      sum16      + threeninethsc *      sum25      + fourninethsc *      sum34      + partI) - alpha2_sum;

	  //  alphax.t = p7 * 0/7  + p6 * 1/7  + p5 *  2/7  + p4 *  3/7  + p3 *  4/7  + p2 *  5/7  + p1 * 6/7  + p0 * 7/7
	  //   betax.t = p7 * 7/7  + p6 * 6/7  + p5 *  5/7  + p4 *  4/7  + p3 *  3/7  + p2 *  2/7  + p1 * 1/7  + p0 * 0/7
	  //   gamma.x = p7 * 7/7  + p6 * 7/7  + p5 *  7/7  + p4 *  7/7  + p3 *  7/7  + p2 *  7/7  + p1 * 7/7  + p0 * 7/7 = alpha.x + beta.x = xsum_wsum;
	  //
	  //  alphax.w =      0/7² +      1/7² +       2/7² +       3.7² +       4/7² +       5/7² +      6/7² +      ?/?²
	  //   betax.w =      ?.?² +      6/7² +       5/7² +       4/7² +       3/7² +       2/7² +      1/7² +      0/7²
	  // alphabeta = p7 * 0/-- + p6 * 6/49 + p5 * 10/49 + p4 * 12/49 + p3 * 12/49 + p2 * 10/49 + p1 * 6/49 + p0 * 0/--

	  // compute the least-squares optimal points
	  Vec4 factor = Reciprocal(NegativeMultiplySubtract(alphabeta_sum, alphabeta_sum, alpha2_sum * beta2_sum));
	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_sum, alphax_sum *  beta2_sum) * factor;
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_sum,  betax_sum * alpha2_sum) * factor;

	  // factor = 1.0 / (alpha.w * beta.w - alphabeta * alphabeta)
	  // a = (alpha.x * beta.w -  beta.x * alphabeta) * factor
	  // b = (beta.x * alpha.w - alpha.x * alphabeta) * factor
	  // a = (alpha.x * beta.w -  beta.x * alphabeta) / (alpha.w * beta.w - alphabeta * alphabeta)
	  // b = (beta.x * alpha.w - alpha.x * alphabeta) / (alpha.w * beta.w - alphabeta * alphabeta)

	  // snap floating-point-values to the integer-lattice
	  a = q.SnapToLattice(a, sb, 1 << SBSTART);
	  b = q.SnapToLattice(b, sb, 1 << SBEND);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b,  betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = Dot(e4, metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (besterror > eS) {
	    besterror = eS;

	    beststart = a;
	    bestend   = b;
	    bestiteration = iterationIndex;

	    besti = i,
	    bestj = j,
	    bestk = k,
	    bestl = l,
	    bestm = m,
	    bestn = n,
	    besto = o;
	  }

      // advance
      if (o == count) break;
      part6 += m_points_weights[set][o]; ++o; }

      // advance
      if (n == count) break;
      part5 += m_points_weights[set][n]; ++n; }

      // advance
      if (m == count) break;
      part4 += m_points_weights[set][m]; ++m; }

      // advance
      if (l == count) break;
      part3 += m_points_weights[set][l]; ++l; }

      // advance
      if (k == count) break;
      part2 += m_points_weights[set][k]; ++k; }

      // advance
      if (j == count) break;
      part1 += m_points_weights[set][j]; ++j; }

      // advance
      part0 += m_points_weights[set][i];
    }
#endif

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;

    // advance if possible
    ++iterationIndex;
    if (iterationIndex == m_iterationCount)
      break;

    // stop if a new iteration is an ordering that has already been tried
    Vec4 axis = KillW(bestend - beststart);
    if (!ConstructOrdering(axis, iterationIndex, set))
      break;
  }

  // assign closest points
  u8 const* order = (u8*)m_order[set] + 16 * bestiteration;

  for (int m =     0; m < besti; ++m)
    closest[set][order[m]] = 0;
  for (int m = besti; m < bestj; ++m)
    closest[set][order[m]] = 1;
  for (int m = bestj; m < bestk; ++m)
    closest[set][order[m]] = 2;
  for (int m = bestk; m < bestl; ++m)
    closest[set][order[m]] = 3;
  for (int m = bestl; m < bestm; ++m)
    closest[set][order[m]] = 4;
  for (int m = bestm; m < bestn; ++m)
    closest[set][order[m]] = 5;
  for (int m = bestn; m < besto; ++m)
    closest[set][order[m]] = 6;
  for (int m = besto; m < count; ++m)
    closest[set][order[m]] = 7;

  // copy rgb into the common start/end definition
  m_start[set] = beststart;
  m_end  [set] = bestend;

  return besterror;
}

Scr4 PaletteClusterFit::ClusterSearch8Constant(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb)
{
  /*
  Vec4 const weight1( 9.0f / 64.0f,  9.0f / 64.0f,  9.0f / 64.0f,   81.0f / 4096.0f);
  Vec4 const weight2(18.0f / 64.0f, 18.0f / 64.0f, 18.0f / 64.0f,  324.0f / 4096.0f);
  Vec4 const weight3(27.0f / 64.0f, 27.0f / 64.0f, 27.0f / 64.0f,  729.0f / 4096.0f);
  Vec4 const weight4(37.0f / 64.0f, 37.0f / 64.0f, 37.0f / 64.0f, 1369.0f / 4096.0f);
  Vec4 const weight5(46.0f / 64.0f, 46.0f / 64.0f, 46.0f / 64.0f, 2116.0f / 4096.0f);
  Vec4 const weight6(55.0f / 64.0f, 55.0f / 64.0f, 55.0f / 64.0f, 3025.0f / 4096.0f);
  Vec4 const twonineths                               = VEC4_CONST(486.0f / 4096.0f);
  */
  Vec4 const weight1(1.0f / 7.0f, 1.0f / 7.0f, 1.0f / 7.0f,  1.0f / 49.0f);
  Vec4 const weight2(2.0f / 7.0f, 2.0f / 7.0f, 2.0f / 7.0f,  4.0f / 49.0f);
  Vec4 const weight3(3.0f / 7.0f, 3.0f / 7.0f, 3.0f / 7.0f,  9.0f / 49.0f);
  Vec4 const weight4(4.0f / 7.0f, 4.0f / 7.0f, 4.0f / 7.0f, 16.0f / 49.0f);
  Vec4 const weight5(5.0f / 7.0f, 5.0f / 7.0f, 5.0f / 7.0f, 25.0f / 49.0f);
  Vec4 const weight6(6.0f / 7.0f, 6.0f / 7.0f, 6.0f / 7.0f, 36.0f / 49.0f);
  Vec4 const twonineths                        = VEC4_CONST( 6.0f / 49.0f);
  Vec4 const threenineths                      = VEC4_CONST(10.0f / 49.0f);
  Vec4 const fournineths                       = VEC4_CONST(12.0f / 49.0f);
  // onenineths   = sqr(1.0 / distance) = sqr(1.0 / 7.0) = (1.0 / 49.0)
  // twonineths   = weight1 * weight6 = (6.0 / 49.0)
  // threenineths = weight2 * weight5 = (10.0 / 49.0)
  // fournineths  = weight3 * weight4 = (12.0 / 49.0)

  // constants if weights == 1

  Vec4 const two  = VEC4_CONST(2.0f);
  Vec4 const half = VEC4_CONST(0.5f);
  Vec4 const zero = VEC4_CONST(0.0f);

  assume((count > 0) && (count <= 16));

  // match each point to the closest code
  int besti = 0, bestj = 0, bestk = 0;
  int bestl = 0, bestm = 0, bestn = 0;
  int besto = 0;

  int bestiteration = 0;
  Vec4 beststart    = VEC4_CONST(0.0f);
  Vec4 bestend      = VEC4_CONST(0.0f);

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle[set], 0, set);

  // check all possible clusters and iterate on the total order
  Scr4 besterror = Scr4(FLT_MAX);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_wsum = m_xsum_wsum[set];

    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

	  Vec4 sum34;
	  Vec4 sum25;
	  Vec4 sum16;

    // second cluster [i,j) is one seventh along
    Vec4 part1 = VEC4_CONST(0.0f);
         sum16 = VEC4_CONST(0.0f);
    for (int j = i;;) {
	  Vec4 sub1  =   part0 + part1;
	  Vec4 add1a = weight6 * part1;
	  Vec4 add1b = weight1 * part1;

    // third cluster [j,k) is one seventh along
    Vec4 part2 = VEC4_CONST(0.0f);
         sum25 = VEC4_CONST(0.0f);
    for (int k = j;;) {
	  Vec4 sub2  =    sub1 + part2;
	  Vec4 add2a = weight5 * part2 + add1a;
	  Vec4 add2b = weight2 * part2 + add1b;

    // fourth cluster [k,l) is one seventh along
    Vec4 part3 = VEC4_CONST(0.0f);
         sum34 = VEC4_CONST(0.0f);
    for (int l = k;;) {
	  Vec4 sub3  =    sub2 + part3;
	  Vec4 add3a = weight4 * part3 + add2a;
	  Vec4 add3b = weight3 * part3 + add2b;

    // fifth cluster [l,m) is one seventh along
    Vec4 part4 = VEC4_CONST(0.0f);
         sum34 = fournineths  * (part3).SplatW();
    for (int m = l;;) {
	  Vec4 sub4  =    sub3 + part4;
	  Vec4 add4a = weight3 * part4 + add3a;
	  Vec4 add4b = weight4 * part4 + add3b;

    // sixths cluster [m,n) is one seventh along
    Vec4 part5 = VEC4_CONST(0.0f);
         sum25 = threenineths * (part2).SplatW();
    for (int n = m;;) {
	  Vec4 sub5  =    sub4 + part5;
	  Vec4 add5a = weight2 * part5 + add4a;
	  Vec4 add5b = weight5 * part5 + add4b;

//	  fprintf(stderr, "/* begin part6 iterations */\n");

    // seventh cluster [n,o) is six seventh along
    Vec4 part6 = (n == 0) ? m_points_weights[set][0] : VEC4_CONST(0.0f);
         sum16 = (n == 0) ? twonineths * (part1.SplatW() + Vec4(1.0f)) :
         	            twonineths * (part1.SplatW());
    int omin = (n == 0) ? 1 : n;
    for (int o = omin;;) {
	  Vec4 sub6  =    sub5 + part6;
	  Vec4 add6a = weight1 * part6 + add5a;
	  Vec4 add6b = weight6 * part6 + add5b;

#if 0
#define	limit 1e-5
//	  assert(fabs(_alpha2_sum.W() - alpha2_sum.X()) < limit);
//	  assert(fabs(_beta2_sum.W() - beta2_sum.X()) < limit);
//	  assert(fabs(_alphabeta_sum.W() - alphabeta_sum.X()) < limit);

	  Vec4 _sum34 = fournineths  * (part3 + part4).SplatW();
	  Vec4 _sum25 = threenineths * (part2 + part5).SplatW();
	  Vec4 _sum16 = twonineths   * (part1 + part6).SplatW();
#endif

	  // last cluster [k,count) is at the end
	  Vec4 part7 = xsum_wsum - sub6;

	  // compute least squares terms directly
	  Vec4 const alphax_sum = add6a + part0;
	  Vec4 const  betax_sum = add6b + part7;
//	  Vec4 const  betac_sum = xsum_wsum - alphax_sum;

	  Vec4 const alpha2_sum = alphax_sum.SplatW();
	  Vec4 const  beta2_sum =  betax_sum.SplatW();

	  // TODO: the inner alphabeta_sum seems always to be the same sequence
	  Vec4 const alphabeta_sum = sum16 + sum25 + sum34;

//	  fprintf(stderr, "  {%.9ff},\n", alphabeta_sum.W());

	  //  alphax.t = p7 * 0/7  + p6 * 1/7  + p5 *  2/7  + p4 *  3/7  + p3 *  4/7  + p2 *  5/7  + p1 * 6/7  + p0 * 7/7
	  //   betax.t = p7 * 7/7  + p6 * 6/7  + p5 *  5/7  + p4 *  4/7  + p3 *  3/7  + p2 *  2/7  + p1 * 1/7  + p0 * 0/7
	  //   gamma.x = p7 * 7/7  + p6 * 7/7  + p5 *  7/7  + p4 *  7/7  + p3 *  7/7  + p2 *  7/7  + p1 * 7/7  + p0 * 7/7 = alpha.x + beta.x = xsum_wsum;
	  //
	  //  alphax.w =      0/7² +      1/7² +       2/7² +       3.7² +       4/7² +       5/7² +      6/7² +      ?/?²
	  //   betax.w =      ?.?² +      6/7² +       5/7² +       4/7² +       3/7² +       2/7² +      1/7² +      0/7²
	  // alphabeta = p7 * 0/-- + p6 * 6/49 + p5 * 10/49 + p4 * 12/49 + p3 * 12/49 + p2 * 10/49 + p1 * 6/49 + p0 * 0/--

	  // compute the least-squares optimal points
	  Vec4 factor = Reciprocal(NegativeMultiplySubtract(alphabeta_sum, alphabeta_sum, alpha2_sum * beta2_sum));
	  Vec4 a = NegativeMultiplySubtract( betax_sum, alphabeta_sum, alphax_sum *  beta2_sum) * factor;
	  Vec4 b = NegativeMultiplySubtract(alphax_sum, alphabeta_sum,  betax_sum * alpha2_sum) * factor;

	  // snap floating-point-values to the integer-lattice
	  a = q.SnapToLattice(a, sb, 1 << SBSTART);
	  b = q.SnapToLattice(b, sb, 1 << SBEND);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b,  betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Scr4 eS = Dot(e4, metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (besterror > eS) {
	    besterror = eS;

	    beststart = a;
	    bestend   = b;
	    bestiteration = iterationIndex;

	    besti = i,
	    bestj = j,
	    bestk = k,
	    bestl = l,
	    bestm = m,
	    bestn = n,
	    besto = o;
	  }

      // advance
      if (o == count) break;
      sum16 += twonineths;
      part6 += m_points_weights[set][o]; ++o; }

      // advance
      if (n == count) break;
      sum25 += threenineths;
      part5 += m_points_weights[set][n]; ++n; }

      // advance
      if (m == count) break;
      sum34 += fournineths;
      part4 += m_points_weights[set][m]; ++m; }

      // advance
      if (l == count) break;
      sum34 += fournineths;
      part3 += m_points_weights[set][l]; ++l; }

      // advance
      if (k == count) break;
      sum25 += threenineths;
      part2 += m_points_weights[set][k]; ++k; }

      // advance
      if (j == count) break;
      sum16 += twonineths;
      part1 += m_points_weights[set][j]; ++j; }

      // advance
      part0 += m_points_weights[set][i];
    }

    // stop if we didn't improve in this iteration
    if (bestiteration != iterationIndex)
      break;

    // advance if possible
    ++iterationIndex;
    if (iterationIndex == m_iterationCount)
      break;

    // stop if a new iteration is an ordering that has already been tried
    Vec4 axis = KillW(bestend - beststart);
    if (!ConstructOrdering(axis, iterationIndex, set))
      break;
  }

  // assign closest points
  u8 const* order = (u8*)m_order[set] + 16 * bestiteration;

  for (int m =     0; m < besti; ++m)
    closest[set][order[m]] = 0;
  for (int m = besti; m < bestj; ++m)
    closest[set][order[m]] = 1;
  for (int m = bestj; m < bestk; ++m)
    closest[set][order[m]] = 2;
  for (int m = bestk; m < bestl; ++m)
    closest[set][order[m]] = 3;
  for (int m = bestl; m < bestm; ++m)
    closest[set][order[m]] = 4;
  for (int m = bestm; m < bestn; ++m)
    closest[set][order[m]] = 5;
  for (int m = bestn; m < besto; ++m)
    closest[set][order[m]] = 6;
  for (int m = besto; m < count; ++m)
    closest[set][order[m]] = 7;

  // copy rgb into the common start/end definition
  m_start[set] = beststart;
  m_end  [set] = bestend;

  return besterror;
}

#ifdef	FEATURE_METRIC_SQUARED
#define CMetric(m)  m * m
#else
#define CMetric(m)  m
#endif

void PaletteClusterFit::CompressS23(void* block, vQuantizer &q, int mode)
{
  int ib = GetIndexBits(mode);
  int jb = ib >> 16; ib = ib & 0xFF;
  int cb = GetPrecisionBits(mode);
  int ab = cb >> 16; cb = cb & 0xFF;
  int zb = GetSharedField();

  q.ChangeShared(cb, cb, cb, ab, zb);

  // match each point to the closest code
  a16 u8 closest[4][16];
  Scr4 error = Scr4(0.0f);

  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;
  bool const trns = m_palette->IsMergedAlpha() && m_palette->IsTransparent();
  u8  const tmask = trns ? 0xFF : 0x00;

  assume((isets >  0) && (isets <= 3));
  assume((asets >= 0) && (asets <= 3));
  assume(((isets    +    asets) <= 3));

  // create a codebook
  Vec4 codes[1 << 4];

  // loop over all multi-channel sets
  for (int s = 0, sb = zb; s < (isets + asets); s++, sb >>= 1) {
    // how big is the codebook for the current set
    int kb = ((s < isets) ^ (!!m_swapindex)) ? ib : jb;

    // declare variables
    int const count = m_palette->GetCount(s);

    // in case of separate alpha the colors of the alpha-set have all been set to alpha
    Vec4 rmetric = m_metric[s < isets ? 0 : 1];
    Vec4 xmetric = m_metric[s < isets ? 0 : 2];

    // we do single entry fit for sparse sets
    if (count == 1) {
      // clear alpha-weight if alpha is disabled
      // in case of separate alpha the colors of the alpha-set have all been set to alpha
      u8 mask = ((s < isets) ? 0x7 : 0x8) | tmask;

      // find the closest code
      Scr4 dist = ComputeEndPoints(s, rmetric, cb, ab, sb, kb, mask);

      // save the index (it's just a single one)
      closest[s][0] = GetIndex();

      // accumulate the error
      error += dist * m_palette->GetWeights(s)[0];
    }
#if 0 // cluster fit is very likely much better than the quick index-fit
    // we do dual entry fit for sparse sets
    else if (count == 2) {
      // find the closest codes (it's just two)
      Scr4 dist = StretchEndPoints(s, rmetric, q, sb, kb, closest[s]);

      // accumulate the error
      error += dist;
    }
#endif
    // we do single channel fit for single component sets
    // (cluster-fit is really really really bad on one channel cases)
    // the separate alpha-channel is splatted into the rgb channels
    else if (IsChannel(s)) {
      // find the closest code
      Scr4 dist = ComputeCodebook(s, rmetric, q, sb, kb, closest[s]);

      // accumulate the error
      error += dist;
    }
    else {
      // metric is squared as well
      Vec4 cmetric = CMetric(xmetric);
      Scr4 cerror = Scr4(0.0f);

#if defined(TRACK_STATISTICS)
      /* there is a clear skew towards unweighted clusterfit (all weights == 1)
       *
       * kb == 2, numset ==
       *  [1]	0x00340950 {4657, 9679}
       *  [2]	0x00340958 {379, 110}
       *  [3]	0x00340960 {1, 0}
       * kb == 3, numset ==
       *  [1]	0x00340950 {11587, 5155}
       *  [2]	0x00340958 {9950, 605}
       *  [3]	0x00340960 {194, 0}
       */
      if (m_palette->GetCount(s) > 1) {
	if (m_optimizable[s])
	  gstat.has_noweightsets[kb][s][0]++;
	else
	  gstat.has_noweightsets[kb][s][1]++;
      }
#endif

      // accumulate the error
      assume(kb >= 2 && kb <= 3);
      switch(kb) {
        case 2:
          if (trns)
            cerror = ClusterSearch4Alpha   (closest, count, s,       cmetric , q, sb);
          else if (m_optimizable[s])
            cerror = ClusterSearch4Constant(closest, count, s, KillW(cmetric), q, sb);
          else
            cerror = ClusterSearch4        (closest, count, s, KillW(cmetric), q, sb);
          break;
        case 3:
          if (m_optimizable[s])
	    cerror = ClusterSearch8Constant(closest, count, s, KillW(cmetric), q, sb);
          else
	    cerror = ClusterSearch8        (closest, count, s, KillW(cmetric), q, sb);
          break;
      }

//    cerror = cerror + Dot(m_xxsum_wwsum[s], cmetric);

      // accumulate the real! error
      error += cerror;
    }

    // kill early if this scheme looses
    if (!(error < m_besterror))
      return;
  }

  // because the original alpha-channel's weight was killed it is completely random and need to be set to 1.0f
  if (!m_palette->IsTransparent()) {
    switch (m_palette->GetRotation()) {
      default: for (int a = 0; a < isets + asets; a++) m_start[a].Set<3>(1.0f), m_end[a].Set<3>(1.0f); break;
      case  1: for (int a = 0; a <         isets; a++) m_start[a].Set<0>(1.0f), m_end[a].Set<0>(1.0f); break;
      case  2: for (int a = 0; a <         isets; a++) m_start[a].Set<1>(1.0f), m_end[a].Set<1>(1.0f); break;
      case  3: for (int a = 0; a <         isets; a++) m_start[a].Set<2>(1.0f), m_end[a].Set<2>(1.0f); break;
    }
  }

  // TODO: cluster-fit is still worst than range-fit currently, why?
#if 1 //ndef NDEBUG
  // kill late if this scheme looses
  error = Scr4(0.0f); SumError(closest, q, mode, error);
  // the error coming back from the palettesinglefit is not entirely exact with OLD_QUANTIZERR
  if (!(error < m_besterror))
    return;
  
#ifndef NDEBUG
  fprintf(stderr, "ClstrFit m: %1d, s: %1d+%1d, n:", mode, isets, asets);
  fprintf(stderr, " %2d", 0 < (isets + asets) ? m_palette->GetCount(0) : 0);
  fprintf(stderr, " %2d", 1 < (isets + asets) ? m_palette->GetCount(1) : 0);
  fprintf(stderr, " %2d", 2 < (isets + asets) ? m_palette->GetCount(2) : 0);
  fprintf(stderr, ", c: %1d, a: %1d", cb, ab);
  if (GetPartitionBits(mode) > 0)
    fprintf(stderr, ", p: %2d",  m_palette->GetPartition());
  else if (GetRotationBits(mode) > 0)
    fprintf(stderr, ", r: %2d", m_palette->GetRotation());
  else
    fprintf(stderr, ", r: --");
  if (GetSelectionBits(mode) > 0)
    fprintf(stderr, ", i: %1d,%1d", !m_swapindex ? ib : jb, !m_swapindex ? jb : ib);
  else
    fprintf(stderr, ", i:  %1d ", ib);
  
  fprintf(stderr, ", e: %.8f (< %.8f)\n", error.X(), m_besterror.X() == FLT_MAX ? -1.0f : m_besterror.X());
#endif
#endif

  // use these shared bit configuration
  m_sharedbits = m_sharedbits;

  for (int s = 0;     s <  isets         ; s++)
    m_palette->RemapIndices(closest[s], m_indices[0], s);
  for (int a = isets; a < (isets + asets); a++) {
    m_palette->RemapIndices(closest[a], m_indices[1], a);

    // copy alpha into the common start/end definition
    m_start[a - isets] = TransferW(m_start[a - isets], m_start[a]);
    m_end  [a - isets] = TransferW(m_end  [a - isets], m_end  [a]);
  }

  typedef u8 (&Itwo)[2][16];
  typedef u8 (&Ione)[1][16];

  typedef Vec4 (&V4thr)[3];
  typedef Vec4 (&V4two)[2];
  typedef Vec4 (&V4one)[1];

  // save the block
  int partition = m_palette->GetPartition();
  int rot = m_palette->GetRotation(), sel = m_swapindex;
  switch (mode) {
    case 0: WritePaletteBlock3_m1(partition, (V4thr)m_start, (V4thr)m_end, m_sharedbits, (Ione)m_indices, block); break;
    case 1: WritePaletteBlock3_m2(partition, (V4two)m_start, (V4two)m_end, m_sharedbits, (Ione)m_indices, block); break;
    case 2: WritePaletteBlock3_m3(partition, (V4thr)m_start, (V4thr)m_end, m_sharedbits, (Ione)m_indices, block); break;
    case 3: WritePaletteBlock3_m4(partition, (V4two)m_start, (V4two)m_end, m_sharedbits, (Ione)m_indices, block); break;
    case 4: WritePaletteBlock4_m5(rot, sel , (V4one)m_start, (V4one)m_end, m_sharedbits, (Itwo)m_indices, block); break;
    case 5: WritePaletteBlock4_m6(rot,       (V4one)m_start, (V4one)m_end, m_sharedbits, (Itwo)m_indices, block); break;
    case 6: WritePaletteBlock4_m7(partition, (V4one)m_start, (V4one)m_end, m_sharedbits, (Ione)m_indices, block); break;
    case 7: WritePaletteBlock4_m8(partition, (V4two)m_start, (V4two)m_end, m_sharedbits, (Ione)m_indices, block); break;
  }

  // save the error
  m_besterror = error;
  m_best = true;
}

void PaletteClusterFit::CompressC2(void* block, vQuantizer &q, int mode) {
  /* 2bit it can be done by CompressS23 as well */
  CompressS23(block, q, mode);
}

void PaletteClusterFit::Compress(void* block, vQuantizer &q, int mode) {
  switch (mode) {
#if (CLUSTERINDICES >= 2)
    case 2: /*2*/  CompressS23(block, q, mode); break;
    case 3: /*2*/  CompressS23(block, q, mode); break;
    case 5: /*22*/ CompressS23(block, q, mode); break;

    case 7: /*2*/  CompressC2 (block, q, mode); break;
#endif

#if (CLUSTERINDICES >= 3)
    case 0: /*3*/  CompressS23(block, q, mode); break;
    case 1: /*3*/  CompressS23(block, q, mode); break;
    case 4: /*23*/ CompressS23(block, q, mode); break;
#endif

#if (CLUSTERINDICES >= 4)
    case 6: /*CompressC4(block, q, mode);*/ break;
#endif
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
void ClusterFit_CCR::AssignSet(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, const int metric, const int fit ) amp_restricted
{
  PaletteFit_CCR::AssignSet(barrier, thread, m_palette, metric, fit);

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
  const int count = m_palette.GetCount();
  point16 values = m_palette.GetPoints();
  weight16 weights = m_palette.GetWeights();

  // get the covariance smatrix
  Sym3x3 covariance = ComputeWeightedCovariance3(barrier, thread, count, values, weights);

  // compute the principle component
  float4 principle = GetPrincipleComponent(barrier, thread, covariance);

  threaded_cse(0) {
    // compute the principle component
    m_principle = principle;
  }
}

bool ClusterFit_CCR::ConstructOrdering(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, float4r axis, int iteration) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // cache some values (AMP: causes a full copy, for indexibility)
  const int count = m_palette.GetCount();
  point16 values = m_palette.GetPoints();

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
  point16 unweighted = m_palette.GetPoints();
  weight16 weights = m_palette.GetWeights();

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

void ClusterFit_CCR::Compress(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, out code64 block, const bool trans,
			        IndexBlockLUT yArr) amp_restricted
{
  /* all or nothing branches, OK, same for all threads */
  bool isBtc1 = (trans);
  if (isBtc1) {
    Compress3(barrier, thread, m_palette, block, yArr);
    if (!m_palette.IsTransparent())
      Compress4(barrier, thread, m_palette, block, yArr);
  }
  else
    Compress4(barrier, thread, m_palette, block, yArr);
}

void ClusterFit_CCR::Compress3(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, out code64 block,
			         IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // declare variables
  const int count = m_palette.GetCount();

  const float4 two = 2.0;
  const float4 half_half2 = float4(0.5f, 0.5f, 0.5f, 0.25f);
  const float4 half = 0.5f;
  const float4 grid = float4(31.0f, 63.0f, 31.0f, 0.0f);
  const float4 gridrcp = float4(1.0f/31.0f, 1.0f/63.0f, 1.0f/31.0f, 0.0f);

  // prepare an ordering using the principle axis
  float4 bestline[CVALS];

  bestline[CSTRT] = 0.0f;
  bestline[CSTOP] = m_principle;

  // check all possible clusters and iterate on the total order
  float besterror = m_besterror;
  int bestiteration = 0;
  int besti = 0, bestj = 0;

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0; iterationIndex < m_iterationCount; iterationIndex++) {
    float4 axis = (m_line[CSTOP] - m_line[CSTRT]);

    // stop if a new iteration is an ordering that has already been tried
    if (!ConstructOrdering(barrier, thread, m_palette, axis, iterationIndex) && iterationIndex)
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
	  bestline[CSTRT] = a.GetVec4();
	  bestline[CSTOP] = b.GetVec4();
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
    PaletteFit_CCR::Compress3(barrier, thread, m_palette, block, yArr);
  }
}

void ClusterFit_CCR::Compress4(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, out code64 block,
			         IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // declare variables
  const int count = m_palette.GetCount();

  const float4 two = 2.0f;
  const float4 weight1 = float4( 1.0f/3.0f, 1.0f/3.0f, 1.0f/3.0f, 1.0f/9.0f );
  const float4 weight2 = float4( 2.0f/3.0f, 2.0f/3.0f, 2.0f/3.0f, 4.0f/9.0f );
  const float4 twonineths = 2.0f/9.0f;
  const float4 half = 0.5f;
  const float4 grid = float4( 31.0f, 63.0f, 31.0f, 0.0f );
  const float4 gridrcp = float4( 1.0f/31.0f, 1.0f/63.0f, 1.0f/31.0f, 0.0f );

  // prepare an ordering using the principle axis
  float4 bestline[CVALS];

  bestline[CSTRT] = 0.0f;
  bestline[CSTOP] = m_principle;

  // check all possible clusters and iterate on the total order
  float besterror = m_besterror;
  int bestiteration = 0;
  int besti = 0, bestj = 0, bestk = 0;

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0; iterationIndex < m_iterationCount; iterationIndex++) {
    float4 axis = (m_line[CSTOP] - m_line[CSTRT]);

    // stop if a new iteration is an ordering that has already been tried
    if (!ConstructOrdering(barrier, thread, m_palette, axis, iterationIndex) && iterationIndex)
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
	  float4 alphax_sum = muladd( part2, weight1, muladd( part1, weight2, part0 ) );
	  float4 alpha2_sum; alpha2_sum = alphax_sum.w;

	  float4 betax_sum = muladd( part1, weight1, muladd( part2, weight2, part3 ) );
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
	    bestline[CSTRT] = a.GetVec4();
	    bestline[CSTOP] = b.GetVec4();
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
    PaletteFit_CCR::Compress4(barrier, thread, m_palette, block, yArr);
  }
}
#endif

} // namespace squish
