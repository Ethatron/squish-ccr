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

#include "paletteclusterfit.h"
#include "paletteset.h"
#include "paletteblock.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(USE_PRE)
PaletteClusterFit::PaletteClusterFit(PaletteSet const* palettes, int flags, int swap)
  : SinglePaletteFit(palettes, flags, swap)
{
  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;

  // set the iteration count
  m_iterationCount = (m_flags & kColourIterativeClusterFits) / kColourClusterFit;

  if (m_iterationCount > kMaxIterations) m_iterationCount = kMaxIterations;
  if (m_iterationCount < kMinIterations) m_iterationCount = kMinIterations;

  // loop over all sets
  for (int s = 0; s < (isets + asets); s++) {
    // cache some values
    int const count = m_palette->GetCount(s);
    Vec4 const* values = m_palette->GetPoints(s);
    float const* weights = m_palette->GetWeights(s);

    // we don't do this for sparse sets
    if (count != 1) {
      // get the covariance matrix
      Sym3x3 covariance = ComputeWeightedCovariance(count, values, weights);

      // compute the principle component
      m_principle[s] = Vec4(ComputePrincipleComponent(covariance), 0.0f);
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
  float const* weights = m_palette->GetWeights(set);

  m_xsum_wsum[set] = VEC4_CONST(0.0f);
  for (int i = 0; i < count; ++i) {
    int j = order[i];

    Vec4 p = TransferW(unweighted[j], Vec4(1.0f));
    Vec4 w(weights[j]);
    Vec4 x = p * w;

    m_points_weights[set][i] = x;
    m_xsum_wsum     [set]   += x;
  }

  return true;
}

Vec4 PaletteClusterFit::ClusterSearch4(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q)
{
  /*
  Vec4 const weight1(21.0f / 64.0f, 21.0f / 64.0f, 21.0f / 64.0f,  441.0f / 4096.0f);
  Vec4 const weight2(43.0f / 64.0f, 43.0f / 64.0f, 43.0f / 64.0f, 1849.0f / 4096.0f);
  Vec4 const twonineths                               = VEC4_CONST(882.0f / 4096.0f);
  */
  Vec4 const weight1(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 9.0f);
  Vec4 const weight2(2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f, 4.0f / 9.0f);
  Vec4 const twonineths                        = VEC4_CONST(2.0f / 9.0f);
  // onenineths = sqr(1.0 / distance) = sqr(1.0 / 3.0) = (1.0 / 9.0)
  // twonineths = weight1 * weight2 = (2.0 / 9.0)

  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const half = VEC4_CONST(0.5f);

  // match each point to the closest code
  int besti = 0, bestj = 0, bestk = 0;
  int bestiteration = 0;
  Vec4 beststart    = VEC4_CONST(0.0f);
  Vec4 bestend      = VEC4_CONST(0.0f);

  // prepare an ordering using the principle axis
  ConstructOrdering(m_principle[set], 0, set);

  // check all possible clusters and iterate on the total order
  Vec4 besterror = Vec4(FLT_MAX);

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

	  // compute least squares terms directly
	  Vec4 const alphax_sum = MultiplyAdd(part2, weight1, MultiplyAdd(part1, weight2, part0));
	  Vec4 const  betax_sum = MultiplyAdd(part1, weight1, MultiplyAdd(part2, weight2, part3));

	  Vec4 const alpha2_sum = alphax_sum.SplatW();
	  Vec4 const  beta2_sum =  betax_sum.SplatW();

	  Vec4 const alphabeta_sum = twonineths * (part1 + part2).SplatW();

	  //   alpha.x = p3 * 0/3  + p2 * 1/3  + p1 * 2/3  + p0 * 3/3
	  //    beta.x = p3 * 3/3  + p2 * 2/3  + p1 * 1/3  + p0 * 0/3
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
	  a = q.SnapToLattice(a);
	  b = q.SnapToLattice(b);

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
	  Vec4 eS = HorizontalAdd(e4 * metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (CompareFirstLessThan(eS, besterror)) {
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

Vec4 PaletteClusterFit::ClusterSearch8(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q)
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

  Vec4 const two = VEC4_CONST(2.0f);
  Vec4 const half = VEC4_CONST(0.5f);

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
  Vec4 besterror = Vec4(FLT_MAX);

  // loop over iterations (we avoid the case that all points in first or last cluster)
  for (int iterationIndex = 0;;) {
    // cache some values
    Vec4 const xsum_wsum = m_xsum_wsum[set];

    // first cluster [0,i) is at the start
    Vec4 part0 = VEC4_CONST(0.0f);
    for (int i = 0; i < count; ++i) {

    // second cluster [i,j) is one seventh along
    Vec4 part1 = VEC4_CONST(0.0f);
    for (int j = i;;) {

    // third cluster [j,k) is one seventh along
    Vec4 part2 = VEC4_CONST(0.0f);
    for (int k = j;;) {

    // fourth cluster [k,l) is one seventh along
    Vec4 part3 = VEC4_CONST(0.0f);
    for (int l = k;;) {

    // fifth cluster [l,m) is one seventh along
    Vec4 part4 = VEC4_CONST(0.0f);
    for (int m = l;;) {

    // sixths cluster [m,n) is one seventh along
    Vec4 part5 = VEC4_CONST(0.0f);
    for (int n = m;;) {

    // seventh cluster [n,o) is six seventh along
    Vec4 part6 = (n == 0) ? m_points_weights[set][0] : VEC4_CONST(0.0f);
    int omin = (n == 0) ? 1 : n;
    for (int o = omin;;) {

	  // last cluster [k,count) is at the end
	  Vec4 part7 = xsum_wsum - part6 - part5 - part4 - part3 - part2 - part1 - part0;

	  // compute least squares terms directly
	  Vec4 const alphax_sum =
	  	MultiplyAdd(part6, weight1,
	  	MultiplyAdd(part5, weight2,
	  	MultiplyAdd(part4, weight3,
	  	MultiplyAdd(part3, weight4,
	  	MultiplyAdd(part2, weight5,
	  	MultiplyAdd(part1, weight6,
	  	            part0))))));

	  Vec4 const betax_sum =
	  	MultiplyAdd(part1, weight1,
	  	MultiplyAdd(part2, weight2,
	  	MultiplyAdd(part3, weight3,
	  	MultiplyAdd(part4, weight4,
	  	MultiplyAdd(part5, weight5,
	  	MultiplyAdd(part6, weight6,
	  	            part7))))));

	  Vec4 const alpha2_sum = alphax_sum.SplatW();
	  Vec4 const  beta2_sum =  betax_sum.SplatW();

	  Vec4 const alphabeta_sum =
	  	twonineths   * (part1 + part6).SplatW() +
	  	threenineths * (part2 + part5).SplatW() +
	  	fournineths  * (part3 + part4).SplatW();

	  //  alphax.t = p7 * 0/7  + p6 * 1/7  + p5 *  2/7  + p4 *  3/7  + p3 *  4/7  + p2 *  5/7  + p1 * 6/7  + p0 * 7/7
	  //   betax.t = p7 * 7/7  + p6 * 6/7  + p5 *  5/7  + p4 *  4/7  + p3 *  3/7  + p2 *  2/7  + p1 * 1/7  + p0 * 0/7
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
	  a = q.SnapToLattice(a);
	  b = q.SnapToLattice(b);

	  // compute the error (we skip the constant xxsum)
	  Vec4 e1 = MultiplyAdd(a * a, alpha2_sum, b * b * beta2_sum);
	  Vec4 e2 = NegativeMultiplySubtract(a, alphax_sum, a * b * alphabeta_sum);
	  Vec4 e3 = NegativeMultiplySubtract(b,  betax_sum, e2);
	  Vec4 e4 = MultiplyAdd(two, e3, e1);

	  // apply the metric to the error term
	  Vec4 eS = HorizontalAdd(e4 * metric);

	  // keep the solution if it wins (error can be negative ...)
	  if (CompareFirstLessThan(eS, besterror)) {
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

void PaletteClusterFit::CompressS23(void* block, int mode)
{
  int ib = GetIndexBits(mode);
  int jb = ib >> 16; ib = ib & 0xFF;
  int cb = GetPrecisionBits(mode);
  int ab = cb >> 16; cb = cb & 0xFF;

//vQuantizer c = vQuantizer(cb, cb, cb, ab);
//vQuantizer a = vQuantizer(ab, ab, ab, ab);

  vQuantizer c = vQuantizer(5, 6, 5, ab);
  vQuantizer a = vQuantizer(ab, ab, ab, ab);

  // match each point to the closest code
  a16 u8 closest[4][16];
  Vec4 error = Vec4(0.0f);

  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;

  // create a codebook
  Vec4 codes[1 << 4];

  // loop over all sets
  for (int s = 0; s < (isets + asets); s++) {
    // how big is the codebook for the current set
    int kb = ((s < isets) ^ (!!m_swapindex)) ? ib : jb;

    // the separate alpha-channel is splatted into the rgb channels
    vQuantizer &q = (s < isets ? c : a);

    // declare variables
    int const count = m_palette->GetCount(s);

    // in case of separate alpha the colors of the alpha-set have all been set to alpha
    Vec4 metric = m_metric[s < isets ? 0 : 2];

    // we do single entry fit for sparse sets
    if (count == 1) {
      // clear alpha-weight if alpha is disabled
      // in case of separate alpha the colors of the alpha-set have all been set to alpha
      u8 mask = (ab ? ((s < isets) ? 0xF : 0x8) : 0x7);

      // find the closest code
      Vec4 dist = ComputeEndPoints(s, metric, q, cb, ab, kb, mask);

      // save the index (it's just a single one)
      closest[s][0] = GetIndex();

      // accumulate the error
      error += dist * 16.0f;
    }
    else {
      // accumulate the error
      assume(kb >= 2 && kb <= 3);
      switch(kb) {
        case 2: error += ClusterSearch4(closest, count, s, KillW(metric), q); break;
        case 3: error += ClusterSearch8(closest, count, s, KillW(metric), q); break;
      }
    }

    // kill early if this scheme looses
    if (CompareFirstGreaterThan(error, m_besterror))
      return;
  }

  // because the original alpha-channel's weight was killed it is completely random and need to be set to 1.0f
  if (!m_palette->IsTransparent()) {
    switch (m_palette->GetRotation()) {
      case 0: for (int a = 0; a < (isets + asets); a++) m_start[a].GetW() = m_end[a].GetW() = 1.0f; break;
      case 1: for (int a = 0; a < (isets + asets); a++) m_start[a].GetX() = m_end[a].GetX() = 1.0f; break;
      case 2: for (int a = 0; a < (isets + asets); a++) m_start[a].GetY() = m_end[a].GetY() = 1.0f; break;
      case 3: for (int a = 0; a < (isets + asets); a++) m_start[a].GetZ() = m_end[a].GetZ() = 1.0f; break;
    }
  }

  // TODO: cluster-fit is still worst than range-fit currently, why?
#if 1 //ndef NDEBUG
  // kill late if this scheme looses
  error = Vec4(0.0f); SumError(closest, mode, error);
  // the error coming back from the singlepalettefit is not entirely exact with OLD_QUANTIZERR
  if (CompareFirstGreaterThan(error, m_besterror))
    return;
#endif

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
    case 0: WritePaletteBlock3_m1(partition, (V4thr)m_start, (V4thr)m_end, (Ione)m_indices, block); break;
    case 1: WritePaletteBlock3_m2(partition, (V4two)m_start, (V4two)m_end, (Ione)m_indices, block); break;
    case 2: WritePaletteBlock3_m3(partition, (V4thr)m_start, (V4thr)m_end, (Ione)m_indices, block); break;
    case 3: WritePaletteBlock3_m4(partition, (V4two)m_start, (V4two)m_end, (Ione)m_indices, block); break;
    case 4: WritePaletteBlock4_m5(rot, sel , (V4one)m_start, (V4one)m_end, (Itwo)m_indices, block); break;
    case 5: WritePaletteBlock4_m6(rot,       (V4one)m_start, (V4one)m_end, (Itwo)m_indices, block); break;
    case 6: WritePaletteBlock4_m7(partition, (V4one)m_start, (V4one)m_end, (Ione)m_indices, block); break;
    case 7: WritePaletteBlock4_m8(partition, (V4two)m_start, (V4two)m_end, (Ione)m_indices, block); break;
  }

#if !defined(NDEBUG) && defined(DEBUG_QUANTIZER)
  // m4/m5 may have swap alpha, get the change back
  for (int a = isets; a < (isets + asets); a++) {
    m_start[a] = TransferW(m_start[a], m_start[a - isets]);
    m_end  [a] = TransferW(m_end  [a], m_end  [a - isets]);
  }
#endif

  // save the error
  m_besterror = error;
  m_best = true;
}

void PaletteClusterFit::Compress(void* block, int mode) {
  // TODO: cluster-fit of 8 points doesn't work at all yet
  switch (mode) {
#if (CLUSTERINDICES >= 2)
    case 2: /*2*/ CompressS23(block, mode); break;
    case 3: /*2*/ CompressS23(block, mode); break;
    case 5: /*22*/ CompressS23(block, mode); break;
#endif

#if (CLUSTERINDICES >= 3)
    case 0: /*3*/ CompressS23(block, mode); break;
    case 1: /*3*/ CompressS23(block, mode); break;
    case 4: /*23*/ CompressS23(block, mode); break;
#endif

    case 6: /*CompressC4(block, mode);*/ break;
    case 7: /*CompressC2(block, mode);*/ break;
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(USE_AMP) || defined(USE_COMPUTE)
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
  Sym3x3 covariance = ComputeWeightedCovariance(barrier, thread, count, values, weights);

  // compute the principle component
  float4 principle = ComputePrincipleComponent(barrier, thread, covariance);

  threaded_cse(0) {
    // compute the principle component
    m_principle = principle;
  }
}

bool ClusterFit_CCR::ConstructOrdering(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, float4r axis, int iteration) amp_restricted
{
#if	!defined(USE_COMPUTE)
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
#if	!defined(USE_COMPUTE)
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
#if	!defined(USE_COMPUTE)
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
