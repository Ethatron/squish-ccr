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

#include "colourrangefit.h"
#include "colourset.h"
#include "colourblock.h"

#include "coloursinglefit.h"
#include "coloursinglesnap.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
ColourRangeFit::ColourRangeFit(ColourSet const* colours, int flags)
  : ColourFit(colours, flags)
{
  // initialize endpoints
  ComputeEndPoints();
}

void ColourRangeFit::ComputeEndPoints()
{
  cQuantizer3<5,6,5> q = cQuantizer3<5,6,5>();

  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* weights = m_colours->GetWeights();

  Sym3x3 covariance;
  Vec3 centroid;
  Vec3 principle;

  // get the covariance matrix
  if (m_colours->IsUnweighted())
    ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric);
  else
    ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric, weights);

  // compute the principle component
  GetPrincipleComponent(covariance, principle);

  // get the min and max range as the codebook endpoints
  Vec3 start(0.0f);
  Vec3 end(0.0f);

  if (count > 0) {
#ifdef	FEATURE_RANGEFIT_PROJECT
    // compute the projection
    GetPrincipleProjection(start, end, principle, centroid, count, values);
#else
    Scr3 min, max;

    // compute the range
    start = end = values[0];
    min = max = Dot(values[0], principle);

    for (int i = 1; i < count; ++i) {
      Scr3 val = Dot(values[i], principle);

      if (val < min) {
	start = values[i];
	min = val;
      }
      else if (val > max) {
	end = values[i];
	max = val;
      }
    }
#endif
  }

  // snap floating-point-values to the integer-lattice and save
  m_start = q.SnapToLattice(start);
  m_end   = q.SnapToLattice(end  );
}

void ColourRangeFit::Compress3b(void* block)
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

void ColourRangeFit::Compress3(void* block)
{
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();

  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[3]; Codebook3(codes, m_metric * m_start, m_metric * m_end);

  // match each point to the closest code
  u8 closest[16];

  Scr3 error = Scr3(DISTANCE_BASE);
  for (int i = 0; i < count; ++i) {
    int idx = 0;

    // find the closest code
    Vec3 value = m_metric * values[i];
    Scr3 dist; MinDistance3<true>(dist, idx, value, codes);

    // save the index
    closest[i] = (u8)idx;

    // accumulate the error
    error += dist * freq[i];
  }

  // save this scheme if it wins
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16]; m_colours->RemapIndices(closest, indices);

    // save the block
    WriteColourBlock3(m_start, m_end, indices, block);
  }
}

void ColourRangeFit::Compress4(void* block)
{
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();

  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[4]; Codebook4(codes, m_metric * m_start, m_metric * m_end);

  // match each point to the closest code
  u8 closest[16];

  Scr3 error = Scr3(DISTANCE_BASE);
  for (int i = 0; i < count; ++i) {
    int idx = 0;

    // find the closest code
    Vec3 value = m_metric * values[i];
    Scr3 dist; MinDistance4<true>(dist, idx, value, codes);

    // accumulate the error
    AddDistance(dist, error, freq[i]);

    // save the index
    closest[i] = (u8)idx;
  }

  // save this scheme if it wins
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16]; m_colours->RemapIndices(closest, indices);

    // save the block
    WriteColourBlock4(m_start, m_end, indices, block);
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#if	defined(SQUISH_USE_COMPUTE)
    tile_static float dist[16];
    tile_static int   mins[8];
    tile_static int   maxs[8];
#endif

void ColourRangeFit_CCR::AssignSet(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, const int metric, const int fit ) amp_restricted
{
  ColourFit_CCR::AssignSet(barrier, thread, m_colours, metric, fit);

#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  threaded_cse(0) {
    // initialize the metric
    if (metric == SQUISH_METRIC_UNIT)
      m_metric = float3(0.5000f, 0.5000f, 0.0000f);
    else if (metric == SQUISH_METRIC_PERCEPTUAL)
      m_metric = float3(0.2126f, 0.7152f, 0.0722f);
    else
      m_metric = float3(0.3333f, 0.3333f, 0.3333f);
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

  // get the min and max range as the codebook endpoints
  float3 cline[CVALS];

  /* start = float3(0.0f);
   * end   = float3(0.0f);
   */
  cline[CSTRT] = 0.0f;
  cline[CSTOP] = 0.0f;

  if (count > 0) {
    // calculate min/max
#if	!defined(SQUISH_USE_COMPUTE)
    tile_static float dist[16];
    tile_static int   mins[8];
    tile_static int   maxs[8];
#endif

    // same for all
    const int dn = (thread << 1) + 0;
    const int up = (thread << 1) + 1;

    /* calculate min/max
     *
     * AMP: reduction, O(count) vs. O(ln(16)) = O(4))
     * COMPUTE: error X3671: loop termination conditions in varying flow
     *          control cannot depend on data read from a UAV
     *          mean we can't do "for (i < count)" anyway
     */
    threaded_for(mm1, 8) {
      // dummy-fill values beyond count with value
      // from 0, this remains neutral in the minmax
      const float3 value0 = values[(dn < count ? dn : 0)];
      const float3 value1 = values[(up < count ? up : 0)];
      float d0 = dot(value0, principle);
      float d1 = dot(value1, principle);

      dist[dn] = d0;
      dist[up] = d1;

      // prefer lower indices (will mask the identical prefixed values)
      mins[mm1] = (d0 <= d1 ? dn : up);
      maxs[mm1] = (d0 >= d1 ? dn : up);
    }

    threaded_for(mm2, 4) {
      // prefer lower indices (will mask the identical prefixed values)
      mins[mm2] = (dist[mins[dn]] <= dist[mins[up]] ? mins[dn] : mins[up]);
      maxs[mm2] = (dist[maxs[dn]] >= dist[maxs[up]] ? maxs[dn] : maxs[up]);
    }

    threaded_for(mm3, 2) {
      // prefer lower indices (will mask the identical prefixed values)
      mins[mm3] = (dist[mins[dn]] <= dist[mins[up]] ? mins[dn] : mins[up]);
      maxs[mm3] = (dist[maxs[dn]] >= dist[maxs[up]] ? maxs[dn] : maxs[up]);
    }

    // make writes to "mins"/"maxs" visible to all
    tile_static_memory_fence(barrier);

    // prefer lower indices (will mask the identical prefixed values)
    int min = (dist[mins[0]] <= dist[mins[1]] ? mins[0] : mins[1]);
    int max = (dist[maxs[0]] >= dist[maxs[1]] ? maxs[0] : maxs[1]);

    /* start = values[min];
     * end   = values[max];
     */
    cline[CSTRT] = values[min];
    cline[CSTOP] = values[max];
  }

  // snap floating-point-values to the integer-lattice and save
  const float3 grid = float3( 31.0f, 63.0f, 31.0f );
  const float3 gridrcp = float3( 1.0f/31.0f, 1.0f/63.0f, 1.0f/31.0f );
  const float3 half = 0.5f;

  /* start = saturate(start);
   * end   = saturate(end  );
   *
   * m_line[CSTRT] = truncate(grid * start + half) * gridrcp;
   * m_line[CSTOP] = truncate(grid * end   + half) * gridrcp;
   */
  threaded_for(cs, CVALS) {
    cline[cs] = saturate(cline[cs]);
    m_line[cs] = truncate(grid * cline[cs] + half) * gridrcp;
  }
}

void ColourRangeFit_CCR::Compress(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block, const bool trans,
			    IndexBlockLUT yArr) amp_restricted
{
  /* all or nothing branches, OK, same for all threads */
  bool isBtc1 = (trans);
  if (isBtc1 && m_colours.IsTransparent())
    Compress3 (barrier, thread, m_colours, block, yArr);
  else if (!isBtc1)
    Compress4 (barrier, thread, m_colours, block, yArr);
  else
    Compress34(barrier, thread, m_colours, block, yArr);
}

#define CASE3	0
#define CASE4	1
#define CASES	2

void ColourRangeFit_CCR::Compress3(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
			     IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // cache some values
  // AMP: causes a full copy, for indexibility
  const int count = m_colours.GetCount();
  point16 values = m_colours.GetPoints();

  float3 codes[3];
  float  dists[3];

  // create a codebook
  codes[0] =        m_line[CSTRT]                       ;
  codes[1] =                               m_line[CSTOP];
  codes[2] = 0.5f * m_line[CSTRT] + 0.5f * m_line[CSTOP];

  // match each point to the closest code
  threaded_for(fcscan, count) {
    // find the least error and corresponding index (AMP: pull from shared to local)
    const float3 value = values[fcscan];

    dists[0] = lengthsquared(m_metric * (value - codes[0]));
    dists[1] = lengthsquared(m_metric * (value - codes[1]));
    dists[2] = lengthsquared(m_metric * (value - codes[2]));

    // find the closest code (AMP: vectorized reduction, cset)
    int idx10, idx21, idx;

    idx10 = 0 + (dists[1] < dists[0] ? 1 : 0);
    idx21 = 1 + (dists[2] < dists[1] ? 1 : 0);

    idx = (dists[idx10] < dists[idx21] ? idx10 : idx21);

    // save the index
    m_matches[1][fcscan] = (ccr8)idx;
  }

  // make writes to "m_matches" visible to all
  tile_static_memory_fence(barrier);

  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    ColourFit_CCR::Compress3(barrier, thread, m_colours, block, yArr);
  }
}

void ColourRangeFit_CCR::Compress4(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
			     IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // cache some values
  // AMP: causes a full copy, for indexibility
  const int count = m_colours.GetCount();
  point16 values = m_colours.GetPoints();

  float3 codes[4];
  float  dists[4];

  // create a codebook
  codes[0] =                 m_line[CSTRT]                                ;
  codes[1] =                                                 m_line[CSTOP];
  codes[2] = (2.0f / 3.0f) * m_line[CSTRT] + (1.0f / 3.0f) * m_line[CSTOP];
  codes[3] = (1.0f / 3.0f) * m_line[CSTRT] + (2.0f / 3.0f) * m_line[CSTOP];

  // match each point to the closest code
  threaded_for(fcscan, count) {
    // find the least error and corresponding index (AMP: pull from shared to local)
    const float3 value = values[fcscan];

    dists[0] = lengthsquared(m_metric * (value - codes[0]));
    dists[1] = lengthsquared(m_metric * (value - codes[1]));
    dists[2] = lengthsquared(m_metric * (value - codes[2]));
    dists[3] = lengthsquared(m_metric * (value - codes[3]));

    // find the closest code (AMP: vectorized reduction, cset)
    int idx10, idx32, idx;

    idx10 = 0 + (dists[1] < dists[0] ? 1 : 0);
    idx32 = 2 + (dists[3] < dists[2] ? 1 : 0);

    idx = (dists[idx10] < dists[idx32] ? idx10 : idx32);

    // save the index
    m_matches[1][fcscan] = (ccr8)idx;
  }

  // make writes to "m_matches" visible to all
  tile_static_memory_fence(barrier);

  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    ColourFit_CCR::Compress4(barrier, thread, m_colours, block, yArr);
  }
}

#if	defined(SQUISH_USE_COMPUTE)
  tile_static float3 codes[8];
  tile_static ccr8 closest[16][CASES];
  tile_static float errors[16][CASES];
#endif

void ColourRangeFit_CCR::Compress34(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
			      IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // cache some values
  // AMP: causes a full copy, for indexibility
  const int count = m_colours.GetCount();
  point16 values = m_colours.GetPoints();

  // match each point to the closest code
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static float3 codes[8];
  tile_static ccr8 closest[16][CASES];
  tile_static float errors[16][CASES];
#endif
  float error[CASES];

  /* create a codebook:
   *
   *  codes[0] = (6.0f / 6.0f) * m_line[CSTRT] + (0.0f / 6.0f) * m_line[CSTOP];
   *  codes[-] = (5.0f / 6.0f) * m_line[CSTRT] + (1.0f / 6.0f) * m_line[CSTOP];
   *  codes[2] = (4.0f / 6.0f) * m_line[CSTRT] + (2.0f / 6.0f) * m_line[CSTOP];
   *  codes[3] = (3.0f / 6.0f) * m_line[CSTRT] + (3.0f / 6.0f) * m_line[CSTOP];
   *  codes[4] = (2.0f / 6.0f) * m_line[CSTRT] + (4.0f / 6.0f) * m_line[CSTOP];
   *  codes[-] = (1.0f / 6.0f) * m_line[CSTRT] + (5.0f / 6.0f) * m_line[CSTOP];
   *  codes[6] = (0.0f / 6.0f) * m_line[CSTRT] + (6.0f / 6.0f) * m_line[CSTOP];
   */
  threaded_for(es, 7) {
    codes[es] = ((6.0f - es) / 6.0f) * m_line[CSTRT] +
	        ((0.0f + es) / 6.0f) * m_line[CSTOP];
  }

  // make it visible
  tile_static_memory_fence(barrier);

  /* error[CASE3] = 0.0f;
   * error[CASE4] = 0.0f;
   */
  threaded_for(fcscan, count) {
    // find the least error and corresponding index (AMP: pull from shared to local)
    const float3 value = values[fcscan];
    float dists[5];

    dists[0] = lengthsquared(m_metric * (value - codes[0]));
    dists[1] = lengthsquared(m_metric * (value - codes[2]));
    dists[2] = lengthsquared(m_metric * (value - codes[3]));
    dists[3] = lengthsquared(m_metric * (value - codes[4]));
    dists[4] = lengthsquared(m_metric * (value - codes[6]));

    // find the closest code (AMP: vectorized reduction, cset)
    int4 idx420; int idx[CASES];

    idx420.x = 2 * (0 + (dists[2] < dists[0] ? 1 : 0));
    idx420.y = 2 * (2 + (dists[4] < dists[2] ? 1 : 0));
    idx420.z = 1 * (0 + (dists[1] < dists[0] ? 1 : 0));
    idx420.w = 1 * (3 + (dists[4] < dists[3] ? 1 : 0));

    idx[CASE3] = (dists[idx420.x] < dists[idx420.y] ? idx420.x : idx420.y);
    idx[CASE4] = (dists[idx420.z] < dists[idx420.w] ? idx420.z : idx420.w);

    // save the index
    closest[fcscan][CASE3] = (ccr8)idx[CASE3];
    closest[fcscan][CASE4] = (ccr8)idx[CASE4];

    // accumulate the error of case 3/4
    errors[fcscan][CASE3] = dists[idx[CASE3]];
    errors[fcscan][CASE4] = dists[idx[CASE4]];
  }

  // same for all
  const int dn = (thread << 1) + 0;
  const int up = (thread << 1) + 1;

  // accumulate the error of case 3/4
  threaded_for(er1, 8) {
    // AMP: prefer 2-wide vectorized op
    errors[er1][CASE3] = (dn < count ? errors[dn][CASE3] : 0.0f) +
    	                 (up < count ? errors[up][CASE3] : 0.0f);
    errors[er1][CASE4] = (dn < count ? errors[dn][CASE4] : 0.0f) +
    	                 (up < count ? errors[up][CASE4] : 0.0f);
  }

  // accumulate the error of case 3/4
  threaded_for(er2, 4) {
    // AMP: prefer 2-wide vectorized op
    errors[er2][CASE3] = errors[dn][CASE3] + errors[up][CASE3];
    errors[er2][CASE4] = errors[dn][CASE4] + errors[up][CASE4];
  }

  // accumulate the error of case 3/4
  threaded_for(er3, 2) {
    // AMP: prefer 2-wide vectorized op
    errors[er3][CASE3] = errors[dn][CASE3] + errors[up][CASE3];
    errors[er3][CASE4] = errors[dn][CASE4] + errors[up][CASE4];
  }

  // make writes to "closest"/"errors" visible to all
  tile_static_memory_fence(barrier);

  // accumulate the error of case 3/4
  error[CASE3] = errors[0][CASE3] + errors[1][CASE3];
  error[CASE4] = errors[0][CASE4] + errors[1][CASE4];

  // AMP: all thread end up here, and all make this compare
  int is4 = (error[CASE4] < error[CASE3] ? CASE4 : CASE3);

  // save this scheme if it wins
  {
    // remap the indices
    threaded_for(mscan, count) {
      m_matches[1][mscan] = closest[mscan][is4];
    }

    // remap the indices
    // save the block
    ColourFit_CCR::Compress34(barrier, thread, m_colours, block, is4, yArr);
  }

#undef	CASE3
#undef	CASE4
#undef	CASES
}
#endif

} // namespace squish
