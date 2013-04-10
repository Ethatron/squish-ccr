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

#include "bitonerangefit.h"
#include "bitoneset.h"
#include "bitoneblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
BitoneRangeFit::BitoneRangeFit(BitoneSet const* bitones, int flags)
  : BitoneFit(bitones, flags)
{
  // initialize the best error
  m_besterror = Scr3(FLT_MAX);

  // cache some values
  int const count = m_bitones->GetCount();
  Vec3 const* values = m_bitones->GetPoints();
  float const* weights = m_bitones->GetWeights();

  Sym3x3 covariance;
  Vec3 centroid;
  Vec3 principle;

  // get the covariance matrix
  if (m_bitones->IsUnweighted())
    ComputeWeightedCovariance3(covariance, centroid, count, values, Vec3(1.0f));
  else
    ComputeWeightedCovariance3(covariance, centroid, count, values, Vec3(1.0f), weights);

  // compute the principle component
  GetPrincipleComponent(covariance, principle);

  // get the min and max range as the codebook endpoints
  Vec3 start(0.0f);
  Vec3 end(0.0f);

  if (count > 0) {
#ifdef	FEATURE_RANGEFIT_PROJECT
    Scr3 div = Reciprocal(Dot(principle, principle));
    Vec3 rec = Reciprocal(    principle            );
    Scr3 len, min, max;
    Vec3 chk;

    // compute the projection
    min = max = Dot(values[0] - centroid, principle);

    for (int i = 1; i < count; ++i) {
      len = Dot(values[i] - centroid, principle);
      min = Min(min, len);
      max = Max(max, len);
    }

    start = centroid + principle * min * div;
    end   = centroid + principle * max * div;

    // intersect negative undershoot with axis-plane(s), clamp to 0.0
    chk = start;
    while (CompareAnyLessThan(chk, Vec3(-1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hin = Select(fct, chk, HorizontalMin(chk));

      start -= principle * hin;
      chk = start;
    }

    // intersect negative undershoot with axis-plane(s), clamp to 0.0
    chk = end;
    while (CompareAnyLessThan(chk, Vec3(-1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hin = Select(fct, chk, HorizontalMin(chk));

      end -= principle * hin;
      chk = end;
    }

    // intersect positive overshoot with axis-plane(s), clamp to 1.0
    chk = start - Vec3(1.0f);
    while (CompareAnyGreaterThan(chk, Vec3(1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hax = Select(fct, chk, HorizontalMax(chk));

      start -= principle * hax;
      chk = start - Vec3(1.0f);
    }

    // intersect positive overshoot with axis-plane(s), clamp to 1.0
    chk = end - Vec3(1.0f);
    while (CompareAnyGreaterThan(chk, Vec3(1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hax = Select(fct, chk, HorizontalMax(chk));

      end -= principle * hax;
      chk = end - Vec3(1.0f);
    }

/*  assert(HorizontalMin(start).X() > -0.0001);
    assert(HorizontalMin(end  ).X() > -0.0001);
    assert(HorizontalMax(start).X() <  1.0001);
    assert(HorizontalMax(end  ).X() <  1.0001);	 */
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
  m_start = Truncate(start * 255.0f + Vec3(0.5f)) * (1.0f / 255.0f);
  m_end   = Truncate(end   * 255.0f + Vec3(0.5f)) * (1.0f / 255.0f);
}

void BitoneRangeFit::Compress4(void* block)
{
  // cache some values
  int const count = m_bitones->GetCount();
  Vec3 const* values = m_bitones->GetPoints();
  u8 const* freq = m_bitones->GetFrequencies();

  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[4]; Codebook4(codes, m_start, m_end);

  // match each point to the closest code
  u8 closest[16];

  Scr3 error = Scr3(0.0f);
  for (int i = 0; i < count; ++i) {
    // find the closest code
    Scr3 dist;
    Vec3 value = values[i];
    int idx = 0;

    {
      Scr3 d0 = LengthSquared(value - codes[0]);
      Scr3 d1 = LengthSquared(value - codes[1]);
      Scr3 d2 = LengthSquared(value - codes[2]);
      Scr3 d3 = LengthSquared(value - codes[3]);

      // encourage OoO
      Scr3 da = Min(d0, d1);
      Scr3 db = Min(d2, d3);
      dist    = Min(da, db);

      // will cause VS to make them all cmovs
      if (d3 == dist) { idx = 3; }
      if (d2 == dist) { idx = 2; }
      if (d1 == dist) { idx = 1; }
      if (d0 == dist) { idx = 0; }
    }

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
    u8 indices[16];

    m_bitones->RemapIndices(closest, indices);

    // save the block
    WriteBitoneBlock4(m_start, m_end, indices, block);
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
