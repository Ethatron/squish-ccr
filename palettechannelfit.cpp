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

#include "palettechannelfit.h"
#include "paletteset.h"
#include "paletteblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
PaletteChannelFit::PaletteChannelFit(PaletteSet const* palette, int flags, int swap, int shared)
  : PaletteFit(palette, flags, swap, shared)
{
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;
  
  assume((isets >  0) && (isets <= 3));
  assume((asets >= 0) && (asets <= 3));
  assume(((isets    +    asets) <= 3));

  for (int a = isets; a < (isets + asets); a++) {
    // cache some values
    int const count = m_palette->GetCount(a);
    Vec4 const* values = m_palette->GetPoints(a);

    // we don't do this for sparse sets
    if (count != 1) {
      // get the min and max range as the codebook endpoints
      Vec4 start(1.0f);
      Vec4 end(1.0f);

      if (count > 0) {
	// compute the range
	start = end = values[0];

	for (int i = 1; i < count; ++i) {
	  start = Min(start, values[i]);
	  end   = Max(end  , values[i]);
	}
      }

      // clamp the output to [0, 1]
      m_start_candidate[a] = Min(start, end);
      m_end_candidate  [a] = Max(start, end);
    }
  }
}

Scr4 PaletteChannelFit::ComputeCodebook(int set, Vec4 const &metric, vQuantizer &q, int sb, int ib, u8 (&closest)[16])
{
  // cache some values
  int const count = m_palette->GetCount(set);
  Vec4 const* values = m_palette->GetPoints(set);
  u8 const* freq = m_palette->GetFrequencies(set);

  // snap floating-point-values to the integer-lattice
  Vec4 cstart = q.SnapToLattice(m_start_candidate[set], sb, 1 << SBSTART);
  Vec4 cend   = q.SnapToLattice(m_end_candidate  [set], sb, 1 << SBEND);

  // the lattice to walk over
  Vec4 gridrcp = Reciprocal(q.grid + Vec4(1.0f));
  Vec4 wstart;
  Vec4 wend;
  Vec4 wdelta = gridrcp; if (sb) wdelta *= Vec4(2.0f);

  // lower the lower by 1 and raise the higher by 1
  // compensates a bit the rounding error of the end-points
  cstart = Max(cstart - wdelta, Vec4(0.0f));
  cend   = Min(cend   + wdelta, Vec4(1.0f));

  // create a codebook
  Vec4 codes[1 << 4];

  Scr4 besterror = Scr4(FLT_MAX);
  Vec4 beststart = cstart;
  Vec4 bestend = cend;

  // Brute force approach, try all the possible endpoints with g0 > g1.
  wend = cstart + wdelta;
  while (!CompareFirstGreaterThan(wend, cend)) {
    wstart = cstart;
    while (!CompareFirstGreaterThan(wstart, wend)) {
      // resolve "metric * (value - code)" to "metric * value - metric * code"
      int ccs = CodebookP(codes, ib, metric * wstart, metric * wend);

      Scr4 error = Scr4(0.0f);
      for (int i = 0; i < count; ++i) {
	Scr4 dist = Scr4(FLT_MAX);
	Vec4 value = metric * values[i];
	
	for (int j = 0; j < ccs; j += 4) {
	  Scr4 d0 = LengthSquared(value - codes[j + 0]);
	  Scr4 d1 = LengthSquared(value - codes[j + 1]);
	  Scr4 d2 = LengthSquared(value - codes[j + 2]);
	  Scr4 d3 = LengthSquared(value - codes[j + 3]);

	  // encourage OoO
	  Scr4 da = Min(d0, d1);
	  Scr4 db = Min(d2, d3);
	  dist = Min(da, dist);
	  dist = Min(db, dist);
	}

	// accumulate the error
	error += dist * freq[i];
      }

      if (besterror > error) {
	besterror = error;
	beststart = wstart;
	bestend   = wend;
      }

      wstart += wdelta;
    }

    wend += wdelta;
  }

  // snap floating-point-values to the integer-lattice with up/down skew
  m_start[set] = beststart;
  m_end  [set] = bestend;
  
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  int ccs = CodebookP(codes, ib, metric * beststart, metric * bestend);

  for (int i = 0; i < count; ++i) {
    // find the closest code
    Scr4 dist = Scr4(FLT_MAX);
    Vec4 value = metric * values[i];
    int idx = 0;
    
    for (int j = 0; j < ccs; j += 0) {
      Scr4 d0 = LengthSquared(value - codes[j + 0]);
      Scr4 d1 = LengthSquared(value - codes[j + 1]);
      Scr4 d2 = LengthSquared(value - codes[j + 2]);
      Scr4 d3 = LengthSquared(value - codes[j + 3]);

      // encourage OoO
      Scr4 da = Min(d0, d1);
      Scr4 db = Min(d2, d3);
      dist = Min(da, dist);
      dist = Min(db, dist);

      // will cause VS to make them all cmovs
      if (d0 == dist) { idx = j; } j++;
      if (d1 == dist) { idx = j; } j++;
      if (d2 == dist) { idx = j; } j++;
      if (d3 == dist) { idx = j; } j++;
    }

    // save the index
    closest[i] = (u8)idx;
  }

  return besterror;
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
