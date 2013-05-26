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
  
#ifdef	FEATURE_TEST_LINES
  for (int s =     0; s < (isets + asets); s++) {
#else
  for (int s = isets; s < (isets + asets); s++) {
#endif
    // cache some values
    int const channel = m_palette->GetChannel(s);
    int const count = m_palette->GetCount(s);
    Vec4 const* values = m_palette->GetPoints(s);

    // we don't do this for sparse sets
#ifdef	FEATURE_TEST_LINES
    if ((count != 1) && (channel >= 0)) {
#else
    if ((count != 1)) {
#endif
      // the codebook endpoints
      Vec4 start(1.0f);
      Vec4 end(0.0f);
      
      // get the min and max range as the codebook endpoints
      // compute the range
      for (int i = 0; i < count; ++i) {
	start = Min(start, values[i]);
	end   = Max(end  , values[i]);
      }

      // clamp the output to [0, 1]
      m_start_candidate[s] = Min(start, end);
      m_end_candidate  [s] = Max(start, end);
    }
    
    // store the channel
    m_channel[s] = channel;
  }
}

Scr4 PaletteChannelFit::ComputeCodebook(int set, Vec4 const &metric, vQuantizer &q, int sb, int ib, u8 (&closest)[16])
{
  // cache some values
  int const count = m_palette->GetCount(set);
  Vec4 const* values = m_palette->GetPoints(set);
  Scr4 const* freq = m_palette->GetWeights(set);

  // snap floating-point-values to the integer-lattice
  Vec4 cstart = q.SnapToLatticeClamped(m_start_candidate[set], sb, 1 << SBSTART);
  Vec4 cend   = q.SnapToLatticeClamped(m_end_candidate  [set], sb, 1 << SBEND  );

  // the lattice to walk over
  Vec4 gridrcp = Reciprocal(q.grid + Vec4(1.0f));
  Vec4 wstart;
  Vec4 wend;
  Vec4 wdelta = gridrcp;
  Vec4 wmask;

  switch (m_channel[set]) {
    // select the only non-constant channel
    case  0: wmask = Vec4(true, false, false, false); break;
    case  1: wmask = Vec4(false, true, false, false); break;
    case  2: wmask = Vec4(false, false, true, false); break;
    case  3: wmask = Vec4(false, false, false, true); break;
    // select all grey non-constant channels
    case  4: wmask = Vec4(false, true , true , true); break;
    case  5: wmask = Vec4(true , false, true , true); break;
    case  6: wmask = Vec4(true , true , false, true); break;
    case  7: wmask = Vec4(true , true , true, false); break;
    // all channels were the same
    case  8: wmask = Vec4(true , true , true , true); break;
    // select all grey non-constant channels
    case  9: wmask = Vec4(false, false, true , true); break;
    case 10: wmask = Vec4(true , false, false, true); break;
    case 11: wmask = Vec4(true , true, false, false); break;
    case 12: wmask = Vec4(false, true , true, false); break;
  }

  if (~sb) wdelta *= Vec4(2.0f);
  wdelta = wdelta & wmask;

  // lower the lower by 1 and raise the higher by 1
  // compensates a bit the rounding error of the end-points
  cstart = Max(cstart - wdelta, Vec4(0.0f));
  cend   = Min(cend   + wdelta, Vec4(1.0f));
  
  // create a codebook
  Vec4 codes[1 << 4];

  Scr4 besterror = Scr4(FLT_MAX);
  Vec4 beststart = q.SnapToLatticeClamped(cstart, sb, 1 << SBSTART);
  Vec4 bestend   = q.SnapToLatticeClamped(cend  , sb, 1 << SBEND  );

  // Brute force approach, try all the possible endpoints with g0 > g1.
  wend = cstart + wdelta;
  while (!CompareAnyGreaterThan(wend, cend)) {
    wstart = cstart;
    while (!CompareAnyGreaterThan(wstart, wend)) {
      Vec4 vstart = q.SnapToLatticeClamped(wstart, sb, 1 << SBSTART);
      Vec4 vend   = q.SnapToLatticeClamped(wend  , sb, 1 << SBEND  );

      // resolve "metric * (value - code)" to "metric * value - metric * code"
      int ccs = CodebookP(codes, ib, metric * vstart, metric * vend);

      Scr4 error = Scr4(DISTANCE_BASE);
      for (int i = 0; i < count; ++i) {
	Scr4 dist = Scr4(FLT_MAX);
	Vec4 value = metric * values[i];
	
	for (int j = 0; j < ccs; j += 4)
	  MinDistance4<false>(dist, i, value, codes, j);

	// accumulate the error
	AddDistance(dist, error, freq[i]);
      }

      if (besterror > error) {
	besterror = error;
	beststart = vstart;
	bestend   = vend;
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
    
    for (int j = 0; j < ccs; j += 0)
      MinDistance4<true>(dist, idx, value, codes, j);

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
