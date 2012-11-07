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

#include <assert.h>
#include "colourset.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
ColourSet::ColourSet(u8 const* rgba, int mask, int flags)
  : m_count(0), m_unweighted(true), m_transparent(false)
{
  const float *rgbLUT = ComputeGammaLUT((flags & kSrgbIn) != 0);

  // check the compression mode for dxt1
  bool const isBtc1        = ((flags & kBtc1                   ) != 0);
  bool const clearAlpha    = ((flags & kExcludeAlphaFromPalette) != 0);
  bool const weightByAlpha = ((flags & kWeightColourByAlpha    ) != 0);

#ifdef	FEATURE_TEST_LINES
  Scr3 angl, epsl = Scr3(1.0f - (1.0f / 256));
  Vec3 chkl, line;
#endif

  // build mapped data
  u8 clra = clearAlpha || !isBtc1 ? 0xFF : 0x00;
  u8 wgta = weightByAlpha         ? 0x00 : 0xFF;

  u8 rgbx[4 * 16];
  u8 ___a[1 * 16];

  for (int i = 0; i < 16; ++i) {
    // clear alpha
    rgbx[4 * i + 0] = rgba[4 * i + 0];
    rgbx[4 * i + 1] = rgba[4 * i + 1];
    rgbx[4 * i + 2] = rgba[4 * i + 2];
    rgbx[4 * i + 3] = 0;

    // threshold alpha
    ___a[1 * i + 0] = (((signed char)rgba[4 * i + 3]) >> 7) | clra;

#ifdef FEATURE_IGNORE_ALPHA0
    // threshold color
    if (!(rgba[4 * i + 3] | wgta))
      ___a[1 * i + 0] = 0x00;
#endif
  }

  // create the minimal set
  for (int i = 0; i < 16; ++i) {
    // check this pixel is enabled
    int bit = 1 << i;
    if ((mask & bit) == 0) {
      m_remap[i] = -1;
      continue;
    }

    /* check for transparent pixels when using dxt1
     * check for blanked out pixels when weighting
     */
    if (!___a[i]) {
      m_remap[i] = -1;
      m_transparent = true;
      continue;
    }

    // ensure there is always non-zero weight even for zero alpha
    u8    w = rgba[4 * i + 3] | wgta;
    float W = (float)(w + 1) / 256.0f;

    // loop over previous points for a match
    u8 *rgbvalue = &rgbx[4 * i + 0];
    for (int j = 0;; ++j) {
      u8 *crgbvalue = &rgbx[4 * j + 0];

      // allocate a new point
      if (j == i) {
	// normalize coordinates to [0,1]
	float r = rgbLUT[rgbvalue[0]];
	float g = rgbLUT[rgbvalue[1]];
	float b = rgbLUT[rgbvalue[2]];

	// add the point
	m_remap[i] = m_count;
	m_points[m_count] = Vec3(r, g, b);
	m_weights[m_count] = W;
	m_unweighted = m_unweighted && !(u8)(~w);
#ifdef	FEATURE_EXACT_ERROR
	m_frequencies[m_count] = 1;
#endif

#ifdef	FEATURE_TEST_LINES
        // straight line test
	if (m_count >= 2) {
	  // (a/n * b/m + c/n * d/m + e/n * f/m)
	  chkl = Normalize(m_points[m_count - 1] - m_points[m_count - 2]);
	  angl = Abs(Dot(line, chkl));

	  m_straight = m_straight && (angl < epsl);
        }
	else if (m_count >= 1)
	  line = Normalize(m_points[m_count - 0] - m_points[m_count - 1]);
#endif

	// advance
	++m_count;
	break;
      }

      // check for a match
      int oldbit = 1 << j;
      bool match = ((mask & oldbit) != 0)
	&& (*((int *)rgbvalue) == *((int *)crgbvalue)) && ___a[j]/*
	&& (rgba[4 * i + 0] == rgba[4 * j + 0])
	&& (rgba[4 * i + 1] == rgba[4 * j + 1])
	&& (rgba[4 * i + 2] == rgba[4 * j + 2])
	&& (rgba[4 * j + 3] >= 128 || !isBtc1 || clearAlpha)*/;

      if (match) {
	// get the index of the match
	int const index = m_remap[j];
	assume (index >= 0 && index < 16);

	// map to this point and increase the weight
	m_remap[i] = index;
	m_weights[index] += W;
	m_unweighted = false;
#ifdef	FEATURE_EXACT_ERROR
	m_frequencies[index] += 1;
#endif
	break;
      }
    }
  }

#ifdef FEATURE_WEIGHTS_ROOTED
  // square root the weights
  for (int i = 0; i < m_count; ++i)
    m_weights[i] = math::sqrt(m_weights[i]);
#endif

  // clear if we're suppose to throw alway alpha
  m_transparent = m_transparent && !clearAlpha;

  // we have tables for this
  m_unweighted = m_unweighted && ((m_count == 16) || isBtc1);
}

void ColourSet::RemapIndices(u8 const* source, u8* target) const
{
  for (int i = 0; i < 16; ++i) {
    u8 t = 3; t = ((m_remap[i] == -1) ? t : source[m_remap[i]]); target[i] = t;
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
void ColourSet_CCR::CountSet(tile_barrier barrier, const int thread, pixel16 rgba, int mask, const bool tresh, const bool trans) amp_restricted
{
  // check the compression mode for dxt1
  const bool isBtc1 = (tresh);
  const bool weightByAlpha = (trans);

  // clear counters
  threaded_cse(0) {
    m_count = 0;
    m_transparent = 0;
  }

  wavefrnt_for(rscan, 16) {
    // clear out all values, so we may use all 16 of
    // them, instead of variable count at any time
    m_points [rscan] = 0.0f;
    m_weights[rscan] = 0;

    // use no code by default
    m_remap  [rscan] = -1;
  }

  // create the minimal set
  for (int i = 0; i < 16; ++i) {
    // check this pixel is enabled
    int bit = 1 << i;
    if ((mask & bit) != 0) {
      // check for transparent pixels when using dxt1
      if (isBtc1 && (rgba[i][0] < 128)) {
        // AMP: all thread end up here (CXCHG won't block concurrent reads)
        int zero = 0; threaded_set(m_transparent, zero, 1);
      }
      else {
        // make writes to "m_remap" visible to all
        tile_static_memory_fence(barrier);

        // loop over previous points for a match (AMP: parallel search)
        threaded_for(j, i) {
          // check for a match
          int oldbit = 1 << j;

          // get the index of the match
          bool match = ((mask & oldbit) != 0)
	    && (rgba[i][3] == rgba[j][3])
	    && (rgba[i][2] == rgba[j][2])
	    && (rgba[i][1] == rgba[j][1])
	    && (!isBtc1 || (rgba[j][0] >= 128));

          // AMP: there is only one thread max. which could find a match
          //      so it's only one thread which writes to index, no
          //      atomic needed
          if (match)
	    m_remap[i] = m_remap[j];
        }

        // make writes to "m_remap" visible to thread 0
        tile_static_memory_fence(barrier);

        // allocate a new point
        threaded_cse(0) {
          int index = m_remap[i];
          if (index < 0) {
	    // get the index of the new point
	    index = m_count++;

	    // add the point
	    m_weights[index] = 0;
	    m_points [index] = float3(
	      (float)rgba[i][3],
	      (float)rgba[i][2],
	      (float)rgba[i][1]
	    );

	    m_remap[i] = index;
          }

          // ensure there is always non-zero weight even for zero alpha
          float w = (float)(rgba[i][0] + 1) / 256.0f;

          // map to this point and increase the weight
          m_weights[index] += (weightByAlpha ? w : 1.0f);
        }
      }
    }
  }

  // square root the weights
#if	!defined(SQUISH_USE_COMPUTE)
  wavefrnt_for(wscan, 16) {
    // normalize coordinates to [0,1]
    m_points [wscan] /= 255.0f;
    m_weights[wscan]  = sqrtf(m_weights[wscan]);
  }
#else
  // error X3695: race condition writing to shared memory detected, consider making this write conditional.
  threaded_for(wscan, m_count) {
    // normalize coordinates to [0,1]
    m_points [wscan] /= 255.0f;
    m_weights[wscan]  = sqrtf(m_weights[wscan]);
  }
#endif
}

int ColourSet_CCR::GetCount() amp_restricted {
  return m_count;
}

point16 ColourSet_CCR::GetPoints() amp_restricted {
  return m_points;
}

weight16 ColourSet_CCR::GetWeights() amp_restricted {
  return m_weights;
}

bool ColourSet_CCR::IsTransparent() amp_restricted {
  return !!m_transparent;
}

void ColourSet_CCR::WriteSet(tile_barrier barrier, const int thread, lineC2 cline, inout index16x2 source, out code64 block, const int is4,
			     IndexBlockLUT yArr) amp_restricted
{
  // build the block if we win
  {
    // remap the indices
    RemapIndices(barrier, thread, source);

    // save the block
    {
      if (!is4)
	WriteColourBlock3(barrier, thread, cline, m_indices, block, yArr);
      else
	WriteColourBlock4(barrier, thread, cline, m_indices, block, yArr);
    }
  }
}

void ColourSet_CCR::RemapIndices(tile_barrier barrier, const int thread, inout index16x2 source) amp_restricted
{
  // make writes to "source"/"m_remap" visible to all
  tile_static_memory_fence(barrier);

#if	!defined(SQUISH_USE_COMPUTE)
  // remap to palette-indices
  wavefrnt_for(rmscan, 16) {
    m_indices[rmscan] = source[1][m_remap[rmscan]];
  }
#else
  // error X3500: array reference cannot be used as an l-value not natively addressable
  threaded_cse(0) {
    m_indices[ 0] = source[1][m_remap[ 0]];
    m_indices[ 1] = source[1][m_remap[ 1]];
    m_indices[ 2] = source[1][m_remap[ 2]];
    m_indices[ 3] = source[1][m_remap[ 3]];
    m_indices[ 4] = source[1][m_remap[ 4]];
    m_indices[ 5] = source[1][m_remap[ 5]];
    m_indices[ 6] = source[1][m_remap[ 6]];
    m_indices[ 7] = source[1][m_remap[ 7]];
    m_indices[ 8] = source[1][m_remap[ 8]];
    m_indices[ 9] = source[1][m_remap[ 9]];
    m_indices[10] = source[1][m_remap[10]];
    m_indices[11] = source[1][m_remap[11]];
    m_indices[12] = source[1][m_remap[12]];
    m_indices[13] = source[1][m_remap[13]];
    m_indices[14] = source[1][m_remap[14]];
    m_indices[15] = source[1][m_remap[15]];
  }
#endif
}
#endif

} // namespace squish
