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
#include "helpers.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
ColourSet::ColourSet(u8 const* rgba, int mask, int flags)
  : m_count(0)
  , m_unweighted(true)
  , m_transparent(false)
{
  const float *rgbLUT = ComputeGammaLUT((flags & kSrgbExternal) != 0);

  // check the compression mode for dxt1
  bool const isBtc1        = ((flags & kBtcp                   ) == kBtc1);
  bool const clearAlpha    = ((flags & kExcludeAlphaFromPalette) != 0);
  bool const weightByAlpha = ((flags & kWeightColourByAlpha    ) != 0);

  // build mapped data
  int clra = clearAlpha || !isBtc1 ? 0xFFFF : 0x0000;
  int wgta = weightByAlpha         ? 0x0000 : 0xFFFF;

  a16 u8 rgbx[4 * 16];
  int amask;

  Col4 m0 = Col4(&rgba[0 * 16]);
  Col4 m1 = Col4(&rgba[1 * 16]);
  Col4 m2 = Col4(&rgba[2 * 16]);
  Col4 m3 = Col4(&rgba[3 * 16]);
  Col4 al = CollapseA(m0, m1, m2, m3);

  // clear alpha
  m0 &= Col4(0x00FFFFFF);
  m1 &= Col4(0x00FFFFFF);
  m2 &= Col4(0x00FFFFFF);
  m3 &= Col4(0x00FFFFFF);
  
  StoreAligned(m0, &rgbx[0 * 16]);
  StoreAligned(m1, &rgbx[1 * 16]);
  StoreAligned(m2, &rgbx[2 * 16]);
  StoreAligned(m3, &rgbx[3 * 16]);
  
  // threshold alpha (signed char)
  amask  =   CompareAllLessThan_M8(al, Col4(0)).GetM8() | (clra);
#ifdef FEATURE_IGNORE_ALPHA0
  // threshold color
  amask &= (~CompareAllEqualTo_M8(al, Col4(0)).GetM8()) | (wgta);
#endif
  // combined mask
  amask &= mask;

  // create the minimal set, O(16*count/2)
  for (int i = 0, imask = amask, index; i < 16; ++i, imask >>= 1) {
    // check this pixel is enabled
    if ((imask & 1) == 0) {
      m_remap[i] = -1;

      /* check for transparent pixels when using dxt1
       * check for blanked out pixels when weighting
       */
      m_transparent = m_transparent | ((mask & (1 << i)) != 0);

      continue;
    }
    
    // calculate point's weights
    Weight<u8> wa(rgba, i, (u8)wgta);

    // loop over previous matches for a match
    u8 *rgbvalue = &rgbx[4 * i + 0];
    for (index = 0; index < m_count; ++index) {
      u8 *crgbvalue = &rgbx[4 * index + 0];
      
      // check for a match
      if (*((int *)rgbvalue) == *((int *)crgbvalue)) {
	// get the index of the match
	assume (index >= 0 && index < 16);

	// map to this point and increase the weight
	m_remap[i] = (char)index;
	m_weights[index] += wa.GetWeights();
	m_unweighted = false;
	break;
      }
    }
    
    // re-use the memory of already processed pixels
    assert(index <= i); {
      u8 *crgbvalue = &rgbx[4 * index + 0];

      // allocate a new point
      if (index == m_count) {
	// get the index of the match and advance
	m_count = index + 1;

	// normalize coordinates to [0,1]
	const float *r = &rgbLUT[rgbvalue[0]];
	const float *g = &rgbLUT[rgbvalue[1]];
	const float *b = &rgbLUT[rgbvalue[2]];

	// add the point
	m_remap[i] = (char)index;
	m_points[index] = Vec3(r, g, b);
	m_weights[index] = wa.GetWeights();
	m_unweighted = m_unweighted & wa.IsOne();
	
	// remember match for successive checks
	*((int *)crgbvalue) = *((int *)rgbvalue);
      }
    }
  }

#ifdef FEATURE_IGNORE_ALPHA0
  if ((clra == 0xFFFF) && (amask != 0xFFFF)) {
    if (!m_count) {
      Vec3 sum = Vec3(0.0f);

      for (int i = 0, imask = amask; i < 16; ++i, imask >>= 1) {
	/* assign blanked out pixels when weighting
	 */
	if ((imask & 1) == 0) {
	  m_remap[i] = 0;

	  u8 *rgbvalue = &rgbx[4 * i + 0];

	  // normalize coordinates to [0,1]
	  const float *r = &rgbLUT[rgbvalue[0]];
	  const float *g = &rgbLUT[rgbvalue[1]];
	  const float *b = &rgbLUT[rgbvalue[2]];

	  sum += Vec3(r, g, b);
	}
      }

      // add the point
      m_count = 1;
      m_points[0] = sum * (1.0f / 16.0f);
      m_weights[0] = Scr3(1.0f);
      m_unweighted = true;
    }
    else if (m_transparent) {
      for (int i = 0, imask = amask, index; i < 16; ++i, imask >>= 1) {
	/* assign blanked out pixels when weighting
	 */
	if ((imask & 1) == 0) {
	  u8 *rgbvalue = &rgbx[4 * i + 0];

	  // normalize coordinates to [0,1]
	  const float *r = &rgbLUT[rgbvalue[0]];
	  const float *g = &rgbLUT[rgbvalue[1]];
	  const float *b = &rgbLUT[rgbvalue[2]];

	  // loop over previous matches for a match
	  Scr3 d = Scr3(FLT_MAX);
	  for (index = 0; index < m_count; ++index) {
	    Vec3 diff = m_points[index] - Vec3(r, g, b);
	    Scr3 dist = Dot(diff, diff);

	    if (d > dist) {
	      d = dist;

	      m_remap[i] = (char)index;
	    }
	  }
	}
      }
    }
  }
#endif

#ifdef FEATURE_WEIGHTS_ROOTED
  // square root the weights
  for (int i = 0; i < m_count; ++i)
    m_weights[i] = math::sqrt(m_weights[i]);
#endif

  // clear if we're suppose to throw alway alpha
  m_transparent = m_transparent & !clearAlpha;
}

ColourSet::ColourSet(u16 const* rgba, int mask, int flags)
  : m_count(0)
  , m_unweighted(true)
  , m_transparent(false)
{
}

ColourSet::ColourSet(f23 const* rgba, int mask, int flags)
  : m_count(0)
  , m_unweighted(true)
  , m_transparent(false)
{
//const float *rgbLUT = ComputeGammaLUT((flags & kSrgbIn) != 0);

  // check the compression mode for dxt1
  bool const isBtc1        = ((flags & kBtcp                   ) == kBtc1);
  bool const clearAlpha    = ((flags & kExcludeAlphaFromPalette) != 0);
  bool const weightByAlpha = ((flags & kWeightColourByAlpha    ) != 0);

  // build mapped data
  Scr4 clra = clearAlpha || !isBtc1 ? Scr4(1.0f) : Scr4(0.0f);
  Scr4 wgta = weightByAlpha         ? Scr4(0.0f) : Scr4(1.0f);
  
  Vec3 rgbx[16];
  int amask = 0;

  for (int i = 0; i < 16; i += 4) {
    Vec4 m0; LoadUnaligned(m0, &rgba[4 * i + 4 * 0]);
    Vec4 m1; LoadUnaligned(m1, &rgba[4 * i + 4 * 1]);
    Vec4 m2; LoadUnaligned(m2, &rgba[4 * i + 4 * 2]);
    Vec4 m3; LoadUnaligned(m3, &rgba[4 * i + 4 * 3]);
    Vec4 al = CollapseW(m0, m1, m2, m3);
    Vec4 ibit;

    // clear alpha
    rgbx[i + 0] = KillW(m0).GetVec3();
    rgbx[i + 1] = KillW(m1).GetVec3();
    rgbx[i + 2] = KillW(m2).GetVec3();
    rgbx[i + 3] = KillW(m3).GetVec3();

    // threshold alpha
    ibit = IsGreaterEqual(Max(al, clra), Vec4(0.5f));
#ifdef FEATURE_IGNORE_ALPHA0
    // threshold color
    ibit &=  IsNotEqualTo(Max(al, wgta), Vec4(0.0f));
#endif

    amask += ibit.GetM4() << i;
  }
  
  // combined mask
  amask &= mask;

  // create the minimal set, O(16*count/2)
  for (int i = 0, imask = amask, index; i < 16; ++i, imask >>= 1) {
    // check this pixel is enabled
    if ((imask & 1) == 0) {
      m_remap[i] = -1;

      /* check for transparent pixels when using dxt1
       * check for blanked out pixels when weighting
       */
      m_transparent = m_transparent | ((mask & (1 << i)) != 0);

      continue;
    }
    
    // calculate point's weights
    Weight<f23> wa(rgba, i, wgta);

    // loop over previous matches for a match
    Vec3 *rgbvalue = &rgbx[i];
    for (index = 0; index < m_count; ++index) {
      Vec3 *crgbvalue = &rgbx[index];
      
      // check for a match
      if (CompareAllEqualTo((*rgbvalue), (*crgbvalue))) {
	// get the index of the match
	assume (index >= 0 && index < 16);

	// map to this point and increase the weight
	m_remap[i] = (char)index;
	m_weights[index] += wa.GetWeights();
	m_unweighted = false;
	break;
      }
    }
    
    // re-use the memory of already processed pixels
    assert(index <= i); {
      Vec3 *crgbvalue = &rgbx[index];

      // allocate a new point
      if (index == m_count) {
	// get the index of the match and advance
	m_count = index + 1;

#if 0
	// normalize coordinates to [0,1]
	const float *r = &rgbLUT[rgbvalue[0]];
	const float *g = &rgbLUT[rgbvalue[1]];
	const float *b = &rgbLUT[rgbvalue[2]];
#endif

	// add the point
	m_remap[i] = (char)index;
	m_points[index] = *(rgbvalue);
	m_weights[index] = wa.GetWeights();
	m_unweighted = m_unweighted & wa.IsOne();
	
	// remember match for successive checks
	*(crgbvalue) = *(rgbvalue);
      }
    }
  }

#ifdef FEATURE_IGNORE_ALPHA0
  if ((clra == Scr4(1.0f)) && (amask != 0xFFFF)) {
    if (!m_count) {
      Vec3 sum = Vec3(0.0f);

      for (int i = 0, imask = amask; i < 16; ++i, imask >>= 1) {
	/* assign blanked out pixels when weighting
	 */
	if ((imask & 1) == 0) {
	  m_remap[i] = 0;

	  Vec3 *rgbvalue = &rgbx[i];
	  
#if 0
	  // normalize coordinates to [0,1]
	  const float *r = &rgbLUT[rgbvalue[0]];
	  const float *g = &rgbLUT[rgbvalue[1]];
	  const float *b = &rgbLUT[rgbvalue[2]];
#endif

	  sum += *(rgbvalue);
	}
      }

      // add the point
      m_count = 1;
      m_points[0] = sum * (1.0f / 16.0f);
      m_weights[0] = Scr3(1.0f);
      m_unweighted = true;
    }
    else if (m_transparent) {
      for (int i = 0, imask = amask, index; i < 16; ++i, imask >>= 1) {
	/* assign blanked out pixels when weighting
	 */
	if ((imask & 1) == 0) {
	  Vec3 *rgbvalue = &rgbx[i];
	  
#if 0
	  // normalize coordinates to [0,1]
	  const float *r = &rgbLUT[rgbvalue[0]];
	  const float *g = &rgbLUT[rgbvalue[1]];
	  const float *b = &rgbLUT[rgbvalue[2]];
#endif

	  // loop over previous matches for a match
	  Scr3 d = Scr3(FLT_MAX);
	  for (index = 0; index < m_count; ++index) {
	    Vec3 diff = m_points[index] - *(rgbvalue);
	    Scr3 dist = Dot(diff, diff);

	    if (d > dist) {
	      d = dist;

	      m_remap[i] = (char)index;
	    }
	  }
	}
      }
    }
  }
#endif

#ifdef FEATURE_WEIGHTS_ROOTED
  // square root the weights
  for (int i = 0; i < m_count; ++i)
    m_weights[i] = math::sqrt(m_weights[i]);
#endif

  // clear if we're suppose to throw alway alpha
  m_transparent = m_transparent & !clearAlpha;
}

bool ColourSet::RemoveBlack(const Vec3 &metric, Scr3 &error)
{
  cQuantizer4<5,6,5,0> q = cQuantizer4<5,6,5,0>();
  bool reduced = false;
  
  while (m_count > 1) {
    Scr3 lowest = LengthSquared(metric * Vec3(32.0f / 255.0f));
    int index = -1;

    for (int i = 0; i < 16; ++i) {
      if (m_remap[i] == -1)
	continue;

      // maps to black
      Vec3 colour = m_points[m_remap[i]];
      /*Vec3 result = q.SnapToLattice(colour);*/
      if (true /*CompareAllEqualTo(result, Vec3(0.0f))*/) {
	Scr3 len = LengthSquared(metric * colour);
	if (len < lowest) {
	  lowest = len;
	  index = m_remap[i];
	}
      }
    }

    if (index >= 0) {
      m_count--;
      m_unweighted = false;

      for (int i = 0; i < 16; ++i) {
	if (m_remap[i] == index)
	  m_remap[i] = -1, m_transparent = true;
	else if (m_remap[i] > index)
	  m_remap[i] = -1 + m_remap[i];
      }

      error += LengthSquared(metric * m_points[index]) * m_weights[index];
      
      if (m_count > index) {
	memcpy(&m_points [index + 0], &m_points [index + 1], sizeof(m_points [0]) * (m_count - index));
	memcpy(&m_weights[index + 0], &m_weights[index + 1], sizeof(m_weights[0]) * (m_count - index));
      }

      reduced = true;
      return reduced;
    }
    else
      break;
  }

  return reduced;
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
