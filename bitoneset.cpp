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
#include "bitoneset.h"
#include "helpers.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
BitoneSet::BitoneSet(u8 const* rgba, int mask, int flags)
  : m_count(0)
  , m_unweighted(true)
{
  const float *rgbLUT = ComputeGammaLUT((flags & kSrgbExternal) != 0);
  
  bool const preserveThird = ((flags & kColourMetricUnit   ) != 0);
  bool const weightByAlpha = ((flags & kWeightColourByAlpha) != 0);

  // build mapped data
  Col4 kill = preserveThird ? Col4(0x00FFFFFF) : Col4(0x0000FFFF);
  int wgta = weightByAlpha ? 0x0000 : 0xFFFF;

  a16 u8 rgbx[4 * 16];
  int amask;

  Col4 m0 = Col4(&rgba[0 * 16]);
  Col4 m1 = Col4(&rgba[1 * 16]);
  Col4 m2 = Col4(&rgba[2 * 16]);
  Col4 m3 = Col4(&rgba[3 * 16]);
  Col4 al = CollapseA(m0, m1, m2, m3);
  
  // clear alpha and maybe b
  m0 &= kill;
  m1 &= kill;
  m2 &= kill;
  m3 &= kill;
  
  StoreAligned(m0, &rgbx[0 * 16]);
  StoreAligned(m1, &rgbx[1 * 16]);
  StoreAligned(m2, &rgbx[2 * 16]);
  StoreAligned(m3, &rgbx[3 * 16]);
  
  // combined mask
  amask  = mask;
#ifdef FEATURE_IGNORE_ALPHA0
  // threshold color
  amask &= (~CompareAllEqualTo_M8(al, Col4(0)).GetM8()) | (wgta);
#endif

  // create the minimal set, O(16*count/2)
  for (int i = 0, imask = amask, index; i < 16; ++i, imask >>= 1) {
    // check this pixel is enabled
    if ((imask & 1) == 0) {
      m_remap[i] = -1;
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
  if (amask != 0xFFFF) {
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
    else {
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
}

BitoneSet::BitoneSet(u16 const* rgba, int mask, int flags)
  : m_count(0)
  , m_unweighted(true)
{
}

BitoneSet::BitoneSet(f23 const* rgba, int mask, int flags)
  : m_count(0)
  , m_unweighted(true)
{
//const float *rgbLUT = ComputeGammaLUT((flags & kSrgbIn) != 0);
  
  bool const preserveThird = ((flags & kColourMetricUnit   ) != 0);
  bool const weightByAlpha = ((flags & kWeightColourByAlpha) != 0);

  // build mapped data
  Vec4 kill = preserveThird ? Vec4(true, true, true, false) : Vec4(true, true, false, false);
  Scr4 wgta = weightByAlpha ? Scr4(0.0f) : Scr4(1.0f);

  Vec3 rgbx[16];
  int amask = 0;
  
  for (int i = 0; i < 16; i += 4) {
    Vec4 m0; LoadUnaligned(m0, &rgba[4 * i + 4 * 0]);
    Vec4 m1; LoadUnaligned(m1, &rgba[4 * i + 4 * 1]);
    Vec4 m2; LoadUnaligned(m2, &rgba[4 * i + 4 * 2]);
    Vec4 m3; LoadUnaligned(m3, &rgba[4 * i + 4 * 3]);
    Vec4 al = CollapseW(m0, m1, m2, m3);
    Vec4 ibit;

    // clear alpha and maybe b
    rgbx[i + 0] = (m0 & kill).GetVec3();
    rgbx[i + 1] = (m1 & kill).GetVec3();
    rgbx[i + 2] = (m2 & kill).GetVec3();
    rgbx[i + 3] = (m3 & kill).GetVec3();

#ifdef FEATURE_IGNORE_ALPHA0
    // threshold color
    ibit = IsNotEqualTo(Max(al, wgta), Vec4(0.0f));
#else
    ibit = 0xF;
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
  if (amask != 0xFFFF) {
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
    else {
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
}

void BitoneSet::RemapIndices(u8 const* source, u8* target) const
{
  for (int i = 0; i < 16; ++i) {
    u8 t = 3; t = ((m_remap[i] == -1) ? t : source[m_remap[i]]); target[i] = t;
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
