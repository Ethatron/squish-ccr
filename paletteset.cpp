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

#include "paletteset.h"

namespace squish {

// Associated to partition -1, 16 * 0 bit
static const u16 partitionmasks_1[1] =
{
  0xFFFF
};

// Associated to partition 0-63, 16 * 1 bit
static const u16 partitionmasks_2[64] =
{
  0xCCCC, 0x8888, 0xEEEE, 0xECC8,
  0xC880, 0xFEEC, 0xFEC8, 0xEC80,
  0xC800, 0xFFEC, 0xFE80, 0xE800,
  0xFFE8, 0xFF00, 0xFFF0, 0xF000,
  0xF710, 0x008E, 0x7100, 0x08CE,
  0x008C, 0x7310, 0x3100, 0x8CCE,
  0x088C, 0x3110, 0x6666, 0x366C,
  0x17E8, 0x0FF0, 0x718E, 0x399C,

  0xaaaa, 0xf0f0, 0x5a5a, 0x33cc,
  0x3c3c, 0x55aa, 0x9696, 0xa55a,
  0x73ce, 0x13c8, 0x324c, 0x3bdc,
  0x6996, 0xc33c, 0x9966, 0x0660,
  0x0272, 0x04e4, 0x4e40, 0x2720,
  0xc936, 0x936c, 0x39c6, 0x639c,
  0x9336, 0x9cc6, 0x817e, 0xe718,
  0xccf0, 0x0fcc, 0x7744, 0xee22,
};

// Associated to partition 64-127, 16 * 2 bit
static const unsigned int partitionmasks_3[64] =
{
  0xf60008cc, 0x73008cc8, 0x3310cc80, 0x00ceec00,
  0xcc003300, 0xcc0000cc, 0x00ccff00, 0x3300cccc,
  0xf0000f00, 0xf0000ff0, 0xff0000f0, 0x88884444,
  0x88886666, 0xcccc2222, 0xec80136c, 0x7310008c,
  0xc80036c8, 0x310008ce, 0xccc03330, 0x0cccf000,
  0xee0000ee, 0x77008888, 0xcc0022c0, 0x33004430,
  0x00cc0c22, 0xfc880344, 0x06606996, 0x66009960,
  0xc88c0330, 0xf9000066, 0x0cc0c22c, 0x73108c00,

  0xec801300, 0x08cec400, 0xec80004c, 0x44442222,
  0x0f0000f0, 0x49242492, 0x42942942, 0x0c30c30c,
  0x03c0c03c, 0xff0000aa, 0x5500aa00, 0xcccc3030,
  0x0c0cc0c0, 0x66669090, 0x0ff0a00a, 0x5550aaa0,
  0xf0000aaa, 0x0e0ee0e0, 0x88887070, 0x99906660,
  0xe00e0ee0, 0x88880770, 0xf0000666, 0x99006600,
  0xff000066, 0xc00c0cc0, 0xcccc0330, 0x90006000,
  0x08088080, 0xeeee1010, 0xfff0000a, 0x731008ce,
};

/* *****************************************************************************
 */
#if	!defined(USE_PRE)
void PaletteSet::GetMasks(int flags, int partition, int (&masks)[4]) {
  unsigned int partmask = 0;
  if (((flags & kVariableCodingModes) == kVariableCodingMode2) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode4) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode8))
    partmask = partitionmasks_2[partition];
  if (((flags & kVariableCodingModes) == kVariableCodingMode1) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode3))
    partmask = partitionmasks_3[partition];

  // determine the number of partitions
  if (((flags & kVariableCodingModes) == kVariableCodingMode5) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode6) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode7))
    masks[0] = ( partmask   & 0xFFFF) & ( 0xFFFFFFFF >> 16),	// color-set
    masks[1] = ( partmask   & 0xFFFF) & ( 0xFFFFFFFF >> 16),	// alpha-set
    masks[2] = 0;
  if (((flags & kVariableCodingModes) == kVariableCodingMode1) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode3))
    masks[0] = (~partmask   & 0xFFFF) & ( 0xFFFFFFFF >> 16),
    masks[1] = ( partmask   & 0xFFFF) & ( 0xFFFFFFFF >> 16),
    masks[2] = 0;
  if (((flags & kVariableCodingModes) == kVariableCodingMode2) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode4) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode8))
    masks[0] = (~partmask   & 0xFFFF) & (~partmask >> 16),
    masks[1] = ( partmask   & 0xFFFF) & (~partmask >> 16),
    masks[2] = ( 0xFFFFFFFF & 0xFFFF) & ( partmask >> 16);
}

PaletteSet::PaletteSet(u8 const* rgba, int mask, int flags, int partition, int rotation)
  : m_numsets(1), m_rotid(0), m_partid(0), m_partmask(0xFFFF), m_seperatealpha(false), m_transparent(false)
{
  m_count[0] = m_count[1] = m_count[2] = m_count[3] = 0;

  /* determine the number of sets and select the partition
  if ((0))
    m_numsets = 1, m_partmask = partitionmasks_1[m_partid = 0]; */
  if (((flags & kVariableCodingModes) == kVariableCodingMode2) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode4) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode8))
    m_numsets = 2, m_partmask = partitionmasks_2[m_partid = partition];
  if (((flags & kVariableCodingModes) == kVariableCodingMode1) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode3))
    m_numsets = 3, m_partmask = partitionmasks_3[m_partid = partition];
  if (((flags & kVariableCodingModes) == kVariableCodingMode5) ||
      ((flags & kVariableCodingModes) == kVariableCodingMode6))
    m_seperatealpha = true, m_rotid = rotation;

  // partition_1 mask is: bit cleared -> set 1, bit set -> set 2
  // partition_2 mask is: bit cleared -> set 1, bit set -> set 2, hi bit set -> set 3
  if (m_numsets == 1)
    m_mask[0] = ( m_partmask & 0xFFFF) & ( 0xFFFFFFFF >> 16),	// color-set
    m_mask[1] = ( m_partmask & 0xFFFF) & ( 0xFFFFFFFF >> 16),	// alpha-set
    m_mask[2] = 0;
  if (m_numsets == 2)
    m_mask[0] = (~m_partmask & 0xFFFF) & ( 0xFFFFFFFF >> 16),
    m_mask[1] = ( m_partmask & 0xFFFF) & ( 0xFFFFFFFF >> 16),
    m_mask[2] = 0;
  if (m_numsets == 3)
    m_mask[0] = (~m_partmask & 0xFFFF) & (~m_partmask >> 16),
    m_mask[1] = ( m_partmask & 0xFFFF) & (~m_partmask >> 16),
    m_mask[2] = ( 0xFFFFFFFF & 0xFFFF) & ( m_partmask >> 16);

  /*
  static const float dw[] = {
    sqrtf(1.5 * 1.5 + 1.5 * 1.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 1.5 * 1.5),

    sqrtf(0.5 * 0.5 + 1.5 * 1.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 1.5 * 1.5),

    sqrtf(0.5 * 0.5 + 1.5 * 1.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 1.5 * 1.5),

    sqrtf(1.5 * 1.5 + 1.5 * 1.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 1.5 * 1.5)
  };
  */

  // check the compression mode for btc
  bool clearAlpha    = ((flags & kExcludeAlphaFromPalette) != 0);
//bool indexAlpha    = ((flags & kExcludeAlphaFromPalette) == 0) && !m_seperatealpha;
  bool seperateAlpha = ((flags & kExcludeAlphaFromPalette) == 0) &&  m_seperatealpha;
  bool weightByAlpha = ((flags & kWeightColourByAlpha    ) != 0);
  const float *rgbLUT = ComputeGammaLUT((flags & kSrgbIn) != 0);
  const float *aLUT   = ComputeGammaLUT(false);

  // swap: ra, ga, ba
  int rot[4] = {0,1,2,3};
  switch (m_rotid) {
    case 1: rot[0] = 3; rot[3] = 0; break;
    case 2: rot[1] = 3; rot[3] = 1; break;
    case 3: rot[2] = 3; rot[3] = 2; break;
  }

  // build mapped data
  u8 mska = !seperateAlpha ? 0xFF : 0x00;
  u8 clra = !clearAlpha    ? 0x00 : 0xFF;
  u8 rgbx[4 * 16];
  u8 ___a[1 * 16];

  for (int i = 0; i < 16; ++i) {
    u8 temp[4];

    // clear alpha
    temp[0] = rgba[4 * i + 0];
    temp[1] = rgba[4 * i + 1];
    temp[2] = rgba[4 * i + 2];
    temp[3] = rgba[4 * i + 3] | clra;

    // clear channel 3
    rgbx[4 * i + 0] = temp[rot[0]];
    rgbx[4 * i + 1] = temp[rot[1]];
    rgbx[4 * i + 2] = temp[rot[2]];
    rgbx[4 * i + 3] = temp[rot[3]] & mska;

    // separate channel 3
    ___a[1 * i + 0] = temp[rot[3]];

    // check for transparency (after blanking out)
    m_transparent = m_transparent || (temp[3] < 255);
  }

  for (int s = 0; s < m_numsets; s++) {
    // combined exclusion and selection mask
    int pmask = mask & m_mask[s];

    // create the minimal set
    for (int i = 0; i < 16; ++i) {
      // check this pixel is enabled
      int bit = 1 << i;
      if ((pmask & bit) == 0) {
	if ((mask & bit) == 0)
	  m_remap[s][i] = -1;
	continue;
      }

      // ensure there is always non-zero weight even for zero alpha
      float W = (float)(rgba[4 * i + 3] + 1) / 256.0f;

      // loop over previous points for a match
      u8 *rgbvalue = &rgbx[4 * i + 0];
      for (int j = 0;; ++j) {
	u8 *crgbvalue = &rgbx[4 * j + 0];

	// allocate a new point
	if (j == i) {
	  // get the index of the match and advance
	  int index = m_count[s]++;

	  // normalize coordinates to [0,1]
	  float x = rgbLUT[rgbvalue[0]];
	  float y = rgbLUT[rgbvalue[1]];
	  float z = rgbLUT[rgbvalue[2]];
	  float w =   aLUT[rgbvalue[3]];

	  // add the point
	  m_points[s][index] = Vec4(x, y, z, w);
	  m_weights[s][index] = (weightByAlpha ? W : 1.0f);
	  m_frequencies[s][index] = 1;
	  m_remap[s][i] = index;
	  break;
	}

	// check for a match
	int oldbit = 1 << j;
	bool match = ((pmask & oldbit) != 0)
	  && (rgbvalue[0] == crgbvalue[0])
	  && (rgbvalue[1] == crgbvalue[1])
	  && (rgbvalue[2] == crgbvalue[2])
	  && (rgbvalue[3] == crgbvalue[3]);

	if (match) {
	  // get the index of the match
	  int index = m_remap[s][j];

	  // map to this point and increase the weight
	  m_weights[s][index] += (weightByAlpha ? W : 1.0f);
	  m_frequencies[s][index] += 1;
	  m_remap[s][i] = index;
	  break;
	}
      }
    }

    // square root the weights
    for (int i = 0; i < m_count[s]; ++i)
      m_weights[s][i] = Vec4::sqrt(m_weights[s][i]);

    // TODO: if not m_transparent this all becomes a constant!
    if (seperateAlpha) {
      // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
      int a = s + m_numsets;

      // create the minimal set
      for (int i = 0; i < 16; ++i) {
	// check this pixel is enabled
	int bit = 1 << i;
	if ((pmask & bit) == 0) {
	  if ((mask & bit) == 0)
	    m_remap[a][i] = -1;
	  continue;
	}

	// ensure there is always non-zero weight even for zero alpha
	float W = (float)(rgba[4 * i + 3] + 1) / 256.0f;

	// loop over previous points for a match
	u8 *avalue = &___a[1 * i + 0];
	for (int j = 0;; ++j) {
	  u8 *cavalue = &___a[1 * j + 0];

	  // allocate a new point
	  if (j == i) {
	    // get the index of the match and advance
	    int index = m_count[a]++;

	    // normalize coordinates to [0,1]
	    float w = aLUT[avalue[0]];

	    // add the point
	    m_points[a][index] = Vec4(w);
	    m_weights[a][index] = (weightByAlpha ? W : 1.0f);
	    m_frequencies[a][index] = 1;
	    m_remap[a][i] = index;
	    break;
	  }

	  // check for a match
	  int oldbit = 1 << j;
	  bool match = ((pmask & oldbit) != 0)
	    && (avalue[0] == cavalue[0]);

	  if (match) {
	    // get the index of the match
	    int index = m_remap[a][j];

	    // map to this point and increase the weight
	    m_weights[a][index] += (weightByAlpha ? W : 1.0f);
	    m_frequencies[a][index] += 1;
	    m_remap[a][i] = index;
	    break;
	  }
	}
      }

      // square root the weights
      for (int i = 0; i < m_count[a]; ++i)
	m_weights[a][i] = Vec4::sqrt(m_weights[a][i]);
    }
  }

  // clear if we're suppose to throw alway alpha
  m_transparent = m_transparent && !clearAlpha;
}

void PaletteSet::RemapIndices(u8 const* source, u8* target, int set) const
{
  const int s = set; {
    // selection mask
    int pmask = m_mask[s];

    for (int i = 0; i < 16; ++i) {
      // check this pixel is enabled
      int bit = 1 << i;
      if ((pmask & bit) == 0)
	continue;

      u8 t = 0; t = ((m_remap[s][i] == -1) ? t : source[m_remap[s][i]]); target[i] = t;
    }
  }
}

void PaletteSet::UnmapIndices(u8 const* source, u8* destination, int set, int *codes, int cmask) const
{
  const int s = set; {
    // selection mask
    int pmask = m_mask[s];

    for (int i = 0; i < 16; ++i) {
      // check this pixel is enabled
      int bit = 1 << i;
      if ((pmask & bit) == 0)
	continue;

      if ((cmask >>  0) & 0xFF) destination[4 * i + 0] = (u8)(codes[source[i]] >>  0);
      if ((cmask >>  8) & 0xFF) destination[4 * i + 1] = (u8)(codes[source[i]] >>  8);
      if ((cmask >> 16) & 0xFF) destination[4 * i + 2] = (u8)(codes[source[i]] >> 16);
      if ((cmask >> 24) & 0xFF) destination[4 * i + 3] = (u8)(codes[source[i]] >> 24);
    }
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(USE_AMP) || defined(USE_COMPUTE)
void PaletteSet_CCR::CountSet(tile_barrier barrier, const int thread, pixel16 rgba, int mask, const bool alpha, const bool trans) amp_restricted
{
  // check the compression mode for dxt1
  const bool indexAlpha = (alpha);
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
      {
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
	    && (rgba[i][0] == rgba[j][0] || !indexAlpha);

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
	      (float)rgba[i][1],
	      (float)rgba[i][0]
	    );

	    m_remap[i] = index;
	    m_transparent = m_transparent || (indexAlpha && rgba[i][0]);
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
#if	!defined(USE_COMPUTE)
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

int PaletteSet_CCR::GetCount() amp_restricted {
  return m_count;
}

point16 PaletteSet_CCR::GetPoints() amp_restricted {
  return m_points;
}

weight16 PaletteSet_CCR::GetWeights() amp_restricted {
  return m_weights;
}

bool PaletteSet_CCR::IsTransparent() amp_restricted {
  return !!m_transparent;
}

void PaletteSet_CCR::WriteSet(tile_barrier barrier, const int thread, lineC2 cline, inout index16x2 source, out code64 block, const int is4,
			     IndexBlockLUT yArr) amp_restricted
{
  // build the block if we win
  {
    // remap the indices
    RemapIndices(barrier, thread, source);

    // save the block
    {
      if (!is4)
	WritePaletteBlock3(barrier, thread, cline, m_indices, block, yArr);
      else
	WritePaletteBlock4(barrier, thread, cline, m_indices, block, yArr);
    }
  }
}

void PaletteSet_CCR::RemapIndices(tile_barrier barrier, const int thread, inout index16x2 source) amp_restricted
{
  // make writes to "source"/"m_remap" visible to all
  tile_static_memory_fence(barrier);

#if	!defined(USE_COMPUTE)
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
