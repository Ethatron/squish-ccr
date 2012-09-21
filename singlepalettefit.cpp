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

#include "singlepalettefit.h"
#include "paletteset.h"
#include "paletteblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(USE_PRE)
struct SourceBlock
{
  u8 start;
  u8 end;
  u8 error;
};

struct SinglePaletteLookup2
{
  SourceBlock sources[2];
};

struct SinglePaletteLookup4
{
  SourceBlock sources[4];
};

struct SinglePaletteLookup8
{
  SourceBlock sources[8];
};

#include "singlepalettelookup.inl"

SinglePaletteFit::SinglePaletteFit(PaletteSet const* palette, int flags, int swap)
  : PaletteFit(palette, flags, swap)
{
}

Vec4 SinglePaletteFit::ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, int cb, int ab, int ib, u8 cmask)
{
  assume(ib >= 2 && ib <= 4);
  switch (ib) {
    case 2: {
      SinglePaletteLookup2 const* cl;
      SinglePaletteLookup2 const* al;

      assume(cb >= 5 && cb <= 8);
      switch (cb) {
	case 5:  cl = lookup_5_4; break;
	case 6:  cl = lookup_6_4; break;
	case 7:  cl = lookup_7_4; break;
	case 8:  cl = lookup_8_4; break;
	default: cl = lookup_8_4; break;
      }

      assume(ab >= 5 && ab <= 8);
      switch (ab) {
	case 5:  al = lookup_5_4; break;
	case 6:  al = lookup_6_4; break;
	case 7:  al = lookup_7_4; break;
	case 8:  al = lookup_8_4; break;
	default: al = lookup_8_4; break;
      }

      SinglePaletteLookup2 const* const lookups[] =
      { cl, cl, cl, al };

      return ComputeEndPoints(set, metric, q, lookups, cmask);
    } break;
    case 3: {
      SinglePaletteLookup4 const* cl;
      SinglePaletteLookup4 const* al;

      assume(cb >= 5 && cb <= 7);
      switch (cb) {
	case 5:  cl = lookup_5_8; break;
	case 6:  cl = lookup_6_8; break;
	case 7:  cl = lookup_7_8; break;
	default: cl = lookup_7_8; break;
      }

      assume(ab >= 5 && ab <= 7);
      switch (ab) {
	case 5:  al = lookup_5_8; break;
	case 6:  al = lookup_6_8; break;
	case 7:  al = lookup_7_8; break;
	default: al = lookup_7_8; break;
      }

      SinglePaletteLookup4 const* const lookups[] =
      { cl, cl, cl, al };

      return ComputeEndPoints(set, metric, q, lookups, cmask);
    } break;
    case 4: {
      SinglePaletteLookup8 const* cl;
      SinglePaletteLookup8 const* al;

      assume(cb >= 8 && cb <= 8);
      switch (cb) {
	case 8:  cl = lookup_8_16; break;
	default: cl = lookup_8_16; break;
      }

      assume(ab >= 8 && ab <= 8);
      switch (ab) {
	case 8:  al = lookup_8_16; break;
	default: al = lookup_8_16; break;
      }

      SinglePaletteLookup8 const* const lookups[] =
      { cl, cl, cl, al };

      return ComputeEndPoints(set, metric, q, lookups, cmask);
    } break;

    default:
      return Vec4(INT_MAX);
      break;
  }
}

Vec4 SinglePaletteFit::ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, SinglePaletteLookup2 const* const* lookups, u8 cmask)
{
  // check each index combination (endpoint or intermediate)
  m_error = Vec4(FLT_MAX);

  // grab the single entry
  Vec4 const* values = m_palette->GetPoints(set);
  const float *eLUT = ComputeGammaLUT(false);

  // TODO: vectorize
  m_entry[set][0] = (u8)FloatToInt(255.0f * values->X(), 255);
  m_entry[set][1] = (u8)FloatToInt(255.0f * values->Y(), 255);
  m_entry[set][2] = (u8)FloatToInt(255.0f * values->Z(), 255);
  m_entry[set][3] = (u8)FloatToInt(255.0f * values->W(), 255);

  for (int index = 0; index < 2; ++index) {
    // check the error for this codebook index
    SourceBlock const* sources[4] = {NULL};
    Vec4 error(0.0f);

    for (int channel = 0; channel < 4; ++channel) {
      // skip if it's completely irrelevant what's in a specific channel
      if (cmask & (1 << channel)) {
	// grab the lookup table and index for this channel
	SinglePaletteLookup2 const* lookup = lookups[channel];
	int target = m_entry[set][channel];

	// store a pointer to the source for this channel
	sources[channel] = lookup[target].sources + index;

	// accumulate the error
	error.GetO(channel) = eLUT[sources[channel]->error];
      }
    }

    // calculate the error
    error = LengthSquared(metric * error);

    // keep it if the error is lower
    if (CompareFirstLessThan(error, m_error)) {
      m_start[set] = Vec4(
	sources[0] ? (float)sources[0]->start : 0.0f,
	sources[1] ? (float)sources[1]->start : 0.0f,
	sources[2] ? (float)sources[2]->start : 0.0f,
	sources[3] ? (float)sources[3]->start : 255.0f
      ) * q.gridrcp;

      m_end[set] = Vec4(
	sources[0] ? (float)sources[0]->end : 0.0f,
	sources[1] ? (float)sources[1]->end : 0.0f,
	sources[2] ? (float)sources[2]->end : 0.0f,
	sources[3] ? (float)sources[3]->end : 255.0f
      ) * q.gridrcp;

      m_index = (u8)(1 * index);
      m_error = error;

      // early out
      if (!CompareFirstGreaterThan(m_error, Vec4(0.0f)))
	return m_error;
    }
  }

  return m_error;
}

Vec4 SinglePaletteFit::ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, SinglePaletteLookup4 const* const* lookups, u8 cmask)
{
  // check each index combination (endpoint or intermediate)
  m_error = Vec4(FLT_MAX);

  // grab the single entry
  Vec4 const* values = m_palette->GetPoints(set);
  const float *eLUT = ComputeGammaLUT(false);

  /// TODO: vectorize
  m_entry[set][0] = (u8)FloatToInt(255.0f * values->X(), 255);
  m_entry[set][1] = (u8)FloatToInt(255.0f * values->Y(), 255);
  m_entry[set][2] = (u8)FloatToInt(255.0f * values->Z(), 255);
  m_entry[set][3] = (u8)FloatToInt(255.0f * values->W(), 255);

  for (int index = 0; index < 4; ++index) {
    // check the error for this codebook index
    SourceBlock const* sources[4] = {NULL};
    Vec4 error(0.0f);

    for (int channel = 0; channel < 4; ++channel) {
      // skip if it's completely irrelevant what's in a specific channel
      if (cmask & (1 << channel)) {
	// grab the lookup table and index for this channel
	SinglePaletteLookup4 const* lookup = lookups[channel];
	int target = m_entry[set][channel];

	// store a pointer to the source for this channel
	sources[channel] = lookup[target].sources + index;

	// accumulate the error
	error.GetO(channel) = eLUT[sources[channel]->error];
      }
    }

    // calculate the error
    error = LengthSquared(metric * error);

    // keep it if the error is lower
    if (CompareFirstLessThan(error, m_error)) {
      m_start[set] = Vec4(
	sources[0] ? (float)sources[0]->start : 0.0f,
	sources[1] ? (float)sources[1]->start : 0.0f,
	sources[2] ? (float)sources[2]->start : 0.0f,
	sources[3] ? (float)sources[3]->start : 255.0f
      ) * q.gridrcp;

      m_end[set] = Vec4(
	sources[0] ? (float)sources[0]->end : 0.0f,
	sources[1] ? (float)sources[1]->end : 0.0f,
	sources[2] ? (float)sources[2]->end : 0.0f,
	sources[3] ? (float)sources[3]->end : 255.0f
      ) * q.gridrcp;

      m_index = (u8)(1 * index);
      m_error = error;

      // early out
      if (!CompareFirstGreaterThan(m_error, Vec4(0.0f)))
	return m_error;
    }
  }

  return m_error;
}

Vec4 SinglePaletteFit::ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, SinglePaletteLookup8 const* const* lookups, u8 cmask)
{
  // check each index combination (endpoint or intermediate)
  m_error = Vec4(FLT_MAX);

  // grab the single entry
  Vec4 const* values = m_palette->GetPoints(set);
  const float *eLUT = ComputeGammaLUT(false);

  /// TODO: vectorize
  m_entry[set][0] = (u8)FloatToInt(255.0f * values->X(), 255);
  m_entry[set][1] = (u8)FloatToInt(255.0f * values->Y(), 255);
  m_entry[set][2] = (u8)FloatToInt(255.0f * values->Z(), 255);
  m_entry[set][3] = (u8)FloatToInt(255.0f * values->W(), 255);

  for (int index = 0; index < 8; ++index) {
    // check the error for this codebook index
    SourceBlock const* sources[4] = {NULL};
    Vec4 error(0.0f);

    for (int channel = 0; channel < 4; ++channel) {
      // skip if it's completely irrelevant what's in a specific channel
      if (cmask & (1 << channel)) {
	// grab the lookup table and index for this channel
	SinglePaletteLookup8 const* lookup = lookups[channel];
	int target = m_entry[set][channel];

	// store a pointer to the source for this channel
	sources[channel] = lookup[target].sources + index;

	// accumulate the error
	error.GetO(channel) = eLUT[sources[channel]->error];
      }
    }

    // calculate the error
    error = LengthSquared(metric * error);

    // keep it if the error is lower
    if (CompareFirstLessThan(error, m_error)) {
      m_start[set] = Vec4(
	sources[0] ? (float)sources[0]->start : 0.0f,
	sources[1] ? (float)sources[1]->start : 0.0f,
	sources[2] ? (float)sources[2]->start : 0.0f,
	sources[3] ? (float)sources[3]->start : 255.0f
      ) * q.gridrcp;

      m_end[set] = Vec4(
	sources[0] ? (float)sources[0]->end : 0.0f,
	sources[1] ? (float)sources[1]->end : 0.0f,
	sources[2] ? (float)sources[2]->end : 0.0f,
	sources[3] ? (float)sources[3]->end : 255.0f
      ) * q.gridrcp;

      m_index = (u8)(1 * index);
      m_error = error;

      // early out
      if (!CompareFirstGreaterThan(m_error, Vec4(0.0f)))
	return m_error;
    }
  }

  return m_error;
}
#endif

} // namespace squish

#if	defined(SBL_FLAT)
#include "singlepalettefit_ccr_flat.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 1)
#include "singlepalettefit_ccr_packed.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 2)
#include "singlepalettefit_ccr_packed_copy.cpp"
#elif	defined(SBL_VECTOR)
#include "singlepalettefit_ccr_vector.cpp"
#endif
