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

#include "coloursinglefit.h"
#include "colourset.h"
#include "colourblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
struct SC_SourceBlock
{
  u8 start;
  u8 end;
  u8 error;
};

struct ColourSingleLookup
{
  SC_SourceBlock sources[2];
};

#define	SCL_ITERATIVE
#include "coloursinglelookup.inl"

ColourSingleFit::ColourSingleFit(ColourSet const* colours, int flags)
  : ColourFit(colours, flags)
{
  // grab the single colour
  Vec3 const* values = m_colours->GetPoints();
  Col3 integers = Min(FloatToInt<true>(*values * Vec3(255.0f)), Col3(255));

  m_colour[0] = (u8)integers.R();
  m_colour[1] = (u8)integers.G();
  m_colour[2] = (u8)integers.B();

  // initialize the best error
  m_besterror = Scr3(FLT_MAX);
}

void ColourSingleFit::Compress3(void* block)
{
  // build the table of lookups
  ColourSingleLookup const* const lookups[] =
  {
    sc_lookup_5_3,
    sc_lookup_6_3,
    sc_lookup_5_3
  };

  // find the best end-points and index
  Scr3 error = Scr3(1.0f / (255.0f * 255.0f)) * ComputeEndPoints(lookups);

  // build the block if we win
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16]; m_colours->RemapIndices(&m_index, indices);

    // save the block
    WriteColourBlock3(m_start, m_end, indices, block);
  }
}

void ColourSingleFit::Compress4(void* block)
{
  // build the table of lookups
  ColourSingleLookup const* const lookups[] =
  {
    sc_lookup_5_4,
    sc_lookup_6_4,
    sc_lookup_5_4
  };

  // find the best end-points and index
  Scr3 error = Scr3(1.0f / (255.0f * 255.0f)) * ComputeEndPoints(lookups);

  // build the block if we win
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16]; m_colours->RemapIndices(&m_index, indices);

    // save the block
    WriteColourBlock4(m_start, m_end, indices, block);
  }
}

int ColourSingleFit::ComputeEndPoints(ColourSingleLookup const* const* lookups)
{
  // check each index combination (endpoint or intermediate)
  int besterror = INT_MAX;

  for (int index = 0; index < 2; ++index) {
    // check the error for this codebook index
    SC_SourceBlock const* sources[3];
    int error = 0;

    for (int channel = 0; channel < 3; ++channel) {
      // grab the lookup table and index for this channel
      ColourSingleLookup const* lookup = lookups[channel];
      int target = m_colour[channel];

      // store a pointer to the source for this channel
      sources[channel] = lookup[target].sources + index;

      // accumulate the error
      int diff = sources[channel]->error;
      error += diff * diff;
    }

    // keep it if the error is lower
    if (error < besterror) {
      // save the error
      besterror = error;

      m_start = Vec3(
	(float)sources[0]->start,
	(float)sources[1]->start,
	(float)sources[2]->start
      );

      m_end = Vec3(
	(float)sources[0]->end,
	(float)sources[1]->end,
	(float)sources[2]->end
      );

      m_start /= Vec3(31.0f, 63.0f, 31.0f);
      m_end   /= Vec3(31.0f, 63.0f, 31.0f);

      m_index = (u8)(2 * index);

      // early out
      if (!besterror)
	return besterror;
    }
  }

  return besterror;
}
#endif

} // namespace squish

#if	defined(SBL_FLAT)
#include "coloursinglefit_ccr_flat.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 1)
#include "coloursinglefit_ccr_packed.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 2)
#include "coloursinglefit_ccr_packed_copy.cpp"
#elif	defined(SBL_VECTOR)
#include "coloursinglefit_ccr_vector.cpp"
#endif
