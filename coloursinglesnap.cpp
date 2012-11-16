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

#include "coloursinglesnap.h"
#include "colourset.h"
#include "colourblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
struct ColourSingleLookup
{
  u8 start;
  u8 end;
};

#undef	SCL_ITERATIVE
#include "coloursinglelookup.inl"

ColourSingleSnap::ColourSingleSnap(ColourSet const* colours, int flags)
  : ColourFit(colours, flags)
{
  // grab the single colour
  Vec3 const* values = m_colours->GetPoints();

  // in the 3/4 codebook case it can be observed
  // that the error of the end-point is always
  // higher than the error of any interpolated value
  // in addition the codebook is symmetric, means
  // index 2 can always be made index 1
  // as a result the lookup can be made in such a way
  // that index 1 always contains the best match
  //
  // this holds even if a metric is used, the metric does
  // linear a rescale of the RGB-cube: scaling preserves
  // relative lengths, in either cube the error-vector with
  // smallest length is the same
  m_colour[0] = (u8)FloatToInt(255.0f * values->X(), 255);
  m_colour[1] = (u8)FloatToInt(255.0f * values->Y(), 255);
  m_colour[2] = (u8)FloatToInt(255.0f * values->Z(), 255);
}

void ColourSingleSnap::Compress3(void* block)
{
  const float *eLUT = ComputeGammaLUT(false);

  // just assign the end-points of index 1
  m_start = Vec3(
    eLUT[sc_lookup_5_3[m_colour[0]].start],
    eLUT[sc_lookup_6_3[m_colour[1]].start],
    eLUT[sc_lookup_5_3[m_colour[2]].start]);
  m_end = Vec3(
    eLUT[sc_lookup_5_3[m_colour[0]].end],
    eLUT[sc_lookup_6_3[m_colour[1]].end],
    eLUT[sc_lookup_5_3[m_colour[2]].end]);

  // build the block
  u8 idx = 1, indices[16];
  m_colours->RemapIndices(&idx, indices);

  // save the block
  WriteColourBlock3(m_start, m_end, indices, block);
}

void ColourSingleSnap::Compress4(void* block)
{
  const float *eLUT = ComputeGammaLUT(false);

  // find the best end-points and index
  m_start = Vec3(
    eLUT[sc_lookup_5_4[m_colour[0]].start],
    eLUT[sc_lookup_6_4[m_colour[1]].start],
    eLUT[sc_lookup_5_4[m_colour[2]].start]);
  m_end = Vec3(
    eLUT[sc_lookup_5_4[m_colour[0]].end],
    eLUT[sc_lookup_6_4[m_colour[1]].end],
    eLUT[sc_lookup_5_4[m_colour[2]].end]);

  // build the block
  u8 idx = 1, indices[16];
  m_colours->RemapIndices(&idx, indices);

  // save the block
  WriteColourBlock4(m_start, m_end, indices, block);
}
#endif

} // namespace squish

#if	defined(SBL_FLAT)
#include "coloursinglesnap_ccr_flat.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 1)
#include "coloursinglesnap_ccr_packed.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 2)
#include "coloursinglesnap_ccr_packed_copy.cpp"
#elif	defined(SBL_VECTOR)
#include "coloursinglesnap_ccr_vector.cpp"
#endif
