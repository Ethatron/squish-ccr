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
  // a linear rescale of the RGB-cube: scaling preserves
  // relative lengths, in either cube the error-vector with
  // smallest length is the same

  // values are directly out of the codebook and
  // natural numbers / 255, no need to round
  PackBytes(FloatToInt<true>((*values) * Vec3(255.0f)), (unsigned int &)(m_colour));

  /*
  assert(m_colour[0] == (u8)FloatToInt<true,false>(255.0f * values->X(), 255));
  assert(m_colour[1] == (u8)FloatToInt<true,false>(255.0f * values->Y(), 255));
  assert(m_colour[2] == (u8)FloatToInt<true,false>(255.0f * values->Z(), 255));
   */
}

void ColourSingleSnap::Compress3b(void* block)
{
  // grab the single colour
  Vec3 const* values = m_colours->GetPoints();
  
  // if it's black, make it index 3
  if (values[0] == Vec3(0.0f)) {
    *((unsigned__int64 *)block) = 0xFFFFFFFF00000000ULL;
  }
}

void ColourSingleSnap::Compress3(void* block)
{
  // grab the single colour
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();

  // just assign the end-points of index 2 (interpolant 1)
  Col3 s = Col3(
    sc_lookup_5_3[m_colour[0]].start,
    sc_lookup_6_3[m_colour[1]].start,
    sc_lookup_5_3[m_colour[2]].start);
  Col3 e = Col3(
    sc_lookup_5_3[m_colour[0]].end,
    sc_lookup_6_3[m_colour[1]].end,
    sc_lookup_5_3[m_colour[2]].end);

  m_start = Vec3(s) * (1.0f / 255.0f);
  m_end   = Vec3(e) * (1.0f / 255.0f);
  
  // created interpolated value and error
  Vec3 code = (m_start + m_end) * 0.5f;
  Scr3 error = LengthSquared(m_metric * (values[0] - code)) * freq[0];

  // build the block if we win
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // build the block
    u8 idx = 2, indices[16];
    m_colours->RemapIndices(&idx, indices);

    // save the block
    WriteColourBlock3(m_start, m_end, indices, block);
  }
}

void ColourSingleSnap::Compress4(void* block)
{
  // grab the single colour
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();

  // just assign the end-points of index 2 (interpolant 1)
  Col3 s = Col3(
    sc_lookup_5_4[m_colour[0]].start,
    sc_lookup_6_4[m_colour[1]].start,
    sc_lookup_5_4[m_colour[2]].start);
  Col3 e = Col3(
    sc_lookup_5_4[m_colour[0]].end,
    sc_lookup_6_4[m_colour[1]].end,
    sc_lookup_5_4[m_colour[2]].end);

  m_start = Vec3(s) * (1.0f / 255.0f);
  m_end   = Vec3(e) * (1.0f / 255.0f);
  
  // created interpolated value and error
  Vec3 code = (2.0f * m_start + m_end) * (1.0f / 3.0f);
  Scr3 error = LengthSquared(m_metric * (values[0] - code)) * freq[0];

  // build the block if we win
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // build the block
    u8 idx = 2, indices[16];
    m_colours->RemapIndices(&idx, indices);

    // save the block
    WriteColourBlock4(m_start, m_end, indices, block);
  }
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
