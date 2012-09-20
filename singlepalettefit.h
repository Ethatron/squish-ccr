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

#ifndef SQUISH_SINGLEPALETTEFIT_H
#define SQUISH_SINGLEPALETTEFIT_H

#include <squish.h>
#include <limits.h>

#include "palettefit.h"

// pull in structure definitions
#if	defined(USE_AMP)
#include "singlepalettelookup_ccr.inl"
#endif

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(USE_PRE)
class PaletteSet;
struct SinglePaletteLookup2;
struct SinglePaletteLookup4;
struct SinglePaletteLookup8;
class SinglePaletteFit : public PaletteFit
{
public:
  SinglePaletteFit(PaletteSet const* colours, int flags, int swap);

private:
  Vec4 ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, SinglePaletteLookup2 const* const* lookups, u8 cmask);
  Vec4 ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, SinglePaletteLookup4 const* const* lookups, u8 cmask);
  Vec4 ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, SinglePaletteLookup8 const* const* lookups, u8 cmask);

  u8 m_entry[4][4];
  u8 m_index;
  Vec4 m_error;

protected:
  Vec4 ComputeEndPoints(int set, Vec4 const &metric, vQuantizer &q, int cb, int ab, int ib, u8 cmask);
  u8 GetIndex() { return m_index; }
};
#endif

// -----------------------------------------------------------------------------
#if	defined(USE_AMP) || defined(USE_COMPUTE)
struct SinglePaletteFit_CCR : inherit_hlsl PaletteFit_CCR
{
public_hlsl
  void AssignSet (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, const int metric, const int fit) amp_restricted;
  void Compress  (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block, const bool trans,
		  IndexBlockLUT yArr, SinglePaletteLUT lArr) amp_restricted;

protected_hlsl
  void Compress3 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
		  IndexBlockLUT yArr, SinglePaletteLUT lArr) amp_restricted;
  void Compress4 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
		  IndexBlockLUT yArr, SinglePaletteLUT lArr) amp_restricted;
  void Compress34(tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
		  IndexBlockLUT yArr, SinglePaletteLUT lArr) amp_restricted;

  void ComputeEndPoints(tile_barrier barrier, const int thread, const int is4,
		        SinglePaletteLUT lArr) amp_restricted;
  int  ComputeEndPoints(tile_barrier barrier, const int thread,
		        SinglePaletteLUT lArr) amp_restricted;

#if	!defined(USE_COMPUTE)
private_hlsl
  int3 m_entry;
  ccr8 m_index;
#endif
};

#if	defined(USE_COMPUTE)
  tile_static int3 m_entry;
  tile_static ccr8 m_index;
#endif
#endif

} // namespace squish

#endif // ndef SQUISH_SINGLEPALETTEFIT_H
