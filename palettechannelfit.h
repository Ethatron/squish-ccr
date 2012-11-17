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

#ifndef SQUISH_PALETTECHANNELFIT_H
#define SQUISH_PALETTECHANNELFIT_H

#include <squish.h>
#include <limits.h>

#include "palettefit.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class PaletteSet;
class PaletteChannelFit : public virtual PaletteFit
{
public:
  PaletteChannelFit(PaletteSet const* colours, int flags, int swap = -1, int shared = 0);

private:
  int  m_channel[4];
  Vec4 m_start_candidate[4];
  Vec4 m_end_candidate  [4];

protected:
  Scr4 ComputeCodebook(int set, Vec4 const &metric, vQuantizer &q, int sb, int ib, u8 (&closest)[16]);

  bool IsChannel(int set) { return m_channel[set] >= 0; }
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
struct PaletteChannelFit_CCR : inherit_hlsl PaletteFit_CCR
{
public_hlsl
  void AssignSet (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, const int metric, const int fit) amp_restricted;
  void Compress  (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block, const bool trans,
		  IndexBlockLUT yArr, PaletteChannelLUT lArr) amp_restricted;

protected_hlsl
  void Compress3 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
		  IndexBlockLUT yArr, PaletteChannelLUT lArr) amp_restricted;
  void Compress4 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
		  IndexBlockLUT yArr, PaletteChannelLUT lArr) amp_restricted;
  void Compress34(tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
		  IndexBlockLUT yArr, PaletteChannelLUT lArr) amp_restricted;

  void ComputeEndPoints(tile_barrier barrier, const int thread, const int is4,
		        PaletteChannelLUT lArr) amp_restricted;
  int  ComputeEndPoints(tile_barrier barrier, const int thread,
		        PaletteChannelLUT lArr) amp_restricted;

#if	!defined(SQUISH_USE_COMPUTE)
private_hlsl
  int3 m_entry;
  ccr8 m_index;
#endif
};

#if	defined(SQUISH_USE_COMPUTE)
  tile_static int3 m_entry;
  tile_static ccr8 m_index;
#endif
#endif

} // namespace squish

#endif // ndef SQUISH_PALETTECHANNELFIT_H
