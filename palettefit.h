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

#ifndef SQUISH_PALETTEFIT_H
#define SQUISH_PALETTEFIT_H

#include <squish.h>
#include "maths.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class PaletteSet;
class PaletteFit
{
public:
  static int GetNumSets(int mode);
  static int GetPartitionBits(int mode);
  static int GetIndexBits(int mode);
  static int GetRotationBits(int mode);
  static int GetSelectionBits(int mode);
  static int GetSharedBits(int mode);
  static int GetPrecisionBits(int mode);

  static const int *GetSharedMap(int mode);
  static int GetSharedSkip(int mode);

  // rotate shared bit definition: 0=0 (-), 1=3 (s), 2=2 (u), 3=1 (u)
  // unique p-bit permutations: upper bit start bit set, lower bit stop bit set
  // shared p-bit permutations: upper & lower bit start & stop bit set
  // makes it easier to loop
#define SBSTART	0
#define SBEND	3
#define SR(s)	(s < 0 ? s : m_sharedmap[s])
#define SBSKIP	-1
#define SK(s)	(!(~s))

public:
  PaletteFit(PaletteSet const* palette, int flags, int swap = -1, int shared = -1);

  // change parameters while iterating
  void ChangeFit(PaletteSet const* palette, int flags, int swap, int shared) { m_palette = palette; m_flags = flags; m_swapindex = swap; m_sharedbits = SR(shared); m_best = false; }
  void ChangePalette(PaletteSet const* palette) { m_palette = palette; }
  void ChangeFlags(int flags) { m_flags = flags; }
  void ChangeSwap(int swap) { m_swapindex = swap; }
  void ChangeMode(int mode) { m_sharedmap = GetSharedMap(m_mode = mode); }
  void ChangeShared(int shared) { m_sharedbits = SR(shared); }

  // query some values
  PaletteSet const* GetPalette() const { return m_palette; }
  int GetFlags() const { return m_flags; }
  int GetSwap() const { return m_swapindex; }
  int GetSharedField() const { return m_sharedbits; }

  // error management
  void SetError(Scr4 &error) { m_besterror = error; m_best = false; }
  Scr4 GetError() { return m_besterror; }

  void Compress(void* block, vQuantizer &q);
  virtual void Compress(void* block, vQuantizer &q, int mode) = 0;

#if 1 //ndef NDEBUG
  void Decompress(u8 *rgba, vQuantizer &q, int mode);
  void SumError(u8 (&closest)[4][16], vQuantizer &q, int mode, Scr4 &error);
#endif

  bool Lossless() { return !(m_besterror > Scr4(0.0f)); }
  bool IsBest() { return m_best; }

protected:
  PaletteSet const* m_palette;
  const int *m_sharedmap;
  int m_flags;
  int m_mode;
  int m_swapindex;
  int m_sharedbits;

  Vec4 m_start[4];
  Vec4 m_end[4];
  a16 u8 m_indices[2][16];
  
  Vec4 m_metric[3];
  Scr4 m_besterror;
  bool m_best;
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
struct PaletteFit_CCR
{
public_hlsl
  void AssignSet (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, int metric, int fit) amp_restricted;
#if 0
  void Compress  (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block, bool trans,
                  IndexBlockLUT yArr) amp_restricted;
#endif

protected_hlsl
  void Compress3 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;
  void Compress4 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;
  void Compress34(tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palette, out code64 block, const int is4,
                  IndexBlockLUT yArr) amp_restricted;

// Start- and end-point positions in arrays, for colors
#define	CSTRT	0
#define	CSTOP	1
#define	CVALS	2

#if	!defined(SQUISH_USE_COMPUTE)
protected_hlsl
  float4 m_line[CVALS];
  ccr8 m_matches[2][16];
#endif
};

#if	defined(SQUISH_USE_COMPUTE)
  tile_static float4 m_line[CVALS];
  tile_static ccr8 m_matches[2][16];
#endif
#endif

} // namespace squish

#endif // ndef SQUISH_PALETTEFIT_H
