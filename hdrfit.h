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

#ifndef SQUISH_HDRFIT_H
#define SQUISH_HDRFIT_H

#include <squish.h>
#include "maths.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class HDRSet;
class HDRFit
{
public:
  static int GetNumSets(int mode);
  static int GetPartitionBits(int mode);
  static int GetIndexBits(int mode);

  static int GetExplicitBits(int mode);
  static int GetSharedBits(int mode);
  static int GetDeltaBits(int mode);
  static int GetTruncationBits(int mode);
  static int GetPrecisionBits(int mode);

  static int GetSharedBits(int mode, int channel);
  static int GetDeltaBits(int mode, int channel);

public:
  HDRFit(HDRSet const* palette, int flags);

  // change parameters while iterating
  void ChangeFit(HDRSet const* palette, int flags) { m_palette = palette; m_flags = flags; m_best = false; }
  void ChangePalette(HDRSet const* palette) { m_palette = palette; }
  void ChangeFlags(int flags) { m_flags = flags; }

  // query some values
  HDRSet const* GetPalette() const { return m_palette; }
  int GetFlags() const { return m_flags; }

  // error management
  void SetError(Scr3 &error) { m_besterror = error; m_best = false; }
  Scr3 GetError() { return m_besterror; }

  void Compress(void* block, fQuantizer &q);
  virtual void Compress(void* block, fQuantizer &q, int mode) = 0;

#if 1 //ndef NDEBUG
  void Decompress(u16 *rgb, fQuantizer &q, int mode);
  void SumError(u8 (&closest)[4][16], fQuantizer &q, int mode, Scr3 &error);
#endif

  bool Lossless() { return !(m_besterror > Scr3(0.0f)); }
  bool IsBest() { return m_best; }

protected:
  HDRSet const* m_palette;
  int m_flags;
  int m_mode;

  Vec3 m_start[2];
  Vec3 m_end[2];
  a16 u8 m_indices[16];

  Vec3 m_metric;
  Scr3 m_besterror;
  bool m_best;
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#endif // ndef SQUISH_HDRFIT_H
