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

#ifndef SQUISH_HDRBLOCK_H
#define SQUISH_HDRBLOCK_H

#include <squish.h>
#include "maths.h"

// pull in structure definitions
#if	defined(USE_AMP)
#include "degeneracy_ccr.inl"
#endif

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(USE_PRE)
  void WriteDynamicBlock3_m1(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m2(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m3(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m4(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m5(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m6(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m7(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m8(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_m9(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_mA(int partition, Vec4 (&start)[2], Vec4 (&end)[2], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_mB(int partition, Vec4 (&start)[1], Vec4 (&end)[1], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_mC(int partition, Vec4 (&start)[1], Vec4 (&end)[1], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_mD(int partition, Vec4 (&start)[1], Vec4 (&end)[1], u16 (&indices)[1][16], void* block);
  void WriteDynamicBlock3_mE(int partition, Vec4 (&start)[1], Vec4 (&end)[1], u16 (&indices)[1][16], void* block);

  void DecompressDynamicBtc(u16* rgba, void const* block);
#endif

// -----------------------------------------------------------------------------
#if	defined(USE_AMP) || defined(USE_COMPUTE)
  void WriteDynamicBlock3(tile_barrier barrier, const int thread,
			  lineC2 cline, inout index16 indices, out code64 block) amp_restricted;
  void WriteDynamicBlock4(tile_barrier barrier, const int thread,
			  lineC2 cline, inout index16 indices, out code64 block) amp_restricted;

/*void DecompressDynamicBtc(tile_barrier barrier, const int thread,
			out pixel16 rgba, bool isBtc1) amp_restricted;*/
#endif

} // namespace squish

#endif // ndef SQUISH_HDRBLOCK_H
