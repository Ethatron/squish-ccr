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
#if	defined(SQUISH_USE_AMP)
#include "degeneracy_ccr.inl"
#endif

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
  void WriteHDRBlock_m1(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m2(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m3(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m4(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m5(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m6(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m7(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m8(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_m9(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_mA(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block);

  void WriteHDRBlock_mB(	       Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_mC(	       Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_mD(	       Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block);
  void WriteHDRBlock_mE(	       Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block);

  void ReadHDRBlock_m1(u16* rgb, void const* block);
  void ReadHDRBlock_m2(u16* rgb, void const* block);
  void ReadHDRBlock_m3(u16* rgb, void const* block);
  void ReadHDRBlock_m4(u16* rgb, void const* block);
  void ReadHDRBlock_m5(u16* rgb, void const* block);
  void ReadHDRBlock_m6(u16* rgb, void const* block);
  void ReadHDRBlock_m7(u16* rgb, void const* block);
  void ReadHDRBlock_m8(u16* rgb, void const* block);
  void ReadHDRBlock_m9(u16* rgb, void const* block);
  void ReadHDRBlock_mA(u16* rgb, void const* block);

  void ReadHDRBlock_mB(u16* rgb, void const* block);
  void ReadHDRBlock_mC(u16* rgb, void const* block);
  void ReadHDRBlock_mD(u16* rgb, void const* block);
  void ReadHDRBlock_mE(u16* rgb, void const* block);

  void DecompressHDRsBtc6u(u16* rgb, void const* block);
  void DecompressHDRsBtc6u(f23* rgb, void const* block);
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#endif // ndef SQUISH_HDRBLOCK_H
