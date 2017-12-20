/* -----------------------------------------------------------------------------

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

#if defined(_MSC_VER) && (_MSC_VER > 1300)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>
#endif

#include "hdrblock.h"

#define SBSTART	0
#define SBEND	3

#include "inlineables.cpp"

#pragma warning(disable: 4100)

namespace squish {

/* *****************************************************************************
 */
extern const u16 weights_u16[5][16];
extern const Vec4 weights_V4[5][16];
extern const Col4 weights_C4[5][16];

/* *****************************************************************************
*/
#if	!defined(SQUISH_USE_PRE)
#ifdef __GNUC__
	__attribute__((__aligned__(16)))
#else
	__declspec(align(16))
#endif

extern const unsigned int blockxor[64][/*6*/5][4];
extern const int shorterindex[64][/*6*/5];

/* -----------------------------------------------------------------------------
 */
template<const int sets, const int ibits, const int begin>
static void passreg WriteHDRBlock(int partition, Col4 (&idx)[1], Col4 &blkl, Col4 &blkh)
{
  Col4 iidx = idx[0];
  Col4 iblk;

  /* none of the cases straddles the lo to hi border, all go into hi
   *
   * WriteHDRBlock<2, 3, 82>(partition, remapped, blkl, blkh);
   * WriteHDRBlock<1, 4, 66>(partition, remapped, blkl, blkh);
   */
  assume(partition >= 0 && partition <= 32);

  // max 16*4 -> 64bits
  iblk = CopyBits<ibits, ibits *  0>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  1>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  2>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  3>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  4>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  5>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  6>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  7>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  8>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  9>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 10>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 11>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 12>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 13>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 14>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 15>(iblk, iidx); //dx = ShiftRight<8>(iidx);

  // one bit is omitted per set, so it's one instruction per set + 1
  switch (sets) {
    case 1: {
      // always index 0
      blkh = CopyBits<ibits *  1 - 1, begin - 64 +         0>(blkh, iblk);
      iblk = ShiftRightHalf<ibits *  1>(iblk);
      blkh = CopyBits<ibits * 15 + 0, begin - 64 + ibits - 1>(blkh, iblk);
//    iblk = ShiftRightHalf<ibits * 15>(iblk);
    } break;
    case 2: {
      // always index 0
      blkh = CopyBits<ibits *  1 - 1, begin - 64 +         0>(blkh, iblk);
      iblk = ShiftRightHalf<ibits *  1>(iblk);
      blkh = CopyBits<ibits * 15 + 0, begin - 64 + ibits - 1>(blkh, iblk);
//    iblk = ShiftRightHalf<ibits * 15>(iblk);

      // if it's not the last bit which is cut (low probability)
      int len, bgn = begin - 64;
      if ((len = shorterindex[partition][0]) < 15) {
        len = (len * ibits) + ibits;
        bgn = (bgn + len - 2);

      	// subtract the already conducted shifts
        iblk = ShiftRightHalf(iblk, len - ibits);

        // remaining length can be anything, length overflow is silent
        blkh = CopyBits(blkh, iblk, ibits * 14, bgn);
//      iblk = ShiftRightHalf<ibits * 15>(iblk);
      }
    } break;
    case 3: {
      // always index 0
      blkh = CopyBits<ibits *  1 - 1, begin - 64 +         0>(blkh, iblk);
      iblk = ShiftRightHalf<ibits *  1>(iblk);
      blkh = CopyBits<ibits * 15 + 0, begin - 64 + ibits - 1>(blkh, iblk);
//    iblk = ShiftRightHalf<ibits * 15>(iblk);

      int bgn = begin - 64;
      int fln, fbn;
      fln = shorterindex[partition][3]; {
        fln = (fln * ibits) + ibits;
        fbn = (bgn + fln - 2);

      	// subtract the already conducted shifts
        iblk = ShiftRightHalf(iblk, fln - ibits);

        // remaining length can be anything, length overflow is silent
        blkh = CopyBits(blkh, iblk, ibits * 14, fbn);
//      iblk = ShiftRightHalf(iblk, ibits * 15);
      }

      // if it's not the last bit which is cut (low probability)
      int sln, sbn;
      if ((sln = shorterindex[partition][4]) < 15) {
        sln = (sln * ibits) + ibits;
        sbn = (bgn + sln - 3);

      	// subtract the already conducted shifts
        iblk = ShiftRightHalf(iblk, sln - fln);

        // remaining length can be anything, length overflow is silent
        blkh = CopyBits(blkh, iblk, ibits * 14, sbn);
//      iblk = ShiftRightHalf(iblk, ibits * 15);
      }
    } break;
  }
}

/* -----------------------------------------------------------------------------
 */
#define	ihibit	(1 << (ibits - 1))
#define	ihimsk	((1 << ibits) - 1)
#define	ihixor	Col4((ihimsk << 24) + (ihimsk << 16) + (ihimsk << 8) + (ihimsk << 0))

template<const int ibits>
static void passreg RemapHDRBlock(int partition, Col4 (&start)[2], Col4 (&end)[2], Col4 (&idxs)[1], u8 const (&indices)[1][16])
{
#if 0
  int masks[4], xors[4] = {0,0,0,0};

  PaletteSet::GetMasks(flags, masks);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapRGBA(end[0]); xors[0] = ihimsk; }
  /* set 2 */
  if (indices[0][shorterindex[partition][0]] & ihibit) {
    start[1].SwapRGBA(end[1]); xors[1] = ihimsk; }

  /* TODO: 16 bytes fit into one SSE2 register */
  for (int i = 0; i < 16; ++i) {
    int set =
      (((masks[1] >> i) & 1) * 1);

    indices[0][i] = indices[0][i] ^ xors[set];
  }
#else
  Col4 xors = Col4(0);

  LoadAligned(idxs[0], indices[0]);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapRGBA(end[0]); xors  = Col4(blockxor[partition][0]); }
  /* set 2 */
  if (indices[0][shorterindex[partition][0]] & ihibit) {
    start[1].SwapRGBA(end[1]); xors |= Col4(blockxor[partition][1]); }

  idxs[0] ^= (xors & ihixor);
#endif
}

template<const int ibits>
static void passreg RemapHDRBlock(int partition, Col4 (&start)[1], Col4 (&end)[1], Col4 (&idxs)[1], u8 const (&indices)[1][16])
{
#if 0
  int masks[4], xors[4] = {0,0,0,0};
  u8 const ihibit = 1 << ibits;

//PaletteSet::GetMasks(flags, masks);

  /* same for all set 1s */
  if (indices[0][                         0 ] & ihibit) {
    start[0].SwapRGBA(end[0]); xors[0] = ihimsk; }

  /* TODO: 16 bytes fit into one SSE2 register */
  for (int i = 0; i < 16; ++i) {
    int set = 0;

    indices[0][i] = indices[0][i] ^ xors[set];
  }
#else
  Col4 xors = Col4(0);

  LoadAligned(idxs[0], indices[0]);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapRGBA(end[0]); xors  = Col4(0xFFFFFFFF); }

  idxs[0] ^= (xors & ihixor);
#endif
}

/* -----------------------------------------------------------------------------
 */

enum EField
{
  NA,  // N/A
   M,  // Mode
   D,  // Shape
  RE2,
  RS1,
  RE1,
  RS2,
  GE2,
  GS1,
  GE1,
  GS2,
  BE2,
  BS1,
  BE1,
  BS2,
};

struct ModeDescriptor
{
  EField m_eField;
  char   m_uBit;
};

void VerifyHDRBlock(int num, ModeDescriptor *descs, int mode, int partition, Col4 *start, Col4 *end, int tb, Col4 &blkl, Col4 &blkh)
{
  for (int i = 0; i < num; i++) {
		int chkf;
		int ibit;
    int obit;

    switch (descs[i].m_eField) {
      case   M: ibit =         mode >>  descs[i].m_uBit           ; break;
      case   D: ibit =    partition >>  descs[i].m_uBit           ; break;

      case RS1: { int &ptf = start[0].GetR(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
      case GS1: { int &ptf = start[0].GetR(); ibit = (ptf) >> (descs[i].m_uBit + 16 + tb); } break;
      case BS1: { int &ptf = start[0].GetG(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
                                                                                                             
      case RS2: { int &ptf = start[1].GetR(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
      case GS2: { int &ptf = start[1].GetR(); ibit = (ptf) >> (descs[i].m_uBit + 16 + tb); } break;
      case BS2: { int &ptf = start[1].GetG(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
                                                                                                            
      case RE1: { int &ptf =   end[0].GetR(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
      case GE1: { int &ptf =   end[0].GetR(); ibit = (ptf) >> (descs[i].m_uBit + 16 + tb); } break;
      case BE1: { int &ptf =   end[0].GetG(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
                                                                                                            
      case RE2: { int &ptf =   end[1].GetR(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
      case GE2: { int &ptf =   end[1].GetR(); ibit = (ptf) >> (descs[i].m_uBit + 16 + tb); } break;
      case BE2: { int &ptf =   end[1].GetG(); ibit = (ptf) >> (descs[i].m_uBit      + tb); } break;
    }

    switch (i >> 5) {
      case 0: { int &blkf = blkl.GetR(); chkf = blkf; obit = (blkf) >> (i & 31); } break;
      case 1: { int &blkf = blkl.GetG(); chkf = blkf; obit = (blkf) >> (i & 31); } break;
      case 2: { int &blkf = blkh.GetR(); chkf = blkf; obit = (blkf) >> (i & 31); } break;
      case 3: { int &blkf = blkh.GetG(); chkf = blkf; obit = (blkf) >> (i & 31); } break;
    }

    ibit &= 1;
    obit &= 1;

    assert(ibit == obit);
  }
}

/* -----------------------------------------------------------------------------
 * Remarks
 *
 * ...
 *
 */

void WriteHDRBlock_m1(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);
   
  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
  s[0] = a[0];
  s[1] = a[1];
  e[0] = b[0];
  e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x00 - 10 5 5 5
        {  M, 0}, {  M, 1}, {GS2, 4}, {BS2, 4}, {BE2, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(0);		// 2 mode bit

  // truncate bottom 6 bits
  a[0] = ShiftRightHalf<6>(a[0]);
  b[0] = ShiftRightHalf<6>(b[0]);
  a[1] = ShiftRightHalf<6>(a[1]);
  b[1] = ShiftRightHalf<6>(b[1]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits< 5, 35 -  0>(blkl,                b[0]);	// 1- 5 set 1 red end
  blkh = InjtBits< 5, 65 - 64>(blkh,                a[1]);	// 1- 5 set 2 red start
  blkh = InjtBits< 5, 71 - 64>(blkh,                b[1]);	// 1- 5 set 2 red end

  // truncate bottom 6 bits
  a[0] = ShiftRightHalf<10 + 6>(a[0]);
  b[0] = ShiftRightHalf<10 + 6>(b[0]);
  a[1] = ShiftRightHalf<10 + 6>(a[1]);
  b[1] = ShiftRightHalf<10 + 6>(b[1]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits< 5, 45 -  0>(blkl,                b[0]);	// 1- 5 set 1 green end
  blkl = InjtBits< 4, 41 -  0>(blkl,                a[1]);	// 1- 4 set 2 green start
  blkl = InjtBits< 4, 51 -  0>(blkl,                b[1]);	// 1- 4 set 2 green end

  a[1] = ShiftRightHalf<4>(a[1]);
  b[1] = ShiftRightHalf<4>(b[1]);

  blkl = InjtBits< 1,  2 -  0>(blkl,                a[1]);	// 5- 5 set 2 green start
  blkl = InjtBits< 1, 40 -  0>(blkl,                b[1]);	// 5- 5 set 2 green end

  // truncate bottom 6 bits
  a[0] = ShiftRightHalf<10 + 6    >(a[0]);
  b[0] = ShiftRightHalf<10 + 6    >(b[0]);
  a[1] = ShiftRightHalf<10 + 6 - 4>(a[1]);
  b[1] = ShiftRightHalf<10 + 6 - 4>(b[1]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 5, 55 -  0>(blkl,                b[0]);	// 1- 5 set 1 blue end
  blkl = InjtBits< 3, 61 -  0>(blkl,                a[1]);	// 1- 3 set 2 blue start
  blkl = InjtBits< 1, 50 -  0>(blkl,                b[1]);	// 1- 1 set 2 blue end

  a[1] = ShiftRightHalf<3>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits< 1, 64 - 64>(blkh,                a[1]);	// 4- 4 set 2 blue start
  blkl = InjtBits< 1, 60 -  0>(blkl,                b[1]);	// 2- 2 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits< 1,  3 -  0>(blkl,                a[1]);	// 5- 5 set 2 blue start
  blkh = InjtBits< 1, 70 - 64>(blkh,                b[1]);	// 3- 3 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits< 1, 76 - 64>(blkh,                b[1]);	// 4- 4 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits< 1,  4 -  0>(blkl,                b[1]);	// 5- 5 set 2 blue end
  
  blkh = CopyBits< 5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x00 - 10 5 5 5
    {  M, 0}, {  M, 1}, {GS2, 4}, {BS2, 4}, {BE2, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 0, partition, s, e, 16 - 10, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m2(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
	s[0] = a[0];
	s[1] = a[1];
	e[0] = b[0];
	e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x01 - 7 6 6 6
        {  M, 0}, {  M, 1}, {GS2, 5}, {GE2, 4}, {GE2, 5}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {BE2, 0}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {BS2, 5}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BE2, 3}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(1);		// 2 mode bit

  // truncate bottom 9 bits
  a[0] = ShiftRightHalf<9>(a[0]);
  b[0] = ShiftRightHalf<9>(b[0]);
  a[1] = ShiftRightHalf<9>(a[1]);
  b[1] = ShiftRightHalf<9>(b[1]);

  blkl = InjtBits<7,  5 -  0>(blkl,                a[0]);	// 1-7 set 1 red start
  blkl = InjtBits<6, 35 -  0>(blkl,                b[0]);	// 1-6 set 1 red end
  blkh = InjtBits<6, 65 - 64>(blkh,                a[1]);	// 1-6 set 2 red start
  blkh = InjtBits<6, 71 - 64>(blkh,                b[1]);	// 1-6 set 2 red end

  // truncate bottom 9 bits
  a[0] = ShiftRightHalf<7 + 9>(a[0]);
  b[0] = ShiftRightHalf<7 + 9>(b[0]);
  a[1] = ShiftRightHalf<7 + 9>(a[1]);
  b[1] = ShiftRightHalf<7 + 9>(b[1]);

  blkl = InjtBits<7, 15 -  0>(blkl,                a[0]);	// 1-7 set 1 green start
  blkl = InjtBits<6, 45 -  0>(blkl,                b[0]);	// 1-6 set 1 green end
  blkl = InjtBits<4, 41 -  0>(blkl,                a[1]);	// 1-4 set 2 green start
  blkl = InjtBits<4, 51 -  0>(blkl,                b[1]);	// 1-4 set 2 green end

  a[1] = ShiftRightHalf<4>(a[1]);
  b[1] = ShiftRightHalf<4>(b[1]);

  blkl = InjtBits<1, 24 -  0>(blkl,                a[1]);	// 5-5 set 2 green start
  blkl = InjtBits<2,  3 -  0>(blkl,                b[1]);	// 5-6 set 2 green end

  a[1] = ShiftRightHalf<1>(a[1]);

  blkl = InjtBits<1,  2 -  0>(blkl,                a[1]);	// 6-6 set 2 green start

  // truncate bottom 9 bits
  a[0] = ShiftRightHalf<7 + 9    >(a[0]);
  b[0] = ShiftRightHalf<7 + 9    >(b[0]);
  a[1] = ShiftRightHalf<7 + 9 - 5>(a[1]);
  b[1] = ShiftRightHalf<7 + 9 - 4>(b[1]);

  blkl = InjtBits<7, 25 -  0>(blkl,                a[0]);	// 1-7 set 1 blue start
  blkl = InjtBits<6, 55 -  0>(blkl,                b[0]);	// 1-6 set 1 blue end
  blkl = InjtBits<3, 61 -  0>(blkl,                a[1]);	// 1-3 set 2 blue start
  blkl = InjtBits<2, 12 -  0>(blkl,                b[1]);	// 1-2 set 2 blue end

  a[1] = ShiftRightHalf<3>(a[1]);
  b[1] = ShiftRightHalf<2>(b[1]);

  blkh = InjtBits<1, 64 - 64>(blkh,                a[1]);	// 4-4 set 2 blue start
  blkl = InjtBits<1, 23 -  0>(blkl,                b[1]);	// 3-3 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 14 -  0>(blkl,                a[1]);	// 5-5 set 2 blue start
  blkl = InjtBits<1, 32 -  0>(blkl,                b[1]);	// 4-4 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 22 -  0>(blkl,                a[1]);	// 6-6 set 2 blue start
  blkl = InjtBits<1, 34 -  0>(blkl,                b[1]);	// 5-5 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 33 -  0>(blkl,                b[1]);	// 6-6 set 2 blue end

  blkh = CopyBits<5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x01 - 7 6 6 6
    {  M, 0}, {  M, 1}, {GS2, 5}, {GE2, 4}, {GE2, 5}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {BE2, 0}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {BS2, 5}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BE2, 3}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 1, partition, s, e, 16 - 7, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m3(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
  s[0] = a[0];
  s[1] = a[1];
  e[0] = b[0];
  e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x02 - 11 5 4 4
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RS1,10}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,10},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,10},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(2);		// 5 mode bit

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<5>(a[0]);
  b[0] = ShiftRightHalf<5>(b[0]);
  a[1] = ShiftRightHalf<5>(a[1]);
  b[1] = ShiftRightHalf<5>(b[1]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits< 5, 35 -  0>(blkl,                b[0]);	// 1- 5 set 1 red end
  blkh = InjtBits< 5, 65 - 64>(blkh,                a[1]);	// 1- 5 set 2 red start
  blkh = InjtBits< 5, 71 - 64>(blkh,                b[1]);	// 1- 5 set 2 red end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 40 -  0>(blkl,                a[0]);	// 11-11 set 1 red start

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);
  a[1] = ShiftRightHalf<11 + 5     >(a[1]);
  b[1] = ShiftRightHalf<11 + 5     >(b[1]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits< 4, 45 -  0>(blkl,                b[0]);	// 1- 4 set 1 green end
  blkl = InjtBits< 4, 41 -  0>(blkl,                a[1]);	// 1- 4 set 2 green start
  blkl = InjtBits< 4, 51 -  0>(blkl,                b[1]);	// 1- 4 set 2 green end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 49 -  0>(blkl,                a[0]);	// 11-11 set 1 green start

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);
  a[1] = ShiftRightHalf<11 + 5     >(a[1]);
  b[1] = ShiftRightHalf<11 + 5     >(b[1]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 4, 55 -  0>(blkl,                b[0]);	// 1- 4 set 1 blue end
  blkl = InjtBits< 3, 61 -  0>(blkl,                a[1]);	// 1- 3 set 2 blue start
  blkl = InjtBits< 1, 50 -  0>(blkl,                b[1]);	// 1- 1 set 2 blue end

  a[0] = ShiftRightHalf<10>(a[0]);
  a[1] = ShiftRightHalf< 3>(a[1]);
  b[1] = ShiftRightHalf< 1>(b[1]);

  blkl = InjtBits< 1, 59 -  0>(blkl,                a[0]);	// 11-11 set 1 blue start
  blkh = InjtBits< 1, 64 - 64>(blkh,                a[1]);	// 4- 4 set 2 blue start
  blkl = InjtBits< 1, 60 -  0>(blkl,                b[1]);	// 2- 2 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits< 1, 70 - 64>(blkh,                b[1]);	// 3- 3 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits< 1, 76 - 64>(blkh,                b[1]);	// 4- 4 set 2 blue end

  blkh = CopyBits< 5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x02 - 11 5 4 4
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {RS1,10}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,10},
    {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,10},
    {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 0x02, partition, s, e, 16 - 11, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m4(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
  s[0] = a[0];
  s[1] = a[1];
  e[0] = b[0];
  e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x06 - 11 4 5 4
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,10},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GS1,10}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,10},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {BE2, 0},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {GS2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(6);		// 5 mode bit

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<5>(a[0]);
  b[0] = ShiftRightHalf<5>(b[0]);
  a[1] = ShiftRightHalf<5>(a[1]);
  b[1] = ShiftRightHalf<5>(b[1]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits< 4, 35 -  0>(blkl,                b[0]);	// 1- 4 set 1 red end
  blkh = InjtBits< 4, 65 - 64>(blkh,                a[1]);	// 1- 4 set 2 red start
  blkh = InjtBits< 4, 71 - 64>(blkh,                b[1]);	// 1- 4 set 2 red end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 39 -  0>(blkl,                a[0]);	// 11-11 set 1 red start

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);
  a[1] = ShiftRightHalf<11 + 5     >(a[1]);
  b[1] = ShiftRightHalf<11 + 5     >(b[1]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits< 5, 45 -  0>(blkl,                b[0]);	// 1- 5 set 1 green end
  blkl = InjtBits< 4, 41 -  0>(blkl,                a[1]);	// 1- 4 set 2 green start
  blkl = InjtBits< 4, 51 -  0>(blkl,                b[1]);	// 1- 4 set 2 green end

  a[0] = ShiftRightHalf<10>(a[0]);
  a[1] = ShiftRightHalf< 4>(a[1]);
  b[1] = ShiftRightHalf< 4>(b[1]);

  blkl = InjtBits< 1, 50 -  0>(blkl,                a[0]);	// 11-11 set 1 green start
  blkh = InjtBits< 1, 75 - 64>(blkh,                a[1]);	// 5-5 set 2 green start
  blkl = InjtBits< 1, 40 -  0>(blkl,                b[1]);	// 5-5 set 2 green end

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);
  a[1] = ShiftRightHalf<11 + 5 -  4>(a[1]);
  b[1] = ShiftRightHalf<11 + 5 -  4>(b[1]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 4, 55 -  0>(blkl,                b[0]);	// 1- 4 set 1 blue end
  blkl = InjtBits< 3, 61 -  0>(blkl,                a[1]);	// 1- 3 set 2 blue start
  blkh = InjtBits< 1, 69 - 64>(blkh,                b[1]);	// 1- 1 set 2 blue end

  a[0] = ShiftRightHalf<10>(a[0]);
  a[1] = ShiftRightHalf< 3>(a[1]);
  b[1] = ShiftRightHalf< 1>(b[1]);

  blkl = InjtBits< 1, 59 -  0>(blkl,                a[0]);	// 11-11 set 1 blue start
  blkh = InjtBits< 1, 64 - 64>(blkh,                a[1]);	// 4- 4 set 2 blue start
  blkl = InjtBits< 1, 60 -  0>(blkl,                b[1]);	// 2- 2 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits< 1, 70 - 64>(blkh,                b[1]);	// 3- 3 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits< 1, 76 - 64>(blkh,                b[1]);	// 4- 4 set 2 blue end

  blkh = CopyBits< 5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
	ModeDescriptor layout[82] = {   // 0x06 - 11 4 5 4
		{  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
		{RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
		{GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
		{BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,10},
		{GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
		{GS1,10}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,10},
		{BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {BE2, 0},
		{BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {GS2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
		{  D, 3}, {  D, 4}
  };

  VerifyHDRBlock(82, layout, 0x06, partition, s, e, 16 - 11, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m5(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
  s[0] = a[0];
  s[1] = a[1];
  e[0] = b[0];
  e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x0a - 11 4 4 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,10},
        {BS2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,10},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BS1,10}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {BE2, 1},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {BE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(10);		// 5 mode bit

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<5>(a[0]);
  b[0] = ShiftRightHalf<5>(b[0]);
  a[1] = ShiftRightHalf<5>(a[1]);
  b[1] = ShiftRightHalf<5>(b[1]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits< 4, 35 -  0>(blkl,                b[0]);	// 1- 4 set 1 red end
  blkh = InjtBits< 4, 65 - 64>(blkh,                a[1]);	// 1- 4 set 2 red start
  blkh = InjtBits< 4, 71 - 64>(blkh,                b[1]);	// 1- 4 set 2 red end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 39 -  0>(blkl,                a[0]);	// 11-11 set 1 red start

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);
  a[1] = ShiftRightHalf<11 + 5     >(a[1]);
  b[1] = ShiftRightHalf<11 + 5     >(b[1]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits< 4, 45 -  0>(blkl,                b[0]);	// 1- 4 set 1 green end
  blkl = InjtBits< 4, 41 -  0>(blkl,                a[1]);	// 1- 4 set 2 green start
  blkl = InjtBits< 4, 51 -  0>(blkl,                b[1]);	// 1- 4 set 2 green end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 49 -  0>(blkl,                a[0]);	// 11-11 set 1 green start

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);
  a[1] = ShiftRightHalf<11 + 5     >(a[1]);
  b[1] = ShiftRightHalf<11 + 5     >(b[1]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 5, 55 -  0>(blkl,                b[0]);	// 1- 5 set 1 blue end
  blkl = InjtBits< 3, 61 -  0>(blkl,                a[1]);	// 1- 3 set 2 blue start
  blkl = InjtBits< 1, 50 -  0>(blkl,                b[1]);	// 1- 1 set 2 blue end

  a[0] = ShiftRightHalf<10>(a[0]);
  a[1] = ShiftRightHalf< 3>(a[1]);
  b[1] = ShiftRightHalf< 1>(b[1]);

  blkl = InjtBits< 1, 60 -  0>(blkl,                a[0]);	// 11-11 set 1 blue start
  blkh = InjtBits< 1, 64 - 64>(blkh,                a[1]);	// 4- 4 set 2 blue start
  blkh = InjtBits< 2, 69 - 64>(blkh,                b[1]);	// 2- 3 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<2>(b[1]);

  blkl = InjtBits< 1, 40 -  0>(blkl,                a[1]);	// 5- 5 set 2 blue start
  blkh = InjtBits< 1, 76 - 64>(blkh,                b[1]);	// 4- 4 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits< 1, 75 - 64>(blkh,                b[1]);	// 5- 5 set 2 blue end

  blkh = CopyBits< 5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x0a - 11 4 4 5
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,10},
    {BS2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,10},
    {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BS1,10}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {BE2, 1},
    {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {BE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4}
  };

  VerifyHDRBlock(82, layout, 0x0a, partition, s, e, 16 - 11, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m6(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
  s[0] = a[0];
  s[1] = a[1];
  e[0] = b[0];
  e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x0e - 9 5 5 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(14);		// 5 mode bit

  // truncate bottom 7 bits
  a[0] = ShiftRightHalf<7>(a[0]);
  b[0] = ShiftRightHalf<7>(b[0]);
  a[1] = ShiftRightHalf<7>(a[1]);
  b[1] = ShiftRightHalf<7>(b[1]);

  blkl = InjtBits<9,  5 -  0>(blkl,                a[0]);	// 1-9 set 1 red start
  blkl = InjtBits<5, 35 -  0>(blkl,                b[0]);	// 1-5 set 1 red end
  blkh = InjtBits<5, 65 - 64>(blkh,                a[1]);	// 1-5 set 2 red start
  blkh = InjtBits<5, 71 - 64>(blkh,                b[1]);	// 1-5 set 2 red end

  // truncate bottom 7 bits
  a[0] = ShiftRightHalf<9 + 7>(a[0]);
  b[0] = ShiftRightHalf<9 + 7>(b[0]);
  a[1] = ShiftRightHalf<9 + 7>(a[1]);
  b[1] = ShiftRightHalf<9 + 7>(b[1]);

  blkl = InjtBits<9, 15 -  0>(blkl,                a[0]);	// 1-9 set 1 green start
  blkl = InjtBits<5, 45 -  0>(blkl,                b[0]);	// 1-5 set 1 green end
  blkl = InjtBits<4, 41 -  0>(blkl,                a[1]);	// 1-4 set 2 green start
  blkl = InjtBits<4, 51 -  0>(blkl,                b[1]);	// 1-4 set 2 green end

  a[1] = ShiftRightHalf<4>(a[1]);
  b[1] = ShiftRightHalf<4>(b[1]);

  blkl = InjtBits<1, 24 -  0>(blkl,                a[1]);	// 5-5 set 2 green start
  blkl = InjtBits<1, 40 -  0>(blkl,                b[1]);	// 5-5 set 2 green start

  // truncate bottom 7 bits
  a[0] = ShiftRightHalf<9 + 7    >(a[0]);
  b[0] = ShiftRightHalf<9 + 7    >(b[0]);
  a[1] = ShiftRightHalf<9 + 7 - 4>(a[1]);
  b[1] = ShiftRightHalf<9 + 7 - 4>(b[1]);

  blkl = InjtBits<9, 25 -  0>(blkl,                a[0]);	// 1-9 set 1 blue start
  blkl = InjtBits<5, 55 -  0>(blkl,                b[0]);	// 1-5 set 1 blue end
  blkl = InjtBits<3, 61 -  0>(blkl,                a[1]);	// 1-3 set 2 blue start
  blkl = InjtBits<1, 50 -  0>(blkl,                b[1]);	// 1-1 set 2 blue end

  a[1] = ShiftRightHalf<3>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits<1, 64 - 64>(blkh,                a[1]);	// 4-4 set 2 blue start
  blkl = InjtBits<1, 60 -  0>(blkl,                b[1]);	// 2-2 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 14 -  0>(blkl,                a[1]);	// 5-5 set 2 blue start
  blkh = InjtBits<1, 70 - 64>(blkh,                b[1]);	// 3-3 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits<1, 76 - 64>(blkh,                b[1]);	// 4-4 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 34 -  0>(blkl,                b[1]);	// 5-5 set 2 blue end

  blkh = CopyBits<5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x0e - 9 5 5 5
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 0x0e, partition, s, e, 16 - 9, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m7(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
	s[0] = a[0];
	s[1] = a[1];
	e[0] = b[0];
	e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x12 - 8 6 5 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {GE2, 4}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BE2, 3}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(18);		// 5 mode bit

  // truncate bottom 8 bits
  a[0] = ShiftRightHalf<8>(a[0]);
  b[0] = ShiftRightHalf<8>(b[0]);
  a[1] = ShiftRightHalf<8>(a[1]);
  b[1] = ShiftRightHalf<8>(b[1]);

  blkl = InjtBits<8,  5 -  0>(blkl,                a[0]);	// 1-8 set 1 red start
  blkl = InjtBits<6, 35 -  0>(blkl,                b[0]);	// 1-6 set 1 red end
  blkh = InjtBits<6, 65 - 64>(blkh,                a[1]);	// 1-6 set 2 red start
  blkh = InjtBits<6, 71 - 64>(blkh,                b[1]);	// 1-6 set 2 red end

  // truncate bottom 8 bits
  a[0] = ShiftRightHalf<8 + 8>(a[0]);
  b[0] = ShiftRightHalf<8 + 8>(b[0]);
  a[1] = ShiftRightHalf<8 + 8>(a[1]);
  b[1] = ShiftRightHalf<8 + 8>(b[1]);

  blkl = InjtBits<8, 15 -  0>(blkl,                a[0]);	// 1-8 set 1 green start
  blkl = InjtBits<5, 45 -  0>(blkl,                b[0]);	// 1-5 set 1 green end
  blkl = InjtBits<4, 41 -  0>(blkl,                a[1]);	// 1-4 set 2 green start
  blkl = InjtBits<4, 51 -  0>(blkl,                b[1]);	// 1-4 set 2 green end

  a[1] = ShiftRightHalf<4>(a[1]);
  b[1] = ShiftRightHalf<4>(b[1]);

  blkl = InjtBits<1, 24 -  0>(blkl,                a[1]);	// 5-5 set 2 green start
  blkl = InjtBits<1, 13 -  0>(blkl,                b[1]);	// 5-5 set 2 green end

  // truncate bottom 8 bits
  a[0] = ShiftRightHalf<8 + 8    >(a[0]);
  b[0] = ShiftRightHalf<8 + 8    >(b[0]);
  a[1] = ShiftRightHalf<8 + 8 - 4>(a[1]);
  b[1] = ShiftRightHalf<8 + 8 - 4>(b[1]);

  blkl = InjtBits<8, 25 -  0>(blkl,                a[0]);	// 1-8 set 1 blue start
  blkl = InjtBits<5, 55 -  0>(blkl,                b[0]);	// 1-5 set 1 blue end
  blkl = InjtBits<3, 61 -  0>(blkl,                a[1]);	// 1-3 set 2 blue start
  blkl = InjtBits<1, 50 -  0>(blkl,                b[1]);	// 1-1 set 2 blue end

  a[1] = ShiftRightHalf<3>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits<1, 64 - 64>(blkh,                a[1]);	// 4-4 set 2 blue start
  blkl = InjtBits<1, 60 -  0>(blkl,                b[1]);	// 2-2 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 14 -  0>(blkl,                a[1]);	// 5-5 set 2 blue start
  blkl = InjtBits<1, 23 -  0>(blkl,                b[1]);	// 3-3 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<2, 33 -  0>(blkl,                b[1]);	// 4-5 set 2 blue end

  blkh = CopyBits<5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x12 - 8 6 5 5
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {GE2, 4}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BE2, 3}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 0x12, partition, s, e, 16 - 8, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m8(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
	s[0] = a[0];
	s[1] = a[1];
	e[0] = b[0];
	e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x16 - 8 5 6 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {BE2, 0}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS2, 5}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {GE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(22);		// 5 mode bit

  // truncate bottom 8 bits
  a[0] = ShiftRightHalf<8>(a[0]);
  b[0] = ShiftRightHalf<8>(b[0]);
  a[1] = ShiftRightHalf<8>(a[1]);
  b[1] = ShiftRightHalf<8>(b[1]);

  blkl = InjtBits<8,  5 -  0>(blkl,                a[0]);	// 1-8 set 1 red start
  blkl = InjtBits<5, 35 -  0>(blkl,                b[0]);	// 1-5 set 1 red end
  blkh = InjtBits<5, 65 - 64>(blkh,                a[1]);	// 1-5 set 2 red start
  blkh = InjtBits<5, 71 - 64>(blkh,                b[1]);	// 1-5 set 2 red end

  // truncate bottom 7 bits
  a[0] = ShiftRightHalf<8 + 8>(a[0]);
  b[0] = ShiftRightHalf<8 + 8>(b[0]);
  a[1] = ShiftRightHalf<8 + 8>(a[1]);
  b[1] = ShiftRightHalf<8 + 8>(b[1]);

  blkl = InjtBits<8, 15 -  0>(blkl,                a[0]);	// 1-8 set 1 green start
  blkl = InjtBits<6, 45 -  0>(blkl,                b[0]);	// 1-6 set 1 green end
  blkl = InjtBits<4, 41 -  0>(blkl,                a[1]);	// 1-4 set 2 green start
  blkl = InjtBits<4, 51 -  0>(blkl,                b[1]);	// 1-4 set 2 green end

  a[1] = ShiftRightHalf<4>(a[1]);
  b[1] = ShiftRightHalf<4>(b[1]);

  blkl = InjtBits<1, 24 -  0>(blkl,                a[1]);	// 5-5 set 2 green start
  blkl = InjtBits<1, 40 -  0>(blkl,                b[1]);	// 5-5 set 2 green end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 23 -  0>(blkl,                a[1]);	// 6-6 set 2 green start
  blkl = InjtBits<1, 33 -  0>(blkl,                b[1]);	// 6-6 set 2 green end

  // truncate bottom 8 bits
  a[0] = ShiftRightHalf<8 + 8    >(a[0]);
  b[0] = ShiftRightHalf<8 + 8    >(b[0]);
  a[1] = ShiftRightHalf<8 + 8 - 5>(a[1]);
  b[1] = ShiftRightHalf<8 + 8 - 5>(b[1]);

  blkl = InjtBits<8, 25 -  0>(blkl,                a[0]);	// 1-8 set 1 blue start
  blkl = InjtBits<5, 55 -  0>(blkl,                b[0]);	// 1-5 set 1 blue end
  blkl = InjtBits<3, 61 -  0>(blkl,                a[1]);	// 1-3 set 2 blue start
  blkl = InjtBits<1, 13 -  0>(blkl,                b[1]);	// 1-1 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 60 -  0>(blkl,                b[1]);	// 2-2 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits<1, 70 - 64>(blkh,                b[1]);	// 3-3 set 2 blue end

  a[1] = ShiftRightHalf<3>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits<1, 64 - 64>(blkh,                a[1]);	// 4-4 set 2 blue start
  blkh = InjtBits<1, 76 - 64>(blkh,                b[1]);	// 4-4 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 14 -  0>(blkl,                a[1]);	// 5-5 set 2 blue start
  blkl = InjtBits<1, 34 -  0>(blkl,                b[1]);	// 5-5 set 2 blue end

  blkh = CopyBits<5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x16 - 8 5 6 5
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {BE2, 0}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS2, 5}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {GE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 0x16, partition, s, e, 16 - 8, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_m9(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  a[1] = s[1] - s[0];
  b[0] = e[0] - s[0];
  b[1] = e[1] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  a[1] = PackSShorts(a[1]);
  b[0] = PackSShorts(b[0]);
  b[1] = PackSShorts(b[1]);

#ifndef NDEBUG
	s[0] = a[0];
	s[1] = a[1];
	e[0] = b[0];
	e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x1a - 8 5 5 6
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {BS2, 5}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(26);		// 5 mode bit

  // truncate bottom 8 bits
  a[0] = ShiftRightHalf<8>(a[0]);
  b[0] = ShiftRightHalf<8>(b[0]);
  a[1] = ShiftRightHalf<8>(a[1]);
  b[1] = ShiftRightHalf<8>(b[1]);

  blkl = InjtBits<8,  5 -  0>(blkl,                a[0]);	// 1-8 set 1 red start
  blkl = InjtBits<5, 35 -  0>(blkl,                b[0]);	// 1-5 set 1 red end
  blkh = InjtBits<5, 65 - 64>(blkh,                a[1]);	// 1-5 set 2 red start
  blkh = InjtBits<5, 71 - 64>(blkh,                b[1]);	// 1-5 set 2 red end

  // truncate bottom 7 bits
  a[0] = ShiftRightHalf<8 + 8>(a[0]);
  b[0] = ShiftRightHalf<8 + 8>(b[0]);
  a[1] = ShiftRightHalf<8 + 8>(a[1]);
  b[1] = ShiftRightHalf<8 + 8>(b[1]);

  blkl = InjtBits<8, 15 -  0>(blkl,                a[0]);	// 1-8 set 1 green start
  blkl = InjtBits<5, 45 -  0>(blkl,                b[0]);	// 1-5 set 1 green end
  blkl = InjtBits<4, 41 -  0>(blkl,                a[1]);	// 1-4 set 2 green start
  blkl = InjtBits<4, 51 -  0>(blkl,                b[1]);	// 1-4 set 2 green end

  a[1] = ShiftRightHalf<4>(a[1]);
  b[1] = ShiftRightHalf<4>(b[1]);

  blkl = InjtBits<1, 24 -  0>(blkl,                a[1]);	// 5-5 set 2 green start
  blkl = InjtBits<1, 40 -  0>(blkl,                b[1]);	// 5-5 set 2 green end

  // truncate bottom 8 bits
  a[0] = ShiftRightHalf<8 + 8    >(a[0]);
  b[0] = ShiftRightHalf<8 + 8    >(b[0]);
  a[1] = ShiftRightHalf<8 + 8 - 4>(a[1]);
  b[1] = ShiftRightHalf<8 + 8 - 4>(b[1]);

  blkl = InjtBits<8, 25 -  0>(blkl,                a[0]);	// 1-8 set 1 blue start
  blkl = InjtBits<6, 55 -  0>(blkl,                b[0]);	// 1-6 set 1 blue end
  blkl = InjtBits<3, 61 -  0>(blkl,                a[1]);	// 1-3 set 2 blue start
  blkl = InjtBits<1, 50 -  0>(blkl,                b[1]);	// 1-1 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 13 -  0>(blkl,                b[1]);	// 2-2 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits<1, 70 - 64>(blkh,                b[1]);	// 3-3 set 2 blue end

  a[1] = ShiftRightHalf<3>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkh = InjtBits<1, 64 - 64>(blkh,                a[1]);	// 4-4 set 2 blue start
  blkh = InjtBits<1, 76 - 64>(blkh,                b[1]);	// 4-4 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 14 -  0>(blkl,                a[1]);	// 5-5 set 2 blue start
  blkl = InjtBits<1, 34 -  0>(blkl,                b[1]);	// 5-5 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 23 -  0>(blkl,                a[1]);	// 6-6 set 2 blue start
  blkl = InjtBits<1, 33 -  0>(blkl,                b[1]);	// 6-6 set 2 blue end

  blkh = CopyBits<5, 77 - 64>(blkh, blkh.SetLong(partition));// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x1a - 8 5 5 6
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {BS2, 5}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 0x1a, partition, s, e, 16 - 8, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_mA(int partition, Col4 const (&start)[2], Col4 const (&end)[2], u8 const (&indices)[1][16], void* block)
{
  Col4 s[2] = {start[0], start[1]};
  Col4 e[2] = {end  [0], end  [1]};
  Col4 a[2];
  Col4 b[2];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapHDRBlock<3>(partition, s, e, idxs, indices);

  // create explicits
  a[0] = s[0];
  a[1] = s[1];
  b[0] = e[0];
  b[1] = e[1];

  // pack explicits
  a[0] = PackUShorts(a[0]);
  b[0] = PackUShorts(b[0]);
  a[1] = PackUShorts(a[1]);
  b[1] = PackUShorts(b[1]);

#ifndef NDEBUG
	s[0] = a[0];
	s[1] = a[1];
	e[0] = b[0];
	e[1] = b[1];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x1e - 6 6 6 6
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {GE2, 4}, {BE2, 0}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS2, 5}, {BS2, 5}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {GE2, 5}, {BE2, 3}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  blkl = blkl.SetLong(30);		// 5 mode bit (01111)

  // truncate bottom 10 bits
  a[0] = ShiftRightHalf<10>(a[0]);
  b[0] = ShiftRightHalf<10>(b[0]);
  a[1] = ShiftRightHalf<10>(a[1]);
  b[1] = ShiftRightHalf<10>(b[1]);

  blkl = InjtBits<6,  5 -  0>(blkl,                a[0]);	// 1-6 set 1 red start
  blkl = InjtBits<6, 35 -  0>(blkl,                b[0]);	// 1-6 set 1 red end
  blkh = InjtBits<6, 65 - 64>(blkh,                a[1]);	// 1-6 set 2 red start
  blkh = InjtBits<6, 71 - 64>(blkh,                b[1]);	// 1-6 set 2 red end

  // truncate bottom 10 bits
  a[0] = ShiftRightHalf<6 + 10>(a[0]);
  b[0] = ShiftRightHalf<6 + 10>(b[0]);
  a[1] = ShiftRightHalf<6 + 10>(a[1]);
  b[1] = ShiftRightHalf<6 + 10>(b[1]);

  blkl = InjtBits<6, 15 -  0>(blkl,                a[0]);	// 1-6 set 1 green start
  blkl = InjtBits<6, 45 -  0>(blkl,                b[0]);	// 1-6 set 1 green end
  blkl = InjtBits<4, 41 -  0>(blkl,                a[1]);	// 1-4 set 2 green start
  blkl = InjtBits<4, 51 -  0>(blkl,                b[1]);	// 1-4 set 2 green end

  a[1] = ShiftRightHalf<4>(a[1]);
  b[1] = ShiftRightHalf<4>(b[1]);

  blkl = InjtBits<1, 24 -  0>(blkl,                a[1]);	// 5-5 set 2 green start
  blkl = InjtBits<1, 11 -  0>(blkl,                b[1]);	// 5-5 set 2 green end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 21 -  0>(blkl,                a[1]);	// 6-6 set 2 green start
  blkl = InjtBits<1, 31 -  0>(blkl,                b[1]);	// 6-6 set 2 green end

  // truncate bottom 10 bits
  a[0] = ShiftRightHalf<6 + 10    >(a[0]);
  b[0] = ShiftRightHalf<6 + 10    >(b[0]);
  a[1] = ShiftRightHalf<6 + 10 - 5>(a[1]);
  b[1] = ShiftRightHalf<6 + 10 - 5>(b[1]);

  blkl = InjtBits<6, 25 -  0>(blkl,                a[0]);	// 1-6 set 1 blue start
  blkl = InjtBits<6, 55 -  0>(blkl,                b[0]);	// 1-6 set 1 blue end
  blkl = InjtBits<3, 61 -  0>(blkl,                a[1]);	// 1-3 set 2 blue start
  blkl = InjtBits<2, 12 -  0>(blkl,                b[1]);	// 1-2 set 2 blue end

  a[1] = ShiftRightHalf<3>(a[1]);
  b[1] = ShiftRightHalf<2>(b[1]);

  blkh = InjtBits<1, 64 - 64>(blkh,                a[1]);	// 4-4 set 2 blue start
  blkl = InjtBits<1, 23 -  0>(blkl,                b[1]);	// 3-3 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 14 -  0>(blkl,                a[1]);	// 5-5 set 2 blue start
  blkl = InjtBits<1, 32 -  0>(blkl,                b[1]);	// 4-4 set 2 blue end

  a[1] = ShiftRightHalf<1>(a[1]);
  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 22 -  0>(blkl,                a[1]);	// 6-6 set 2 blue start
  blkl = InjtBits<1, 34 -  0>(blkl,                b[1]);	// 5-5 set 2 blue end

  b[1] = ShiftRightHalf<1>(b[1]);

  blkl = InjtBits<1, 33 -  0>(blkl,                b[1]);	// 6-6 set 2 blue end

  blkh = CopyBits<5, 77 - 64>(blkh, blkh.SetLong(partition));	// 5 bit partition

#ifndef NDEBUG
  ModeDescriptor layout[82] = {   // 0x1e - 6 6 6 6
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {GE2, 4}, {BE2, 0}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS2, 5}, {BS2, 5}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {GE2, 5}, {BE2, 3}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
    {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
    {  D, 3}, {  D, 4},
  };

  VerifyHDRBlock(82, layout, 0x1e, partition, s, e, 16 - 6, blkl, blkh);
#endif // NDEBUG

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WriteHDRBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_mB(               Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block)
{
  Col4 s[1] = {start[0]};
  Col4 e[1] = {end  [0]};
  Col4 a[1];
  Col4 b[1];
  Col4 blkl, blkh;
  Col4 idxs[1];
  int partition = 0;

  // remap the indices
  RemapHDRBlock<4>(partition, s, e, idxs, indices);

  // create explicits
  a[0] = s[0];
  b[0] = e[0];

  // pack explicits
  a[0] = PackUShorts(a[0]);
  b[0] = PackUShorts(b[0]);

#ifndef NDEBUG
	s[0] = a[0];
	e[0] = b[0];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x03 - 10 10
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {RE1, 6}, {RE1, 7}, {RE1, 8}, {RE1, 9}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE1, 6}, {GE1, 7}, {GE1, 8}, {GE1, 9}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BE1, 6}, {BE1, 7}, {BE1, 8}, {BE1, 9}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0},
    },
   */

  blkl = blkl.SetLong(3);		// 5 mode bit (11000), 0 bit partition

  // truncate bottom 6 bits
  a[0] = ShiftRightHalf<6>(a[0]);
  b[0] = ShiftRightHalf<6>(b[0]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits<10, 35 -  0>(blkl,                b[0]);	// 1-10 set 1 red end

  // truncate bottom 6 bits
  a[0] = ShiftRightHalf<10 + 6>(a[0]);
  b[0] = ShiftRightHalf<10 + 6>(b[0]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits<10, 45 -  0>(blkl,                b[0]);	// 1-10 set 1 green end

  // truncate bottom 6 bits
  a[0] = ShiftRightHalf<10 + 6>(a[0]);
  b[0] = ShiftRightHalf<10 + 6>(b[0]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 9, 55 -  0>(blkl,                b[0]);	// 1-9 set 1 blue end

    b[0] = ShiftRightHalf<9>(b[0]);

    blkh = InjtBits< 1, 64 - 64>(blkh,                b[0]);	// 10-10 set 1 blue end

#ifndef NDEBUG
  ModeDescriptor layout[65] = {   // 0x03 - 10 10
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {RE1, 5}, {RE1, 6}, {RE1, 7}, {RE1, 8}, {RE1, 9}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {GE1, 5}, {GE1, 6}, {GE1, 7}, {GE1, 8}, {GE1, 9}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE1, 5}, {BE1, 6}, {BE1, 7}, {BE1, 8}, {BE1, 9}
  };

  VerifyHDRBlock(65, layout, 0x03, partition, s, e, 16 - 10, blkl, blkh);
#endif // NDEBUG

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  WriteHDRBlock<1, 4, 65>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_mC(               Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block)
{
  Col4 s[1] = {start[0]};
  Col4 e[1] = {end  [0]};
  Col4 a[1];
  Col4 b[1];
  Col4 blkl, blkh;
  Col4 idxs[1];
  int partition = 0;

  // remap the indices
  RemapHDRBlock<4>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  b[0] = e[0] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  b[0] = PackSShorts(b[0]);

#ifndef NDEBUG
	s[0] = a[0];
	e[0] = b[0];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x07 - 11 9
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {RE1, 6}, {RE1, 7}, {RE1, 8}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE1, 6}, {GE1, 7}, {GE1, 8}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BE1, 6}, {BE1, 7}, {BE1, 8}, {BS1,10}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0},
    },
   */

  blkl = blkl.SetLong(7);		// 5 mode bit (11100), 0 bit partition

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<5>(a[0]);
  b[0] = ShiftRightHalf<5>(b[0]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits< 9, 35 -  0>(blkl,                b[0]);	// 1-9 set 1 red end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 44 -  0>(blkl,                a[0]);	// 11-11 set 1 red start

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits< 9, 45 -  0>(blkl,                b[0]);	// 1-9 set 1 green end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 54 -  0>(blkl,                a[0]);	// 11-11 set 1 green start

  // truncate bottom 5 bits
  a[0] = ShiftRightHalf<11 + 5 - 10>(a[0]);
  b[0] = ShiftRightHalf<11 + 5     >(b[0]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 9, 55 -  0>(blkl,                b[0]);	// 1-9 set 1 blue end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkh = InjtBits< 1, 64 - 64>(blkh,                a[0]);	// 11-11 set 1 blue start

#ifndef NDEBUG
  ModeDescriptor layout[65] = {   // 0x07 - 11 9
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {RE1, 5}, {RE1, 6}, {RE1, 7}, {RE1, 8}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {GE1, 5}, {GE1, 6}, {GE1, 7}, {GE1, 8}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE1, 5}, {BE1, 6}, {BE1, 7}, {BE1, 8}, {BS1,10}
  };

  VerifyHDRBlock(65, layout, 0x07, partition, s, e, 16 - 11, blkl, blkh);
#endif // NDEBUG

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  WriteHDRBlock<1, 4, 65>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_mD(               Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block)
{
  Col4 s[1] = {start[0]};
  Col4 e[1] = {end  [0]};
  Col4 a[1];
  Col4 b[1];
  Col4 blkl, blkh;
  Col4 idxs[1];
  int partition = 0;

  // remap the indices
  RemapHDRBlock<4>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  b[0] = e[0] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  b[0] = PackSShorts(b[0]);

#ifndef NDEBUG
	s[0] = a[0];
	e[0] = b[0];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x0b - 12 8
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {RE1, 6}, {RE1, 7}, {RS1,11}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE1, 6}, {GE1, 7}, {GS1,11}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BE1, 6}, {BE1, 7}, {BS1,11}, {BS1,10}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0},
    },
   */

  blkl = blkl.SetLong(11);		// 5 mode bit (11010), 0 bit partition

  // truncate bottom 4 bits
  a[0] = ShiftRightHalf<4>(a[0]);
  b[0] = ShiftRightHalf<4>(b[0]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits< 8, 35 -  0>(blkl,                b[0]);	// 1-8 set 1 red end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 44 -  0>(blkl,                a[0]);	// 11-11 set 1 red start

  a[0] = ShiftRightHalf< 1>(a[0]);

  blkl = InjtBits< 1, 43 -  0>(blkl,                a[0]);	// 12-12 set 1 red start

  // truncate bottom 4 bits
  a[0] = ShiftRightHalf<12 + 4 - 11>(a[0]);
  b[0] = ShiftRightHalf<12 + 4     >(b[0]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits< 8, 45 -  0>(blkl,                b[0]);	// 1-8 set 1 green end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkl = InjtBits< 1, 54 -  0>(blkl,                a[0]);	// 11-11 set 1 green start

  a[0] = ShiftRightHalf< 1>(a[0]);

  blkl = InjtBits< 1, 53 -  0>(blkl,                a[0]);	// 12-12 set 1 green start

  // truncate bottom 4 bits
  a[0] = ShiftRightHalf<12 + 4 - 11>(a[0]);
  b[0] = ShiftRightHalf<12 + 4     >(b[0]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 8, 55 -  0>(blkl,                b[0]);	// 1-8 set 1 blue end

  a[0] = ShiftRightHalf<10>(a[0]);

  blkh = InjtBits< 1, 64 - 64>(blkh,                a[0]);	// 11-11 set 1 blue start

  a[0] = ShiftRightHalf< 1>(a[0]);

  blkl = InjtBits< 1, 63 -  0>(blkl,                a[0]);	// 12-12 set 1 blue start

#ifndef NDEBUG
  ModeDescriptor layout[65] = {   // 0x0b - 12 8
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
    {RE1, 5}, {RE1, 6}, {RE1, 7}, {RS1,11}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
    {GE1, 5}, {GE1, 6}, {GE1, 7}, {GS1,11}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
    {BE1, 5}, {BE1, 6}, {BE1, 7}, {BS1,11}, {BS1,10}
  };

  VerifyHDRBlock(65, layout, 0x0b, partition, s, e, 16 - 12, blkl, blkh);
#endif // NDEBUG

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  WriteHDRBlock<1, 4, 65>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WriteHDRBlock_mE(               Col4 const (&start)[1], Col4 const (&end)[1], u8 const (&indices)[1][16], void* block)
{
  Col4 s[1] = {start[0]};
  Col4 e[1] = {end  [0]};
  Col4 a[1];
  Col4 b[1];
  Col4 t[1];
  Col4 blkl, blkh;
  Col4 idxs[1];
  int partition = 0;

  // remap the indices
  RemapHDRBlock<4>(partition, s, e, idxs, indices);

  // create deltas
  a[0] =        s[0];
  b[0] = e[0] - s[0];

  // pack deltas
  a[0] = PackUShorts(a[0]);
  b[0] = PackSShorts(b[0]);

#ifndef NDEBUG
	s[0] = a[0];
	e[0] = b[0];
#endif

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
    {   // 0x0f - 16 4
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,15},
        {RS1,14}, {RS1,13}, {RS1,12}, {RS1,11}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,15},
        {GS1,14}, {GS1,13}, {GS1,12}, {GS1,11}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,15},
        {BS1,14}, {BS1,13}, {BS1,12}, {BS1,11}, {BS1,10}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0}, { NA, 0},
        { NA, 0}, { NA, 0},
    },
   */

  blkl = blkl.SetLong(15);		// 5 mode bit (11110), 0 bit partition

  // truncate bottom 0 bits
//a[0] = ShiftRightHalf<0>(a[0]);
//b[0] = ShiftRightHalf<0>(b[0]);

  blkl = InjtBits<10,  5 -  0>(blkl,                a[0]);	// 1-10 set 1 red start
  blkl = InjtBits< 4, 35 -  0>(blkl,                b[0]);	// 1-4 set 1 red end

  a[0] = ShiftRightHalf<10 - 2>(a[0]);
  t[0] = RevsBits(a[0]);

  blkl = InjtBits< 6, 39 -  0>(blkl,                t[0]);	// 16-11 set 1 red start

  // truncate bottom 0 bits
  a[0] = ShiftRightHalf<16 - 10 + 2>(a[0]);
  b[0] = ShiftRightHalf<16         >(b[0]);

  blkl = InjtBits<10, 15 -  0>(blkl,                a[0]);	// 1-10 set 1 green start
  blkl = InjtBits< 4, 45 -  0>(blkl,                b[0]);	// 1-4 set 1 green end

  a[0] = ShiftRightHalf<10 - 2>(a[0]);
  t[0] = RevsBits(a[0]);

  blkl = InjtBits< 6, 49 -  0>(blkl,                t[0]);	// 16-11 set 1 green start

  // truncate bottom 0 bits
  a[0] = ShiftRightHalf<16 - 10 + 2>(a[0]);
  b[0] = ShiftRightHalf<16         >(b[0]);

  blkl = InjtBits<10, 25 -  0>(blkl,                a[0]);	// 1-10 set 1 blue start
  blkl = InjtBits< 4, 55 -  0>(blkl,                b[0]);	// 1-4 set 1 blue end

  a[0] = ShiftRightHalf<10 - 2>(a[0]);
  t[0] = RevsBits(a[0]);

  blkl = InjtBits< 5, 59 -  0>(blkl,                t[0]);	// 16-12 set 1 blue start

  t[0] = ShiftRightHalf<5>(t[0]);

  blkh = InjtBits< 1, 64 - 64>(blkh,                t[0]);	// 11-11 set 1 blue start

#ifndef NDEBUG
  ModeDescriptor layout[65] = {   // 0x0f - 16 4
    {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
    {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
    {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
    {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,15},
    {RS1,14}, {RS1,13}, {RS1,12}, {RS1,11}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,15},
    {GS1,14}, {GS1,13}, {GS1,12}, {GS1,11}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,15},
    {BS1,14}, {BS1,13}, {BS1,12}, {BS1,11}, {BS1,10}
  };

  VerifyHDRBlock(65, layout, 0x0f, partition, s, e, 16 - 16, blkl, blkh);
#endif // NDEBUG

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  WriteHDRBlock<1, 4, 65>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

/* *****************************************************************************
 */
extern const u8 whichsetinpartition[64][/*3*/2][16];

/* -----------------------------------------------------------------------------
 */
template<const int sets, const int ibits, const int begin>
void ReadHDRBlock(int partition, __int64 *codes, Col4 &blkl, Col4 &blkh, __int64 *out)
{
  Col4 iblk;

  /* none of the cases straddles the lo to hi border, all go into hi
   *
   * WriteHDRBlock<2, 3, 82>(partition, remapped, blkl, blkh);
   * WriteHDRBlock<1, 4, 65>(partition, remapped, blkl, blkh);
   */

  // one bit is omitted per set, so it's one instruction per set + 1
  assume(sets >= 1 && sets <= 3);
  switch(sets) {
    case 1: {
      // always index 0
      iblk = ShiftLeftHalf<ibits>(
        ExtrBits<ibits * 15 - 0, begin - 64 + ibits *  1 - 1>(blkh)) |
        ExtrBits<ibits *  1 - 1, begin - 64 +              0>(blkh)
      ;
    } break;
    case 2: {
      // always index 0
      iblk = ShiftLeftHalf<ibits>(
        ExtrBits<ibits * 15 - 0, begin - 64 + ibits *  1 - 1>(blkh)) |
        ExtrBits<ibits *  1 - 1, begin - 64 +              0>(blkh)
      ;

      // if it's not the last bit which is cut (low probability)
      int len, bgn = begin - 64;
      if ((len = shorterindex[partition][0]) < 15) {
        len = (len * ibits) + ibits;
        bgn = (bgn + len - 2);

        /* no obvious mask-pattern, set of 9 distinct lens
         *
         * bgn=82-64 len={2,6,8}*3={6,18,24}+b={ 9,21,27}-1={ 8,20,26};
         * bgn=65-64 len={2,6,8}*4={8,24,32}+b={12,32,40}-1={11,31,39};
         */

      	// remaining length can be anything, length overflow is silent
        blkh = ExtrBits(blkh, ibits * 14, bgn);

        // do the whole shift
	iblk = ShiftLeftHalf(blkh, len) |
	  MaskBits(iblk, len - 1, 0)
	;
      }
    } break;
    case 3: {
      iblk = ShiftLeftHalf<ibits>(
        ExtrBits<ibits * 15 - 0, begin - 64 + ibits *  1 - 1>(blkh)) |
        ExtrBits<ibits *  1 - 1, begin - 64 +              0>(blkh)
      ;

      int bgn = begin - 64;
      int fln, fbn;
      fln = shorterindex[partition][3]; {
        fln = (fln * ibits) + ibits;
        fbn = (bgn + fln - 2);

      	// remaining length can be anything, length overflow is silent
        blkh = ExtrBits(blkh, ibits * 14, fbn);

        // do the whole shift
	iblk = ShiftLeftHalf(blkh, fln) |
	  MaskBits(iblk, fln - 1, 0)
	;
      }

      // if it's not the last bit which is cut (low probability)
      int sln, sbn;
      if ((sln = shorterindex[partition][4]) < 15) {
        sln = (sln * ibits) + ibits;
        sbn = (bgn + sln - 3);

      	// subtract the already conducted shifts
        // remaining length can be anything, length overflow is silent
        blkh = ExtrBits(blkh, ibits * 14, sbn - fbn);

        // do the whole shift
	iblk = ShiftLeftHalf(blkh, sln) |
	  MaskBits(iblk, sln - 1, 0)
	;
      }
    } break;
  }

  // store out the palettes
  u8 const *ios = (sets >= 2 ? whichsetinpartition[partition][sets - 2] : NULL);

  // max 16*4 -> 64bits
  // no scattered writes on SSE
  // movd: 3clk 1/1, extrq: 2clk 1/2
  Col4 c0, c1, c2, c3,
       c4, c5, c6, c7;

  // m1:
  //  c00 = 0 -> [0][0] -> 95 cd e1 ff, c08 = 6 -> [0][6] -> 80 c3 d0 ff
  //  c01 = 4 -> [0][4] -> 87 c6 d5 ff, c09 = 4 -> [1][4] -> 79 c3 cb ff
  //  c02 = 7 -> [0][7] -> 7c c1 cd ff, c10 = 2 -> [1][2] -> 80 c3 cb ff
  //  c03 = 7 -> [0][7] -> 7c c1 cd ff, c11 = 0 -> [1][0] -> 87 c3 cb ff
  //  c04 = 2 -> [0][2] -> 8e ca db ff, c12 = 3 -> [1][3] -> 7c c3 cb ff
  //  c05 = 7 -> [0][7] -> 7c c1 cd ff, c13 = 4 -> [1][4] -> 79 c3 cb ff
  //  c06 = 6 -> [0][6] -> 80 c3 d0 ff, c14 = 3 -> [1][3] -> 7c c3 cb ff
  //  c07 = 0 -> [1][0] -> 87 c3 cb ff, c15 = 2 -> [1][2] -> 80 c3 cb ff

  // m3:
  //  c00 = 0 -> [0][0] -> 92 c8 dc ff, c08 = 3 -> [1][3] -> 8c c8 dc ff
  //  c01 = 2 -> [1][2] -> 8b c6 d9 ff, c09 = 3 -> [1][3] -> 8c c8 dc ff
  //  c02 = 0 -> [1][0] -> 8a c2 d2 ff, c10 = 2 -> [0][2] -> 87 c4 d4 ff
  //  c03 = 1 -> [1][1] -> 8b c4 d5 ff, c11 = 3 -> [0][3] -> 82 c2 d0 ff
  //  c04 = 1 -> [0][1] -> 8d c6 d8 ff, c12 = 3 -> [1][3] -> 8c c8 dc ff
  //  c05 = 1 -> [0][1] -> 8d c6 d8 ff, c13 = 3 -> [1][3] -> 8c c8 dc ff
  //  c06 = 1 -> [1][1] -> 8b c4 d5 ff, c14 = 3 -> [1][3] -> 8c c8 dc ff
  //  c07 = 0 -> [1][0] -> 8a c2 d2 ff, c15 = 3 -> [0][3] -> 82 c2 d0 ff

  // m7:
  //  c00 = 6 -> [-][6] -> cb f1 fa ff, c08 = b -> [-][b] -> d4 f4 fb fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c01 = 5 -> [-][5] -> c9 f0 fa ff, c09 = 8 -> [-][8] -> ce f2 fb fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c02 = 0 -> [-][0] -> bf ed f9 ff, c10 = 3 -> [-][3] -> c5 ef fa ff
  //  c03 = 1 -> [-][1] -> c1 ee f9 ff, c11 = 0 -> [-][0] -> bf ed f9 ff
  //  c04 = 7 -> [-][7] -> cd f1 fa ff, c12 = e -> [-][e] -> da f5 fc fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c05 = 5 -> [-][5] -> c9 f0 fa ff, c13 = b -> [-][b] -> d4 f4 fb fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c06 = 4 -> [-][4] -> c7 ef fa ff, c14 = 6 -> [-][6] -> cb f1 fa ff
  //  c07 = 2 -> [-][2] -> c3 ee f9 ff, c15 = 0 -> [-][0] -> bf ed f9 ff

  // m0:
  //  c00 = 3 -> [0][3] -> 8c d7 e1 ff, c08 = 5 -> [0][5] -> 8c d3 d7 ff
  //  c01 = 3 -> [0][3] -> 8c d7 e1 ff, c09 = 6 -> [0][6] -> 8c d0 d3 ff
  //  c02 = 4 -> [1][4] -> 84 c6 d9 ff, c10 = 0 -> [1][0] -> 84 c6 c6 ff
  //  c03 = 3 -> [2][3] -> 95 c9 d0 ff, c11 = 0 -> [2][0] -> 84 c6 c6 ff
  //  c04 = 4 -> [0][4] -> 8c d5 dc ff, c12 = 6 -> [0][6] -> 8c d0 d3 ff
  //  c05 = 4 -> [0][4] -> 8c d5 dc ff, c13 = 7 -> [0][7] -> 8c ce ce ff
  //  c06 = 2 -> [1][2] -> 84 c6 cf ff, c14 = 0 -> [1][0] -> 84 c6 c6 ff
  //  c07 = 0 -> [2][0] -> 84 c6 c6 ff, c15 = 1 -> [2][1] -> 8a c7 c9 ff

  // m3:
  //  c00 = 0 -> [0][0] -> 9c b5 94 ff, c08 = 1 -> [0][1] -> a4 b8 91 ff
  //  c01 = 1 -> [0][1] -> a4 b8 91 ff, c09 = 3 -> [0][3] -> b5 bd 8c ff
  //  c02 = 3 -> [2][3] -> ad b5 8c ff, c10 = 3 -> [2][3] -> ad b5 8c ff
  //  c03 = 1 -> [2][1] -> b2 c0 9d ff, c11 = 2 -> [2][2] -> b0 bb 94 ff
  //  c04 = 0 -> [0][0] -> 9c b5 94 ff, c12 = 2 -> [1][2] -> aa b2 94 ff
  //  c05 = 1 -> [0][1] -> a4 b8 91 ff, c13 = 2 -> [1][2] -> aa b2 94 ff
  //  c06 = 2 -> [2][2] -> b0 bb 94 ff, c14 = 0 -> [1][0] -> 94 ad 94 ff
  //  c07 = 2 -> [2][2] -> b0 bb 94 ff, c15 = 1 -> [1][1] -> 9f b0 94 ff

  // m8:
  //  c00 = 1 -> [0][1] -> c3 3c 10 ee, c08 = 3 -> [0][3] -> 51 00 00 d3
  //  c01 = 0 -> [0][0] -> fb 59 18 fb, c09 = 2 -> [1][2] -> 7b 28 0d aa
  //  c02 = 0 -> [0][0] -> fb 59 18 fb, c10 = 1 -> [1][1] -> b5 4d 17 cf
  //  c03 = 1 -> [0][1] -> c3 3c 10 ee, c11 = 0 -> [1][0] -> eb 71 20 f3
  //  c04 = 2 -> [0][2] -> 89 1d 08 e0, c12 = 3 -> [1][3] -> 45 04 04 86
  //  c05 = 1 -> [0][1] -> c3 3c 10 ee, c13 = 3 -> [1][3] -> 45 04 04 86
  //  c06 = 0 -> [0][0] -> fb 59 18 fb, c14 = 3 -> [1][3] -> 45 04 04 86
  //  c07 = 0 -> [1][0] -> eb 71 20 f3, c15 = 1 -> [1][1] -> b5 4d 17 cf

  c0 = ExtrBits<ibits, ibits *  0>(iblk);
  c1 = ExtrBits<ibits, ibits *  1>(iblk);
  c2 = ExtrBits<ibits, ibits *  2>(iblk);
  c3 = ExtrBits<ibits, ibits *  3>(iblk);
  c4 = ExtrBits<ibits, ibits *  4>(iblk);
  c5 = ExtrBits<ibits, ibits *  5>(iblk);
  c6 = ExtrBits<ibits, ibits *  6>(iblk);
  c7 = ExtrBits<ibits, ibits *  7>(iblk);

  out[ 0] = codes[(sets >= 2 ? ios[ 0] << ibits : 0) + c0.GetLong()];
  out[ 1] = codes[(sets >= 2 ? ios[ 1] << ibits : 0) + c1.GetLong()];
  out[ 2] = codes[(sets >= 2 ? ios[ 2] << ibits : 0) + c2.GetLong()];
  out[ 3] = codes[(sets >= 2 ? ios[ 3] << ibits : 0) + c3.GetLong()];
  out[ 4] = codes[(sets >= 2 ? ios[ 4] << ibits : 0) + c4.GetLong()];
  out[ 5] = codes[(sets >= 2 ? ios[ 5] << ibits : 0) + c5.GetLong()];
  out[ 6] = codes[(sets >= 2 ? ios[ 6] << ibits : 0) + c6.GetLong()];
  out[ 7] = codes[(sets >= 2 ? ios[ 7] << ibits : 0) + c7.GetLong()];

  c0 = ExtrBits<ibits, ibits *  8>(iblk);
  c1 = ExtrBits<ibits, ibits *  9>(iblk);
  c2 = ExtrBits<ibits, ibits * 10>(iblk);
  c3 = ExtrBits<ibits, ibits * 11>(iblk);
  c4 = ExtrBits<ibits, ibits * 12>(iblk);
  c5 = ExtrBits<ibits, ibits * 13>(iblk);
  c6 = ExtrBits<ibits, ibits * 14>(iblk);
  c7 = ExtrBits<ibits, ibits * 15>(iblk);

  out[ 8] = codes[(sets >= 2 ? ios[ 8] << ibits : 0) + c0.GetLong()];
  out[ 9] = codes[(sets >= 2 ? ios[ 9] << ibits : 0) + c1.GetLong()];
  out[10] = codes[(sets >= 2 ? ios[10] << ibits : 0) + c2.GetLong()];
  out[11] = codes[(sets >= 2 ? ios[11] << ibits : 0) + c3.GetLong()];
  out[12] = codes[(sets >= 2 ? ios[12] << ibits : 0) + c4.GetLong()];
  out[13] = codes[(sets >= 2 ? ios[13] << ibits : 0) + c5.GetLong()];
  out[14] = codes[(sets >= 2 ? ios[14] << ibits : 0) + c6.GetLong()];
  out[15] = codes[(sets >= 2 ? ios[15] << ibits : 0) + c7.GetLong()];
}

/* -----------------------------------------------------------------------------
 */
void ReadHDRBlock_m1(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 2 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   *   5 bits of shared bitfield          4 bits of shared exponent          5 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                      1 bit  deltacoded exponent
   *   5 bits of deltacoded bitfield      4 bits of deltacoded mantissa      5 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   6 bits of lost bitfield            6 bits of lost mantissa            6 bits of lost mantissa
   *

    {   // 0x00 - 10 5 5 5
        {  M, 0}, {  M, 1}, {GS2, 4}, {BS2, 4}, {BE2, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits<10, 25 -  0>(blkl, a[0]);		// 1-10 set 1 blue start
  ExtrBits< 5, 55 -  0>(blkl, b[0]);		// 1-5 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 50 -  0>(blkl, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 60 -  0>(blkl, t[1]);		// 2-2 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 1;

    ExtrBits< 1,  3 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 70 - 64>(blkh, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 2;

    ExtrBits< 1, 76 - 64>(blkh, t[0]);		// 4-4 set 2 blue end
    ExtrBits< 1,  4 -  0>(blkl, t[1]);		// 5-5 set 2 blue end

    b[1] |= t[0] << 3;
    b[1] |= t[1] << 4;

  ConcBits<10, 15 -  0>(blkl, a[0]);		// 1-10 set 1 green start
  ConcBits< 5, 45 -  0>(blkl, b[0]);		// 1-5 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1,  2 -  0>(blkl, t[0]);		// 5-5 set 2 green start
    ExtrBits< 1, 40 -  0>(blkl, t[1]);		// 5-5 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

  ConcBits<10,  5 -  0>(blkl, a[0]);		// 1-10 set 1 red start
  ConcBits< 5, 35 -  0>(blkl, b[0]);		// 1-5 set 1 red end
  ConcBits< 5, 65 - 64>(blkh, a[1]);		// 1-5 set 2 red start
  ConcBits< 5, 71 - 64>(blkh, b[1]);		// 1-5 set 2 red end

  // zero-fill bottom 6 bits
  a[0] = (a[0] << (16 - 10));
  a[1] = (a[1] << (16 - 10));
  b[0] = (b[0] << (16 - 10));
  b[1] = (b[1] << (16 - 10));

  // undo deltas
  s[0] = a[0];
  s[1] = a[0] + FillSign<5 + 16>(a[1]);
  e[0] = a[0] + FillSign<5 + 16>(b[0]);
  e[1] = a[0] + FillSign<5 + 16>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 10);
	s[1] += Col4(0x8000 >> 10);
	e[0] += Col4(0x8000 >> 10);
	e[1] += Col4(0x8000 >> 10);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m2(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 2 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   *   1 bit  of shared bitfield                                             1 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                      5 bit  deltacoded exponent         4 bits deltacoded exponent
   *   6 bits of deltacoded bitfield      1 bits of deltacoded mantissa      2 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   9 bits of lost bitfield            9 bits of lost mantissa            9 bits of lost mantissa
   *

    {   // 0x01 - 7 6 6 6
        {  M, 0}, {  M, 1}, {GS2, 5}, {GE2, 4}, {GE2, 5}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {BE2, 0}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {BS2, 5}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BE2, 3}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits< 7, 25 -  0>(blkl, a[0]);		// 1-7 set 1 blue start
  ExtrBits< 6, 55 -  0>(blkl, b[0]);		// 1-6 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 2, 12 -  0>(blkl, b[1]);		// 1-2 set 2 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 23 -  0>(blkl, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 2;

    ExtrBits< 1, 14 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 32 -  0>(blkl, t[1]);		// 4-4 set 2 blue end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 3;

    ExtrBits< 1, 22 -  0>(blkl, t[0]);		// 6-6 set 2 blue start
    ExtrBits< 1, 34 -  0>(blkl, t[1]);		// 5-5 set 2 blue end

    a[1] |= t[0] << 5;
    b[1] |= t[1] << 4;

    ExtrBits< 1, 33 -  0>(blkl, t[0]);		// 6-6 set 2 blue end

    b[1] |= t[0] << 5;

  ConcBits< 7, 15 -  0>(blkl, a[0]);		// 1-7 set 1 green start
  ConcBits< 6, 45 -  0>(blkl, b[0]);		// 1-6 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 24 -  0>(blkl, t[0]);		// 5-5 set 2 green start
    ExtrBits< 2,  3 -  0>(blkl, t[1]);		// 5-6 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

    ExtrBits< 1,  2 -  0>(blkl, t[0]);		// 6-6 set 2 green start

    a[1] |= t[0] << 5;

  ConcBits< 7,  5 -  0>(blkl, a[0]);		// 1-7 set 1 red start
  ConcBits< 6, 35 -  0>(blkl, b[0]);		// 1-6 set 1 red end
  ConcBits< 6, 65 - 64>(blkh, a[1]);		// 1-6 set 2 red start
  ConcBits< 6, 71 - 64>(blkh, b[1]);		// 1-6 set 2 red end

  // zero-fill bottom 9 bits
  a[0] = (a[0] << (16 - 7));
  a[1] = (a[1] << (16 - 7));
  b[0] = (b[0] << (16 - 7));
  b[1] = (b[1] << (16 - 7));

  // undo deltas
  s[0] = a[0];
  s[1] = a[0] + FillSign<1 + 16>(a[1]);
  e[0] = a[0] + FillSign<1 + 16>(b[0]);
  e[1] = a[0] + FillSign<1 + 16>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 7);
	s[1] += Col4(0x8000 >> 7);
	e[0] += Col4(0x8000 >> 7);
	e[1] += Col4(0x8000 >> 7);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m3(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   * 6/7 bits of shared bitfield          5 bits of shared exponent          5 bits of shared exponent
   *                                    0/1 bits of shared mantissa        1/2 bits of shared mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   * 5/4 bits of deltacoded bitfield    5/4 bits of deltacoded mantissa    5/4 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   5 bits of lost bitfield            5 bits of lost mantissa            5 bits of lost mantissa
   *

    {   // 0x02 - 11 5 4 4
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RS1,10}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,10},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,10},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
    */

  ExtrBits<10, 25 -  0>(blkl, a[0]);		// 1-10 set 1 blue start
  ExtrBits< 4, 55 -  0>(blkl, b[0]);		// 1-4 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 50 -  0>(blkl, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 59 -  0>(blkl, t[0]);		// 11-11 set 1 blue start
    ExtrBits< 1, 60 -  0>(blkl, t[1]);		// 2-2 set 2 blue end

    a[0] |= t[0] << 10;
    b[1] |= t[1] << 1;

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 70 - 64>(blkh, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 2;

    ExtrBits< 1, 76 - 64>(blkh, t[0]);		// 4-4 set 2 blue end

    b[1] |= t[0] << 3;

  ConcBits<10, 15 -  0>(blkl, a[0]);		// 1-10 set 1 green start
  ConcBits< 4, 45 -  0>(blkl, b[0]);		// 1-4 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 49 -  0>(blkl, t[0]);		// 11-11 set 1 green start

    a[0] |= t[0] << 10;

  ConcBits<10,  5 -  0>(blkl, a[0]);		// 1-10 set 1 red start
  ConcBits< 5, 35 -  0>(blkl, b[0]);		// 1-5 set 1 red end
  ConcBits< 5, 65 - 64>(blkh, a[1]);		// 1-5 set 2 red start
  ConcBits< 5, 71 - 64>(blkh, b[1]);		// 1-5 set 2 red end

    ExtrBits< 1, 40 -  0>(blkl, t[0]);		// 11-11 set 1 red start

    a[0] |= t[0] << 10;

  // zero-fill bottom 5 bits
  a[0] = (a[0] << (16 - 11));
  a[1] = (a[1] << (16 - 11));
  b[0] = (b[0] << (16 - 11));
  b[1] = (b[1] << (16 - 11));

  // variable shift up by [6,7,7,0]
  a[1] = (a[1] * Col4((1 << 16) + 64, (1 << 16) + 128, (1 << 16) + 128, 0));
  b[0] = (b[0] * Col4((1 << 16) + 64, (1 << 16) + 128, (1 << 16) + 128, 0));
  b[1] = (b[1] * Col4((1 << 16) + 64, (1 << 16) + 128, (1 << 16) + 128, 0));

  // sign extend
  a[1] = ExtendSign<6 + 16>(a[1] << 16);
  b[0] = ExtendSign<6 + 16>(b[0] << 16);
  b[1] = ExtendSign<6 + 16>(b[1] << 16);

  // variable shift down by [0,1,1,0]
  a[1] = (a[1] * Col4((1 << 16) + 2, (1 << 16) + 1, (1 << 16) + 1, 0));
  b[0] = (b[0] * Col4((1 << 16) + 2, (1 << 16) + 1, (1 << 16) + 1, 0));
  b[1] = (b[1] * Col4((1 << 16) + 2, (1 << 16) + 1, (1 << 16) + 1, 0));

  // undo deltas
  s[0] = a[0];
  s[1] = a[0] + ExtendSign<1>(a[1]);
  e[0] = a[0] + ExtendSign<1>(b[0]);
  e[1] = a[0] + ExtendSign<1>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 11);
	s[1] += Col4(0x8000 >> 11);
	e[0] += Col4(0x8000 >> 11);
	e[1] += Col4(0x8000 >> 11);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m4(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   * 6/7 bits of shared bitfield          5 bits of shared exponent          5 bits of shared exponent
   *                                    0/1 bits of shared mantissa        1/2 bits of shared mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   * 5/4 bits of deltacoded bitfield    5/4 bits of deltacoded mantissa    5/4 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   5 bits of lost bitfield            5 bits of lost mantissa            5 bits of lost mantissa
   *

    {   // 0x06 - 11 4 5 4
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,10},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GS1,10}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,10},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {BE2, 0},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {GS2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits<10, 25 -  0>(blkl, a[0]);		// 1-10 set 1 blue start
  ExtrBits< 4, 55 -  0>(blkl, b[0]);		// 1-4 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 69 - 64>(blkh, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 59 -  0>(blkl, t[0]);		// 11-11 set 1 blue start
    ExtrBits< 1, 60 -  0>(blkl, t[1]);		// 2-2 set 2 blue end

    a[0] |= t[0] << 10;
    b[1] |= t[1] << 1;

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 70 - 64>(blkh, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 2;

    ExtrBits< 1, 76 - 64>(blkh, t[0]);		// 4-4 set 2 blue end

    b[1] |= t[0] << 3;

  ConcBits<10, 15 -  0>(blkl, a[0]);		// 1-10 set 1 green start
  ConcBits< 5, 45 -  0>(blkl, b[0]);		// 1-5 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 50 -  0>(blkl, t[1]);		// 11-11 set 1 green start

    a[0] |= t[1] << 10;

    ExtrBits< 1, 75 - 64>(blkh, t[0]);		// 5-5 set 2 green start
    ExtrBits< 1, 40 -  0>(blkl, t[1]);		// 5-5 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

  ConcBits<10,  5 -  0>(blkl, a[0]);		// 1-10 set 1 red start
  ConcBits< 4, 35 -  0>(blkl, b[0]);		// 1-4 set 1 red end
  ConcBits< 4, 65 - 64>(blkh, a[1]);		// 1-4 set 2 red start
  ConcBits< 4, 71 - 64>(blkh, b[1]);		// 1-4 set 2 red end

    ExtrBits< 1, 39 -  0>(blkl, t[0]);		// 11-11 set 1 red start

    a[0] |= t[0] << 10;

  // zero-fill bottom 5 bits
  a[0] = (a[0] << (16 - 11));
  a[1] = (a[1] << (16 - 11));
  b[0] = (b[0] << (16 - 11));
  b[1] = (b[1] << (16 - 11));

  // variable shift up by [7,6,7,0]
  a[1] = (a[1] * Col4((1 << 16) + 128, (1 << 16) + 64, (1 << 16) + 128, 0));
  b[0] = (b[0] * Col4((1 << 16) + 128, (1 << 16) + 64, (1 << 16) + 128, 0));
  b[1] = (b[1] * Col4((1 << 16) + 128, (1 << 16) + 64, (1 << 16) + 128, 0));

  // sign extend
  a[1] = ExtendSign<6 + 16>(a[1] << 16);
  b[0] = ExtendSign<6 + 16>(b[0] << 16);
  b[1] = ExtendSign<6 + 16>(b[1] << 16);

  // variable shift down by [1,0,1,0]
  a[1] = (a[1] * Col4((1 << 16) + 1, (1 << 16) + 2, (1 << 16) + 1, 0));
  b[0] = (b[0] * Col4((1 << 16) + 1, (1 << 16) + 2, (1 << 16) + 1, 0));
  b[1] = (b[1] * Col4((1 << 16) + 1, (1 << 16) + 2, (1 << 16) + 1, 0));

  // undo deltas
  s[0] = a[0];
	s[1] = a[0] + ExtendSign<1>(a[1]);
	e[0] = a[0] + ExtendSign<1>(b[0]);
	e[1] = a[0] + ExtendSign<1>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 11);
	s[1] += Col4(0x8000 >> 11);
	e[0] += Col4(0x8000 >> 11);
	e[1] += Col4(0x8000 >> 11);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m5(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   * 6/7 bits of shared bitfield          5 bits of shared exponent          5 bits of shared exponent
   *                                    0/1 bits of shared mantissa        1/2 bits of shared mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   * 5/4 bits of deltacoded bitfield    5/4 bits of deltacoded mantissa    5/4 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   5 bits of lost bitfield            5 bits of lost mantissa            5 bits of lost mantissa
   *

    {   // 0x0a - 11 4 4 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,10},
        {BS2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,10},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BS1,10}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {BE2, 1},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {BE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits<10, 25 -  0>(blkl, a[0]);		// 1-10 set 1 blue start
  ExtrBits< 5, 55 -  0>(blkl, b[0]);		// 1-5 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 50 -  0>(blkl, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 60 -  0>(blkl, t[0]);		// 11-11 set 1 blue start
    ExtrBits< 2, 69 - 64>(blkh, t[1]);		// 2-3 set 2 blue end

    a[0] |= t[0] << 10;
    b[1] |= t[1] << 1;

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 76 - 64>(blkh, t[1]);		// 4-4 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 3;

    ExtrBits< 1, 40 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 75 - 64>(blkh, t[1]);		// 5-5 set 2 blue end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

  ConcBits<10, 15 -  0>(blkl, a[0]);		// 1-10 set 1 green start
  ConcBits< 4, 45 -  0>(blkl, b[0]);		// 1-4 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 49 -  0>(blkl, t[1]);		// 11-11 set 1 green start

    a[0] |= t[1] << 10;

  ConcBits<10,  5 -  0>(blkl, a[0]);		// 1-10 set 1 red start
  ConcBits< 4, 35 -  0>(blkl, b[0]);		// 1-4 set 1 red end
  ConcBits< 4, 65 - 64>(blkh, a[1]);		// 1-4 set 2 red start
  ConcBits< 4, 71 - 64>(blkh, b[1]);		// 1-4 set 2 red end

    ExtrBits< 1, 39 -  0>(blkl, t[0]);		// 11-11 set 1 red start

    a[0] |= t[0] << 10;

  // zero-fill bottom 5 bits
  a[0] = (a[0] << (16 - 11));
  a[1] = (a[1] << (16 - 11));
  b[0] = (b[0] << (16 - 11));
  b[1] = (b[1] << (16 - 11));

  // variable shift up by [7,7,6,0]
  a[1] = (a[1] * Col4((1 << 16) + 128, (1 << 16) + 128, (1 << 16) + 64, 0));
  b[0] = (b[0] * Col4((1 << 16) + 128, (1 << 16) + 128, (1 << 16) + 64, 0));
  b[1] = (b[1] * Col4((1 << 16) + 128, (1 << 16) + 128, (1 << 16) + 64, 0));

	// sign extend
	a[1] = ExtendSign<6 + 16>(a[1] << 16);
	b[0] = ExtendSign<6 + 16>(b[0] << 16);
	b[1] = ExtendSign<6 + 16>(b[1] << 16);

  // variable shift down by [1,1,0,0]
  a[1] = (a[1] * Col4((1 << 16) + 1, (1 << 16) + 1, (1 << 16) + 2, 0));
  b[0] = (b[0] * Col4((1 << 16) + 1, (1 << 16) + 1, (1 << 16) + 2, 0));
  b[1] = (b[1] * Col4((1 << 16) + 1, (1 << 16) + 1, (1 << 16) + 2, 0));

  // undo deltas
	s[0] = a[0];
	s[1] = a[0] + ExtendSign<1>(a[1]);
	e[0] = a[0] + ExtendSign<1>(b[0]);
	e[1] = a[0] + ExtendSign<1>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 11);
	s[1] += Col4(0x8000 >> 11);
	e[0] += Col4(0x8000 >> 11);
	e[1] += Col4(0x8000 >> 11);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m6(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   *   4 bits of shared bitfield          3 bits of shared exponent          4 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                      2 bits of deltacoded exponent      1 bit  of deltacoded exponent
   *   5 bits of deltacoded bitfield      3 bits of deltacoded mantissa      4 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   7 bits of lost bitfield            7 bits of lost mantissa            7 bits of lost mantissa
   *

    {   // 0x0e - 9 5 5 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits< 9, 25 -  0>(blkl, a[0]);		// 1-9 set 1 blue start
  ExtrBits< 5, 55 -  0>(blkl, b[0]);		// 1-5 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 50 -  0>(blkl, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 60 -  0>(blkl, t[1]);		// 2-2 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 1;

    ExtrBits< 1, 14 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 70 - 64>(blkh, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 2;

    ExtrBits< 1, 76 - 64>(blkh, t[0]);		// 4-4 set 2 blue end
    ExtrBits< 1, 34 -  0>(blkl, t[1]);		// 5-5 set 2 blue end

    b[1] |= t[0] << 3;
    b[1] |= t[1] << 4;

  ConcBits< 9, 15 -  0>(blkl, a[0]);		// 1-9 set 1 green start
  ConcBits< 6, 45 -  0>(blkl, b[0]);		// 1-5 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 24 -  0>(blkl, t[0]);		// 5-5 set 2 green start
    ExtrBits< 2, 40 -  0>(blkl, t[1]);		// 5-5 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

  ConcBits< 9,  5 -  0>(blkl, a[0]);		// 1-9 set 1 red start
  ConcBits< 5, 35 -  0>(blkl, b[0]);		// 1-5 set 1 red end
  ConcBits< 5, 65 - 64>(blkh, a[1]);		// 1-5 set 2 red start
  ConcBits< 5, 71 - 64>(blkh, b[1]);		// 1-5 set 2 red end

  // zero-fill bottom 7 bits
  a[0] = (a[0] << (16 - 9));
  a[1] = (a[1] << (16 - 9));
  b[0] = (b[0] << (16 - 9));
  b[1] = (b[1] << (16 - 9));

  // undo deltas
  s[0] = a[0];
  s[1] = a[0] + FillSign<4 + 16>(a[1]);
  e[0] = a[0] + FillSign<4 + 16>(b[0]);
  e[1] = a[0] + FillSign<4 + 16>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 9);
	s[1] += Col4(0x8000 >> 9);
	e[0] += Col4(0x8000 >> 9);
	e[1] += Col4(0x8000 >> 9);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m7(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   * 2/3 bits of shared bitfield        1/2 bits of shared exponent        2/3 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                    4/3 bits of deltacoded exponent    3/2 bit  of deltacoded exponent
   * 6/5 bits of deltacoded bitfield    2/2 bits of deltacoded mantissa    3/3 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   8 bits of lost bitfield            8 bits of lost mantissa            8 bits of lost mantissa
   *

    {   // 0x12 - 8 6 5 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {GE2, 4}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BE2, 3}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits< 8, 25 -  0>(blkl, a[0]);		// 1-8 set 1 blue start
  ExtrBits< 5, 55 -  0>(blkl, b[0]);		// 1-5 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 50 -  0>(blkl, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 60 -  0>(blkl, t[1]);		// 2-2 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 1;

    ExtrBits< 1, 14 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 23 -  0>(blkh, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 2;

    ExtrBits< 2, 33 -  0>(blkl, t[0]);		// 4-5 set 2 blue end

    b[1] |= t[0] << 3;

  ConcBits< 8, 15 -  0>(blkl, a[0]);		// 1-8 set 1 green start
  ConcBits< 5, 45 -  0>(blkl, b[0]);		// 1-5 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 24 -  0>(blkl, t[0]);		// 5-5 set 2 green start
    ExtrBits< 1, 13 -  0>(blkl, t[1]);		// 5-5 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

  ConcBits< 8,  5 -  0>(blkl, a[0]);		// 1-8 set 1 red start
  ConcBits< 6, 35 -  0>(blkl, b[0]);		// 1-6 set 1 red end
  ConcBits< 6, 65 - 64>(blkh, a[1]);		// 1-6 set 2 red start
  ConcBits< 6, 71 - 64>(blkh, b[1]);		// 1-6 set 2 red end

  // zero-fill bottom 8 bits
  a[0] = (a[0] << (16 - 8));
  a[1] = (a[1] << (16 - 8));
  b[0] = (b[0] << (16 - 8));
  b[1] = (b[1] << (16 - 8));

  // variable shift up by [2,3,3,0]
  a[1] = (a[1] * Col4((1 << 16) + 4, (1 << 16) + 8, (1 << 16) + 8, 0));
  b[0] = (b[0] * Col4((1 << 16) + 4, (1 << 16) + 8, (1 << 16) + 8, 0));
  b[1] = (b[1] * Col4((1 << 16) + 4, (1 << 16) + 8, (1 << 16) + 8, 0));

  // sign extend
  a[1] = ExtendSign<2 + 16>(a[1] << 16);
  b[0] = ExtendSign<2 + 16>(b[0] << 16);
  b[1] = ExtendSign<2 + 16>(b[1] << 16);

  // variable shift down by [0,1,1,0]
  a[1] = (a[1] * Col4((1 << 16) + 2, (1 << 16) + 1, (1 << 16) + 1, 0));
  b[0] = (b[0] * Col4((1 << 16) + 2, (1 << 16) + 1, (1 << 16) + 1, 0));
  b[1] = (b[1] * Col4((1 << 16) + 2, (1 << 16) + 1, (1 << 16) + 1, 0));

  // undo deltas
  s[0] = a[0];
	s[1] = a[0] + ExtendSign<1>(a[1]);
	e[0] = a[0] + ExtendSign<1>(b[0]);
	e[1] = a[0] + ExtendSign<1>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 8);
	s[1] += Col4(0x8000 >> 8);
	e[0] += Col4(0x8000 >> 8);
	e[1] += Col4(0x8000 >> 8);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m8(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   * 2/3 bits of shared bitfield        1/2 bits of shared exponent        2/3 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                    4/3 bits of deltacoded exponent    3/2 bit  of deltacoded exponent
   * 6/5 bits of deltacoded bitfield    2/2 bits of deltacoded mantissa    3/3 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   8 bits of lost bitfield            8 bits of lost mantissa            8 bits of lost mantissa
   *

    {   // 0x16 - 8 5 6 5
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {BE2, 0}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS2, 5}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {GE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE2, 1}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits< 8, 25 -  0>(blkl, a[0]);		// 1-8 set 1 blue start
  ExtrBits< 5, 55 -  0>(blkl, b[0]);		// 1-5 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 13 -  0>(blkl, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 60 -  0>(blkl, t[1]);		// 2-2 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 1;

    ExtrBits< 1, 14 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 70 - 64>(blkh, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 2;

    ExtrBits< 1, 76 - 64>(blkh, t[0]);		// 4-4 set 2 blue end
    ExtrBits< 1, 34 -  0>(blkl, t[1]);		// 5-5 set 2 blue end

    b[1] |= t[0] << 3;
    b[1] |= t[1] << 4;

  ConcBits< 8, 15 -  0>(blkl, a[0]);		// 1-8 set 1 green start
  ConcBits< 6, 45 -  0>(blkl, b[0]);		// 1-6 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 24 -  0>(blkl, t[0]);		// 5-5 set 2 green start
    ExtrBits< 1, 40 -  0>(blkl, t[1]);		// 5-5 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

    ExtrBits< 1, 23 -  0>(blkl, t[0]);		// 6-6 set 2 green start
    ExtrBits< 1, 33 -  0>(blkl, t[1]);		// 6-6 set 2 green end

    a[1] |= t[0] << 5;
    b[1] |= t[1] << 5;

  ConcBits< 8,  5 -  0>(blkl, a[0]);		// 1-8 set 1 red start
  ConcBits< 5, 35 -  0>(blkl, b[0]);		// 1-5 set 1 red end
  ConcBits< 5, 65 - 64>(blkh, a[1]);		// 1-5 set 2 red start
  ConcBits< 5, 71 - 64>(blkh, b[1]);		// 1-5 set 2 red end

  // zero-fill bottom 8 bits
  a[0] = (a[0] << (16 - 8));
  a[1] = (a[1] << (16 - 8));
  b[0] = (b[0] << (16 - 8));
  b[1] = (b[1] << (16 - 8));

  // variable shift up by [3,2,3,0]
  a[1] = (a[1] * Col4((1 << 16) + 8, (1 << 16) + 4, (1 << 16) + 8, 0));
  b[0] = (b[0] * Col4((1 << 16) + 8, (1 << 16) + 4, (1 << 16) + 8, 0));
  b[1] = (b[1] * Col4((1 << 16) + 8, (1 << 16) + 4, (1 << 16) + 8, 0));

	// sign extend
	a[1] = ExtendSign<2 + 16>(a[1] << 16);
	b[0] = ExtendSign<2 + 16>(b[0] << 16);
	b[1] = ExtendSign<2 + 16>(b[1] << 16);

  // variable shift down by [1,0,1,0]
  a[1] = (a[1] * Col4((1 << 16) + 1, (1 << 16) + 2, (1 << 16) + 1, 0));
  b[0] = (b[0] * Col4((1 << 16) + 1, (1 << 16) + 2, (1 << 16) + 1, 0));
  b[1] = (b[1] * Col4((1 << 16) + 1, (1 << 16) + 2, (1 << 16) + 1, 0));

  // undo deltas
	s[0] = a[0];
	s[1] = a[0] + ExtendSign<1>(a[1]);
	e[0] = a[0] + ExtendSign<1>(b[0]);
	e[1] = a[0] + ExtendSign<1>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 8);
	s[1] += Col4(0x8000 >> 8);
	e[0] += Col4(0x8000 >> 8);
	e[1] += Col4(0x8000 >> 8);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_m9(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   * 2/3 bits of shared bitfield        1/2 bits of shared exponent        2/3 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                    4/3 bits of deltacoded exponent    3/2 bit  of deltacoded exponent
   * 6/5 bits of deltacoded bitfield    2/2 bits of deltacoded mantissa    3/3 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   8 bits of lost bitfield            8 bits of lost mantissa            8 bits of lost mantissa
   *

    {   // 0x1a - 8 5 5 6
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {BS2, 5}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {GE2, 4}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {BE2, 0}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {BE2, 2}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {BE2, 3}, {  D, 0}, {  D, 1}, {  D, 2},
        {  D, 3}, {  D, 4},
    },
   */

  ExtrBits< 8, 25 -  0>(blkl, a[0]);		// 1-8 set 1 blue start
  ExtrBits< 6, 55 -  0>(blkl, b[0]);		// 1-6 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 1, 50 -  0>(blkl, b[1]);		// 1-1 set 2 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 13 -  0>(blkl, t[1]);		// 2-2 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 1;

    ExtrBits< 1, 70 - 64>(blkh, t[0]);		// 3-3 set 2 blue end
    ExtrBits< 1, 76 - 64>(blkh, t[1]);		// 4-4 set 2 blue end

    b[1] |= t[0] << 2;
    b[1] |= t[1] << 3;

    ExtrBits< 1, 14 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 34 -  0>(blkl, t[1]);		// 5-5 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 4;

    ExtrBits< 1, 23 -  0>(blkl, t[0]);		// 6-6 set 2 blue start
    ExtrBits< 1, 33 -  0>(blkl, t[1]);		// 6-6 set 2 blue end

    a[1] |= t[0] << 5;
    b[1] |= t[1] << 5;

  ConcBits< 8, 15 -  0>(blkl, a[0]);		// 1-8 set 1 green start
  ConcBits< 5, 45 -  0>(blkl, b[0]);		// 1-5 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 24 -  0>(blkl, t[0]);		// 5-5 set 2 green start
    ExtrBits< 2, 40 -  0>(blkl, t[1]);		// 5-5 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

  ConcBits< 8,  5 -  0>(blkl, a[0]);		// 1-8 set 1 red start
  ConcBits< 5, 35 -  0>(blkl, b[0]);		// 1-5 set 1 red end
  ConcBits< 5, 65 - 64>(blkh, a[1]);		// 1-5 set 2 red start
  ConcBits< 5, 71 - 64>(blkh, b[1]);		// 1-5 set 2 red end

  // zero-fill bottom 8 bits
  a[0] = (a[0] << (16 - 8));
  a[1] = (a[1] << (16 - 8));
  b[0] = (b[0] << (16 - 8));
  b[1] = (b[1] << (16 - 8));

  // variable shift up by [3,3,2,0]
  a[1] = (a[1] * Col4((1 << 16) + 8, (1 << 16) + 8, (1 << 16) + 4, 0));
  b[0] = (b[0] * Col4((1 << 16) + 8, (1 << 16) + 8, (1 << 16) + 4, 0));
  b[1] = (b[1] * Col4((1 << 16) + 8, (1 << 16) + 8, (1 << 16) + 4, 0));

	// sign extend
	a[1] = ExtendSign<2 + 16>(a[1] << 16);
	b[0] = ExtendSign<2 + 16>(b[0] << 16);
	b[1] = ExtendSign<2 + 16>(b[1] << 16);

  // variable shift down by [1,1,0,0]
  a[1] = (a[1] * Col4((1 << 16) + 1, (1 << 16) + 1, (1 << 16) + 2, 0));
  b[0] = (b[0] * Col4((1 << 16) + 1, (1 << 16) + 1, (1 << 16) + 2, 0));
  b[1] = (b[1] * Col4((1 << 16) + 1, (1 << 16) + 1, (1 << 16) + 2, 0));

  // undo deltas
	s[0] = a[0];
	s[1] = a[0] + ExtendSign<1>(a[1]);
	e[0] = a[0] + ExtendSign<1>(b[0]);
	e[1] = a[0] + ExtendSign<1>(b[1]);

	// add half a unit
	s[0] += Col4(0x8000 >> 8);
	s[1] += Col4(0x8000 >> 8);
	e[0] += Col4(0x8000 >> 8);
	e[1] += Col4(0x8000 >> 8);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_mA(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[2];
  Col4 e[2];
  Col4 a[2];
  Col4 b[2];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 2 mode bits, 5 partition bits
  partition = ExtrBits< 5, 77 - 64>(blkh).GetLong();
  
  /* ushort:				half-float:
   *
   *                                      1 bit  of explicit sign
   *   6 bits of explicit bitfield        5 bits of explicit exponent        5 bits of explicit exponent
   *                                                                         1 bit  of explicit mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *  10 bits of lost bitfield           10 bits of lost mantissa           10 bits of lost mantissa
   *

    {   // 0x1e - 6 6 6 6
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {GE2, 4}, {BE2, 0}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS2, 5}, {BS2, 5}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {GE2, 5}, {BE2, 3}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5},
    },
   */

  ExtrBits< 6, 25 -  0>(blkl, a[0]);		// 1-6 set 1 blue start
  ExtrBits< 6, 55 -  0>(blkl, b[0]);		// 1-6 set 1 blue end
  ExtrBits< 3, 61 -  0>(blkl, a[1]);		// 1-3 set 2 blue start
  ExtrBits< 2, 12 -  0>(blkl, b[1]);		// 1-2 set 2 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);		// 4-4 set 2 blue start
    ExtrBits< 1, 23 -  0>(blkl, t[1]);		// 3-3 set 2 blue end

    a[1] |= t[0] << 3;
    b[1] |= t[1] << 2;

    ExtrBits< 1, 14 -  0>(blkl, t[0]);		// 5-5 set 2 blue start
    ExtrBits< 1, 32 -  0>(blkl, t[1]);		// 4-4 set 2 blue end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 3;

    ExtrBits< 1, 22 -  0>(blkl, t[0]);		// 6-6 set 2 blue start
    ExtrBits< 1, 34 -  0>(blkl, t[1]);		// 5-5 set 2 blue end

    a[1] |= t[0] << 5;
    b[1] |= t[1] << 4;

    ExtrBits< 1, 33 -  0>(blkl, t[0]);		// 6-6 set 2 blue end

    b[1] |= t[0] << 5;

  ConcBits< 6, 15 -  0>(blkl, a[0]);		// 1-6 set 1 green start
  ConcBits< 6, 45 -  0>(blkl, b[0]);		// 1-6 set 1 green end
  ConcBits< 4, 41 -  0>(blkl, a[1]);		// 1-4 set 2 green start
  ConcBits< 4, 51 -  0>(blkl, b[1]);		// 1-4 set 2 green end

    ExtrBits< 1, 24 -  0>(blkl, t[0]);		// 5-5 set 2 green start
    ExtrBits< 1, 11 -  0>(blkl, t[1]);		// 5-5 set 2 green end

    a[1] |= t[0] << 4;
    b[1] |= t[1] << 4;

    ExtrBits< 1, 21 -  0>(blkl, t[0]);		// 6-6 set 2 green start
    ExtrBits< 1, 31 -  0>(blkl, t[1]);		// 6-6 set 2 green end

    a[1] |= t[0] << 5;
    b[1] |= t[1] << 5;

  ConcBits< 6,  5 -  0>(blkl, a[0]);		// 1-6 set 1 red start
  ConcBits< 6, 35 -  0>(blkl, b[0]);		// 1-6 set 1 red end
  ConcBits< 6, 65 - 64>(blkh, a[1]);		// 1-6 set 2 red start
  ConcBits< 6, 71 - 64>(blkh, b[1]);		// 1-6 set 2 red end

  // zero-fill bottom 10 bits
  a[0] = (a[0] << (16 - 6));
  a[1] = (a[1] << (16 - 6));
  b[0] = (b[0] << (16 - 6));
  b[1] = (b[1] << (16 - 6));

  // undo explicits
  s[0] = a[0];
  s[1] = a[1];
  e[0] = b[0];
  e[1] = b[1];

	// add half a unit
	s[0] += Col4(0x8000 >> 6);
	s[1] += Col4(0x8000 >> 6);
	e[0] += Col4(0x8000 >> 6);
	e[1] += Col4(0x8000 >> 6);

  // generate the midpoints
  CodebookP<3>(codes[0], s[0].GetCol3(), e[0].GetCol3());
  CodebookP<3>(codes[1], s[1].GetCol3(), e[1].GetCol3());

  // 128 - 82 -> 46 index bits + 2 bits from 2 set start/end order -> 16 * 3bit
  ReadHDRBlock<2, 3, 82>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_mB(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[1];
  Col4 e[1];
  Col4 a[1];
  Col4 b[1];
  Col4 t[1];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[1][1 << 4];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 0 partition bits
  partition = 0/*ExtrBits< 5, 77 - 64>(blkh).GetLong()*/;
  
  /* ushort:				half-float:
   *
   *                                      1 bit  of explicit sign
   *   6 bits of explicit bitfield        5 bits of explicit exponent        5 bits of explicit exponent
   *                                                                         1 bit  of explicit mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *  10 bits of lost bitfield           10 bits of lost mantissa           10 bits of lost mantissa
   *

    {   // 0x1e - 6 6 6 6
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {GE2, 4}, {BE2, 0}, {BE2, 1}, {BS2, 4}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS2, 5}, {BS2, 5}, {BE2, 2}, {GS2, 4}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {GE2, 5}, {BE2, 3}, {BE2, 5}, {BE2, 4}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {GS2, 0}, {GS2, 1}, {GS2, 2}, {GS2, 3}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE2, 0}, {GE2, 1}, {GE2, 2}, {GE2, 3}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BS2, 0}, {BS2, 1}, {BS2, 2}, {BS2, 3}, {RS2, 0}, {RS2, 1}, {RS2, 2}, {RS2, 3}, {RS2, 4},
        {RS2, 5}, {RE2, 0}, {RE2, 1}, {RE2, 2}, {RE2, 3}, {RE2, 4}, {RE2, 5},
    },
   */

  ExtrBits<10, 25 -  0>(blkl, a[0]);	// 1-10 set 1 blue start
  ConcBits<10, 15 -  0>(blkl, a[0]);	// 1-10 set 1 green start
  ConcBits<10,  5 -  0>(blkl, a[0]);	// 1-10 set 1 red start

  ExtrBits< 9, 55 -  0>(blkl, b[0]);	// 1-9 set 1 blue end

    ExtrBits< 1, 64 - 64>(blkh, t[0]);	// 10-10 set 1 blue end

    b[0] |= t[0] << 9;

  ConcBits<10, 45 -  0>(blkl, b[0]);	// 1-10 set 1 green end
  ConcBits<10, 35 -  0>(blkl, b[0]);	// 1-10 set 1 red end

  // zero-fill bottom 6 bits
  a[0] = (a[0] << (16 - 10));
  b[0] = (b[0] << (16 - 10));

  // undo explicits
  s[0] = a[0];
  e[0] = b[0];

	// add half a unit
	s[0] += Col4(0x8000 >> 10);
	e[0] += Col4(0x8000 >> 10);

  // generate the midpoints
  CodebookP<4>(codes[0], s[0].GetCol3(), e[0].GetCol3());

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  ReadHDRBlock<1, 4, 65>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_mC(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[1];
  Col4 e[1];
  Col4 a[1];
  Col4 b[1];
  Col4 t[1];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[1][1 << 4];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 0 partition bits
  partition = 0/*ExtrBits< 5, 77 - 64>(blkh).GetLong()*/;

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   *   2 bits of shared bitfield          1 bits of shared exponent          2 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                      4 bits of deltacoded exponent      3 bit  of deltacoded exponent
   *   9 bits of deltacoded bitfield      5 bits of deltacoded mantissa      6 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   5 bits of lost bitfield            5 bits of lost mantissa            5 bits of lost mantissa
   *

    {   // 0x07 - 11 9
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {RE1, 6}, {RE1, 7}, {RE1, 8}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE1, 6}, {GE1, 7}, {GE1, 8}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BE1, 6}, {BE1, 7}, {BE1, 8}, {BS1,10},
    },
    */

  ExtrBits<10, 25 -  0>(blkl, a[0]);	// 1-10 set 1 blue start

    ExtrBits< 1, 64 - 64>(blkh, t[0]);	// 11-11 set 1 blue start

    a[0] |= t[0] << 10;

  ConcBits<10, 15 -  0>(blkl, a[0]);	// 1-10 set 1 green start

    ExtrBits< 1, 54 -  0>(blkl, t[0]);	// 11-11 set 1 green start

    a[0] |= t[0] << 10;

  ConcBits<10,  5 -  0>(blkl, a[0]);	// 1-10 set 1 red start

    ExtrBits< 1, 44 -  0>(blkl, t[0]);	// 11-11 set 1 red start

    a[0] |= t[0] << 10;

  ExtrBits< 9, 55 -  0>(blkl, b[0]);	// 1-9 set 1 blue end
  ConcBits< 9, 45 -  0>(blkl, b[0]);	// 1-9 set 1 green end
  ConcBits< 9, 35 -  0>(blkl, b[0]);	// 1-9 set 1 red end

  // zero-fill bottom 5 bits
  a[0] = (a[0] << (16 - 11));
  b[0] = (b[0] << (16 - 11));

  // undo deltas
  s[0] = a[0];
  e[0] = a[0] + FillSign<2 + 16>(b[0]);

	// add half a unit
	s[0] += Col4(0x8000 >> 11);
	e[0] += Col4(0x8000 >> 11);

  // generate the midpoints
  CodebookP<4>(codes[0], s[0].GetCol3(), e[0].GetCol3());

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  ReadHDRBlock<1, 4, 65>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_mD(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[1];
  Col4 e[1];
  Col4 a[1];
  Col4 b[1];
  Col4 t[2];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[1][1 << 4];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 0 partition bits
  partition = 0/*ExtrBits< 5, 77 - 64>(blkh).GetLong()*/;

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   *   4 bits of shared bitfield          3 bits of shared exponent          4 bits of shared exponent
   * ---------------------------------- ---------------------------------- ----------------------------------
   *                                      2 bits of deltacoded exponent      1 bit  of deltacoded exponent
   *   8 bits of deltacoded bitfield      6 bits of deltacoded mantissa      7 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   4 bits of lost bitfield            4 bits of lost mantissa            4 bits of lost mantissa
   *

    {   // 0x0b - 12 8
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RE1, 4},
        {RE1, 5}, {RE1, 6}, {RE1, 7}, {RS1,11}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GE1, 4},
        {GE1, 5}, {GE1, 6}, {GE1, 7}, {GS1,11}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BE1, 4},
        {BE1, 5}, {BE1, 6}, {BE1, 7}, {BS1,11}, {BS1,10},
    },
    */

  ExtrBits<10, 25 -  0>(blkl, a[0]);	// 1-10 set 1 blue start

    ExtrBits< 1, 64 - 64>(blkh, t[0]);	// 11-11 set 1 blue start
    ExtrBits< 1, 63 -  0>(blkl, t[1]);	// 12-12 set 1 blue start

    a[0] |= t[0] << 10;
    a[0] |= t[1] << 11;

  ConcBits<10, 15 -  0>(blkl, a[0]);	// 1-10 set 1 green start

    ExtrBits< 1, 54 -  0>(blkl, t[0]);	// 11-11 set 1 green start
    ExtrBits< 1, 53 -  0>(blkl, t[1]);	// 12-12 set 1 green start

    a[0] |= t[0] << 10;
    a[0] |= t[1] << 11;

  ConcBits<10,  5 -  0>(blkl, a[0]);	// 1-10 set 1 red start

    ExtrBits< 1, 44 -  0>(blkl, t[0]);	// 11-11 set 1 red start
    ExtrBits< 1, 43 -  0>(blkl, t[1]);	// 12-12 set 1 red start

    a[0] |= t[0] << 10;
    a[0] |= t[1] << 11;

  ExtrBits< 8, 55 -  0>(blkl, b[0]);	// 1-8 set 1 blue end
  ConcBits< 8, 45 -  0>(blkl, b[0]);	// 1-8 set 1 green end
  ConcBits< 8, 35 -  0>(blkl, b[0]);	// 1-8 set 1 red end

  // zero-fill bottom 4 bits
  a[0] = (a[0] << (16 - 12));
  b[0] = (b[0] << (16 - 12));

  // undo deltas
  s[0] = a[0];
  e[0] = a[0] + FillSign<4 + 16>(b[0]);

	// add half a unit
	s[0] += Col4(0x8000 >> 12);
	e[0] += Col4(0x8000 >> 12);

  // generate the midpoints
  CodebookP<4>(codes[0], s[0].GetCol3(), e[0].GetCol3());

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  ReadHDRBlock<1, 4, 65>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void ReadHDRBlock_mE(u16* rgb, void const* block) {
  // get the packed values
  Col4 s[1];
  Col4 e[1];
  Col4 a[1];
  Col4 b[1];
  Col4 t[3];
  Col4 blkl, blkh;

  // remap the indices
  unsigned__int64 codes[1][1 << 4];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bits, 0 partition bits
  partition = 0/*ExtrBits< 5, 77 - 64>(blkh).GetLong()*/;

  /* ushort:				half-float:
   *
   *                                      1 bit  of shared sign
   *   8 bits of shared bitfield          5 bits of shared exponent          5 bits of shared exponent
   *                                      2 bits of shared mantissa          3 bits of shared mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   4 bits of deltacoded bitfield      4 bits of deltacoded mantissa      4 bits of deltacoded mantissa
   * ---------------------------------- ---------------------------------- ----------------------------------
   *   4 bits of lost bitfield            4 bits of lost mantissa            4 bits of lost mantissa
   *

    {   // 0x0f - 16 4
        {  M, 0}, {  M, 1}, {  M, 2}, {  M, 3}, {  M, 4}, {RS1, 0}, {RS1, 1}, {RS1, 2}, {RS1, 3}, {RS1, 4},
        {RS1, 5}, {RS1, 6}, {RS1, 7}, {RS1, 8}, {RS1, 9}, {GS1, 0}, {GS1, 1}, {GS1, 2}, {GS1, 3}, {GS1, 4},
        {GS1, 5}, {GS1, 6}, {GS1, 7}, {GS1, 8}, {GS1, 9}, {BS1, 0}, {BS1, 1}, {BS1, 2}, {BS1, 3}, {BS1, 4},
        {BS1, 5}, {BS1, 6}, {BS1, 7}, {BS1, 8}, {BS1, 9}, {RE1, 0}, {RE1, 1}, {RE1, 2}, {RE1, 3}, {RS1,15},
        {RS1,14}, {RS1,13}, {RS1,12}, {RS1,11}, {RS1,10}, {GE1, 0}, {GE1, 1}, {GE1, 2}, {GE1, 3}, {GS1,15},
        {GS1,14}, {GS1,13}, {GS1,12}, {GS1,11}, {GS1,10}, {BE1, 0}, {BE1, 1}, {BE1, 2}, {BE1, 3}, {BS1,15},
        {BS1,14}, {BS1,13}, {BS1,12}, {BS1,11}, {BS1,10}
    },
   */

  ExtrBits<10, 25 -  0>(blkl, a[0]);	// 1-10 set 1 blue start

    ExtrBits< 1, 64 - 64>(blkh, t[0]);	// 11-11 set 1 blue start
    ExtrBits< 1, 63 -  0>(blkl, t[1]);	// 12-12 set 1 blue start
    ExtrBits< 1, 62 -  0>(blkl, t[2]);	// 13-13 set 1 blue start

    a[0] |= t[0] << 10;
    a[0] |= t[1] << 11;
    a[0] |= t[2] << 12;

    ExtrBits< 1, 61 -  0>(blkl, t[0]);	// 14-14 set 1 blue start
    ExtrBits< 1, 60 -  0>(blkl, t[1]);	// 15-15 set 1 blue start
    ExtrBits< 1, 59 -  0>(blkl, t[2]);	// 16-16 set 1 blue start

    a[0] |= t[0] << 13;
    a[0] |= t[1] << 14;
    a[0] |= t[2] << 15;

  ConcBits<10, 15 -  0>(blkl, a[0]);	// 1-10 set 1 green start

    ExtrBits< 1, 54 -  0>(blkl, t[0]);	// 11-11 set 1 green start
    ExtrBits< 1, 53 -  0>(blkl, t[1]);	// 12-12 set 1 green start
    ExtrBits< 1, 52 -  0>(blkl, t[2]);	// 13-13 set 1 green start

    a[0] |= t[0] << 10;
    a[0] |= t[1] << 11;
    a[0] |= t[2] << 12;

    ExtrBits< 1, 51 -  0>(blkl, t[0]);	// 14-14 set 1 green start
    ExtrBits< 1, 50 -  0>(blkl, t[1]);	// 15-15 set 1 green start
    ExtrBits< 1, 49 -  0>(blkl, t[2]);	// 16-16 set 1 green start

    a[0] |= t[0] << 13;
    a[0] |= t[1] << 14;
    a[0] |= t[2] << 15;

  ConcBits<10,  5 -  0>(blkl, a[0]);	// 1-10 set 1 red start

    ExtrBits< 1, 44 -  0>(blkl, t[0]);	// 11-11 set 1 red start
    ExtrBits< 1, 43 -  0>(blkl, t[1]);	// 12-12 set 1 red start
    ExtrBits< 1, 42 -  0>(blkl, t[2]);	// 13-13 set 1 red start

    a[0] |= t[0] << 10;
    a[0] |= t[1] << 11;
    a[0] |= t[2] << 12;

    ExtrBits< 1, 41 -  0>(blkl, t[0]);	// 14-14 set 1 red start
    ExtrBits< 1, 40 -  0>(blkl, t[1]);	// 15-15 set 1 red start
    ExtrBits< 1, 39 -  0>(blkl, t[2]);	// 16-16 set 1 red start

    a[0] |= t[0] << 13;
    a[0] |= t[1] << 14;
    a[0] |= t[2] << 15;

  ExtrBits< 4, 55 -  0>(blkl, b[0]);	// 1-4 set 1 blue end
  ConcBits< 4, 45 -  0>(blkl, b[0]);	// 1-4 set 1 green end
  ConcBits< 4, 35 -  0>(blkl, b[0]);	// 1-4 set 1 red end

  // zero-fill bottom 0 bits
//a[0] = (a[0] << (16 - 16));
//b[0] = (b[0] << (16 - 16));

  // undo deltas
  s[0] = a[0];
	e[0] = a[0] + FillSign<12 + 16>(b[0]);

	// add half a unit
//s[0] += Col4(0x8000 >> 16);
//e[0] += Col4(0x8000 >> 16);

  // generate the midpoints
  CodebookP<4>(codes[0], s[0].GetCol3(), e[0].GetCol3());

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  ReadHDRBlock<1, 4, 65>(partition, (__int64 *)codes, blkl, blkh, (__int64 *)rgb);
}

void DecompressHDRsBtc6u(u16* rgb, void const* block)
{
  // get the block bytes
  u8 const* bytes = reinterpret_cast< u8 const* >(block);

  // get the mode-selector
  unsigned long selector = *((int *)bytes);

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308952%28v=vs.85%29.aspx
   *
   * 1	46 bits	5 bits	75 bits (10.555, 10.555, 10.555)	2 bits (    0-0)
   * 2	46 bits	5 bits	75 bits (7666, 7666, 7666)		2 bits (    0-1)
   * 3	46 bits	5 bits	72 bits (11.555, 11.444, 11.444)	5 bits (000-1-0)
   * 4	46 bits	5 bits	72 bits (11.444, 11.555, 11.444)	5 bits (001-1-0)
   * 5	46 bits	5 bits	72 bits (11.444, 11.444, 11.555)	5 bits (010-1-0)
   * 6	46 bits	5 bits	72 bits (9555, 9555, 9555)		5 bits (011-1-0)
   * 7 	46 bits	5 bits	72 bits (8666, 8555, 8555)		5 bits (100-1-0)
   * 8 	46 bits	5 bits	72 bits (8555, 8666, 8555)		5 bits (101-1-0)
   * 9 	46 bits	5 bits	72 bits (8555, 8555, 8666)		5 bits (110-1-0)
   * 10 46 bits	5 bits	72 bits (6666, 6666, 6666)		5 bits (111-1-0)
   * 11	63 bits	0 bits	60 bits (10.10, 10.10, 10.10)		5 bits (000-1-1)
   * 12	63 bits	0 bits	60 bits (11.9, 11.9, 11.9)		5 bits (001-1-1)
   * 13	63 bits	0 bits	60 bits (12.8, 12.8, 12.8)		5 bits (010-1-1)
   * 14	63 bits	0 bits	60 bits (16.4, 16.4, 16.4)		5 bits (011-1-1)
   */

  switch (selector & 3) {
    case 0: ReadHDRBlock_m1(rgb, block); break;
    case 1: ReadHDRBlock_m2(rgb, block); break;
    case 2:
      selector = (selector >> 2) & 0x7;
      switch (selector) {
        case 0: ReadHDRBlock_m3(rgb, block); break;
        case 1: ReadHDRBlock_m4(rgb, block); break;
        case 2: ReadHDRBlock_m5(rgb, block); break;
        case 3: ReadHDRBlock_m6(rgb, block); break;
        case 4: ReadHDRBlock_m7(rgb, block); break;
        case 5: ReadHDRBlock_m8(rgb, block); break;
        case 6: ReadHDRBlock_m9(rgb, block); break;
        case 7: ReadHDRBlock_mA(rgb, block); break;
      }
      break;
    case 3:
      selector = (selector >> 2) & 0x7;
      switch (selector) {
        case 0: ReadHDRBlock_mB(rgb, block); break;
        case 1: ReadHDRBlock_mC(rgb, block); break;
        case 2: ReadHDRBlock_mD(rgb, block); break;
        case 3: ReadHDRBlock_mE(rgb, block); break;
//      case 4: ReadHDRBlock_mF(rgb, block); break;
//      case 5: ReadHDRBlock_mG(rgb, block); break;
//      case 6: ReadHDRBlock_mH(rgb, block); break;
//      case 7: ReadHDRBlock_mI(rgb, block); break;
      }
      break;
  }
}

void DecompressHDRsBtc6u(f23* rgb, void const* block)
{
	/* TODO: direct conversion to float without going over memory */
  u16 xyz[4 * 16];

  DecompressHDRsBtc6u(xyz, block);

	/* TODO: SIMD-acceleration of UHalfToFloat/SHalfToFloat */
	for (int i = 0; i < 16; i++) {
		// scale the magnitude by 31/64
		rgb[(i * 4) + 0] = SHalfToFloat((xyz[(i * 4) + 0] * 31) >> 6);
		rgb[(i * 4) + 1] = SHalfToFloat((xyz[(i * 4) + 1] * 31) >> 6);
		rgb[(i * 4) + 2] = SHalfToFloat((xyz[(i * 4) + 2] * 31) >> 6);
		rgb[(i * 4) + 3] = 1.0f;
	}
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
