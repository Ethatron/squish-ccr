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

#include "colourblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
static void WriteColourBlock(int a, int b, u8* indices, void* block)
{
  // get the block as bytes
  u8* bytes = ( u8* )block;

  // write the endpoints
  bytes[0] = (u8)(a & 0xff);
  bytes[1] = (u8)(a >> 8);
  bytes[2] = (u8)(b & 0xff);
  bytes[3] = (u8)(b >> 8);

  // write the indices
  for (int i = 0; i < 4; ++i) {
    u8 const* ind = indices + 4 * i;

    bytes[4 + i] = 
      (ind[0] << 0) | 
      (ind[1] << 2) | 
      (ind[2] << 4) | 
      (ind[3] << 6);
  }
}

void WriteColourBlock3(Vec3::Arg start, Vec3::Arg end, u8 const* indices, void* block)
{
  // get the packed values
  int a = FloatTo565(start);
  int b = FloatTo565(end);

  // remap the indices
  u8 remapped[16];

  if (a <= b) {
    // use the indices directly
    for (int i = 0; i < 16; ++i)
      remapped[i] = indices[i];
  }
  else {
    // swap a and b
    std::swap(a, b);
    for (int i = 0; i < 16; ++i) {
      if (indices[i] == 0)
	remapped[i] = 1;
      else if (indices[i] == 1)
	remapped[i] = 0;
      else
	remapped[i] = indices[i];
    }
  }

  // write the block
  WriteColourBlock(a, b, remapped, block);
}

void WriteColourBlock4(Vec3::Arg start, Vec3::Arg end, u8 const* indices, void* block)
{
  // get the packed values
  int a = FloatTo565(start);
  int b = FloatTo565(end);

  // remap the indices
  u8 remapped[16];

  if (a < b) {
    // swap a and b
    std::swap(a, b);
    for (int i = 0; i < 16; ++i)
      remapped[i] = (indices[i] ^ 0x1) & 0x3;
  }
  else if (a == b) {
    // use index 0
    for (int i = 0; i < 16; ++i)
      remapped[i] = 0;
  }
  else {
    // use the indices directly
    for (int i = 0; i < 16; ++i)
      remapped[i] = indices[i];
  }

  // write the block
  WriteColourBlock(a, b, remapped, block);
}

void DecompressColourBtc(u8* rgba, void const* block, bool isBtc1)
{
  // get the block bytes
  u8 const* bytes = reinterpret_cast< u8 const* >(block);

  // unpack the endpoints
  u8 codes[16];

  int a = Unpack565(bytes + 0, codes + 0);
  int b = Unpack565(bytes + 2, codes + 4);

  // generate the midpoints
  if (isBtc1)
    Codebook3(codes, a <= b);
  else
    Codebook4(codes);

  // unpack the indices
  u8 indices[16];
  for (int i = 0; i < 4; ++i) {
    u8* ind = indices + 4 * i;
    u8 packed = bytes[4 + i];

    ind[0] = (packed >> 0) & 0x3;
    ind[1] = (packed >> 2) & 0x3;
    ind[2] = (packed >> 4) & 0x3;
    ind[3] = (packed >> 6) & 0x3;
  }

  // store out the colours
  for (int i = 0; i < 16; ++i) {
    u8 offset = 4 * indices[i];
    for (int j = 0; j < 4; ++j)
      rgba[4 * i + j] = codes[offset + j];
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
static void WriteColourBlock(tile_barrier barrier, const int thread, lineI2 colour, index16 indices, out code64 block) amp_restricted
{
  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  threaded_cse(0) {
    // write the endpoints
    block[0] =
      ((colour[CSTRT] /*& 0xFFFF*/) <<  0) |
      ((colour[CSTOP] /*& 0xFFFF*/) << 16);

    block[1] =
      (((indices[ 0] << 0) | (indices[ 1] << 2) | (indices[ 2] << 4) | (indices[ 3] << 6)) <<  0) |
      (((indices[ 4] << 0) | (indices[ 5] << 2) | (indices[ 6] << 4) | (indices[ 7] << 6)) <<  8) |
      (((indices[ 8] << 0) | (indices[ 9] << 2) | (indices[10] << 4) | (indices[11] << 6)) << 16) |
      (((indices[12] << 0) | (indices[13] << 2) | (indices[14] << 4) | (indices[15] << 6)) << 24);
  }
}

void WriteColourBlock3(tile_barrier barrier, const int thread, lineC2 cline, inout index16 indices, out code64 block,
		       IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 2, 3, 4, 5, 6, 7 };	// (a >  b)

  // degenerate case "start > stop" (AMP: local register, not enough for group-shared)
  int colour[CVALS];
  int sorted[CVALS];

  // get the packed values (AMP: cline is guaranteed to be valid/written)
  FloatTo565(cline, colour);

  // colour[CSTRT] > colour[CSTOP] := sorted[CSTRT] == colour[CSTRT]
  sorted[CSTRT] = min(colour[CSTRT], colour[CSTOP]);
  sorted[CSTOP] = max(colour[CSTRT], colour[CSTOP]);

  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  // remap the indices
  wavefrnt_for(i, 16) {
#if	defined(SQUISH_USE_COMPUTE)
    indices[i] = (sorted[CSTRT] == colour[CSTOP] ? indices[i] : yArr[IBL_COLOR3][indices[i]].mapped);
#else
    indices[i] = (sorted[CSTRT] == colour[CSTOP] ? indices[i] : yArr(IBL_COLOR3, indices[i]).mapped);
#endif
  }

  // write the block
  WriteColourBlock(barrier, thread, sorted, indices, block);
}

void WriteColourBlock4(tile_barrier barrier, const int thread, lineC2 cline, inout index16 indices, out code64 block,
		       IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 3, 2, 5, 4, 7, 6 };	// (a <  b)
//static ccr8 slut[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };	// (a == b)

  // degenerate case "start > stop"
  int colour[CVALS];
  int sorted[CVALS];

  // get the packed values (AMP: cline is guaranteed to be valid/written)
  FloatTo565(cline, colour);

  // colour[CSTRT] < colour[CSTOP] := sorted[CSTRT] == colour[CSTOP]
  sorted[CSTRT] = max(colour[CSTRT], colour[CSTOP]);
  sorted[CSTOP] = min(colour[CSTRT], colour[CSTOP]);

  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  // remap the indices
  wavefrnt_for(i, 16) {
    indices[i] = (indices[i] ^ (sorted[CSTRT] == colour[CSTOP] ? 1 : 0)) & (sorted[CSTRT] == sorted[CSTOP] ? 0 : 3);
  }

  // write the block
  WriteColourBlock(barrier, thread, colour, indices, block);
}
#endif

} // namespace squish
