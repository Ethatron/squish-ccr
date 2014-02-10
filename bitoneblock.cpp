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

#include "bitoneblock.h"

#include "inlineables.cpp"

namespace squish {
  
/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
static void WriteBitoneBlock(int a, int b, u8 const* indices, void* block)
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

    // [3-0] [7-4] [11-8] [15-12] big endian dword
    // [15-12] [11-8] [7-4] [3-0] little endian dword
    bytes[4 + i] =
      (ind[0] << 0) +
      (ind[1] << 2) +
      (ind[2] << 4) +
      (ind[3] << 6);
  }
}

static void WriteBitoneBlock(int a, int b, Col4 const& indices, void* block)
{
  // get the block as ints
  unsigned int* ints = (unsigned int*)block;

  Col4 reindexed =
    (indices      ) +
    (indices >>  6) +
    (indices >> 12) +
    (indices >> 18);
  
  // write the endpoints
  ints[0] = (b << 16) + (a);

  // write the indices
  // [3-0] [7-4] [11-8] [15-12] big endian dword
  // [15-12] [11-8] [7-4] [3-0] little endian dword
  StoreUnaligned(reindexed & Col4(0x000000FF), (u8*)&ints[1]);
}

void WriteBitoneBlock4(Vec3::Arg start, Vec3::Arg end, u8 const* indices, void* block)
{
  // get the packed values
  int a = FloatTo88(start);
  int b = FloatTo88(end);

  // remap the indices
  Col4 remapped = Col4(indices);

  if (a < b) {
    // swap a and b
    std::swap(a, b);
    // swap index 0 and 1, 2 and 3
    remapped ^= Col4(0x01010101);
  }
  else if (a == b) {
    // use index 0
    remapped  = Col4(0x00000000);
  }
  
  // write the block
  WriteBitoneBlock(a, b, remapped, block);
}

void ReadBitoneBlock(
  u8 (&codes  )[16],
  u8 (&indices)[16],
  void const* block)
{
  // get the block bytes
  u8 const* bytes = reinterpret_cast< u8 const* >(block);

  // unpack the endpoints
  Unpack88(bytes + 0, codes + 0);
  Unpack88(bytes + 2, codes + 4);

  // generate the midpoints
  Codebook3or4(codes, false);

  // unpack the indices
  for (int i = 0; i < 4; ++i) {
    u8* ind = indices + 4 * i;
    u8 packed = bytes[4 + i];

    ind[0] = (packed >> 0) & 0x3;
    ind[1] = (packed >> 2) & 0x3;
    ind[2] = (packed >> 4) & 0x3;
    ind[3] = (packed >> 6) & 0x3;
  }
}

void DecompressBitonesCtx1u(u8* rgba, void const* block)
{
  u8 codes[16];
  u8 indices[16];

  ReadBitoneBlock(codes, indices, block);

  // store out the bitones
  for (int i = 0; i < 16; ++i) {
    u8 offset = 4 * indices[i];

    rgba[4 * i + 0] = codes[offset + 0] * (255 / 255);
    rgba[4 * i + 1] = codes[offset + 1] * (255 / 255);
    rgba[4 * i + 2] = codes[offset + 2] * (255 / 255);
    rgba[4 * i + 3] = codes[offset + 3] * (255 / 255);
  }
}

void DecompressBitonesCtx1u(u16* rgba, void const* block)
{
  u8 codes[16];
  u8 indices[16];

  ReadBitoneBlock(codes, indices, block);

  // store out the bitones
  for (int i = 0; i < 16; ++i) {
    u8 offset = 4 * indices[i];

    rgba[4 * i + 0] = codes[offset + 0] * (65535 / 255);
    rgba[4 * i + 1] = codes[offset + 1] * (65535 / 255);
    rgba[4 * i + 2] = codes[offset + 2] * (65535 / 255);
    rgba[4 * i + 3] = codes[offset + 3] * (65535 / 255);
  }
}

void DecompressBitonesCtx1u(f23* rgba, void const* block)
{
  u8 codes[16];
  u8 indices[16];

  ReadBitoneBlock(codes, indices, block);

  // store out the bitones
  for (int i = 0; i < 16; ++i) {
    u8 offset = 4 * indices[i];

    rgba[4 * i + 0] = codes[offset + 0] * (1.0f / 255.0f);
    rgba[4 * i + 1] = codes[offset + 1] * (1.0f / 255.0f);
    rgba[4 * i + 2] = codes[offset + 2] * (1.0f / 255.0f);
    rgba[4 * i + 3] = codes[offset + 3] * (1.0f / 255.0f);
  }
}

void DecompressNormalsCtx1u(u8* xyzd, void const* block)
{
  const Vec3 scale  = Vec3( 1.0f / 127.5f);
  const Vec3 offset = Vec3(-1.0f * 127.5f);
  const Vec3 scalei = Vec3( 1.0f * 127.5f);
  
  u8 codes[16];
  u8 indices[16];

  ReadBitoneBlock(codes, indices, block);

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i) {
    Col3 _xyz0  = Col3(codes[indices[i] * 4 + 0], codes[indices[i] * 4 + 1]);
    Vec3 cxyz0  = scale * (offset + _xyz0);
         cxyz0  = Complement<DISARM>(cxyz0);
	 cxyz0  = (scalei * cxyz0) - offset;
	 _xyz0  = FloatToInt<true>(cxyz0);
	 _xyz0 |= Col4(0,0,0,0xFF).GetCol3();

    StoreUnaligned(_xyz0, &xyzd[4 * i]);
  }
}

void DecompressNormalsCtx1u(u16* xyzd, void const* block)
{
  const Vec3 scale   = Vec3( 1.0f / 127.5f);
  const Vec3 offset  = Vec3(-1.0f * 127.5f);
  const Vec3 scalei  = Vec3( 1.0f * 32767.5f);
  const Vec3 offseti = Vec3(-1.0f * 32767.5f);
  
  u8 codes[16];
  u8 indices[16];

  ReadBitoneBlock(codes, indices, block);

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i) {
    Col3 _xyz0  = Col3(codes[indices[i] * 4 + 0], codes[indices[i] * 4 + 1]);
    Vec3 cxyz0  = scale * (offset + _xyz0);
         cxyz0  = Complement<DISARM>(cxyz0);
	 cxyz0  = (scalei * cxyz0) - offseti;
	 _xyz0  = FloatToInt<true>(cxyz0);
	 _xyz0 |= Col4(0,0,0,0xFF).GetCol3();

    StoreUnaligned(_xyz0, &xyzd[4 * i]);
  }
}

void DecompressNormalsCtx1u(f23* xyzd, void const* block)
{
  const Vec3 scale  = Vec3( 1.0f / 127.5f);
  const Vec3 offset = Vec3(-1.0f * 127.5f);
  
  u8 codes[16];
  u8 indices[16];

  ReadBitoneBlock(codes, indices, block);

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i) {
    Col3 _xyz0  = Col3(codes[indices[i] * 4 + 0], codes[indices[i] * 4 + 1]);
    Vec3 cxyz0  = scale * (offset + _xyz0);
         cxyz0  = Complement<DISARM>(cxyz0);
    Vec4 dxyz0  = TransferW(cxyz0, Vec4(1.0f));

    StoreUnaligned(dxyz0, &xyzd[4 * i]);
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
