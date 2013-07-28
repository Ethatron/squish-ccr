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

#ifndef SQUISH_ALPHA_H
#define SQUISH_ALPHA_H

#include <squish.h>
#include <limits.h>

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
  void CompressAlphaBtc2u(u8  const* rgba, int mask, void* block);
  void CompressAlphaBtc3u(u8  const* rgba, int mask, void* block, int flags);
  void CompressAlphaBtc3s(s8  const* rgba, int mask, void* block, int flags);
  void CompressDepthBtc4u(u8  const* rgba, int mask, void* block, int flags);
  void CompressDepthBtc4s(s8  const* rgba, int mask, void* block, int flags);
  
  void CompressAlphaBtc2u(u16 const* rgba, int mask, void* block);
  void CompressAlphaBtc3u(u16 const* rgba, int mask, void* block, int flags);
  void CompressAlphaBtc3s(s16 const* rgba, int mask, void* block, int flags);
  void CompressDepthBtc4u(u16 const* rgba, int mask, void* block, int flags);
  void CompressDepthBtc4s(s16 const* rgba, int mask, void* block, int flags);

  void CompressAlphaBtc2u(f23 const* rgba, int mask, void* block);
  void CompressAlphaBtc3u(f23 const* rgba, int mask, void* block, int flags);
  void CompressAlphaBtc3s(f23 const* rgba, int mask, void* block, int flags);
  void CompressDepthBtc4u(f23 const* rgba, int mask, void* block, int flags);
  void CompressDepthBtc4s(f23 const* rgba, int mask, void* block, int flags);

  void DecompressAlphaBtc2u(u8 * rgba, void const* block);
  void DecompressAlphaBtc3u(u8 * rgba, void const* block, int flags);
  void DecompressAlphaBtc3s(s8 * rgba, void const* block, int flags);
  void DecompressDepthBtc4u(u8 * rgba, void const* block, int flags);
  void DecompressDepthBtc4s(s8 * rgba, void const* block, int flags);
  
  void DecompressAlphaBtc2u(u16* rgba, void const* block);
  void DecompressAlphaBtc3u(u16* rgba, void const* block, int flags);
  void DecompressAlphaBtc3s(s16* rgba, void const* block, int flags);
  void DecompressDepthBtc4u(u16* rgba, void const* block, int flags);
  void DecompressDepthBtc4s(s16* rgba, void const* block, int flags);

  void DecompressAlphaBtc2u(f23* rgba, void const* block);
  void DecompressAlphaBtc3u(f23* rgba, void const* block, int flags);
  void DecompressAlphaBtc3s(f23* rgba, void const* block, int flags);
  void DecompressDepthBtc4u(f23* rgba, void const* block, int flags);
  void DecompressDepthBtc4s(f23* rgba, void const* block, int flags);
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
  void CompressAlphaBtc2u(tile_barrier barrier, const int thread,
			 pixel16 rgba, int mask, out code64 block,
			 IndexBlockLUT yArr) amp_restricted;
  void CompressAlphaBtc3u(tile_barrier barrier, const int thread,
			 pixel16 rgba, int mask, out code64 block,
			 IndexBlockLUT yArr) amp_restricted;

/*void DecompressAlphaBtc2(tile_barrier barrier, const int thread,
			   out pixel16 rgba, code64 block) amp_restricted;
  void DecompressAlphaBtc3(tile_barrier barrier, const int thread,
			   out pixel16 rgba, code64 block) amp_restricted;*/
#endif

} // namespace squish

#endif // ndef SQUISH_ALPHA_H
