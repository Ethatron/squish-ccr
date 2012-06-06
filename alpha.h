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

namespace squish {

#if	!defined(USE_PRE)
  void CompressAlphaDxt3( u8 const* rgba, int mask, void* block );
  void CompressAlphaDxt5( u8 const* rgba, int mask, void* block );

  void DecompressAlphaDxt3( u8* rgba, void const* block );
  void DecompressAlphaDxt5( u8* rgba, void const* block );
#endif

#if	defined(USE_AMP) || defined(USE_COMPUTE)
  void CompressAlphaDxt3(tile_barrier barrier, const int thread, pixel16 rgba, int mask, out code64 block,
			 IndexBlockLUT yArr) amp_restricted;
  void CompressAlphaDxt5(tile_barrier barrier, const int thread, pixel16 rgba, int mask, out code64 block,
			 IndexBlockLUT yArr) amp_restricted;
#endif

} // namespace squish

#endif // ndef SQUISH_ALPHA_H
