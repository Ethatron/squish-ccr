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

#ifndef SQUISH_COLOURFIT_H
#define SQUISH_COLOURFIT_H

#include <squish.h>
#include "maths.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class ColourSet;
class ColourFit
{
public:
  ColourFit( ColourSet const* colours, int flags );

  void Compress( void* block );

protected:
  virtual void Compress3b(void* block) = 0;
  virtual void Compress3(void* block) = 0;
  virtual void Compress4(void* block) = 0;

  ColourSet const* m_colours;
  int m_flags;

  Vec3 m_metric;
  Scr3 m_besterror;
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
struct ColourFit_CCR
{
public_hlsl
  void AssignSet (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, int metric, int fit) amp_restricted;
#if 0
  void Compress  (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, out code64 block, bool trans,
                  IndexBlockLUT yArr) amp_restricted;
#endif

protected_hlsl
  void Compress3 (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;
  void Compress4 (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;
  void Compress34(tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, out code64 block, const int is4,
                  IndexBlockLUT yArr) amp_restricted;

// Start- and end-point positions in arrays, for colors
#define	CSTRT	0
#define	CSTOP	1
#define	CVALS	2

#if	!defined(SQUISH_USE_COMPUTE)
protected_hlsl
  float3 m_line[CVALS];
  ccr8 m_matches[2][16];
#endif
};

#if	defined(SQUISH_USE_COMPUTE)
  tile_static float3 m_line[CVALS];
  tile_static ccr8 m_matches[2][16];
#endif
#endif

} // namespace squish

#endif // ndef SQUISH_COLOURFIT_H
