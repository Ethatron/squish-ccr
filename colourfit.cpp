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

#include "colourfit.h"
#include "colourset.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
extern Vec4 g_metric[8];

ColourFit::ColourFit( ColourSet const* colours, int flags )
  : m_colours(colours), m_flags(flags)
{
  // initialize the metric
  m_metric = KillW(g_metric[(flags & kColourMetrics) >> 4]);

  // initialize the best error
  m_besterror = Scr4(FLT_MAX);
}

void ColourFit::Compress( void* block )
{
  const bool isBtc1f = ((m_flags & kBtcp) == kBtc1);
  const bool isBtc1b = ((m_flags & kExcludeAlphaFromPalette) != 0);

  if (isBtc1f) {
    Compress3(block);
    if (!m_colours->IsTransparent())
      Compress4(block);
    if (isBtc1b)
      {/*Compress3b(block)*/;}
  }
  else
    Compress4(block);
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
void ColourFit_CCR::AssignSet(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, int metric, int fit) amp_restricted
{
  // -1 remap goes to source[0][15]
  wavefrnt_for(miscan, 16) {
    m_matches[0][miscan] = 3;
//  m_matches[1][miscan] = 0;
  }
}

#if 0
void ColourFit_CCR::Compress(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block, const bool trans,
			     IndexBlockLUT yArr) amp_restricted
{
  /* all or nothing branches, OK, same for all threads */
  bool isBtc1 = (trans);
  if (isBtc1 && m_colours.IsTransparent())
    Compress3 (barrier, thread, m_colours, block, yArr);
  else if (!isBtc1)
    Compress4 (barrier, thread, m_colours, block, yArr);
  else
    Compress34(barrier, thread, m_colours, block, yArr);
}
#endif

void ColourFit_CCR::Compress3(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
			      IndexBlockLUT yArr) amp_restricted
{
  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    m_colours.WriteSet(barrier, thread, m_line, m_matches, block, 0, yArr);
  }
}

void ColourFit_CCR::Compress4(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
			      IndexBlockLUT yArr) amp_restricted
{
  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    m_colours.WriteSet(barrier, thread, m_line, m_matches, block, 1, yArr);
  }
}

void ColourFit_CCR::Compress34(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block, const int is4,
			       IndexBlockLUT yArr) amp_restricted
{
  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    m_colours.WriteSet(barrier, thread, m_line, m_matches, block, is4, yArr);
  }
}
#endif

} // namespace squish
