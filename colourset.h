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

#ifndef SQUISH_COLOURSET_H
#define SQUISH_COLOURSET_H

#include <squish.h>
#include "maths.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
/*! @brief Represents a set of block colours
 */
class ColourSet
{
public:
  ColourSet(u8  const* rgba, int mask, int flags);
  ColourSet(u16 const* rgba, int mask, int flags);
  ColourSet(f23 const* rgba, int mask, int flags);

  bool IsTransparent() const { return m_transparent; }
  bool IsUnweighted() const { return m_unweighted; }

  int GetCount() const { return m_count; }
  Vec3 const* GetPoints() const { return m_points; }
  Scr3 const* GetWeights() const { return m_weights; }

  bool RemoveBlack(const Vec3 &metric, Scr3 &error);
  void RemapIndices(u8 const* source, u8* target) const;

private:
  bool  m_transparent, m_unweighted;
  int   m_count;
  Vec3  m_points[16];
  Scr3  m_weights[16];
  char  m_remap[16];
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
struct ColourSet_CCR
{
public_hlsl
  int GetCount() amp_restricted;
  point16 GetPoints() amp_restricted;
  weight16 GetWeights() amp_restricted;
  bool IsTransparent() amp_restricted;

  void CountSet    (tile_barrier barrier, const int thread,
		    pixel16 rgba, int mask, const bool tresh, const bool trans) amp_restricted;
  void WriteSet    (tile_barrier barrier, const int thread,
		    lineC2 cline, inout index16x2 source, out code64 block, const int is4,
		    IndexBlockLUT yArr) amp_restricted;

protected_hlsl
  void RemapIndices(tile_barrier barrier, const int thread,
		    inout index16x2 source) amp_restricted;

#if	!defined(SQUISH_USE_COMPUTE)
private_hlsl
  int m_transparent;
  int     m_count;
  float3 m_points[16];
  float m_weights[16];
  int     m_remap[16];
  ccr8  m_indices[16];
#endif
};

#if	defined(SQUISH_USE_COMPUTE)
  tile_static int m_transparent;
  tile_static int     m_count;
  tile_static float3 m_points[16];
  tile_static float m_weights[16];
  tile_static int     m_remap[16];
  tile_static ccr8  m_indices[16];
#endif

#if	!defined(SQUISH_USE_COMPUTE)
  typedef ColourSet_CCR &ColourSet_CCRr;
#else
  typedef ColourSet_CCR  ColourSet_CCRr;
#endif

#endif

} // namespace sqish

#endif // ndef SQUISH_COLOURSET_H