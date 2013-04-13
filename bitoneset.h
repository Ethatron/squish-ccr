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

#ifndef SQUISH_BITONESET_H
#define SQUISH_BITONESET_H

#include <squish.h>
#include "maths.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
/*! @brief Represents a set of block values
 */
class BitoneSet
{
public:
  BitoneSet(u8  const* rgba, int mask, int flags);
  BitoneSet(u16 const* rgba, int mask, int flags);
  BitoneSet(f23 const* rgba, int mask, int flags);

  bool IsUnweighted() const { return m_unweighted; }

  int GetCount() const { return m_count; }
  Vec3 const* GetPoints() const { return m_points; }
  Scr3 const* GetWeights() const { return m_weights; }

  void RemapIndices(u8 const* source, u8* target) const;

private:
  bool  m_unweighted;
  int   m_count;
  Vec3  m_points[16];
  Scr3  m_weights[16];
  char  m_remap[16];
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace sqish

#endif // ndef SQUISH_BITONESET_H