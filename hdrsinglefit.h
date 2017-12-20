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

#ifndef SQUISH_HDRSINGLEFIT_H
#define SQUISH_HDRSINGLEFIT_H

#include <squish.h>
#include <limits.h>

#include "hdrfit.h"

// pull in structure definitions
#if	defined(SQUISH_USE_AMP)
#include "hdrsinglelookup_ccr.inl"
#endif

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
struct HDRSingleLookup2;
struct HDRSingleLookup4;
struct HDRSingleLookup8;

class HDRSet;
class HDRSingleFit : public virtual HDRFit
{
public:
  HDRSingleFit(HDRSet const* colours, int flags);

private:
  Scr4 ComputeEndPoints(int set, Vec4 const &metric, HDRSingleLookup2 const* const* lookups, u8 cmask);
  Scr4 ComputeEndPoints(int set, Vec4 const &metric, HDRSingleLookup4 const* const* lookups, u8 cmask);
  Scr4 ComputeEndPoints(int set, Vec4 const &metric, HDRSingleLookup8 const* const* lookups, u8 cmask);

  u8 m_entry[4][4];
  u8 m_index;

protected:
  Scr4 ComputeEndPoints(int set, Vec4 const &metric, int tb, int db, int ib, u8 cmask);
  u8 GetIndex() { return m_index; }
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#endif // ndef SQUISH_HDRSINGLEFIT_H
