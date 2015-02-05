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

#ifndef SQUISH_COLOURSINGLESNAP_H
#define SQUISH_COLOURSINGLESNAP_H

#include <squish.h>
#include <limits.h>
#include "colourfit.h"

// pull in structure definitions
#if	defined(SQUISH_USE_AMP)
#include "coloursinglelookup_ccr.inl"
#endif

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class ColourSet;
struct ColourSingleLookup;
class ColourSingleSnap : public ColourFit
{
public:
  ColourSingleSnap(ColourSet const* colours, int flags);
  
  // error management
  void SetError(Scr4 &error) { m_besterror = error; }
  void SetError(Scr3 &error) { m_besterror = error; }
  Scr3 GetError() { return m_besterror; }

private:
  virtual void Compress3b(void* block);
  virtual void Compress3(void* block);
  virtual void Compress4(void* block);

  u8   m_colour[4];
  Vec3 m_start;
  Vec3 m_end;
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#endif // ndef SQUISH_COLOURSINGLESNAP_H
