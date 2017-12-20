/* -----------------------------------------------------------------------------

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

#ifndef SQUISH_HDRINDEXFIT_H
#define SQUISH_HDRINDEXFIT_H

#include <squish.h>
#include <limits.h>

#include "hdrfit.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class HDRSet;
class HDRIndexFit : public virtual HDRFit
{
public:
  HDRIndexFit(HDRSet const* colours, int flags);
  
#ifdef	FEATURE_INDEXFIT_INLINED
private:
  Col3 m_qstart;
  Col3 m_qend;
  Vec3 m_qerror;
  Scr3 m_berror;

  doinline
  void ErrorEndPoints(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
		      Vec3 const* values, Scr3 const* freq,
		      int ib, int idxs);
  
  doinline
  Scr3 ErrorInterpolants(Vec3 const &metric, fQuantizer &q, Vec3 const* values, Scr3 const* freq,
			 int ib, int idxs, Vec3 &value0, Vec3 &value1, int closest0, int closest1);
  doinline
  Scr3 ErrorInterpolantsS(Vec3 const &metric, fQuantizer &q, Vec3 const* values, Scr3 const* freq,
			  int ib, int idxs, Vec3 &value0, int closest0);
  doinline
  Scr3 ErrorInterpolantsE(Vec3 const &metric, fQuantizer &q, Vec3 const* values, Scr3 const* freq,
			  int ib, int idxs, Vec3 &value1, int closest1);
  
  doinline
  void BetterInterpolants(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
			  Vec3 const* values, Scr3 const* freq,
			  int ib, int idxs, Vec3 &value0, Vec3 &value1, int closest0, int closest1);
  doinline
  void BetterInterpolantsS(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
			   Vec3 const* values, Scr3 const* freq,
			   int ib, int idxs, Vec3 &value0, int closest0);
  doinline
  void BetterInterpolantsE(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
			   Vec3 const* values, Scr3 const* freq,
			   int ib, int idxs, Vec3 &value1, int closest1);
#endif

protected:
  Scr3 StretchEndPoints(int set, Vec3 const &metric, fQuantizer &q, int ib, u8 (&closest)[16]);
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#endif // ndef SQUISH_PALETTEINDEXFIT_H
