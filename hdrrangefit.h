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

#ifndef SQUISH_HDRRANGEFIT_H
#define SQUISH_HDRRANGEFIT_H

#include <squish.h>
#include "maths.h"

#include "hdrsinglefit.h"
#include "hdrsinglesnap.h"
#include "hdrindexfit.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class HDRSet;
class HDRRangeFit : public HDRSingleMatch, public HDRIndexFit
{
public:
  HDRRangeFit(HDRSet const* palette, int flags);
  
  virtual void Compress(void* block, fQuantizer &q, int mode);

private:
  Vec3 m_centroid[2];
  Col3 m_qstart[2];
  Col3 m_qend  [2];

#ifdef	FEATURE_ELIMINATE_FLATBOOKS
  Vec4 m_start_candidate[4];
  Vec4 m_end_candidate[4];
#endif
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // squish

#endif // ndef SQUISH_HDRRANGEFIT_H
