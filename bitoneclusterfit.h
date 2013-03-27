/* -----------------------------------------------------------------------------

	Copyright (c) 2006 Simon Brown                          si@sjbrown.co.uk
	Copyright (c) 2007 Ignacio Castano                   icastano@nvidia.com
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

#ifndef SQUISH_BITONECLUSTERFIT_H
#define SQUISH_BITONECLUSTERFIT_H

#include <squish.h>
#include "maths.h"
#include "simd.h"
#include "bitonefit.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class BitoneClusterFit : public BitoneFit
{
public:
  BitoneClusterFit(BitoneSet const* bitones, int flags);

public:
  enum {
    kMinIterations = 1,
    kMaxIterations = 15
  };

  static int SanitizeFlags(int flags) {
    if (flags > (kColourClusterFit * kMaxIterations))
      return (kColourClusterFit * kMaxIterations);
    if (flags < (kColourClusterFit * kMinIterations))
      return (kColourClusterFit * kMinIterations);

    return flags;
  }

private:
  bool ConstructOrdering(Vec3 const& axis, int iteration);

  void ClusterFit4Constant(void* block);
  void ClusterFit4(void* block);

  virtual void Compress4(void* block);

  int  m_iterationCount;
  Vec3 m_principle;
  Scr4 m_besterror;
  Vec4 m_xsum_wsum;
  Vec4 m_points_weights[16];
  a16 u8 m_order[16 * kMaxIterations];

  bool m_optimizable;
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#endif // ndef SQUISH_BITONECLUSTERFIT_H
