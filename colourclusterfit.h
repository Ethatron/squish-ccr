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

#ifndef SQUISH_COLOURCLUSTERFIT_H
#define SQUISH_COLOURCLUSTERFIT_H

#include <squish.h>
#include "maths.h"
#include "simd.h"
#include "colourfit.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class ColourClusterFit : public ColourFit
{
public:
  ColourClusterFit(ColourSet const* colours, int flags);

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
  void SumError3(u8 (&closest)[16], Vec4 &beststart, Vec4 &bestend, Scr4 &besterror);
  void SumError4(u8 (&closest)[16], Vec4 &beststart, Vec4 &bestend, Scr4 &besterror);

  void ComputeEndPoints();
  bool ConstructOrdering(Vec3 const& axis, int iteration);

  void ClusterFit3Constant(void* block);
  void ClusterFit4Constant(void* block);

  void ClusterFit3(void* block);
  void ClusterFit4(void* block);

  virtual void Compress3b(void* block);
  virtual void Compress3(void* block);
  virtual void Compress4(void* block);

  int  m_iterationCount;
  Vec3 m_principle;
  Vec4 m_xsum_wsum;
  Vec4 m_points_weights[16];
  a16 u8 m_order[16 * kMaxIterations];

  bool m_optimizable;
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
struct ClusterFit_CCR : inherit_hlsl ColourFit_CCR
{
public_hlsl
  void AssignSet (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, const int metric, const int fit) amp_restricted;
  void Compress  (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, out code64 block, const bool trans,
                  IndexBlockLUT yArr) amp_restricted;

protected_hlsl
  void Compress3 (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;
  void Compress4 (tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;

private_hlsl
  bool ConstructOrdering(tile_barrier barrier, const int thread,
                  ColourSet_CCRr m_colours, float3r axis, int iteration) amp_restricted;

#define	kMaxIterations	8

#if	!defined(SQUISH_USE_COMPUTE)
  int m_iterationCount;
  float3 m_principle;
  ccr8 m_order[kMaxIterations][16];
  float4 m_points_weights[16];
  float4 m_xsum_wsum;
  float4 m_metric4;
  float m_besterror;
#endif
};

#if	defined(SQUISH_USE_COMPUTE)
  tile_static int m_iterationCount;
  tile_static float3 m_principle;
  tile_static ccr8 m_order[kMaxIterations][16];
  tile_static float4 m_points_weights[16];
  tile_static float4 m_xsum_wsum;
  tile_static float4 m_metric4;
  tile_static float m_besterror;
#endif
#endif
} // namespace squish

#endif // ndef SQUISH_COLOURCLUSTERFIT_H
