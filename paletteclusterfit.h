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

#ifndef SQUISH_PALETTECLUSTERFIT_H
#define SQUISH_PALETTECLUSTERFIT_H

#include <squish.h>
#include "maths.h"
#include "simd.h"

#include "palettefit.h"
#include "palettesinglefit.h"
#include "palettesinglesnap.h"
#include "palettechannelfit.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(SQUISH_USE_PRE)
class PaletteClusterFit : public PaletteSingleMatch, public PaletteChannelFit
{
public:
  PaletteClusterFit(PaletteSet const* palettes, int flags, int swap = -1, int shared = -1);

  virtual void Compress(void* block, vQuantizer &q, int mode);

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
#define CLUSTERINDICES	3
  // separate components, 4/8 colors, 4/8 alphas
  void CompressS23(void* block, vQuantizer &q, int mode);
  // combined components, 4/16 colors+alphas
  void CompressC2(void* block, vQuantizer &q, int mode);
  void CompressC4(void* block, vQuantizer &q, int mode);

  bool ConstructOrdering(Vec4 const& axis, int iteration, int set);

  Scr4 ClusterSearch4Alpha(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb);

  Scr4 ClusterSearch4Constant(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb);
  Scr4 ClusterSearch8Constant(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb);

  Scr4 ClusterSearch4(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb);
  Scr4 ClusterSearch8(u8 (&closest)[4][16], int count, int set, Vec4 const &metric, vQuantizer &q, int sb);

  int  m_iterationCount;
  Vec4 m_principle[4];

  Vec4 m_xsum_wsum[4 * 2];
//Vec4 m_xxsum_wwsum[4];
  Vec4 m_points_weights[4][16 * 2];

  a16 u8 m_order[4][16 * kMaxIterations];

  bool m_optimizable[4];
};
#endif

// -----------------------------------------------------------------------------
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
struct ClusterFit_CCR : inherit_hlsl PaletteFit_CCR
{
public_hlsl
  void AssignSet (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palettes, const int metric, const int fit) amp_restricted;
  void Compress  (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palettes, out code64 block, const bool trans,
                  IndexBlockLUT yArr) amp_restricted;

protected_hlsl
  void Compress3 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palettes, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;
  void Compress4 (tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palettes, out code64 block,
                  IndexBlockLUT yArr) amp_restricted;

private_hlsl
  bool ConstructOrdering(tile_barrier barrier, const int thread,
                  PaletteSet_CCRr m_palettes, float4r axis, int iteration) amp_restricted;

#define	kMaxIterations	8

#if	!defined(SQUISH_USE_COMPUTE)
  int m_iterationCount;
  float4 m_principle;
  ccr8 m_order[kMaxIterations][16];
  float4 m_points_weights[16];
  float4 m_xsum_wsum;
  float4 m_metric4;
  float m_besterror;
#endif
};

#if	defined(SQUISH_USE_COMPUTE)
  tile_static int m_iterationCount;
  tile_static float4 m_principle;
  tile_static ccr8 m_order[kMaxIterations][16];
  tile_static float4 m_points_weights[16];
  tile_static float4 m_xsum_wsum;
  tile_static float4 m_metric4;
  tile_static float m_besterror;
#endif
#endif
} // namespace squish

#endif // ndef SQUISH_PALETTECLUSTERFIT_H
