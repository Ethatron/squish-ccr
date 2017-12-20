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

#include <squish.h>
#include <assert.h>
#include <memory.h>

#include "alpha.h"

#include "bitoneset.h"
#include "colourset.h"
#include "paletteset.h"
#include "hdrset.h"

#include "maths.h"

// Btc2/Btc3/Btc4/Btc5
#include "alphanormalfit.h"

// Ctx1
#include "bitonenormalfit.h"
#include "bitonerangefit.h"
#include "bitoneclusterfit.h"
#include "bitoneblock.h"

// Btc1/Btc2/Btc3
#include "colournormalfit.h"
#include "colourrangefit.h"
#include "colourclusterfit.h"
#include "colourblock.h"

// Btc7
#include "palettenormalfit.h"
#include "paletterangefit.h"
#include "paletteclusterfit.h"
#include "paletteblock.h"

// Btc6
#include "hdrrangefit.h"
//nclude "hdrclusterfit.h"
#include "hdrblock.h"

#include "coloursinglefit.h"
#include "coloursinglesnap.h"
#include "palettesinglefit.h"
#include "palettesinglesnap.h"
#include "hdrsinglefit.h"
#include "hdrsinglesnap.h"

namespace squish {

#pragma warning(disable: 4482)

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
int SanitizeFlags(int flags)
{
  // grab the flag bits
  int method = flags & (kBtcp);
  int fit    = flags & (kColourRangeFit | kAlphaIterativeFit | kColourIterativeClusterFits);
  int metric = flags & (kColourMetrics);
  int extra  = flags & (kWeightColourByAlpha);
  int mode   = flags & (kVariableCodingModes);
  int map    = flags & (kSrgbExternal | kSrgbInternal | kSignedExternal | kSignedInternal);

  // set defaults
  if (!method || ((method > kCtx1)))
    method = kBtc1;
  if (!metric || ((metric > kColourMetricUnit) && (metric != kColourMetricCustom)))
    metric = kColourMetricPerceptual;

  if (!fit)
    fit = (kColourClusterFit * 1);
  if (fit & kColourIterativeClusterFits) {
    if (method <= kBtc3)
      fit = ColourClusterFit::SanitizeFlags(fit);
    if (method == kBtc7)
      fit = PaletteClusterFit::SanitizeFlags(fit);
  }

  if ((method == kBtc6) && (mode > kVariableCodingMode14))
    mode = 0;
  if ((method == kBtc7) && (mode > kVariableCodingMode8))
    mode = 0;

  // done
  return method + fit + metric + extra + mode + map;
}

/* *****************************************************************************
 */
Vec4 g_metric[8] =
{
#ifdef FEATURE_METRIC_ROOTED
  // sum squared is 2.0f
  Vec4(0.5773f, 0.5773f, 0.5773f, 1.0f),
  Vec4(0.4611f, 0.8456f, 0.2687f, 1.0f),	// kColourMetricPerceptual
  Vec4(0.5773f, 0.5773f, 0.5773f, 1.0f),	// kColourMetricUniform
  Vec4(0.7071f, 0.7071f, 0.0000f, 1.0f),	// kColourMetricUnit
  Vec4(0.5773f, 0.5773f, 0.5773f, 1.0f),	// kColourMetricGray
  Vec4(0.5773f, 0.5773f, 0.5773f, 1.0f),
  Vec4(0.5773f, 0.5773f, 0.5773f, 1.0f),
  Vec4(0.5773f, 0.5773f, 0.5773f, 1.0f)		// kColourMetricCustom
#else
  // sum is 2.0f
  Vec4(0.3333f, 0.3333f, 0.3333f, 1.0f),
  Vec4(0.2126f, 0.7152f, 0.0722f, 1.0f),	// kColourMetricPerceptual
  Vec4(0.3333f, 0.3333f, 0.3333f, 1.0f),	// kColourMetricUniform
  Vec4(0.5000f, 0.5000f, 0.0000f, 1.0f),	// kColourMetricUnit
  Vec4(0.3333f, 0.3333f, 0.3333f, 1.0f),	// kColourMetricGray
  Vec4(0.3333f, 0.3333f, 0.3333f, 1.0f),
  Vec4(0.3333f, 0.3333f, 0.3333f, 1.0f),
  Vec4(0.3333f, 0.3333f, 0.3333f, 1.0f)		// kColourMetricCustom
#endif
};

void SetWeights(int flags, const f23* rgba)
{
  // initialize the metric
  const bool custom = ((flags & kColourMetrics) == kColourMetricCustom);

  if (custom)
  {
    g_metric[7] = Vec4(rgba[0], rgba[1], rgba[2], 1.0f);
    g_metric[7] /= Vec4(HorizontalAdd(g_metric[7].GetVec3()), 1.0f);

#ifdef FEATURE_METRIC_ROOTED
    g_metric[7] = Sqrt(g_metric[7]);
#endif
  }
}

/* *****************************************************************************
 */
template<typename dtyp>
void CompressBitoneCtx1u(dtyp const* rgba, int mask, void* block, int flags)
{
  // create the minimal point set
  BitoneSet colours(rgba, mask, flags);

  if (((flags & kColourRangeFit) != 0) || (colours.GetCount() == 0)) {
    // do a range fit
    BitoneRangeFit fit(&colours, flags);
    fit.Compress(block);
  }
  else {
    // default to a cluster fit (could be iterative or not)
    BitoneClusterFit fit(&colours, flags);
    fit.Compress(block);
  }
}

template<typename dtyp>
void CompressNormalCtx1u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // create the minimal point set
  BitoneSet bitones(xyzd, mask, flags);

  // check the compression type and compress normals
  {
    // do a normal fit
    BitoneNormalFit fit(&bitones, flags);
    fit.Compress(block);
  }
}

template<typename dtyp>
void CompressColourBtc1u(dtyp const* rgba, int mask, void* block, int flags)
{
  // create the minimal point set
  ColourSet colours(rgba, mask, flags);

  // check the compression type and compress colour
  if (colours.GetCount() == 1) {
    // always do a single colour fit
    ColourSingleMatch fit(&colours, flags);
    fit.Compress(block);
  }
  else if (((flags & kColourRangeFit) != 0) || (colours.GetCount() == 0)) {
    // do a range fit
    ColourRangeFit fit(&colours, flags);
    fit.Compress(block);
  }
  else {
    // default to a cluster fit (could be iterative or not)
    ColourClusterFit fit(&colours, flags);
    fit.Compress(block);
  }
}

template<typename dtyp>
void CompressNormalBtc1u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // create the minimal point set
  ColourSet normals(xyzd, mask, flags);

  // check the compression type and compress normals
  if (normals.GetCount() == 1) {
    // always do a single colour fit
    ColourSingleMatch fit(&normals, flags);
    fit.Compress(block);
  }
  else {
    // do a range fit
    ColourNormalFit fit(&normals, flags);
    fit.Compress(block);
  }
}

#if defined(TRACK_STATISTICS)
struct statistics gstat = {0};
#endif

template<typename dtyp, class PaletteTypeFit>
Scr4 CompressPaletteBtc7uV1(dtyp const* rgba, int mask, void* block, int flags)
{
#if !defined(NDEBUG) && defined(DEBUG_SETTING)
#define DEBUG_MODE	kVariableCodingMode1
#define DEBUG_FIT	kColourRangeFit   //kColourClusterFit * 15

  flags = (flags & (~kVariableCodingModes)) | (DEBUG_MODE);
  flags = (flags & (~kColourIterativeClusterFit)) | (DEBUG_FIT);
#endif

  /* we start with 1 set so we get some statistics about the color-
   * palette, based on that we decide if we need to search into higher
   * number of sets
   *
   * observations:
   * - if there is 1-2 color(s), we only need 1 set:
   *   the index-precision doesn't matter in that case and we choose
   *   the coding with the highest start/end-point precision
   *   nevertheless the 2 colors are not necessarily also start/end
   *   interpolated colors may achieve superior precision
   * - if there is 3-4 color(s), we only need 2 sets:
   *   the available partitions may not correspond exactly to the
   *   distribution of the 3-4 colors, so for maximum quality we need
   *   to do the whole search regardless (including 3 sets)
   *   if we've found a 2 set partition with 1-2 colors in each we can
   *   abort immediately
   *
   * rangefit searches for the best configuration (partition/swap/rotation)
   * optionally clusterfit makes the best of that partition
   */
  static const int modeorder[3][8] = {
    {
#define MODEORDER_EXPL	    0
#define MODEORDER_EXPL_MIN  0
#define MODEORDER_EXPL_MAX  7
      // order: mode (lo to hi)
      kVariableCodingMode1, //{ 3, 4, 0, 0,  4, 0, 1,  0,  3, 0 },
      kVariableCodingMode2, //{ 2, 6, 0, 0,  6, 0, 0,  1,  3, 0 },
      kVariableCodingMode3, //{ 3, 6, 0, 0,  5, 0, 0,  0,  2, 0 },
      kVariableCodingMode4, //{ 2, 6, 0, 0,  7, 0, 1,  0,  2, 0 },
      kVariableCodingMode5, //{ 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
      kVariableCodingMode6, //{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },
      kVariableCodingMode7, //{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
      kVariableCodingMode8, //{ 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 },
    },
    {
#define MODEORDER_OPAQ	    1
#define MODEORDER_OPAQ_MIN  0
#define MODEORDER_OPAQ_MAX  6
      // order: sets (lo to hi), ibs (hi to lo), prc (hi to lo)
      kVariableCodingMode7, //{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
      kVariableCodingMode5, //{ 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
      kVariableCodingMode6, //{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },

      kVariableCodingMode2, //{ 2, 6, 0, 0,  6, 0, 0,  1,  3, 0 },  // non-alpha variant of mode 8
      kVariableCodingMode4, //{ 2, 6, 0, 0,  7, 0, 1,  0,  2, 0 },  // non-alpha variant of mode 8

      kVariableCodingMode1, //{ 3, 4, 0, 0,  4, 0, 1,  0,  3, 0 },
      kVariableCodingMode3, //{ 3, 6, 0, 0,  5, 0, 0,  0,  2, 0 },

      0,
    },
    {
#define MODEORDER_TRNS	    2
#define MODEORDER_TRNS_MIN  0
#define MODEORDER_TRNS_MAX  3
      // order: sets (lo to hi), ibs (hi to lo), prc (hi to lo)
      kVariableCodingMode7, //{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
      kVariableCodingMode5, //{ 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
      kVariableCodingMode6, //{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },

      kVariableCodingMode8, //{ 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 },  // alpha variant of mode 2/4

      0, 0, 0, 0
    }
  };

  int numm = flags &  ( kVariableCodingModes),
        sm = (numm == 0 ? MODEORDER_OPAQ_MIN : (numm >> 24) - 1),
	em = (numm == 0 ? MODEORDER_OPAQ_MAX :               sm),
	om = (numm == 0 ? MODEORDER_OPAQ     :   MODEORDER_EXPL);
             flags &= (~kVariableCodingModes);

#undef MODEORDER_EXPL
#undef MODEORDER_EXPL_MIN
#undef MODEORDER_EXPL_MAX

  // limits sets to 3 and choose the partition freely
  int lmts =  3;
  int lmtp = -1;

  // use the same data-structure all the time
  PaletteSet bestpal;
  int bestmde = -1;
  int bestswp = -1;
  int bestbit = -1;
  int besttyp = -1;

  Scr4 error(FLT_MAX);

  for (int m = sm; m <= em; m++) {
    int mode = modeorder[om][m];
    int mnum = (mode >> 24) - 1;

    // a mode has a specific number of sets, and variable rotations and partitions
    int nums = PaletteFit::GetNumSets      (mnum);
    int numr = PaletteFit::GetRotationBits (mnum);
    int nump = PaletteFit::GetPartitionBits(mnum);
    int numx = PaletteFit::GetSelectionBits(mnum);
    int numb = PaletteFit::GetSharedBits   (mnum);
    int numi = PaletteFit::GetIndexBits    (mnum);

    // stop if set-limit reached
    if (nums > lmts)
      break;

    // lock on the perfect partition
    int sp = (lmtp == -1 ?               0 : lmtp),
	ep = (lmtp == -1 ? (1 << nump) - 1 : lmtp);
    // search through rotations
    int sr =                             0,
	er =               (1 << numr) - 1;
    // search through index-swaps
    int sx =                             0,
	ex =               (1 << numx) - 1;
    // search through shared bits
#ifdef FEATURE_SHAREDBITS_TRIALS
    int sb = (numb > 0   ?               0 : SBSKIP),
	eb = (numb > 0   ?       numb      : SBSKIP);
#else
    int sb = (numb > 0   ?          SBSKIP : SBSKIP),
	eb = (numb > 0   ?          SBSKIP : SBSKIP);
#endif

#if !defined(NDEBUG) && defined(DEBUG_SETTING)
#define DEBUG_PARTITION	0
#define DEBUG_ROTATION	0
#define DEBUG_SELECTION	0
#define DEBUG_SHAREDBIT	(numb > 0 ? 0 : -1)

//  sp = ep = DEBUG_PARTITION;
    sr = er = DEBUG_ROTATION;
    sx = ex = DEBUG_SELECTION;
//  sb = eb = DEBUG_SHAREDBIT;
#endif

    int cb = PaletteFit::GetPrecisionBits(mnum);
    int ab = cb >> 16; cb = cb & 0xFF;

    // create the initial point set and quantizer
    PaletteSet initial(rgba, mask, flags + mode);
    vQuantizer qnt(cb, cb, cb, ab);

    // signal if we do we have anything better this iteration of the search
    bool better = false;
    // check if we can do a cascade with the cluster-fit (merged alpha 4 bit is the only exception)
    bool cluster = PaletteTypeFit::IsClusterable(flags) && (((numi >>  0) & 0xFF) <= CLUSTERINDICES)
                                                        && (((numi >> 16) & 0xFF) <= CLUSTERINDICES);

    // if we see we have transparent values, back up from trying to test non-alpha only modes
    // this will affect only successive trials, if an explicit mode is requested it's a NOP
    if (initial.IsTransparent())
      om = MODEORDER_TRNS, em = MODEORDER_TRNS_MAX;

    // if we see we have no transparent values, don't try non-rotated palettes (alpha is constant for all)
    if (!initial.IsTransparent() && initial.IsSeperateAlpha())
      sr = 1;
#if	defined(FEATURE_SHAREDBITS_TRIALS)
    // if we see we have no transparent values, force all shared bits to 1, or non-opaque codebook-entries occur
    // the all transparent case isn't so crucial, when we use IGNORE_ALPHA0 it's redundant to force 0 anyway
    if (!initial.IsTransparent() && initial.IsMergedAlpha())
      sb = eb;
    // otherwise just use the most occurring bit (parity) for all other cases
    // otherwise just use the most occurring bit (parity) for all non-alpha cases
    else if (((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_ALPHAONLYOPAQUE)) ||
	     ((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_ALPHAONLY) && (mode != kVariableCodingMode7) && !initial.IsTransparent()) ||
	     ((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_LOWPRC) && (mode < kVariableCodingMode5) && (mode != kVariableCodingMode1)))
      sb = eb = SBSKIP;
#endif

    // choose rotation or partition, they're mutually exclusive
    int spr = (er ? sr : sp),
	epr = (er ? er : ep);

    // search for the best partition/rotation
    for (int pr = spr; pr <= epr; pr++) {
      // create the minimal point set
      PaletteSet palette(initial, mask, flags + mode, pr);

#if 0
      // if we see we have less colors than sets, 
      // then back up from trying to test with more sets
      // or even other partitions
      if (palette.GetOptimal())
	lmtp = pr, lmts = nums;
#endif

#if defined(TRACK_STATISTICS)
      for (int xu = 0; xu < nums; xu++) {
	int cnt = palette.GetCount(xu);
	gstat.num_counts[mnum][p][xu][cnt]++;
#ifdef	FEATURE_TEST_LINES
	if (cnt > 2) {
	  int chn = palette.GetChannel(xu) + 1;
	  gstat.num_channels[mnum][p][xu][chn]++;
	}
#endif
      }

      if (palette.GetCount() <= nums)
	gstat.has_countsets[nums]++;
#endif

      // do a range fit (which uses single palette fit if appropriate)
      PaletteTypeFit fit(&palette, flags + mode);
      
      // TODO: swap & shared are mutual exclusive

      // search for the best swap
      for (int x = sx; x <= ex; x++) {
	fit.ChangeSwap(x);
	// search for the best shared bit
	for (int b = sb; b <= eb; b++) {
	  fit.ChangeShared(b);
	  
	  // update with old best error (reset IsBest)
	  fit.SetError(error);

	  // we could code it lossless, no point in trying any further at all
	  fit.PaletteTypeFit::Compress(block, qnt, mnum);
	  if (fit.IsBest()) {
	    if (fit.Lossless())
	      return Scr4(0.0f);

#if   !defined(TRACK_STATISTICS) && !defined(VERIFY_QUANTIZER)
	    if (cluster)
#endif
	    {
	      bestmde = mode,
	      bestpal = palette,
	      bestswp = x,
	      bestbit = b,
	      besttyp = 0,
	      better  = true;
	    }
	  }
      
	  // update with new best error
	  error = fit.GetError();
	}
      }
    }

    // check the compression type and compress palette of the chosen partition even better
    if (better && cluster) {
      int degree = (flags & kColourIterativeClusterFits);

      // default to a cluster fit (could be iterative or not)
      PaletteClusterFit fit(&bestpal, flags + mode);
      
      // we want the whole shebang, this takes looong!
      if (degree < (kColourClusterFit * 15))
	sb = eb = bestbit;
      if (degree < (kColourClusterFit * 14))
	sx = ex = bestswp;

      // search for the best swap
      for (int x = sx; x <= ex; x++) {
	fit.ChangeSwap(x);
	// search for the best shared bit
	for (int b = sb; b <= eb; b++) {
	  fit.ChangeShared(b);
	  
	  // update with old best error (reset IsBest)
	  fit.SetError(error);

	  // we could code it lossless, no point in trying any further at all
	  fit.PaletteClusterFit::Compress(block, qnt, mnum);
	  if (fit.IsBest()) {
	    if (fit.Lossless())
	      return Scr4(0.0f);

	    if (cluster || 1)
	      besttyp = 1;
	  }

#if defined(TRACK_STATISTICS)
	  gstat.btr_cluster[mnum][fit.IsBest() ? 1 : 0]++;
#endif
      
	  // update with new best error
	  error = fit.GetError();
	}
      }
    }

#if defined(TRACK_STATISTICS)
    gstat.win_partition[mnum][bestpal.GetPartition()]++;
    gstat.win_rotation [mnum][bestpal.GetRotation ()]++;
    gstat.win_swap     [mnum][bestpal.GetRotation ()][bestswp]++;
#endif
  }

#if defined(TRACK_STATISTICS)
  gstat.win_mode[(bestmde >> 24) - 1]++;
  gstat.win_cluster[(bestmde >> 24) - 1][besttyp]++;
#endif

#if defined(VERIFY_QUANTIZER)
  int cb = PaletteFit::GetPrecisionBits((bestmde >> 24) - 1);
  int ab = cb >> 16; cb = cb & 0xFF;

  // create the initial point set and quantizer
  vQuantizer qnt(cb, cb, cb, ab);

  if (!besttyp) {
    // do a range fit (which uses single palette fit if appropriate)
    PaletteRangeFit fit(&bestpal, flags + bestmde, bestswp, bestbit);

    fit.Compress(block, qnt, (bestmde >> 24) - 1);
    fit.Decompress((u8*)rgba, qnt, (bestmde >> 24) - 1);
  }
  else {
    // default to a cluster fit (could be iterative or not)
    PaletteClusterFit fit(&bestpal, flags + bestmde, bestswp, bestbit);

    fit.Compress(block, qnt, (bestmde >> 24) - 1);
    fit.Decompress((u8*)rgba, qnt, (bestmde >> 24) - 1);
  }
#endif

#if defined(VERIFY_ENCODER)
  DecompressColoursBtc7u((u8*)rgba, block);
#endif

  return error;
}

template<typename dtyp, class PaletteTypeFit>
Scr4 CompressPaletteBtc7uV2(dtyp const* rgba, int mask, void* block, int flags)
{
  vQuantizer q7778(7, 7, 7, 8);
  vQuantizer q5556(5, 5, 5, 6);
  vQuantizer q8888(8, 8, 8, 8);
  vQuantizer q6666(6, 6, 6, 6);
  vQuantizer q8880(8, 8, 8, 0);
  vQuantizer q7770(7, 7, 7, 0);
  vQuantizer q5550(5, 5, 5, 0);
  
  const struct {
    vQuantizer *qnt;
  } caseqnt[8] = {
    { &q7778 },//{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 }, // 7+0,8+0  7778
    { &q5556 },//{ 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 }, // 5+0,6+0  5556
    
    { &q8888 },//{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 }, // 7+1,7+1  8888
    { &q6666 },//{ 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 }, // 5+1,5+1  6666
    
    { &q8880 },//{ 2, 6, 0, 0,  7, 0, 1,  0,  2, 0 }, // 7+1,0+0  8880
    { &q7770 },//{ 2, 6, 0, 0,  6, 0, 0,  1,  3, 0 }, // 6+1,0+0  7770

    { &q5550 },//{ 3, 6, 0, 0,  5, 0, 0,  0,  2, 0 }, // 5+0,0+0  5550
    { &q5550 },//{ 3, 4, 0, 0,  4, 0, 1,  0,  3, 0 }, // 4+1,0+0  5550
  };
  
  /* we start with 1 set so we get some statistics about the color-
   * palette, based on that we decide if we need to search into higher
   * number of sets
   *
   * observations:
   * - if there is 1-2 color(s), we only need 1 set:
   *   the index-precision doesn't matter in that case and we choose
   *   the coding with the highest start/end-point precision
   *   nevertheless the 2 colors are not necessarily also start/end
   *   interpolated colors may achieve superior precision
   * - if there is 3-4 color(s), we only need 2 sets:
   *   the available partitions may not correspond exactly to the
   *   distribution of the 3-4 colors, so for maximum quality we need
   *   to do the whole search regardless (including 3 sets)
   *   if we've found a 2 set partition with 1-2 colors in each we can
   *   abort immediately
   *
   * rangefit searches for the best configuration (partition/swap/rotation)
   * optionally clusterfit makes the best of that partition
   */
  static const struct {
    int mode, mnum, casen, groupn;
  } caseorder[8] = {
    { kVariableCodingMode6, 5, 0, 0 },//{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },
    { kVariableCodingMode5, 4, 0, 0 },//{ 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
    
    { kVariableCodingMode7, 6, 1, 1 },//{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
    { kVariableCodingMode8, 7, 1, 0 },//{ 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 },  // alpha variant of mode 2/4
    
    { kVariableCodingMode4, 3, 2, 0 },//{ 2, 6, 0, 0,  7, 0, 1,  0,  2, 0 },  // non-alpha variant of mode 8
    { kVariableCodingMode2, 1, 2, 0 },//{ 2, 6, 0, 0,  6, 0, 0,  1,  3, 0 },  // non-alpha variant of mode 8

    { kVariableCodingMode3, 2, 2, 1 },//{ 3, 6, 0, 0,  5, 0, 0,  0,  2, 0 },
    { kVariableCodingMode1, 0, 2, 1 },//{ 3, 4, 0, 0,  4, 0, 1,  0,  3, 0 },
  };
  
#define MODECASE_MIN  0
#define MODECASE_MAX  2
  static const struct {
    int num[2], size;
  } casegroups[3] = {
    { {1,1}, 2 },// 1x2 for both transparent and opaque
    { {1,2}, 1 },// 2x1 for transparent, 1x1 for opaque
    { {2,0}, 2 },// 0x2 for transparent, 2x2 for opaque
  };
  
  static int caselimit[2] = {
    2,
    1
  };

  // use the same data-structure all the time
  PaletteSet bestpal[2];
  Vec4 bestblock[2];
  int bestqnt[2] = {-1,-1};
  int bestmde[2] = {-1,-1};
  int bestswp[2] = {-1,-1};
  int bestbit[2] = {-1,-1};
  int besttyp[2] = {-1,-1};

  Scr4 error[2]; 
  
  error[0] = Scr4(FLT_MAX);
  error[1] = Scr4(FLT_MAX);

  // number of cases to walk through
  int numm = flags &  ( kVariableCodingModes),
        sc = (numm == 0 ? MODECASE_MIN : caseorder[(numm >> 24) - 1].casen),
	ec = (numm == 0 ? MODECASE_MAX :                                sc);
             flags &= (~kVariableCodingModes);

  // cases: separate (2x), merged alpha (2x), and no alpha (4x)
  for (int mc = sc; mc <= ec; mc++) {
    // offset of the current case
    int go = mc * 2;

    // create the initial point set
    PaletteSet initial(rgba, mask, flags + caseorder[go].mode);
    
    // if we see we have transparent values, back up from trying to test non-alpha only modes
    // this will affect only successive trials, if an explicit mode is requested it's a NOP
    ec = caselimit[initial.IsTransparent()];
    
    // number of groups per case, modes per group (1x2, 2x1, 2x2)
    int ng = casegroups[mc].num[initial.IsTransparent()];
    int gm = casegroups[mc].size;
    
    int sg = 0;
    int eg = ng - 1;
    
    for (int mg = sg; mg <= eg; mg++) {
      // offset of the current group's start and end
      int sm = go + (mg * gm);
      int em = sm + (gm - 1);
      
      // a mode has a specific number of sets, and variable rotations and partitions
      int numr = PaletteFit::GetRotationBits (caseorder[sm].mnum);
      int nump = PaletteFit::GetPartitionBits(caseorder[sm].mnum);
      
      // search through partitions
      int sp =               0,
	  ep = (1 << nump) - 1;
      // search through rotations
      int sr =               0,
	  er = (1 << numr) - 1;
      
      // if we see we have no transparent values, don't try non-rotated palettes (alpha is constant for all)
      if (!initial.IsTransparent() && initial.IsSeperateAlpha())
	sr = 1;

      // choose rotation or partition, they're mutually exclusive
      int spr = (er ? sr : sp),
	  epr = (er ? er : ep);
      
      // signal if we do we have anything better this iteration of the search
      bool better[2] = {false,false};

      for (int pr = spr; pr <= epr; pr++) {
	// create the minimal point set
	PaletteSet palette(initial, mask, flags + caseorder[sm].mode, pr);

	// do a range fit (which uses single palette fit if appropriate)
	PaletteTypeFit fit(&palette, flags + caseorder[sm].mode);
	
	// exclude mode 1 from the upper partitions
	if ((em == 7) && (pr >= (1 << 4)))
	  em = em - 1;

	for (int m = sm; m <= em; m++) {
	  int mode = caseorder[m].mode;
	  int mnum = caseorder[m].mnum;
	  int mofs = (flags & kColourRangeFit ? 0 : m - sm);

	  // a mode has a specific number of sets, and variable rotations and partitions
	  int numx = PaletteFit::GetSelectionBits(mnum);
	  int numb = PaletteFit::GetSharedBits   (mnum);

	  // search through index-swaps
	  int sx =                             0,
	      ex =               (1 << numx) - 1;
	  // search through shared bits
#ifdef FEATURE_SHAREDBITS_TRIALS
	  int sb = (numb > 0   ?               0 : SBSKIP),
	      eb = (numb > 0   ?       numb      : SBSKIP);
#else
	  int sb = (numb > 0   ?          SBSKIP : SBSKIP),
	      eb = (numb > 0   ?          SBSKIP : SBSKIP);
#endif
	  
#if	defined(FEATURE_SHAREDBITS_TRIALS)
	  // if we see we have no transparent values, force all shared bits to 1, or non-opaque codebook-entries occur
	  // the all transparent case isn't so crucial, when we use IGNORE_ALPHA0 it's redundant to force 0 anyway
	  if (!initial.IsTransparent() && initial.IsMergedAlpha())
	    sb = eb;
	  // otherwise just use the most occurring bit (parity) for all other cases
	  // otherwise just use the most occurring bit (parity) for all non-alpha cases
	  else if (((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_ALPHAONLYOPAQUE)) ||
		   ((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_ALPHAONLY) && (mode != kVariableCodingMode7) && !initial.IsTransparent()) ||
		   ((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_LOWPRC) && (mode < kVariableCodingMode5) && (mode != kVariableCodingMode1)))
	    sb = eb = SBSKIP;
#endif

	  // TODO: swap & shared are mutual exclusive
	  fit.ChangeMode(mnum);
	  // search for the best swap
	  for (int x = sx; x <= ex; x++) {
	    fit.ChangeSwap(x);
	    // search for the best shared bit
	    for (int b = sb; b <= eb; b++) {
	      fit.ChangeShared(b);
	      
	      // update with old best error (reset IsBest)
	      fit.SetError(error[mofs]);

	      // we could code it lossless, no point in trying any further at all
	      fit.PaletteTypeFit::Compress(block, *caseqnt[m].qnt, mnum);
	      if (fit.IsBest()) {
		if (fit.Lossless())
		  return Scr4(0.0f);

#if   !defined(TRACK_STATISTICS) && !defined(VERIFY_QUANTIZER)
		if (PaletteTypeFit::IsClusterable(flags))
#endif
		{
		  bestqnt[mofs] = m,
		  bestmde[mofs] = mode,
		  bestpal[mofs] = palette,
		  bestswp[mofs] = x,
		  bestbit[mofs] = b,
		  besttyp[mofs] = 0,
		  better [mofs]  = true;

		  LoadUnaligned(bestblock[mofs], block);
		}
	      
		// update with new best error
		error[mofs] = fit.GetError();
	      }
	    }
	  }
	}

#if 0
	// if we see we have less colors than sets, 
	// then back up from trying to test with more sets
	// or even other partitions
	if (palette.GetOptimal())
	  // modes are ordered by decreasing precision
	  // and increasing number of sets, this guarantees
	  // that no other following mode can possibly have
	  // more precise end-point than the current one
	  return;
#endif
      }
      
      if (!PaletteTypeFit::IsClusterable(flags))
	continue;

      Scr4 tobeat = Min(error[0], error[1]);
      if (better[0] | better[1]) {
	// ambiguous result, choose the better
	if (better[0] & better[1])
	  StoreUnaligned(bestblock[error[1] < error[0]], block);

	for (int m = 0; m <= 1; m++) {
	  if (better[m]) {
	    // a mode has a specific number of sets, and variable rotations and partitions
	    int mode = bestmde[m];
	    int mnum = (mode >> 24) - 1;
	    int numi = PaletteFit::GetIndexBits(mnum);

	    // check if we can do a cascade with the cluster-fit (merged alpha 4 bit is the only exception)
	    bool cluster = (((numi >>  0) & 0xFF) <= CLUSTERINDICES)
			&& (((numi >> 16) & 0xFF) <= CLUSTERINDICES);

	    // check the compression type and compress palette of the chosen partition even better
	    if (cluster) {
	      int degree = (flags & kColourIterativeClusterFits);

	      // default to a cluster fit (could be iterative or not)
	      PaletteClusterFit fit(&bestpal[m], flags + mode);

	      // a mode has a specific number of sets, and variable rotations and partitions
	      int numx = PaletteFit::GetSelectionBits(mnum);
	      int numb = PaletteFit::GetSharedBits   (mnum);

	      // search through index-swaps
	      int sx =                             0,
		  ex =               (1 << numx) - 1;
	      // search through shared bits
#ifdef FEATURE_SHAREDBITS_TRIALS
	      int sb = (numb > 0   ?               0 : SBSKIP),
		  eb = (numb > 0   ?       numb      : SBSKIP);
#else
	      int sb = (numb > 0   ?          SBSKIP : SBSKIP),
		  eb = (numb > 0   ?          SBSKIP : SBSKIP);
#endif
	  
#if	defined(FEATURE_SHAREDBITS_TRIALS)
	      // if we see we have no transparent values, force all shared bits to 1, or non-opaque codebook-entries occur
	      // the all transparent case isn't so crucial, when we use IGNORE_ALPHA0 it's redundant to force 0 anyway
	      if (!bestpal[m].IsTransparent() && bestpal[m].IsMergedAlpha())
		sb = eb;
	      // otherwise just use the most occurring bit (parity) for all other cases
	      // otherwise just use the most occurring bit (parity) for all non-alpha cases
	      else if (((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_ALPHAONLYOPAQUE)) ||
			((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_ALPHAONLY) && (mode != kVariableCodingMode7) && !bestpal[m].IsTransparent()) ||
			((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_LOWPRC) && (mode < kVariableCodingMode5) && (mode != kVariableCodingMode1)))
		sb = eb = SBSKIP;
#endif

	      // we want the whole shebang, this takes looong!
	      if (degree < (kColourClusterFit * 15))
		sb = eb = bestbit[m];
	      if (degree < (kColourClusterFit * 14))
		sx = ex = bestswp[m];
    
	      // TODO: swap & shared are mutual exclusive
	      fit.ChangeMode(mnum);
	      // search for the best swap
	      for (int x = sx; x <= ex; x++) {
		fit.ChangeSwap(x);
		// search for the best shared bit
		for (int b = sb; b <= eb; b++) {
		  fit.ChangeShared(b);
	    
		  // update with old best error (reset IsBest)
		  fit.SetError(tobeat);

		  // we could code it lossless, no point in trying any further at all
		  fit.PaletteClusterFit::Compress(block, *caseqnt[bestqnt[m]].qnt, mnum);
		  if (fit.IsBest()) {
		    if (fit.Lossless())
		      return Scr4(0.0f);

		    if (cluster || 1)
		      besttyp[m] = 1;
      
		    // update with new best error
		    tobeat = fit.GetError();
		  }

#if defined(TRACK_STATISTICS)
		  gstat.btr_cluster[mnum][fit.IsBest() ? 1 : 0]++;
#endif
		}
	      }
	    }
	  }
	}
      }

      error[0] = error[1] = tobeat;
    }
  }

#if defined(TRACK_STATISTICS)
  gstat.win_partition[(bestmde >> 24) - 1][bestpal.GetPartition()]++;
  gstat.win_rotation [(bestmde >> 24) - 1][bestpal.GetRotation ()]++;
  gstat.win_swap     [(bestmde >> 24) - 1][bestpal.GetRotation ()][bestswp]++;

  gstat.win_mode     [(bestmde >> 24) - 1]++;
  gstat.win_cluster  [(bestmde >> 24) - 1][besttyp]++;
#endif
  
#if defined(VERIFY_ENCODER)
  DecompressColoursBtc7u((u8*)rgba, block);
#endif

  return error[0];
}

template<typename dtyp>
void CompressColourBtc6u(dtyp const* rgb, int mask, void* block, int flags)
{
  static const int modeorder[1][14] = {
    {
#define MODEORDER_EXPL	    0
#define MODEORDER_EXPL_MIN  0
#define MODEORDER_EXPL_MAX  13
      // order: mode (lo to hi)
      kVariableCodingMode1, //{ 2, 5,   0, { 5, 5, 5}, {5,5,5},  6, 3 },
      kVariableCodingMode2, //{ 2, 5,   0, { 1, 1, 1}, {6,6,6},  9, 3 },
      kVariableCodingMode3, //{ 2, 5,   0, { 6, 7, 7}, {5,4,4},  5, 3 },
      kVariableCodingMode4, //{ 2, 5,   0, { 7, 6, 7}, {4,5,4},  5, 3 },
      kVariableCodingMode5, //{ 2, 5,   0, { 7, 7, 6}, {4,4,5},  5, 3 },
      kVariableCodingMode6, //{ 2, 5,   0, { 4, 4, 4}, {5,5,5},  7, 3 },
      kVariableCodingMode7, //{ 2, 5,   0, { 2, 3, 3}, {6,5,5},  8, 3 },
      kVariableCodingMode8, //{ 2, 5,   0, { 3, 2, 3}, {5,6,5},  8, 3 },
      kVariableCodingMode9, //{ 2, 5,   0, { 3, 3, 2}, {5,5,5},  8, 3 },
      kVariableCodingMode10,//{ 2, 5,   6, { 0, 0, 0}, {0,0,0}, 10, 3 },
      kVariableCodingMode11,//{ 1, 0,  10, { 0, 0, 0}, {0,0,0},  6, 4 },
      kVariableCodingMode12,//{ 1, 0,   0, { 2, 2, 2}, {9,9,9},  5, 4 },
      kVariableCodingMode13,//{ 1, 0,   0, { 4, 4, 4}, {8,8,8},  4, 4 },
      kVariableCodingMode14,//{ 1, 0,   0, {12,12,12}, {4,4,4},  0, 4 } 
    }
  };

  int numm = flags &  ( kVariableCodingModes),
    sm = (numm == 0 ? MODEORDER_EXPL_MIN : (numm >> 24) - 1),
    em = (numm == 0 ? MODEORDER_EXPL_MAX :               sm),
    om = (numm == 0 ? MODEORDER_EXPL     :   MODEORDER_EXPL);
	     flags &= (~kVariableCodingModes);

#undef MODEORDER_EXPL
#undef MODEORDER_EXPL_MIN
#undef MODEORDER_EXPL_MAX

  // limits sets to 2 and choose the partition freely
  int lmts =  2;
  int lmtp = -1;

  // use the same data-structure all the time
  HDRSet bestpal;
  int bestmde = -1;
  int besttyp = -1;

  Scr3 error(FLT_MAX);

  for (int m = sm; m <= em; m++) {
    int mode = modeorder[om][m];
    int mnum = (mode >> 24) - 1;

    // a mode has a specific number of sets, and variable partitions
    int nums = HDRFit::GetNumSets      (mnum);
    int nump = HDRFit::GetPartitionBits(mnum);
    int numi = HDRFit::GetIndexBits    (mnum); numi = numi;

    // stop if set-limit reached
    if (nums > lmts)
      break;

    // lock on the perfect partition
    int sp = (lmtp == -1 ?               0 : lmtp),
	ep = (lmtp == -1 ? (1 << nump) - 1 : lmtp);

    int tb = HDRFit::GetTruncationBits(mnum);
    int db = HDRFit::GetDeltaBits(mnum);

    // create the initial point set and quantizer
    HDRSet initial(rgb, mask, flags + mode);
    fQuantizer qnt(tb, db);

    // signal if we do we have anything better this iteration of the search
    bool better = false;

    // search for the best partition
    for (int p = sp; p <= ep; p++) {
      // create the minimal point set
      HDRSet palette(initial, mask, flags + mode, p);

#if 0
			// if we see we have less colors than sets, 
			// then back up from trying to test with more sets
			// or even other partitions
      if (palette.GetCount() <= nums)
	lmtp = p, lmts = nums;
#endif

      // do a range fit (which uses single palette fit if appropriate)
      HDRRangeFit fit(&palette, flags + mode);

      // update with old best error (reset IsBest)
      fit.SetError(error);
      fit.Compress(block, qnt, mnum);

      // we could code it lossless, no point in trying any further at all
      if (fit.IsBest()) {
	if (fit.Lossless())
	  return;

	error = fit.GetError();
	if (1)
	  bestmde = mode,
	  bestpal = palette,
	  besttyp = 0,
	  better  = true;
      }
    }
  }

#if defined(VERIFY_ENCODER)
  DecompressHDRsBtc6u((f23*)rgb, block);
#endif
}

/* *****************************************************************************
 */
template<typename dtyp>
void CompressMaskedBitoneCtx1u(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = block;

  // compress color separately if necessary
  CompressBitoneCtx1u(rgba, mask, colourBlock, flags);
}

template<typename dtyp>
void CompressMaskedColourBtc1u(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = block;

  // compress color separately if necessary
  CompressColourBtc1u(rgba, mask, colourBlock, flags);
}

template<typename dtyp>
void CompressMaskedColourBtc2u(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressColourBtc1u(rgba, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc2u(rgba, mask, alphaBlock);
}

template<typename dtyp>
void CompressMaskedColourBtc3u(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressColourBtc1u(rgba, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc3u(rgba, mask, alphaBlock, flags);
}

template<typename dtyp>
void CompressMaskedNormalCtx1u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* normalBlock = block;

  // compress color separately if necessary
  CompressNormalCtx1u(xyzd, mask, normalBlock, flags);
}

template<typename dtyp>
void CompressMaskedNormalBtc1u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* normalBlock = block;

  // compress color separately if necessary
  CompressNormalBtc1u(xyzd, mask, normalBlock, flags);
}

template<typename dtyp>
void CompressMaskedNormalBtc2u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressNormalBtc1u(xyzd, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc2u(xyzd, mask, alphaBlock);
}

template<typename dtyp>
void CompressMaskedNormalBtc3u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressNormalBtc1u(xyzd, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc3u(xyzd, mask, alphaBlock, flags);
}

template<typename dtyp>
void CompressMaskedAlphaBtc4u(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = block;

  // compress a into plane 1
  CompressDepthBtc4u(rgba - 3, mask, plane1Block, flags);
}

template<typename dtyp>
void CompressMaskedAlphaBtc4s(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = block;

  // compress a into plane 1
  CompressDepthBtc4s(rgba - 3, mask, plane1Block, flags);
}

template<typename dtyp>
void CompressMaskedAlphaBtc5u(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = reinterpret_cast<u8*>(block) + 8;
  void* plane2Block = block;

  // compress a into plane 1
  CompressDepthBtc4u(rgba - 3, mask, plane1Block, flags);
  // compress b into plane 2
  CompressDepthBtc4u(rgba - 2, mask, plane2Block, flags);
}

template<typename dtyp>
void CompressMaskedAlphaBtc5s(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = reinterpret_cast<u8*>(block) + 8;
  void* plane2Block = block;

  // compress a into plane 1
  CompressDepthBtc4s(rgba - 3, mask, plane1Block, flags);
  // compress b into plane 2
  CompressDepthBtc4s(rgba - 2, mask, plane2Block, flags);
}

template<typename dtyp>
void CompressMaskedNormalBtc5u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = reinterpret_cast<u8*>(block) + 8;
  void* plane2Block = block;

  // compress xy into plane 1/2
  CompressNormalsBtc5u(xyzd, mask, plane1Block, plane2Block, flags);
}

template<typename dtyp>
void CompressMaskedNormalBtc5s(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = reinterpret_cast<u8*>(block) + 8;
  void* plane2Block = block;

  // compress xy into plane 1/2
  CompressNormalsBtc5s(xyzd, mask, plane1Block, plane2Block, flags);
}

template<typename dtyp>
void CompressMaskedColourBtc6u(dtyp const* rgb, int mask, void* block, int flags)
{
  // get the block locations
  void* mixedBlock = block;

  // compress color and alpha merged if necessary
  CompressColourBtc6u(rgb, mask, mixedBlock, flags);
}

template<typename dtyp>
void CompressMaskedColourBtc7u(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* mixedBlock = block;
  
#ifdef DEBUG_DETAILS
  // compress color and alpha merged if necessary
  fprintf(stderr, "CompressPaletteBtc7uV1\n");
  Scr4 errora = CompressPaletteBtc7uV1<dtyp,PaletteRangeFit>(rgba, mask, mixedBlock, flags);
  fprintf(stderr, "CompressPaletteBtc7uV2\n");
  Scr4 errorb = CompressPaletteBtc7uV2<dtyp,PaletteRangeFit>(rgba, mask, mixedBlock, flags);

  if (errorb > errora) {
    bool damn = true; damn = false;

    fprintf(stderr, "CompressPaletteBtc7uV2\n");
    errorb = CompressPaletteBtc7uV2<dtyp,PaletteRangeFit>(rgba, mask, mixedBlock, flags);
  }
  else if (errorb < errora) {
    bool cool = true; cool = false;
  }
#else
  CompressPaletteBtc7uV2<dtyp,PaletteRangeFit>(rgba, mask, mixedBlock, flags);
#endif
}

template<typename dtyp>
void CompressMaskedNormalBtc7u(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* mixedBlock = block;
  
  CompressPaletteBtc7uV2<dtyp,PaletteNormalFit>(xyzd, mask, mixedBlock, flags);
}

void CompressMasked(u8 const* rgba, int mask, void* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
    CompressMaskedNormalCtx1u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
    CompressMaskedNormalBtc1u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
    CompressMaskedNormalBtc2u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
    CompressMaskedNormalBtc3u(rgba, mask, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
    CompressMaskedNormalBtc5u(rgba, mask, block, flags);
  // BTC-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
    CompressMaskedNormalBtc7u(rgba, mask, block, flags);

  // DXT-type compression
  else if ((flags & kBtcp) == (kBtc1))
    CompressMaskedColourBtc1u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc2))
    CompressMaskedColourBtc2u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc3))
    CompressMaskedColourBtc3u(rgba, mask, block, flags);
  // ATI-type compression
  else if ((flags & kBtcp) == (kBtc4))
    CompressMaskedAlphaBtc4u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc5))
    CompressMaskedAlphaBtc5u(rgba, mask, block, flags);
  // BTC-type compression
  else if ((flags & kBtcp) == (kBtc7))
    CompressMaskedColourBtc7u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc6))
    {}// while this is possible (up-cast), should we support it?
}

void CompressMasked(u16 const* rgba, int mask, void* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
    CompressMaskedNormalCtx1u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
    CompressMaskedNormalBtc1u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
    CompressMaskedNormalBtc2u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
    CompressMaskedNormalBtc3u(rgba, mask, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
    CompressMaskedNormalBtc5u(rgba, mask, block, flags);
  // BTC-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
    CompressMaskedNormalBtc7u(rgba, mask, block, flags);

  // DXT-type compression
  else if ((flags & kBtcp) == (kBtc1))
    CompressMaskedColourBtc1u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc2))
    CompressMaskedColourBtc2u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc3))
    CompressMaskedColourBtc3u(rgba, mask, block, flags);
  // ATI-type compression
  else if ((flags & kBtcp) == (kBtc4))
    CompressMaskedAlphaBtc4u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc5))
    CompressMaskedAlphaBtc5u(rgba, mask, block, flags);
  // BTC-type compression
  else if ((flags & kBtcp) == (kBtc7))
    CompressMaskedColourBtc7u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc6))
    CompressMaskedColourBtc6u(rgba, mask, block, flags);
}

void CompressMasked(f23 const* rgba, int mask, void* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
    CompressMaskedNormalCtx1u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
    CompressMaskedNormalBtc1u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
    CompressMaskedNormalBtc2u(rgba, mask, block, flags);
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
    CompressMaskedNormalBtc1u(rgba, mask, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
    CompressMaskedNormalBtc5u(rgba, mask, block, flags);
  // BTC-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
    CompressMaskedNormalBtc7u(rgba, mask, block, flags);

  // DXT-type compression
  else if ((flags & kBtcp) == (kBtc1))
    CompressMaskedColourBtc1u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc2))
    CompressMaskedColourBtc2u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc3))
    CompressMaskedColourBtc3u(rgba, mask, block, flags);
  // ATI-type compression
  else if ((flags & kBtcp) == (kBtc4))
    CompressMaskedAlphaBtc4u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc5))
    CompressMaskedAlphaBtc5u(rgba, mask, block, flags);
  // BTC-type compression
  else if ((flags & kBtcp) == (kBtc7))
    CompressMaskedColourBtc7u(rgba, mask, block, flags);
  else if ((flags & kBtcp) == (kBtc6))
    CompressMaskedColourBtc6u(rgba, mask, block, flags);
}

void Compress(u8 const* rgba, void* block, int flags)
{
  // compress with full mask
  CompressMasked(rgba, -1, block, flags);
}

void Compress(u16 const* rgb, void* block, int flags)
{
  // compress with full mask
  CompressMasked(rgb, -1, block, flags);
}

void Compress(f23 const* rgba, void* block, int flags)
{
  // compress with full mask
  CompressMasked(rgba, -1, block, flags);
}

/* *****************************************************************************
 */
template<typename dtyp>
void DecompressNormalCtx1u(dtyp* xyzd, void const* block, int flags)
{
  // get the block locations
  void const* normalBlock = block;

  // decompress normals
  DecompressNormalsCtx1u(xyzd, normalBlock);
}

template<typename dtyp>
void DecompressBitoneCtx1u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* colourBlock = block;

  // decompress colour
  DecompressBitonesCtx1u(rgba, colourBlock);
}

template<typename dtyp>
void DecompressColourBtc1u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* colourBlock = block;

  // decompress colour
  DecompressColoursBtc1u(rgba, colourBlock, true);
}

template<typename dtyp>
void DecompressColourBtc2u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* colourBlock = reinterpret_cast<u8 const* >(block) + 8;
  void const*  alphaBlock = block;

  // decompress colour
  DecompressColoursBtc1u(rgba, colourBlock, false);
  // decompress alpha separately if necessary
  DecompressAlphaBtc2u(rgba, alphaBlock);
}

template<typename dtyp>
void DecompressColourBtc3u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* colourBlock = reinterpret_cast<u8 const* >(block) + 8;
  void const*  alphaBlock = block;

  // decompress colour
  DecompressColoursBtc1u(rgba, colourBlock, false);
  // decompress alpha separately if necessary
  DecompressAlphaBtc3u(rgba, alphaBlock, flags);
}

template<typename dtyp>
void DecompressAlphaBtc4u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = block;

  // decompress plane 1 into a
  DecompressDepthBtc4u(rgba - 3, plane1Block, flags);
}

template<typename dtyp>
void DecompressAlphaBtc4s(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = block;

  // decompress plane 1 into a
  DecompressDepthBtc4s(rgba - 3, plane1Block, flags);
}

template<typename dtyp>
void DecompressAlphaBtc5u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = reinterpret_cast<u8 const* >(block) + 8;
  void const* plane2Block = block;

  // decompress plane 1 into a
  DecompressDepthBtc4u(rgba - 3, plane1Block, flags);
  // decompress plane 2 into b
  DecompressDepthBtc4u(rgba - 2, plane2Block, flags);
}

template<typename dtyp>
void DecompressAlphaBtc5s(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = reinterpret_cast<u8 const* >(block) + 8;
  void const* plane2Block = block;

  // decompress plane 1 into a
  DecompressDepthBtc4s(rgba - 3, plane1Block, flags);
  // decompress plane 2 into b
  DecompressDepthBtc4s(rgba - 2, plane2Block, flags);
}

template<typename dtyp>
void DecompressNormalBtc5u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = reinterpret_cast<u8 const* >(block) + 8;
  void const* plane2Block = block;

  // compress xy into plane 1/2
  DecompressNormalsBtc5u(rgba, plane1Block, plane2Block);
}

template<typename dtyp>
void DecompressNormalBtc5s(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = reinterpret_cast<u8 const* >(block) + 8;
  void const* plane2Block = block;

  // compress xy into plane 1/2
  DecompressNormalsBtc5s(rgba, plane1Block, plane2Block);
}

template<typename dtyp>
void DecompressColourBtc6u(dtyp* rgb, void const* block, int flags)
{
  // get the block locations
  void const* mixedBlock = block;

  // decompress color and alpha merged if necessary
	DecompressHDRsBtc6u(rgb, mixedBlock);
}

template<typename dtyp>
void DecompressColourBtc7u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* mixedBlock = block;

  // decompress color and alpha merged if necessary
  DecompressColoursBtc7u(rgba, mixedBlock);
}

template<typename dtyp>
void DecompressNormalBtc7u(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* mixedBlock = block;
  
  // decompress normals
  DecompressNormalsBtc7u(rgba, mixedBlock);
}

void Decompress(u8* rgba, void const* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
    DecompressNormalCtx1u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
//  DecompressNormalBtc1u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
//  DecompressNormalBtc2u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
//  DecompressNormalBtc3u(rgba, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
    DecompressNormalBtc5u(rgba, block, flags);
  // BTC-type compression
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
//  DecompressNormalBtc7u(rgba, block, flags);

  // DXT-type compression
  else if ((flags & kBtcp) == (kBtc1))
    DecompressColourBtc1u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc2))
    DecompressColourBtc2u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc3))
    DecompressColourBtc3u(rgba, block, flags);
  // ATI-type compression
  else if ((flags & kBtcp) == (kBtc4))
    DecompressAlphaBtc4u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc5))
    DecompressAlphaBtc5u(rgba, block, flags);
  // BTC-type compression
  else if ((flags & kBtcp) == (kBtc7))
    DecompressColourBtc7u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc6))
    {}// while this is possible (down-cast), should we support it?
}

void Decompress(u16* rgba, void const* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
    DecompressNormalCtx1u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
//  DecompressNormalBtc1u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
//  DecompressNormalBtc2u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
//  DecompressNormalBtc3u(rgba, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
    DecompressNormalBtc5u(rgba, block, flags);
  // BTC-type compression
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
//  DecompressNormalBtc7u(rgba, block, flags);

  // DXT-type compression
  else if ((flags & kBtcp) == (kBtc1))
    DecompressColourBtc1u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc2))
    DecompressColourBtc2u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc3))
    DecompressColourBtc3u(rgba, block, flags);
  // ATI-type compression
  else if ((flags & kBtcp) == (kBtc4))
    DecompressAlphaBtc4u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc5))
    DecompressAlphaBtc5u(rgba, block, flags);
  // BTC-type compression
  else if ((flags & kBtcp) == (kBtc7))
    DecompressColourBtc7u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc6))
    DecompressColourBtc6u(rgba, block, flags);
}

void Decompress(f23* rgba, void const* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
    DecompressNormalCtx1u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
//  DecompressNormalBtc1u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
//  DecompressNormalBtc2u(rgba, block, flags);
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
//  DecompressNormalBtc3u(rgba, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
    DecompressNormalBtc5u(rgba, block, flags);
  // BTC-type compression
//else if ((flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
//  DecompressNormalBtc7u(rgba, block, flags);

  // DXT-type compression
  else if ((flags & kBtcp) == (kBtc1))
    DecompressColourBtc1u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc2))
    DecompressColourBtc2u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc3))
    DecompressColourBtc3u(rgba, block, flags);
  // ATI-type compression
  else if ((flags & kBtcp) == (kBtc4))
    DecompressAlphaBtc4u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc5))
    DecompressAlphaBtc5u(rgba, block, flags);
  // BTC-type compression
  else if ((flags & kBtcp) == (kBtc7))
    DecompressColourBtc7u(rgba, block, flags);
  else if ((flags & kBtcp) == (kBtc6))
    DecompressColourBtc6u(rgba, block, flags);
}

/* *****************************************************************************
 */
int GetStorageRequirements(int width, int height, int flags)
{
  // compute the storage requirements
  int blockcount = ((width + 3) / 4) * ((height + 3) / 4);
  int blocksize  = 16;

  /**/ if ((flags & kBtcp) <= kBtc3)
    blocksize = ((flags & kBtcp) == kBtc1) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc5)
    blocksize = ((flags & kBtcp) == kBtc4) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc7)
    blocksize =                                  16;
  else if ((flags & kBtcp) == kCtx1)
    blocksize =                              8     ;

  return blockcount * blocksize;
}

/* *****************************************************************************
 */
struct sqio GetSquishIO(int width, int height, sqio::dtp datatype, int flags)
{
  struct sqio s;
  
  s.datatype = datatype;
  s.flags = SanitizeFlags(flags);

  // compute the storage requirements
  s.blockcount = ((width + 3) / 4) * ((height + 3) / 4);
  s.blocksize  = 16;
  
  /**/ if ((flags & kBtcp) <= kBtc3)
    s.blocksize = ((flags & kBtcp) == kBtc1) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc5)
    s.blocksize = ((flags & kBtcp) == kBtc4) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc7)
    s.blocksize =                                  16;
  else if ((flags & kBtcp) == kCtx1)
    s.blocksize =                              8     ;

  s.compressedsize = s.blockcount * s.blocksize;
  s.decompressedsize = sizeof(u8) * 4 * (width * height);

  if (datatype == sqio::DT_U8) {
    // 3Dc-type compression
    /**/ if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc5 | kSignedExternal | kSignedInternal | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5s<s8>,
      s.decoder = (sqio::dec)DecompressNormalBtc5s<s8>;
    // ATI-type compression
    else if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc4 | kSignedExternal | kSignedInternal | kColourMetricUniform))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4s<s8>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4s<s8>;
    else if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc5 | kSignedExternal | kSignedInternal | kColourMetricUniform))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5s<s8>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5s<s8>;

    // DXT-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalCtx1u<u8>,
      s.decoder = (sqio::dec)DecompressNormalCtx1u<u8>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc1u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc1u<u8>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc2u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc2u<u8>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc3u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc3u<u8>;
    // 3Dc-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5u<u8>,
      s.decoder = (sqio::dec)DecompressNormalBtc5u<u8>;
    // BTC-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc7u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc7u<u8>;
    
    // DXT-type compression
    else if ((s.flags & kBtcp) == (kCtx1))
      s.encoder = (sqio::enc)CompressMaskedBitoneCtx1u<u8>,
      s.decoder = (sqio::dec)DecompressBitoneCtx1u<u8>;
    else if ((s.flags & kBtcp) == (kBtc1))
      s.encoder = (sqio::enc)CompressMaskedColourBtc1u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc1u<u8>;
    else if ((s.flags & kBtcp) == (kBtc2))
      s.encoder = (sqio::enc)CompressMaskedColourBtc2u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc2u<u8>;
    else if ((s.flags & kBtcp) == (kBtc3))
      s.encoder = (sqio::enc)CompressMaskedColourBtc3u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc3u<u8>;
    // ATI-type compression
    else if ((s.flags & kBtcp) == (kBtc4))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4u<u8>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4u<u8>;
    else if ((s.flags & kBtcp) == (kBtc5))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5u<u8>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5u<u8>;
    // BTC-type compression
    else if ((s.flags & kBtcp) == (kBtc7))
      s.encoder = (sqio::enc)CompressMaskedColourBtc7u<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc7u<u8>;
    else if ((s.flags & kBtcp) == (kBtc6))
      {}// while this is possible (down-cast), should we support it?
  }
  else if (datatype == sqio::DT_U16) {
    // 3Dc-type compression
    /**/ if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc5 | kSignedExternal | kSignedInternal | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5s<s16>,
      s.decoder = (sqio::dec)DecompressNormalBtc5s<s16>;
    // ATI-type compression
    else if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc4 | kSignedExternal | kSignedInternal | kColourMetricUniform))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4s<s16>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4s<s16>;
    else if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc5 | kSignedExternal | kSignedInternal | kColourMetricUniform))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5s<s16>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5s<s16>;

    // DXT-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalCtx1u<u16>,
      s.decoder = (sqio::dec)DecompressNormalCtx1u<u16>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc1u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc1u<u16>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc2u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc2u<u16>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc3u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc3u<u16>;
    // 3Dc-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5u<u16>,
      s.decoder = (sqio::dec)DecompressNormalBtc5u<u16>;
    // BTC-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc7u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc7u<u16>;

    // DXT-type compression
    else if ((s.flags & kBtcp) == (kCtx1))
      s.encoder = (sqio::enc)CompressMaskedBitoneCtx1u<u16>,
      s.decoder = (sqio::dec)DecompressBitoneCtx1u<u16>;
    else if ((s.flags & kBtcp) == (kBtc1))
      s.encoder = (sqio::enc)CompressMaskedColourBtc1u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc1u<u16>;
    else if ((s.flags & kBtcp) == (kBtc2))
      s.encoder = (sqio::enc)CompressMaskedColourBtc2u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc2u<u16>;
    else if ((s.flags & kBtcp) == (kBtc3))
      s.encoder = (sqio::enc)CompressMaskedColourBtc3u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc3u<u16>;
    // ATI-type compression
    else if ((s.flags & kBtcp) == (kBtc4))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4u<u16>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4u<u16>;
    else if ((s.flags & kBtcp) == (kBtc5))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5u<u16>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5u<u16>;
    // BTC-type compression
    else if ((s.flags & kBtcp) == (kBtc7))
      s.encoder = (sqio::enc)CompressMaskedColourBtc7u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc7u<u16>;
    else if ((s.flags & kBtcp) == (kBtc6))
      s.encoder = (sqio::enc)CompressMaskedColourBtc6u<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc6u<u16>;
  }
  else if (datatype == sqio::DT_F23) {
    // 3Dc-type compression
    /**/ if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc5 | kSignedExternal | kSignedInternal | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5s<f23>,
      s.decoder = (sqio::dec)DecompressNormalBtc5s<f23>;
    // ATI-type compression
    else if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc4 | kSignedExternal | kSignedInternal | kColourMetricUniform))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4s<f23>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4s<f23>;
    else if ((s.flags & (kBtcp | kSignedness | kColourMetrics)) == (kBtc5 | kSignedExternal | kSignedInternal | kColourMetricUniform))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5s<f23>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5s<f23>;

    // DXT-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kCtx1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalCtx1u<f23>,
      s.decoder = (sqio::dec)DecompressNormalCtx1u<f23>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc1u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc1u<f23>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc2 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc2u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc2u<f23>;
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc3 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc3u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc3u<f23>;
    // 3Dc-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc5 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5u<f23>,
      s.decoder = (sqio::dec)DecompressNormalBtc5u<f23>;
    // BTC-type compression
    else if ((s.flags & (kBtcp | kColourMetrics)) == (kBtc7 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc7u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc7u<f23>;

    // DXT-type compression
    else if ((s.flags & kBtcp) == (kCtx1))
      s.encoder = (sqio::enc)CompressMaskedBitoneCtx1u<f23>,
      s.decoder = (sqio::dec)DecompressBitoneCtx1u<f23>;
    else if ((s.flags & kBtcp) == (kBtc1))
      s.encoder = (sqio::enc)CompressMaskedColourBtc1u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc1u<f23>;
    else if ((s.flags & kBtcp) == (kBtc2))
      s.encoder = (sqio::enc)CompressMaskedColourBtc2u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc2u<f23>;
    else if ((s.flags & kBtcp) == (kBtc3))
      s.encoder = (sqio::enc)CompressMaskedColourBtc3u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc3u<f23>;
    // ATI-type compression
    else if ((s.flags & kBtcp) == (kBtc4))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4u<f23>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4u<f23>;
    else if ((s.flags & kBtcp) == (kBtc5))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5u<f23>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5u<f23>;
    // BTC-type compression
    else if ((s.flags & kBtcp) == (kBtc7))
      s.encoder = (sqio::enc)CompressMaskedColourBtc7u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc7u<f23>;
    else if ((s.flags & kBtcp) == (kBtc6))
      s.encoder = (sqio::enc)CompressMaskedColourBtc6u<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc6u<f23>;
  }

  return s;
}

/* *****************************************************************************
 */
void CompressImage(u8 const* rgba, int width, int height, void* blocks, int flags)
{
  // fix any bad flags
  flags = SanitizeFlags(flags);

  // initialize the block output
  unsigned char* targetBlock = reinterpret_cast< unsigned char* >(blocks);
  int bytesPerBlock = 16;
  
  /**/ if ((flags & kBtcp) <= kBtc3)
    bytesPerBlock = ((flags & kBtcp) == kBtc1) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc5)
    bytesPerBlock = ((flags & kBtcp) == kBtc4) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc7)
    bytesPerBlock =                                  16;
  else if ((flags & kBtcp) == kCtx1)
    bytesPerBlock =                              8     ;

  // loop over blocks
  for (int y = 0; y < height; y += 4) {
    for (int x = 0; x < width; x += 4) {
      // build the 4x4 block of pixels
      u8 sourceRgba[16 * 4];
      u8* targetPixel = sourceRgba;

      int mask = 0;
      for (int py = 0; py < 4; ++py) {
	for (int px = 0; px < 4; ++px) {
	  // get the source pixel in the image
	  int sx = x + px;
	  int sy = y + py;

	  // enable if we're in the image
	  if (sx < width && sy < height) {
	    // copy the rgba value
	    u8 const* sourcePixel = rgba + 4 * (width * sy + sx);
	    for (int i = 0; i < 4; ++i)
	      *targetPixel++ = *sourcePixel++;

	    // enable this pixel
	    mask |= (1 << (4 * py + px));
	  }
	  else {
	    // skip this pixel as its outside the image
	    targetPixel += 4;
	  }
	}
      }

      // compress it into the output
      CompressMasked(sourceRgba, mask, targetBlock, flags);

      // advance
      targetBlock += bytesPerBlock;
    }
  }
}

void DecompressImage(u8* rgba, int width, int height, void const* blocks, int flags)
{
  // fix any bad flags
  flags = SanitizeFlags(flags);

  // initialize the block input
  unsigned char const* sourceBlock = reinterpret_cast< unsigned char const* >(blocks);
  int bytesPerBlock = 16;
  
  /**/ if ((flags & kBtcp) <= kBtc3)
    bytesPerBlock = ((flags & kBtcp) == kBtc1) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc5)
    bytesPerBlock = ((flags & kBtcp) == kBtc4) ? 8 : 16;
  else if ((flags & kBtcp) <= kBtc7)
    bytesPerBlock =                                  16;
  else if ((flags & kBtcp) == kCtx1)
    bytesPerBlock =                              8     ;

  // loop over blocks
  for (int y = 0; y < height; y += 4) {
    for (int x = 0; x < width; x += 4) {
      // decompress the block
      u8 targetRgba[4 * 16];

      Decompress(targetRgba, sourceBlock, flags);

      // write the decompressed pixels to the correct image locations
      u8 const* sourcePixel = targetRgba;
      for (int py = 0; py < 4; ++py) {
	for (int px = 0; px < 4; ++px) {
	  // get the target location
	  int sx = x + px;
	  int sy = y + py;

	  if (sx < width && sy < height) {
	    u8* targetPixel = rgba + 4 * (width * sy + sx);

	    // copy the rgba value
	    for (int i = 0; i < 4; ++i)
	      *targetPixel++ = *sourcePixel++;
	  }
	  else {
	    // skip this pixel as its outside the image
	    sourcePixel += 4;
	  }
	}
      }

      // advance
      sourceBlock += bytesPerBlock;
    }
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#if	defined(SQUISH_USE_COMPUTE)
    tile_static ColourSingleFit_CCR sfit;
    tile_static ColourRangeFit_CCR rfit;
    tile_static ClusterFit_CCR cfit;
#endif

void CompressColorBtc ( tile_barrier barrier, const int thread,
			pixel16 rgba, ColourSet_CCRr colours, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted {
  // all of these conditions are identical over the entire
  // thread-group, so all threads take one of 0,1,2, not
  // distinct branches

  // check the compression type and compress colour
  if (colours.GetCount() == 1) {
    // always do a single colour fit
#if	!defined(SQUISH_USE_COMPUTE)
    tile_static ColourSingleFit_CCR sfit;
#endif

    sfit.AssignSet(barrier, thread, colours, metric, fit);
    sfit.Compress (barrier, thread, colours, block, trans, yArr, lArr);
  }
  else if ((fit == SQUISH_FIT_RANGE) || (colours.GetCount() == 0)) {
    // do a range fit
#if	!defined(SQUISH_USE_COMPUTE)
    tile_static ColourRangeFit_CCR rfit;
#endif

    rfit.AssignSet(barrier, thread, colours, metric, fit);
    rfit.Compress (barrier, thread, colours, block, trans, yArr);
  }
  else {
    // default to a cluster fit (could be iterative or not)
#if	!defined(SQUISH_USE_COMPUTE)
    tile_static ClusterFit_CCR cfit;
#endif

    cfit.AssignSet(barrier, thread, colours, metric, fit);
    cfit.Compress (barrier, thread, colours, block, trans, yArr);
  }
}

#if	defined(SQUISH_USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

void CompressColourBtc1u( tile_barrier barrier, const int thread,
			pixel16 rgba, int mask, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted {
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

  // create the minimal point set
  colours.CountSet(barrier, thread, rgba, mask, true, trans);

  // compress colour
  CompressColorBtc(barrier, thread, rgba, colours, block, metric, trans, fit, yArr, lArr);
}

void CompressColorBtc2( tile_barrier barrier, const int thread,
			pixel16 rgba, int mask, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted {
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

  // create the minimal point set
  colours.CountSet(barrier, thread, rgba, mask, false, trans);

  // compress colour
  CompressColorBtc(barrier, thread, rgba, colours, block, metric, trans, fit, yArr, lArr);
}

void CompressColorBtc3( tile_barrier barrier, const int thread,
			pixel16 rgba, int mask, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted {
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

  // create the minimal point set
  colours.CountSet(barrier, thread, rgba, mask, false, trans);

  // compress colour
  CompressColorBtc(barrier, thread, rgba, colours, block, metric, trans, fit, yArr, lArr);
}
#endif

} // namespace squish
