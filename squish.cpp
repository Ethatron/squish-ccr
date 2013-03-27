/* -----------------------------------------------------------------------------

	Copyright (c) 2006 Simon Brown                          si@sjbrown.co.uk
	Copyright (c) 2012 Niels Fr�hling              niels@paradice-insight.us

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

#include "colourset.h"
#include "paletteset.h"
//nclude "hdrset.h"

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
//nclude "hdrrangefit.h"
//nclude "hdrclusterfit.h"
#include "hdrblock.h"

#include "coloursinglefit.h"
#include "coloursinglesnap.h"
#include "palettesinglefit.h"
#include "palettesinglesnap.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
int SanitizeFlags(int flags)
{
  // grab the flag bits
  int method = flags & (kBtc1 | kBtc2 | kBtc3 | kBtc4 | kBtc5 | kBtc6 | kBtc7);
  int fit    = flags & (kColourRangeFit | kAlphaIterativeFit | kColourIterativeClusterFits);
  int metric = flags & (kColourMetricUniform | kColourMetricPerceptual | kColourMetricUnit);
  int extra  = flags & (kWeightColourByAlpha);
  int mode   = flags & (kVariableCodingModes);

  // set defaults
  if (!method)
    method = kBtc1;

  if (!fit)
    fit = (kColourClusterFit * 1);
  if (fit & kColourIterativeClusterFits) {
    if (method & (kBtc1 | kBtc2 | kBtc3))
      fit = ColourClusterFit::SanitizeFlags(fit);
    if (method & (kBtc7))
      fit = PaletteClusterFit::SanitizeFlags(fit);
  }

  if (!metric)
    metric = kColourMetricPerceptual;

  if ((method & kBtc6) && (mode > kVariableCodingMode14))
    mode = 0;
  if ((method & kBtc7) && (mode > kVariableCodingMode8))
    mode = 0;

  // done
  return method + fit + metric + extra + mode;
}

/* *****************************************************************************
 */
template<typename dtyp>
void CompressColorBtc1(dtyp const* rgba, int mask, void* block, int flags)
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
void CompressNormalBtc1(dtyp const* xyzd, int mask, void* block, int flags)
{
  // create the minimal point set
  ColourSet colours(xyzd, mask, flags);

  // check the compression type and compress colour
  if (colours.GetCount() == 1) {
    // always do a single colour fit
    ColourSingleMatch fit(&colours, flags);
    fit.Compress(block);
  }
  else {
    // do a range fit
    ColourNormalFit fit(&colours, flags);
    fit.Compress(block);
  }
}

#if defined(TRACK_STATISTICS)
struct statistics gstat = {0};
#endif

template<typename dtyp>
void CompressPaletteBtc7(dtyp const* rgba, int mask, void* block, int flags)
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
    bool cluster = ((flags & kColourRangeFit) == 0) && (((numi >>  0) & 0xFF) <= CLUSTERINDICES)
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

    // TODO: partition & rotation are mutual exclusive

    // search for the best partition
    for (int p = sp; p <= ep; p++) {
      // search for the best rotation
      for (int r = sr; r <= er; r++) {
	// create the minimal point set
	PaletteSet palette(initial, mask, flags + mode, p, r);

	// if we see we have less colors than sets, 
	// then back up from trying to test with more sets
	// or even other partitions
	if (palette.GetCount() <= nums)
	  lmtp = p, lmts = nums;

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
	PaletteRangeFit fit(&palette, flags + mode);

	// TODO: swap & shared are mutual exclusive

	// search for the best swap
	for (int x = sx; x <= ex; x++) {
	  fit.ChangeSwap(x);
	  // search for the best shared bit
	  for (int b = sb; b <= eb; b++) {
	    fit.ChangeShared(b);

	    // update with old best error (reset IsBest)
	    fit.SetError(error);
	    fit.Compress(block, qnt, mnum);

	    // we could code it lossless, no point in trying any further at all
	    if (fit.IsBest()) {
	      if (fit.Lossless())
		return;

	      error = fit.GetError();
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
	  }
	}
      }
    }

    // check the compression type and compress palette of the chosen partition even better
    if (better && cluster) {
      int degree = flags & (kColourClusterFit * 15);

      // default to a cluster fit (could be iterative or not)
      PaletteClusterFit fit(&bestpal, flags + mode);

      // we want the whole shebang, this takes looong!
      if (degree < (kColourClusterFit * 15))
	sb = eb = bestbit;
      if (degree < (kColourClusterFit * 14))
	sb = eb = bestbit, sx = ex = bestswp;

      // search for the best swap
      for (int x = sx; x <= ex; x++) {
	fit.ChangeSwap(x);
	// search for the best shared bit
	for (int b = sb; b <= eb; b++) {
	  fit.ChangeShared(b);

	  // update with old best error (reset IsBest)
	  fit.SetError(error);
	  fit.Compress(block, qnt, mnum);

#if defined(TRACK_STATISTICS)
	  gstat.btr_cluster[mnum][fit.IsBest() ? 1 : 0]++;
#endif

	  // we could code it lossless, no point in trying any further at all
	  if (fit.IsBest()) {
	    if (fit.Lossless())
	      return;

	    error = fit.GetError();
	    if (cluster || 1)
	      bestmde = mode,
	      besttyp = 1;
	  }
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
  DecompressPaletteBtc((u8*)rgba, block);
#endif
}

template<typename dtyp>
void CompressDynamicBtc6(dtyp const* rgb, int mask, void* block, int flags)
{
  // ...
  flags = flags;
  block = block;
  mask  = mask;
  rgb   = rgb;
}

/* *****************************************************************************
 */
template<typename dtyp>
void CompressMaskedColourBtc1(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = block;

  // compress color separately if necessary
  CompressColorBtc1(rgba, mask, colourBlock, flags);
}

template<typename dtyp>
void CompressMaskedColourBtc2(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressColorBtc1(rgba, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc2(rgba, mask, alphaBlock);
}

template<typename dtyp>
void CompressMaskedColourBtc3(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressColorBtc1(rgba, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc3(rgba, mask, alphaBlock, flags);
}

template<typename dtyp>
void CompressMaskedNormalBtc1(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* normalBlock = block;

  // compress color separately if necessary
  CompressNormalBtc1(xyzd, mask, normalBlock, flags);
}

template<typename dtyp>
void CompressMaskedNormalBtc2(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressNormalBtc1(xyzd, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc2(xyzd, mask, alphaBlock);
}

template<typename dtyp>
void CompressMaskedNormalBtc3(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* colourBlock = reinterpret_cast<u8*>(block) + 8;
  void*  alphaBlock = block;

  // compress color separately if necessary
  CompressNormalBtc1(xyzd, mask, colourBlock, flags);
  // compress alpha separately if necessary
  CompressAlphaBtc3(xyzd, mask, alphaBlock, flags);
}

template<typename dtyp>
void CompressMaskedAlphaBtc4(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = block;

  // compress a into plane 1
  CompressAlphaBtc3(rgba - 3, mask, plane1Block, flags);
}

template<typename dtyp>
void CompressMaskedAlphaBtc5(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = block;
  void* plane2Block = reinterpret_cast<u8*>(block) + 8;

  // compress a into plane 1
  CompressAlphaBtc3(rgba - 3, mask, plane1Block, flags);
  // compress b into plane 2
  CompressAlphaBtc3(rgba - 2, mask, plane2Block, flags);
}

template<typename dtyp>
void CompressMaskedNormalBtc5(dtyp const* xyzd, int mask, void* block, int flags)
{
  // get the block locations
  void* plane1Block = block;
  void* plane2Block = reinterpret_cast<u8*>(block) + 8;

  // compress xy into plane 1/2
  CompressNormalsBtc5(xyzd, mask, plane1Block, plane2Block, flags);
}

template<typename dtyp>
void CompressMaskedDynamicBtc6(dtyp const* rgb, int mask, void* block, int flags)
{
  // get the block locations
  void* mixedBlock = block;

  // compress color and alpha merged if necessary
  CompressDynamicBtc6(rgb, mask, mixedBlock, flags);
}

template<typename dtyp>
void CompressMaskedPaletteBtc7(dtyp const* rgba, int mask, void* block, int flags)
{
  // get the block locations
  void* mixedBlock = block;

  // compress color and alpha merged if necessary
  CompressPaletteBtc7(rgba, mask, mixedBlock, flags);
}

void CompressMasked(u8 const* rgba, int mask, void* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
    CompressMaskedNormalBtc1(rgba, mask, block, flags);
  else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
    CompressMaskedNormalBtc2(rgba, mask, block, flags);
  else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
    CompressMaskedNormalBtc3(rgba, mask, block, flags);
  // DXT-type compression
  else if (flags & (kBtc1))
    CompressMaskedColourBtc1(rgba, mask, block, flags);
  else if (flags & (kBtc2))
    CompressMaskedColourBtc2(rgba, mask, block, flags);
  else if (flags & (kBtc3))
    CompressMaskedColourBtc3(rgba, mask, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
    CompressMaskedNormalBtc5(rgba, mask, block, flags);
  // ATI-type compression
  else if (flags & (kBtc4))
    CompressMaskedAlphaBtc4(rgba, mask, block, flags);
  else if (flags & (kBtc5))
    CompressMaskedAlphaBtc5(rgba, mask, block, flags);
  // BTC-type compression
  else if (flags & (kBtc7))
    CompressMaskedPaletteBtc7(rgba, mask, block, flags);
  else if (flags & (kBtc6))
    {}// while this is possible (up-cast), should we support it?
}

void CompressMasked(u16 const* rgba, int mask, void* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
    CompressMaskedNormalBtc1(rgba, mask, block, flags);
  else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
    CompressMaskedNormalBtc2(rgba, mask, block, flags);
  else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
    CompressMaskedNormalBtc3(rgba, mask, block, flags);
  // DXT-type compression
  else if (flags & (kBtc1))
    CompressMaskedColourBtc1(rgba, mask, block, flags);
  else if (flags & (kBtc2))
    CompressMaskedColourBtc2(rgba, mask, block, flags);
  else if (flags & (kBtc3))
    CompressMaskedColourBtc3(rgba, mask, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
    CompressMaskedNormalBtc5(rgba, mask, block, flags);
  // ATI-type compression
  else if (flags & (kBtc4))
    CompressMaskedAlphaBtc4(rgba, mask, block, flags);
  else if (flags & (kBtc5))
    CompressMaskedAlphaBtc5(rgba, mask, block, flags);
  // BTC-type compression
  else if (flags & (kBtc7))
    CompressMaskedPaletteBtc7(rgba, mask, block, flags);
  else if (flags & (kBtc6))
    CompressMaskedDynamicBtc6(rgba, mask, block, flags);
}

void CompressMasked(f23 const* rgba, int mask, void* block, int flags)
{
  // DXT-type compression
  /**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
    CompressMaskedNormalBtc1(rgba, mask, block, flags);
  else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
    CompressMaskedNormalBtc2(rgba, mask, block, flags);
  else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
    CompressMaskedNormalBtc1(rgba, mask, block, flags);
  // DXT-type compression
  else if (flags & (kBtc1))
    CompressMaskedColourBtc1(rgba, mask, block, flags);
  else if (flags & (kBtc2))
    CompressMaskedColourBtc2(rgba, mask, block, flags);
  else if (flags & (kBtc3))
    CompressMaskedColourBtc3(rgba, mask, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
    CompressMaskedNormalBtc5(rgba, mask, block, flags);
  // ATI-type compression
  else if (flags & (kBtc4))
    CompressMaskedAlphaBtc4(rgba, mask, block, flags);
  else if (flags & (kBtc5))
    CompressMaskedAlphaBtc5(rgba, mask, block, flags);
  // BTC-type compression
  else if (flags & (kBtc7))
    CompressMaskedPaletteBtc7(rgba, mask, block, flags);
  else if (flags & (kBtc6))
    CompressMaskedDynamicBtc6(rgba, mask, block, flags);
}

void Compress(u8 const* rgba, void* block, int flags)
{
  // compress with full mask
  CompressMasked(rgba, 0xFFFF, block, flags);
}

void Compress(u16 const* rgb, void* block, int flags)
{
  // compress with full mask
  CompressMasked(rgb, 0xFFFF, block, flags);
}

void Compress(f23 const* rgba, void* block, int flags)
{
  // compress with full mask
  CompressMasked(rgba, 0xFFFF, block, flags);
}

/* *****************************************************************************
 */
void DecompressDynamicBtc(u16* rgb, void const* block)
{
  // ...
  block = block;
  rgb = rgb;
}

template<typename dtyp>
void DecompressColourBtc1(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* colourBlock = block;

  // decompress colour
  DecompressColoursBtc1(rgba, colourBlock, true);
}

template<typename dtyp>
void DecompressColourBtc2(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* colourBlock = reinterpret_cast<u8 const* >(block) + 8;
  void const*  alphaBlock = block;

  // decompress colour
  DecompressColoursBtc1(rgba, colourBlock, false);
  // decompress alpha separately if necessary
  DecompressAlphaBtc2(rgba, alphaBlock);
}

template<typename dtyp>
void DecompressColourBtc3(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* colourBlock = reinterpret_cast<u8 const* >(block) + 8;
  void const*  alphaBlock = block;

  // decompress colour
  DecompressColoursBtc1(rgba, colourBlock, false);
  // decompress alpha separately if necessary
  DecompressAlphaBtc3(rgba, alphaBlock);
}

template<typename dtyp>
void DecompressAlphaBtc4(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = block;

  // decompress plane 1 into a
  DecompressAlphaBtc3(rgba - 3, plane1Block);
}

template<typename dtyp>
void DecompressAlphaBtc5(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = block;
  void const* plane2Block = reinterpret_cast<u8 const* >(block) + 8;

  // decompress plane 1 into a
  DecompressAlphaBtc3(rgba - 3, plane1Block);
  // decompress plane 2 into b
  DecompressAlphaBtc3(rgba - 2, plane2Block);
}

template<typename dtyp>
void DecompressNormalBtc5(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* plane1Block = block;
  void const* plane2Block = reinterpret_cast<u8 const* >(block) + 8;

  // compress xy into plane 1/2
  DecompressNormalsBtc5(rgba, plane1Block, plane2Block);
}

template<typename dtyp>
void DecompressDynamicBtc6(dtyp* rgb, void const* block, int flags)
{
  // get the block locations
  void const* mixedBlock = block;

  // decompress color and alpha merged if necessary
//DecompressDynamicsBtc6(rgb, mixedBlock);
  
  rgb = rgb;
  mixedBlock = mixedBlock;
  flags = flags;
}

template<typename dtyp>
void DecompressPaletteBtc7(dtyp* rgba, void const* block, int flags)
{
  // get the block locations
  void const* mixedBlock = block;

  // decompress color and alpha merged if necessary
  DecompressPalettesBtc7(rgba, mixedBlock);
}

void Decompress(u8* rgba, void const* block, int flags)
{
  // DXT-type compression
///**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
//  DecompressNormalBtc1(rgba, block, flags);
//else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
//  DecompressNormalBtc2(rgba, block, flags);
//else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
//  DecompressNormalBtc3(rgba, block, flags);

  // DXT-type compression
  if (flags & (kBtc1))
    DecompressColourBtc1(rgba, block, flags);
  else if (flags & (kBtc2))
    DecompressColourBtc2(rgba, block, flags);
  else if (flags & (kBtc3))
    DecompressColourBtc3(rgba, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
    DecompressNormalBtc5(rgba, block, flags);
  // ATI-type compression
  else if (flags & (kBtc4))
    DecompressAlphaBtc4(rgba, block, flags);
  else if (flags & (kBtc5))
    DecompressAlphaBtc5(rgba, block, flags);
  // BTC-type compression
  else if (flags & (kBtc7))
    DecompressPaletteBtc7(rgba, block, flags);
  else if (flags & (kBtc6))
    {}// while this is possible (down-cast), should we support it?
}

void Decompress(u16* rgba, void const* block, int flags)
{
  // DXT-type compression
///**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
//  DecompressNormalBtc1(rgba, block, flags);
//else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
//  DecompressNormalBtc2(rgba, block, flags);
//else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
//  DecompressNormalBtc3(rgba, block, flags);

  // DXT-type compression
  if (flags & (kBtc1))
    DecompressColourBtc1(rgba, block, flags);
  else if (flags & (kBtc2))
    DecompressColourBtc2(rgba, block, flags);
  else if (flags & (kBtc3))
    DecompressColourBtc3(rgba, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
    DecompressNormalBtc5(rgba, block, flags);
  // ATI-type compression
  else if (flags & (kBtc4))
    DecompressAlphaBtc4(rgba, block, flags);
  else if (flags & (kBtc5))
    DecompressAlphaBtc5(rgba, block, flags);
  // BTC-type compression
  else if (flags & (kBtc7))
    DecompressPaletteBtc7(rgba, block, flags);
  else if (flags & (kBtc6))
    DecompressDynamicBtc6(rgba, block, flags);
}

void Decompress(f23* rgba, void const* block, int flags)
{
  // DXT-type compression
///**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
//  DecompressNormalBtc1(rgba, block, flags);
//else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
//  DecompressNormalBtc2(rgba, block, flags);
//else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
//  DecompressNormalBtc3(rgba, block, flags);

  // DXT-type compression
  if (flags & (kBtc1))
    DecompressColourBtc1(rgba, block, flags);
  else if (flags & (kBtc2))
    DecompressColourBtc2(rgba, block, flags);
  else if (flags & (kBtc3))
    DecompressColourBtc3(rgba, block, flags);
  // 3Dc-type compression
  else if ((flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
    DecompressNormalBtc5(rgba, block, flags);
  // ATI-type compression
  else if (flags & (kBtc4))
    DecompressAlphaBtc4(rgba, block, flags);
  else if (flags & (kBtc5))
    DecompressAlphaBtc5(rgba, block, flags);
  // BTC-type compression
  else if (flags & (kBtc7))
    DecompressPaletteBtc7(rgba, block, flags);
  else if (flags & (kBtc6))
    DecompressDynamicBtc6(rgba, block, flags);
}

/* *****************************************************************************
 */
int GetStorageRequirements(int width, int height, int flags)
{
  // compute the storage requirements
  int blockcount = ((width + 3) / 4) * ((height + 3) / 4);
  int blocksize  = 16;

  if (flags & (kBtc1 | kBtc2 | kBtc3))
    blocksize = ((flags & kBtc1) != 0) ? 8 : 16;
  else if (flags & (kBtc4 | kBtc5))
    blocksize = ((flags & kBtc4) != 0) ? 8 : 16;
  else if (flags & (kBtc6 | kBtc7))
    blocksize =                              16;

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

  if (s.flags & (kBtc1 | kBtc2 | kBtc3))
    s.blocksize = ((flags & kBtc1) != 0) ? 8 : 16;
  else if (s.flags & (kBtc4 | kBtc5))
    s.blocksize = ((flags & kBtc4) != 0) ? 8 : 16;
  else if (s.flags & (kBtc6 | kBtc7))
    s.blocksize =                              16;

  s.compressedsize = s.blockcount * s.blocksize;
  s.decompressedsize = sizeof(u8) * 4 * (width * height);

  if (datatype == sqio::dtp::DT_U8) {
    // DXT-type compression
    /**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc1<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc1<u8>;
    else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc2<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc2<u8>;
    else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc3<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc3<u8>;
    // DXT-type compression
    else if (s.flags & (kBtc1))
      s.encoder = (sqio::enc)CompressMaskedColourBtc1<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc1<u8>;
    else if (s.flags & (kBtc2))
      s.encoder = (sqio::enc)CompressMaskedColourBtc2<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc2<u8>;
    else if (s.flags & (kBtc3))
      s.encoder = (sqio::enc)CompressMaskedColourBtc3<u8>,
      s.decoder = (sqio::dec)DecompressColourBtc3<u8>;
    // 3Dc-type compression
    else if ((s.flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5<u8>,
      s.decoder = (sqio::dec)DecompressNormalBtc5<u8>;
    // ATI-type compression
    else if (s.flags & (kBtc4))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4<u8>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4<u8>;
    else if (s.flags & (kBtc5))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5<u8>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5<u8>;
    // BTC-type compression
    else if (s.flags & (kBtc7))
      s.encoder = (sqio::enc)CompressMaskedPaletteBtc7<u8>,
      s.decoder = (sqio::dec)DecompressPaletteBtc7<u8>;
    else if (s.flags & (kBtc6))
      {}// while this is possible (down-cast), should we support it?
  }
  else if (datatype == sqio::dtp::DT_U16) {
    // DXT-type compression
    /**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc1<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc1<u16>;
    else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc2<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc2<u16>;
    else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc3<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc3<u16>;
    // DXT-type compression
    else if (s.flags & (kBtc1))
      s.encoder = (sqio::enc)CompressMaskedColourBtc1<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc1<u16>;
    else if (s.flags & (kBtc2))
      s.encoder = (sqio::enc)CompressMaskedColourBtc2<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc2<u16>;
    else if (s.flags & (kBtc3))
      s.encoder = (sqio::enc)CompressMaskedColourBtc3<u16>,
      s.decoder = (sqio::dec)DecompressColourBtc3<u16>;
    // 3Dc-type compression
    else if ((s.flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5<u16>,
      s.decoder = (sqio::dec)DecompressNormalBtc5<u16>;
    // ATI-type compression
    else if (s.flags & (kBtc4))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4<u16>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4<u16>;
    else if (s.flags & (kBtc5))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5<u16>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5<u16>;
    // BTC-type compression
    else if (s.flags & (kBtc7))
      s.encoder = (sqio::enc)CompressMaskedPaletteBtc7<u16>,
      s.decoder = (sqio::dec)DecompressPaletteBtc7<u16>;
    else if (s.flags & (kBtc6))
      s.encoder = (sqio::enc)CompressMaskedDynamicBtc6<u16>,
      s.decoder = (sqio::dec)DecompressDynamicBtc6<u16>;
  }
  else if (datatype == sqio::dtp::DT_F23) {
    // DXT-type compression
    /**/ if ((flags & (kBtc1 | kColourMetricUnit)) == (kBtc1 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc1<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc1<f23>;
    else if ((flags & (kBtc2 | kColourMetricUnit)) == (kBtc2 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc2<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc2<f23>;
    else if ((flags & (kBtc3 | kColourMetricUnit)) == (kBtc3 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc3<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc3<f23>;
    // DXT-type compression
    else if (s.flags & (kBtc1))
      s.encoder = (sqio::enc)CompressMaskedColourBtc1<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc1<f23>;
    else if (s.flags & (kBtc2))
      s.encoder = (sqio::enc)CompressMaskedColourBtc2<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc2<f23>;
    else if (s.flags & (kBtc3))
      s.encoder = (sqio::enc)CompressMaskedColourBtc3<f23>,
      s.decoder = (sqio::dec)DecompressColourBtc3<f23>;
    // 3Dc-type compression
    else if ((s.flags & (kBtc5 | kColourMetricUnit)) == (kBtc5 | kColourMetricUnit))
      s.encoder = (sqio::enc)CompressMaskedNormalBtc5<f23>,
      s.decoder = (sqio::dec)DecompressNormalBtc5<f23>;
    // ATI-type compression
    else if (s.flags & (kBtc4))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc4<f23>,
      s.decoder = (sqio::dec)DecompressAlphaBtc4<f23>;
    else if (s.flags & (kBtc5))
      s.encoder = (sqio::enc)CompressMaskedAlphaBtc5<f23>,
      s.decoder = (sqio::dec)DecompressAlphaBtc5<f23>;
    // BTC-type compression
    else if (s.flags & (kBtc7))
      s.encoder = (sqio::enc)CompressMaskedPaletteBtc7<f23>,
      s.decoder = (sqio::dec)DecompressPaletteBtc7<f23>;
    else if (s.flags & (kBtc6))
      s.encoder = (sqio::enc)CompressMaskedDynamicBtc6<f23>,
      s.decoder = (sqio::dec)DecompressDynamicBtc6<f23>;
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
  u8* targetBlock = reinterpret_cast< u8* >(blocks);
  int bytesPerBlock = 16;

  if (flags & (kBtc1 | kBtc2 | kBtc3))
    bytesPerBlock = ((flags & kBtc1) != 0) ? 8 : 16;
  else if (flags & (kBtc4 | kBtc5))
    bytesPerBlock = ((flags & kBtc4) != 0) ? 8 : 16;
  else if (flags & (kBtc6 | kBtc7))
    bytesPerBlock =                              16;

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
  u8 const* sourceBlock = reinterpret_cast< u8 const* >(blocks);
  int bytesPerBlock = 16;

  if (flags & (kBtc1 | kBtc2 | kBtc3))
    bytesPerBlock = ((flags & kBtc1) != 0) ? 8 : 16;
  else if (flags & (kBtc4 | kBtc5))
    bytesPerBlock = ((flags & kBtc4) != 0) ? 8 : 16;
  else if (flags & (kBtc6 | kBtc7))
    bytesPerBlock =                              16;

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

void CompressColorBtc1( tile_barrier barrier, const int thread,
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
