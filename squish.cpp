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
#include <memory.h>

#include "colourset.h"
#include "paletteset.h"

#include "maths.h"

#include "colourrangefit.h"
#include "colourclusterfit.h"
#include "colourblock.h"

#include "paletterangefit.h"
#include "paletteclusterfit.h"
#include "paletteblock.h"

#include "alpha.h"
#include "singlecolourfit.h"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(USE_PRE)
static int FixFlags( int flags )
{
  // grab the flag bits
  int method = flags & (kBtc1 | kBtc2 | kBtc3 | kBtc4 | kBtc5 | kBtc6 | kBtc7);
  int fit    = flags & (kColourRangeFit | kColourIterativeClusterFits);
  int metric = flags & (kColourMetricUniform | kColourMetricPerceptual | kColourMetricUnit);
  int extra  = flags & (kWeightColourByAlpha);
  int mode   = flags & (kVariableCodingModes);

  // set defaults
  if (!method)
    method = kBtc1;

  if (!fit)
    fit = kColourClusterFit;
  if (fit > (kColourClusterFit * 15))
    fit = (kColourClusterFit * 15);

  if (!metric)
    metric = kColourMetricPerceptual;

  if ((method == kBtc6) && (mode > kVariableCodingMode14))
    mode = 0;
  if ((method == kBtc7) && (mode > kVariableCodingMode8))
    mode = 0;

  // done
  return method | fit | metric | extra | mode;
}

void CompressColorBtc(u8 const* rgba, int mask, void* block, int flags)
{
  // create the minimal point set
  ColourSet colours(rgba, mask, flags);

  // check the compression type and compress colour
  if (colours.GetCount() == 1) {
    // always do a single colour fit
    SingleColourFit fit(&colours, flags);
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

void CompressMixedBtc(u8 const* rgba, int mask, void* block, int flags)
{
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
#if !defined(NDEBUG) && defined(DEBUG_SETTING)
#define DEBUG_MODE	kVariableCodingMode1
#define DEBUG_FIT	kColourClusterFit * 15

  flags = (flags & (~kVariableCodingModes)) | (DEBUG_MODE);
  flags = (flags & (~kColourIterativeClusterFit)) | (DEBUG_FIT);
#endif

  int numm = flags &  ( kVariableCodingModes),
        sm = (numm == 0 ? 1 : numm >> 24),
	em = (numm == 0 ? 7 : numm >> 24),
	om = (numm == 0 ? 1 :          0);
             flags &= (~kVariableCodingModes);

  static const int modeorder[3][8] = {
    {
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
      // order: sets (lo to hi), ibs (hi to lo), prc (hi to lo)
      kVariableCodingMode7, //{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
      kVariableCodingMode5, //{ 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
      kVariableCodingMode6, //{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },

      kVariableCodingMode8, //{ 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 },  // alpha variant of mode 2/4

      0, 0, 0, 0
    }
  };

  // limits sets to 3 and choose the partition freely
  int lmts =  3;
  int lmtp = -1;

  // use the same data-structure all the time
  PaletteSet bestpal;
  int bestmde = -1;
  int bestswp = -1;

  Vec4 error(FLT_MAX);

  for (int m = sm; m <= em; m++) {
    int mode = modeorder[om][m - 1];
    int mnum = (mode >> 24) - 1;

    // a mode has a specific number of sets, and variable rotations and partitions
    int nums = PaletteFit::GetNumSets      (mnum);
    int numr = PaletteFit::GetRotationBits (mnum);
    int nump = PaletteFit::GetPartitionBits(mnum);
    int numx = PaletteFit::GetSelectionBits(mnum);
    int numi = PaletteFit::GetIndexBits    (mnum);

    // stop if set-limit reached
    if (nums > lmts)
      break;

    // lock on the perfect partition
    int sp = (lmtp == -1 ?               0 : lmtp),
	ep = (lmtp == -1 ? (1 << nump) - 1 : lmtp);
    int sr =                             0,
	er =               (1 << numr) - 1;
    int sx =                             0,
	ex =               (1 << numx) - 1;

#if !defined(NDEBUG) && defined(DEBUG_SETTING)
#define DEBUG_PARTITION	0
#define DEBUG_ROTATION	0
#define DEBUG_SELECTION	0

    sp = ep = DEBUG_PARTITION;
    sr = er = DEBUG_ROTATION;
    sx = ex = DEBUG_SELECTION;
#endif

    // signal if we do we have anything better this iteration of the search
    bool better = false;
    // TODO: fix the cluster-fit for 1-channel cases, and for 3bit
    bool cluster = ((flags & kColourRangeFit) == 0) && (((numi >> 0) & 0xFF) < 3) && (((numi >> 16) & 0xFF) < 3) && (mode != kVariableCodingMode8);

    // search for the best partition
    for (int p = sp; p <= ep; p++) {
      // search for the best rotation
      for (int r = sr; r <= er; r++) {
	// create the minimal point set
	PaletteSet palette(rgba, mask, flags | mode, p, r);

	for (int x = sx; x <= ex; x++) {
	  // do a range fit (which uses single palette fit if appropriate)
	  PaletteRangeFit fit(&palette, flags | mode, x);

	  fit.SetError(error);
	  fit.PaletteFit::Compress(block);

	  // we could code it lossless, no point in trying any further at all
	  if (fit.Lossless())
	    return;
	  if (fit.IsBest()) {
	    error = fit.GetError();

#if !defined(NDEBUG) && defined(DEBUG_SETTING) && (defined(DEBUG_ENCODER) || defined(DEBUG_QUANTIZER))
	    if (cluster || 1)
#else
	    if (cluster)
#endif
	      bestmde = mode,
	      bestpal = palette,
	      bestswp = x,
	      better = true;
	  }
	}

	// if we see we have transparent values, back up from trying to test non-alpha only modes
	if (palette.IsTransparent() && (om == 1))
	  om = 2, em = 4;
	// if we see we have less colors than sets, back up from trying to test with more sets or even other partitions
	if (palette.GetCount() <= nums)
	  lmtp = p, lmts = nums;
      }
    }

    // check the compression type and compress palette of the chosen partition even better
    if (better) {
      // default to a cluster fit (could be iterative or not)
      PaletteClusterFit fit(&bestpal, flags | mode, bestswp);

      fit.SetError(error);
      fit.PaletteFit::Compress(block);

#if !defined(NDEBUG) && defined(DEBUG_QUANTIZER)
      // take the cluster-fit
      if (fit.IsBest())
	fit.PaletteFit::Decompress((u8*)rgba, mnum);

      // or re-create the range-fit
      else {
	PaletteRangeFit fit(&bestpal, flags | mode, bestswp);

	fit.SetError(error);
	fit.PaletteFit::Compress(block);
	fit.PaletteFit::Decompress((u8*)rgba, mnum);
      }
#endif

#if !defined(NDEBUG) && defined(DEBUG_ENCODER)
      DecompressPaletteBtc((u8*)rgba, block);
#endif
    }
  }
}

void CompressMasked(u8 const* rgba, int mask, void* block, int flags)
{
  // fix any bad flags
  flags = FixFlags(flags);

  // DXT-type compression
  if (flags & (kBtc1 | kBtc2 | kBtc3)) {
    // get the block locations
    void* colourBlock = block;
    void*  alphaBlock = block;
    if ((flags & (kBtc2 | kBtc3)) != 0)
      colourBlock = reinterpret_cast<u8*>(block) + 8;

    // compress color separately if necessary
    CompressColorBtc(rgba, mask, colourBlock, flags);

    // compress alpha separately if necessary
    if ((flags & kBtc2) != 0)
      CompressAlphaBtc2(rgba, mask, alphaBlock);
    else if ((flags & kBtc3) != 0)
      CompressAlphaBtc3(rgba, mask, alphaBlock);
  }
  // ATI-type compression
  else if (flags & (kBtc4 | kBtc5)) {
    // get the block locations
    void* plane1Block = block;
    void* plane2Block = block;
    if ((flags & (kBtc5)) != 0)
      plane2Block = reinterpret_cast<u8*>(block) + 8;

    // compress r into plane 1
    CompressAlphaBtc3(rgba - 3, mask, plane1Block);
    // compress g into plane 2 if necessary
    if ((flags & (kBtc5)) != 0)
      CompressAlphaBtc3(rgba - 2, mask, plane2Block);
  }
  // BTC-type compression
  else if (flags & (kBtc6 | kBtc7)) {
    // get the block locations
    void* mixedBlock = block;

    // compress color and alpha merged if necessary
    CompressMixedBtc(rgba, mask, mixedBlock, flags);
  }
}

void Compress(u8 const* rgba, void* block, int flags)
{
  // compress with full mask
  CompressMasked(rgba, 0xFFFF, block, flags);
}

void Decompress(u8* rgba, void const* block, int flags)
{
  // fix any bad flags
  flags = FixFlags(flags);

  // DXT-type compression
  if (flags & (kBtc1 | kBtc2 | kBtc3)) {
    // get the block locations
    void const* colourBlock = block;
    void const*  alphaBlock = block;
    if ((flags & (kBtc2 | kBtc3)) != 0)
      colourBlock = reinterpret_cast<u8 const* >(block) + 8;

    // decompress colour
    DecompressColourBtc(rgba, colourBlock, (flags & kBtc1) != 0);

    // decompress alpha separately if necessary
    if ((flags & kBtc2) != 0)
      DecompressAlphaBtc2(rgba, alphaBlock);
    else if ((flags & kBtc3) != 0)
      DecompressAlphaBtc3(rgba, alphaBlock);
  }
  // ATI-type compression
  else if (flags & (kBtc4 | kBtc5)) {
    // get the block locations
    void const* plane1Block = block;
    void const* plane2Block = block;
    if ((flags & (kBtc5)) != 0)
      plane2Block = reinterpret_cast<u8 const* >(block) + 8;

    // decompress plane 1 into r
    DecompressAlphaBtc3(rgba - 3, plane1Block);
    // decompress plane 2 into g if necessary
    if ((flags & (kBtc5)) != 0)
      DecompressAlphaBtc3(rgba - 2, plane2Block);
  }
  // BTC-type compression
  else if (flags & (kBtc6 | kBtc7)) {
    // get the block locations
    void const* mixedBlock = block;

    // decompress color and alpha merged if necessary
    DecompressPaletteBtc(rgba, mixedBlock);
  }
}

int GetStorageRequirements(int width, int height, int flags)
{
  // fix any bad flags
  flags = FixFlags(flags);

  // compute the storage requirements
  int blockcount = ((width + 3) / 4) * ((height + 3) / 4);
  int blocksize  = 16;

  if (flags & (kBtc1 | kBtc2 | kBtc3))
    blocksize = ((flags & kBtc1) != 0) ? 8 : 16;
  else if (flags & (kBtc4 | kBtc5))
    blocksize = ((flags & kBtc4) != 0) ? 8 : 16;
  else if (flags & (kBtc6 | kBtc7))
    blocksize = 16;

  return blockcount * blocksize;
}

void CompressImage(u8 const* rgba, int width, int height, void* blocks, int flags)
{
  // fix any bad flags
  flags = FixFlags(flags);

  // initialize the block output
  u8* targetBlock = reinterpret_cast< u8* >(blocks);
  int bytesPerBlock = 16;

  if (flags & (kBtc1 | kBtc2 | kBtc3))
    bytesPerBlock = ((flags & kBtc1) != 0) ? 8 : 16;
  else if (flags & (kBtc4 | kBtc5))
    bytesPerBlock = ((flags & kBtc4) != 0) ? 8 : 16;
  else if (flags & (kBtc6 | kBtc7))
    bytesPerBlock = 16;

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
  flags = FixFlags(flags);

  // initialize the block input
  u8 const* sourceBlock = reinterpret_cast< u8 const* >(blocks);
  int bytesPerBlock = 16;

  if (flags & (kBtc1 | kBtc2 | kBtc3))
    bytesPerBlock = ((flags & kBtc1) != 0) ? 8 : 16;
  else if (flags & (kBtc4 | kBtc5))
    bytesPerBlock = ((flags & kBtc4) != 0) ? 8 : 16;
  else if (flags & (kBtc6 | kBtc7))
    bytesPerBlock = 16;

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
#if	defined(USE_AMP) || defined(USE_COMPUTE)
#if	defined(USE_COMPUTE)
    tile_static SingleColourFit_CCR sfit;
    tile_static ColourRangeFit_CCR rfit;
    tile_static ClusterFit_CCR cfit;
#endif

void CompressColorBtc ( tile_barrier barrier, const int thread,
			pixel16 rgba, ColourSet_CCRr colours, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted {
  // all of these conditions are identical over the entire
  // thread-group, so all threads take one of 0,1,2, not
  // distinct branches

  // check the compression type and compress colour
  if (colours.GetCount() == 1) {
    // always do a single colour fit
#if	!defined(USE_COMPUTE)
    tile_static SingleColourFit_CCR sfit;
#endif

    sfit.AssignSet(barrier, thread, colours, metric, fit);
    sfit.Compress (barrier, thread, colours, block, trans, yArr, lArr);
  }
  else if ((fit == SQUISH_FIT_RANGE) || (colours.GetCount() == 0)) {
    // do a range fit
#if	!defined(USE_COMPUTE)
    tile_static ColourRangeFit_CCR rfit;
#endif

    rfit.AssignSet(barrier, thread, colours, metric, fit);
    rfit.Compress (barrier, thread, colours, block, trans, yArr);
  }
  else {
    // default to a cluster fit (could be iterative or not)
#if	!defined(USE_COMPUTE)
    tile_static ClusterFit_CCR cfit;
#endif

    cfit.AssignSet(barrier, thread, colours, metric, fit);
    cfit.Compress (barrier, thread, colours, block, trans, yArr);
  }
}

#if	defined(USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

void CompressColorBtc1( tile_barrier barrier, const int thread,
			pixel16 rgba, int mask, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted {
#if	!defined(USE_COMPUTE)
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
			IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted {
#if	!defined(USE_COMPUTE)
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
			IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted {
#if	!defined(USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

  // create the minimal point set
  colours.CountSet(barrier, thread, rgba, mask, false, trans);

  // compress colour
  CompressColorBtc(barrier, thread, rgba, colours, block, metric, trans, fit, yArr, lArr);
}
#endif

} // namespace squish
