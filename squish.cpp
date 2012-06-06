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
#include "colourset.h"
#include "maths.h"
#include "rangefit.h"
#include "clusterfit.h"
#include "colourblock.h"
#include "alpha.h"
#include "singlecolourfit.h"

namespace squish {

#if	!defined(USE_PRE)
static int FixFlags( int flags )
{
  // grab the flag bits
  int method = flags & ( kDxt1 | kDxt3 | kDxt5 );
  int fit = flags & ( kColourIterativeClusterFit | kColourClusterFit | kColourRangeFit );
  int metric = flags & ( kColourMetricPerceptual | kColourMetricUniform | kColourMetricUnit );
  int extra = flags & ( kWeightColourByAlpha );

  // set defaults
  if( method != kDxt3 && method != kDxt5 )
    method = kDxt1;
  if( fit != kColourRangeFit )
    fit = kColourClusterFit;
  if( metric != kColourMetricUniform )
    if( metric != kColourMetricUnit )
      metric = kColourMetricPerceptual;

  // done
  return method | fit | metric | extra;
}

void CompressColorDxt( u8 const* rgba, int mask, void* block, int flags )
{
  // create the minimal point set
  ColourSet colours( rgba, mask, flags );

  // check the compression type and compress colour
  if( colours.GetCount() == 1 )
  {
    // always do a single colour fit
    SingleColourFit fit( &colours, flags );
    fit.Compress( block );
  }
  else if( ( flags & kColourRangeFit ) != 0 || colours.GetCount() == 0 )
  {
    // do a range fit
    RangeFit fit( &colours, flags );
    fit.Compress( block );
  }
  else
  {
    // default to a cluster fit (could be iterative or not)
    ClusterFit fit( &colours, flags );
    fit.Compress( block );
  }
}

void CompressMasked( u8 const* rgba, int mask, void* block, int flags )
{
  // fix any bad flags
  flags = FixFlags( flags );

  // get the block locations
  void* colourBlock = block;
  void* alphaBock = block;
  if( ( flags & ( kDxt3 | kDxt5 ) ) != 0 )
    colourBlock = reinterpret_cast< u8* >( block ) + 8;

  // compress color separately if necessary
  CompressColorDxt( rgba, mask, colourBlock, flags );

  // compress alpha separately if necessary
  if( ( flags & kDxt3 ) != 0 )
    CompressAlphaDxt3( rgba, mask, alphaBock );
  else if( ( flags & kDxt5 ) != 0 )
    CompressAlphaDxt5( rgba, mask, alphaBock );
}

void Compress( u8 const* rgba, void* block, int flags )
{
  // compress with full mask
  CompressMasked( rgba, 0xffff, block, flags );
}

void Decompress( u8* rgba, void const* block, int flags )
{
  // fix any bad flags
  flags = FixFlags( flags );

  // get the block locations
  void const* colourBlock = block;
  void const* alphaBock = block;
  if( ( flags & ( kDxt3 | kDxt5 ) ) != 0 )
    colourBlock = reinterpret_cast< u8 const* >( block ) + 8;

  // decompress colour
  DecompressColour( rgba, colourBlock, ( flags & kDxt1 ) != 0 );

  // decompress alpha separately if necessary
  if( ( flags & kDxt3 ) != 0 )
    DecompressAlphaDxt3( rgba, alphaBock );
  else if( ( flags & kDxt5 ) != 0 )
    DecompressAlphaDxt5( rgba, alphaBock );
}

int GetStorageRequirements( int width, int height, int flags )
{
  // fix any bad flags
  flags = FixFlags( flags );

  // compute the storage requirements
  int blockcount = ( ( width + 3 )/4 ) * ( ( height + 3 )/4 );
  int blocksize = ( ( flags & kDxt1 ) != 0 ) ? 8 : 16;
  return blockcount*blocksize;
}

void CompressImage( u8 const* rgba, int width, int height, void* blocks, int flags )
{
  // fix any bad flags
  flags = FixFlags( flags );

  // initialise the block output
  u8* targetBlock = reinterpret_cast< u8* >( blocks );
  int bytesPerBlock = ( ( flags & kDxt1 ) != 0 ) ? 8 : 16;

  // loop over blocks
  for( int y = 0; y < height; y += 4 )
  {
    for( int x = 0; x < width; x += 4 )
    {
      // build the 4x4 block of pixels
      u8 sourceRgba[16*4];
      u8* targetPixel = sourceRgba;
      int mask = 0;
      for( int py = 0; py < 4; ++py )
      {
	for( int px = 0; px < 4; ++px )
	{
	  // get the source pixel in the image
	  int sx = x + px;
	  int sy = y + py;

	  // enable if we're in the image
	  if( sx < width && sy < height )
	  {
	    // copy the rgba value
	    u8 const* sourcePixel = rgba + 4*( width*sy + sx );
	    for( int i = 0; i < 4; ++i )
	      *targetPixel++ = *sourcePixel++;

	    // enable this pixel
	    mask |= ( 1 << ( 4*py + px ) );
	  }
	  else
	  {
	    // skip this pixel as its outside the image
	    targetPixel += 4;
	  }
	}
      }

      // compress it into the output
      CompressMasked( sourceRgba, mask, targetBlock, flags );

      // advance
      targetBlock += bytesPerBlock;
    }
  }
}

void DecompressImage( u8* rgba, int width, int height, void const* blocks, int flags )
{
  // fix any bad flags
  flags = FixFlags( flags );

  // initialise the block input
  u8 const* sourceBlock = reinterpret_cast< u8 const* >( blocks );
  int bytesPerBlock = ( ( flags & kDxt1 ) != 0 ) ? 8 : 16;

  // loop over blocks
  for( int y = 0; y < height; y += 4 )
  {
    for( int x = 0; x < width; x += 4 )
    {
      // decompress the block
      u8 targetRgba[4*16];
      Decompress( targetRgba, sourceBlock, flags );

      // write the decompressed pixels to the correct image locations
      u8 const* sourcePixel = targetRgba;
      for( int py = 0; py < 4; ++py )
      {
	for( int px = 0; px < 4; ++px )
	{
	  // get the target location
	  int sx = x + px;
	  int sy = y + py;
	  if( sx < width && sy < height )
	  {
	    u8* targetPixel = rgba + 4*( width*sy + sx );

	    // copy the rgba value
	    for( int i = 0; i < 4; ++i )
	      *targetPixel++ = *sourcePixel++;
	  }
	  else
	  {
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

#if	defined(USE_AMP) || defined(USE_COMPUTE)
#if	defined(USE_COMPUTE)
    tile_static SingleColourFit_CCR sfit;
    tile_static RangeFit_CCR rfit;
    tile_static ClusterFit_CCR cfit;
#endif

void CompressColorDxt ( tile_barrier barrier, const int thread,
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
    tile_static RangeFit_CCR rfit;
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

void CompressColorDxt1( tile_barrier barrier, const int thread,
			pixel16 rgba, int mask, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted {
#if	!defined(USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

  // create the minimal point set
  colours.CountSet(barrier, thread, rgba, mask, true, trans);

  // compress colour
  CompressColorDxt(barrier, thread, rgba, colours, block, metric, trans, fit, yArr, lArr);
}

void CompressColorDxt3( tile_barrier barrier, const int thread,
			pixel16 rgba, int mask, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted {
#if	!defined(USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

  // create the minimal point set
  colours.CountSet(barrier, thread, rgba, mask, false, trans);

  // compress colour
  CompressColorDxt(barrier, thread, rgba, colours, block, metric, trans, fit, yArr, lArr);
}

void CompressColorDxt5( tile_barrier barrier, const int thread,
			pixel16 rgba, int mask, out code64 block,
			int metric, bool trans, int fit,
			IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted {
#if	!defined(USE_COMPUTE)
  tile_static ColourSet_CCR colours;
#endif

  // create the minimal point set
  colours.CountSet(barrier, thread, rgba, mask, false, trans);

  // compress colour
  CompressColorDxt(barrier, thread, rgba, colours, block, metric, trans, fit, yArr, lArr);
}
#endif

} // namespace squish
