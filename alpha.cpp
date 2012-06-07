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

#include "alpha.h"
#include "maths.h"
#include "simd.h"

#if	!defined(USE_COMPUTE)
#include <cmath>
#include <algorithm>
#endif

namespace squish {

#ifndef FLOATTOINT
#define FLOATTOINT
  static int FloatToInt( float a, int limit ) ccr_restricted
  {
    // use ANSI round-to-zero behaviour to get round-to-nearest
    int i = ( int )( a + 0.5f );

    // clamp to the limit
    if( i < 0 )
      i = 0;
    else if( i > limit )
      i = limit;

    // done
    return i;
  }

#if	defined(USE_AMP) || defined(USE_COMPUTE)
  static int3 FloatToInt( float3 a, int3 limit ) amp_restricted
  {
#if	!defined(USE_COMPUTE)
    using namespace Concurrency::vector_math;
#endif

    // use ANSI round-to-zero behaviour to get round-to-nearest
    int3 i = (int3)(a + 0.5f);

    // clamp to the limit
    return minimax(i, 0, limit);
  }

  static int4 FloatToInt( float4 a, int4 limit ) amp_restricted
  {
#if	!defined(USE_COMPUTE)
    using namespace Concurrency::vector_math;
#endif

    // use ANSI round-to-zero behaviour to get round-to-nearest
    int4 i = (int4)(a + 0.5f);

    // clamp to the limit
    return minimax(i, 0, limit);
  }
#endif
#endif

#if	!defined(USE_PRE)
void CompressAlphaDxt3( u8 const* rgba, int mask, void* block )
{
  u8* bytes = reinterpret_cast< u8* >( block );

  // quantise and pack the alpha values pairwise
  for( int i = 0; i < 8; ++i )
  {
    // quantise down to 4 bits
    float alpha1 = ( float )rgba[8*i + 3] * ( 15.0f/255.0f );
    float alpha2 = ( float )rgba[8*i + 7] * ( 15.0f/255.0f );
    int quant1 = FloatToInt( alpha1, 15 );
    int quant2 = FloatToInt( alpha2, 15 );

    // set alpha to zero where masked
    int bit1 = 1 << ( 2*i );
    int bit2 = 1 << ( 2*i + 1 );
    if( ( mask & bit1 ) == 0 )
      quant1 = 0;
    if( ( mask & bit2 ) == 0 )
      quant2 = 0;

    // pack into the byte
    bytes[i] = ( u8 )( quant1 | ( quant2 << 4 ) );
  }
}

void DecompressAlphaDxt3( u8* rgba, void const* block )
{
  u8 const* bytes = reinterpret_cast< u8 const* >( block );

  // unpack the alpha values pairwise
  for( int i = 0; i < 8; ++i )
  {
    // quantise down to 4 bits
    u8 quant = bytes[i];

    // unpack the values
    u8 lo = quant & 0x0f;
    u8 hi = quant & 0xf0;

    // convert back up to bytes
    rgba[8*i + 3] = lo | ( lo << 4 );
    rgba[8*i + 7] = hi | ( hi >> 4 );
  }
}

static void FixRange( int& min, int& max, int steps )
{
  if( max - min < steps )
    max = std::min<int>( min + steps, 255 );
  if( max - min < steps )
    min = std::max<int>( 0, max - steps );
}

static int FitCodes( u8 const* rgba, int mask, u8 const* codes, u8* indices )
{
  // fit each alpha value to the codebook
  int err = 0;
  for( int i = 0; i < 16; ++i )
  {
    // check this pixel is valid
    int bit = 1 << i;
    if( ( mask & bit ) == 0 )
    {
      // use the first code
      indices[i] = 0;
      continue;
    }

    // find the least error and corresponding index
    int value = rgba[4*i + 3];
    int least = INT_MAX;
    int index = 0;
    for( int j = 0; j < 8; ++j )
    {
      // get the squared error from this code
      int dist = ( int )value - ( int )codes[j];
      dist *= dist;

      // compare with the best so far
      if( dist < least )
      {
	least = dist;
	index = j;
      }
    }

    // save this index and accumulate the error
    indices[i] = ( u8 )index;
    err += least;
  }

  // return the total error
  return err;
}

static void WriteAlphaBlock( int alpha0, int alpha1, u8 const* indices, void* block )
{
  u8* bytes = reinterpret_cast< u8* >( block );

  // write the first two bytes
  bytes[0] = ( u8 )alpha0;
  bytes[1] = ( u8 )alpha1;

  // pack the indices with 3 bits each
  u8* dest = bytes + 2;
  u8 const* src = indices;
  for( int i = 0; i < 2; ++i )
  {
    // pack 8 3-bit values
    int value = 0;
    for( int j = 0; j < 8; ++j )
    {
      int index = *src++;
      value |= ( index << 3*j );
    }

    // store in 3 bytes
    for( int j = 0; j < 3; ++j )
    {
      int byte = ( value >> 8*j ) & 0xff;
      *dest++ = ( u8 )byte;
    }
  }
}

static void WriteAlphaBlock5( int alpha0, int alpha1, u8 const* indices, void* block )
{
  // check the relative values of the endpoints
  if( alpha0 > alpha1 )
  {
    // swap the indices
    u8 swapped[16];
    for( int i = 0; i < 16; ++i )
    {
      u8 index = indices[i];
      if( index == 0 )
	swapped[i] = 1;
      else if( index == 1 )
	swapped[i] = 0;
      else if( index <= 5 )
	swapped[i] = 7 - index;
      else
	swapped[i] = index;
    }

    // write the block
    WriteAlphaBlock( alpha1, alpha0, swapped, block );
  }
  else
  {
    // write the block
    WriteAlphaBlock( alpha0, alpha1, indices, block );
  }
}

static void WriteAlphaBlock7( int alpha0, int alpha1, u8 const* indices, void* block )
{
  // check the relative values of the endpoints
  if( alpha0 < alpha1 )
  {
    // swap the indices
    u8 swapped[16];
    for( int i = 0; i < 16; ++i )
    {
      u8 index = indices[i];
      if( index == 0 )
	swapped[i] = 1;
      else if( index == 1 )
	swapped[i] = 0;
      else
	swapped[i] = 9 - index;
    }

    // write the block
    WriteAlphaBlock( alpha1, alpha0, swapped, block );
  }
  else
  {
    // write the block
    WriteAlphaBlock( alpha0, alpha1, indices, block );
  }
}

void CompressAlphaDxt5( u8 const* rgba, int mask, void* block )
{
  // get the range for 5-alpha and 7-alpha interpolation
  int min5 = 255;
  int max5 = 0;
  int min7 = 255;
  int max7 = 0;
  for( int i = 0; i < 16; ++i )
  {
    // check this pixel is valid
    int bit = 1 << i;
    if( ( mask & bit ) == 0 )
      continue;

    // incorporate into the min/max
    int value = rgba[4*i + 3];
    if( value < min7 )
      min7 = value;
    if( value > max7 )
      max7 = value;
    if( value != 0 && value < min5 )
      min5 = value;
    if( value != 255 && value > max5 )
      max5 = value;
  }

  // handle the case that no valid range was found
  if( min5 > max5 )
    min5 = max5;
  if( min7 > max7 )
    min7 = max7;

  // fix the range to be the minimum in each case
  FixRange( min5, max5, 5 );
  FixRange( min7, max7, 7 );

  // set up the 5-alpha code book
  u8 codes5[8];
  codes5[0] = ( u8 )min5;
  codes5[1] = ( u8 )max5;
  for( int i = 1; i < 5; ++i )
    codes5[1 + i] = ( u8 )( ( ( 5 - i )*min5 + i*max5 )/5 );
  codes5[6] = 0;
  codes5[7] = 255;

  // set up the 7-alpha code book
  u8 codes7[8];
  codes7[0] = ( u8 )min7;
  codes7[1] = ( u8 )max7;
  for( int i = 1; i < 7; ++i )
    codes7[1 + i] = ( u8 )( ( ( 7 - i )*min7 + i*max7 )/7 );

  // fit the data to both code books
  u8 indices5[16];
  u8 indices7[16];
  int err5 = FitCodes( rgba, mask, codes5, indices5 );
  int err7 = FitCodes( rgba, mask, codes7, indices7 );

  // save the block with least error
  if( err5 <= err7 )
    WriteAlphaBlock5( min5, max5, indices5, block );
  else
    WriteAlphaBlock7( min7, max7, indices7, block );
}

void DecompressAlphaDxt5( u8* rgba, void const* block )
{
  // get the two alpha values
  u8 const* bytes = reinterpret_cast< u8 const* >( block );
  int alpha0 = bytes[0];
  int alpha1 = bytes[1];

  // compare the values to build the codebook
  u8 codes[8];
  codes[0] = ( u8 )alpha0;
  codes[1] = ( u8 )alpha1;
  if( alpha0 <= alpha1 )
  {
    // use 5-alpha codebook
    for( int i = 1; i < 5; ++i )
      codes[1 + i] = ( u8 )( ( ( 5 - i )*alpha0 + i*alpha1 )/5 );
    codes[6] = 0;
    codes[7] = 255;
  }
  else
  {
    // use 7-alpha codebook
    for( int i = 1; i < 7; ++i )
      codes[1 + i] = ( u8 )( ( ( 7 - i )*alpha0 + i*alpha1 )/7 );
  }

  // decode the indices
  u8 indices[16];
  u8 const* src = bytes + 2;
  u8* dest = indices;
  for( int i = 0; i < 2; ++i )
  {
    // grab 3 bytes
    int value = 0;
    for( int j = 0; j < 3; ++j )
    {
      int byte = *src++;
      value |= ( byte << 8*j );
    }

    // unpack 8 3-bit values from it
    for( int j = 0; j < 8; ++j )
    {
      int index = ( value >> 3*j ) & 0x7;
      *dest++ = ( u8 )index;
    }
  }

  // write out the indexed codebook values
  for( int i = 0; i < 16; ++i )
    rgba[4*i + 3] = codes[indices[i]];
}
#endif

#if	defined(USE_AMP) || defined(USE_COMPUTE)
/* C++ AMP version */
static void FixRange(out lineA2 aline, int steps) amp_restricted
{
  if ((aline[ASTOP] - aline[ASTRT]) < steps)
    aline[ASTOP] = ((aline[ASTRT] + steps) < 255 ? (aline[ASTRT] + steps) : 255);
  if ((aline[ASTOP] - aline[ASTRT]) < steps)
    aline[ASTRT] = (  0 > (aline[ASTOP] - steps) ? (aline[ASTOP] - steps) :   0);
}

#if	defined(USE_COMPUTE)
  tile_static int cerror;
#endif

/* C++ AMP version */
static int FitCodes(tile_barrier barrier, const int thread, pixel16 rgba, int mask, index8 codes, out index16 matched) amp_restricted
{
  // fit each alpha value to the codebook (AMP: prefer vectorization over parallelism)
#if	!defined(USE_COMPUTE)
  tile_static int cerror;
#endif

  threaded_cse(0) {
    cerror = 0; }
  wavefrnt_for(fcscan, 16) {
    // use the first code by default
    matched[fcscan] = 0;

    // check this pixel is valid
    int bit = 1 << fcscan;
    if ((mask & bit) != 0) {
      // find the least error and corresponding index
      const int value = rgba[fcscan][0];
      int dists[8];

      dists[0] = value - (int)codes[0];
      dists[1] = value - (int)codes[1];
      dists[2] = value - (int)codes[2];
      dists[3] = value - (int)codes[3];

      dists[4] = value - (int)codes[4];
      dists[5] = value - (int)codes[5];
      dists[6] = value - (int)codes[6];
      dists[7] = value - (int)codes[7];

      // find the closest code (AMP: vectorized reduction, cset)
      int4 idx07, idx;

      idx07.x = 1 * (0 + (dists[1] < dists[0] ? 1 : 0));
      idx07.y = 1 * (2 + (dists[3] < dists[2] ? 1 : 0));
      idx07.z = 1 * (4 + (dists[5] < dists[4] ? 1 : 0));
      idx07.w = 1 * (6 + (dists[7] < dists[6] ? 1 : 0));

      idx.x = (dists[idx07.x] < dists[idx07.y] ? idx07.x : idx07.y);
      idx.y = (dists[idx07.z] < dists[idx07.w] ? idx07.z : idx07.w);

      idx.z = (dists[idx.x] < dists[idx.y] ? idx.x : idx.y);

      // save the index
      matched[fcscan] = (ccr8)idx.z;

      // accumulate the error
      threaded_add(cerror, dists[idx.z]);
    }
  }

  // return the total error
  return cerror;
}

/* C++ AMP version */
static void WriteAlphaBlock(tile_barrier barrier, const int thread, lineA2 alpha, index16 indices, out code64 block) amp_restricted
{
  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  threaded_cse(0) {
    // write the first two bytes
    block[0] =
      (alpha  [ASTRT] <<  0) |
      (alpha  [ASTOP] <<  8) |
      (indices[ 0] << 16) | (indices[ 1] << 19) |
      (indices[ 2] << 22) | (indices[ 3] << 25) |
      (indices[ 4] << 28) | (indices[ 5] << 31);

    // pack the indices with 3 bits each
    block[1] =
      (indices[ 5] >>  1) |
      (indices[ 6] <<  2) | (indices[ 7] <<  5) |
      (indices[ 8] <<  8) | (indices[ 9] << 11) |
      (indices[10] << 14) | (indices[11] << 17) |
      (indices[12] << 20) | (indices[13] << 23) |
      (indices[14] << 26) | (indices[15] << 29);
  }
}

/* C++ AMP version */
static void WriteAlphaBlock5(tile_barrier barrier, const int thread, lineA2 alpha, inout index16 indices, out code64 block,
			     IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 5, 4, 3, 2, 6, 7 };	// (a >  b)

  // degenerate case "start > stop" (AMP: local register, not enough for group-shared)
  int sorted[AVALS];

  // alpha[ASTRT] > alpha[ASTOP] := sorted[ASTRT] == alpha[ASTOP]
  sorted[ASTRT] = min(alpha[ASTRT], alpha[ASTOP]);
  sorted[ASTOP] = max(alpha[ASTRT], alpha[ASTOP]);

  // swap the indices
  wavefrnt_for(i, 16) {
#if	defined(USE_COMPUTE)
    indices[i] = (sorted[ASTRT] == alpha[ASTOP] ? yArr[IBL_ALPHA5][indices[i]].mapped : indices[i]);
#else
    indices[i] = (sorted[ASTRT] == alpha[ASTOP] ? yArr(IBL_ALPHA5, indices[i]).mapped : indices[i]);
#endif
  }

  // write the block
  WriteAlphaBlock(barrier, thread, sorted, indices, block);
}

/* C++ AMP version */
static void WriteAlphaBlock7(tile_barrier barrier, const int thread, lineA2 alpha, inout index16 indices, out code64 block,
			     IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 7, 6, 5, 4, 3, 2 };	// (a <  b)

  // degenerate case "start < stop" (AMP: local register, not enough for group-shared)
  int sorted[AVALS];

  // alpha[ASTRT] < alpha[ASTOP] := sorted[ASTRT] == alpha[ASTRT]
  sorted[ASTRT] = max(alpha[ASTRT], alpha[ASTOP]);
  sorted[ASTOP] = min(alpha[ASTRT], alpha[ASTOP]);

  // swap the indices
  wavefrnt_for(i, 16) {
#if	defined(USE_COMPUTE)
    indices[i] = (sorted[ASTRT] == alpha[ASTRT] ? yArr[IBL_ALPHA7][indices[i]].mapped : indices[i]);
#else
    indices[i] = (sorted[ASTRT] == alpha[ASTRT] ? yArr(IBL_ALPHA7, indices[i]).mapped : indices[i]);
#endif
  }

  // write the block
  WriteAlphaBlock(barrier, thread, sorted, indices, block);
}

#define	ERR5	0
#define	ERR7	1
#define	ERRS	2
#define CNUMS	4
#define CBITS	2	// number of bits
#define CMASK	3	// mask

#if	defined(USE_COMPUTE)
  tile_static ccr8 acodes [ERRS][ 8];
  tile_static ccr8 matched[ERRS][16];
  tile_static int aline[ERRS][AVALS];
//tile_static int error[ERRS];
  tile_static int error[4];
#endif

/* C++ AMP version */
void CompressAlphaDxt5(tile_barrier barrier, const int thread, pixel16 rgba, int mask, out code64 block,
		       IndexBlockLUT yArr) amp_restricted
{
  // get the range for 5-alpha and 7-alpha interpolation
#if	!defined(USE_COMPUTE)
  tile_static ccr8 acodes [ERRS][ 8];
  tile_static ccr8 matched[ERRS][16];
  tile_static int aline[ERRS][AVALS];
  tile_static int error[ERRS];
#endif

  /* initialize:
   *
   *  aline[ERR7][ASTRT] = 255;
   *  aline[ERR7][ASTOP] = 0;
   *  aline[ERR5][ASTRT] = 255;
   *  aline[ERR5][ASTOP] = 0;
   */
  threaded_for(li, CNUMS) {
    const int ewch = li & 1;
    const int cpos = li >> 1;

    aline[ewch][cpos] = (cpos ? 0 : 255);
  }

  wavefrnt_for(mmscan, 16) {
    // check this pixel is valid
    int bit = 1 << mmscan;
    if ((mask & bit) != 0) {
      // incorporate into the min/max
      const int value = rgba[mmscan][0];

      threaded_min(aline[ERR7][ASTRT], value);
      threaded_max(aline[ERR7][ASTOP], value);
      if ((value != 0) && (value != 255)) {
	threaded_min(aline[ERR5][ASTRT], value);
	threaded_max(aline[ERR5][ASTOP], value);
      }
    }
  }

  /* fix the range to be the minimum in each case
   *
   *  FixRange(aline[ERR5], 5);
   *  FixRange(aline[ERR7], 7);
   */
  threaded_for(ec, ERRS) {
    const int ewch = ec;
    const int step = 5 + (ec * 2);

    /* handle the case that no valid range was found
     *
     *  threaded_min(aline[ERR5][ASTRT], aline[ERR5][ASTOP]);
     *  threaded_min(aline[ERR7][ASTRT], aline[ERR7][ASTOP]);
     */
    threaded_min(aline[ewch][ASTRT], aline[ewch][ASTOP]);

    FixRange(aline[ewch], step);
  }

  /* set up the 5-alpha code book
   * set up the 7-alpha code book
   *
   *  acodes[ERR5][0] = (ccr8)aline[ERR5][ASTRT];
   *  acodes[ERR5][1] = (ccr8)aline[ERR5][ASTOP];
   *
   *  acodes[ERR7][0] = (ccr8)aline[ERR7][ASTRT];
   *  acodes[ERR7][1] = (ccr8)aline[ERR7][ASTOP];
   */
  threaded_for(ei, CNUMS) {
    const int ewch = ei & 1;
    const int cpos = ei >> 1;

    acodes[ewch][cpos] = (ccr8)aline[ewch][cpos];
  }

  /* set up the 5-alpha code book
   * set up the 7-alpha code book
   *
   *  acodes[ERR5][0] = (ccr8)aline[ERR5][ASTRT];
   *  acodes[ERR5][1] = (ccr8)aline[ERR5][ASTOP];
   *
   *  acodes[ERR7][0] = (ccr8)aline[ERR7][ASTRT];
   *  acodes[ERR7][1] = (ccr8)aline[ERR7][ASTOP];
   *
   *  acodes[ERR5][6] = 0;
   *  acodes[ERR5][7] = 255;
   */
  threaded_for(es, 6) {
//  const int ewch = (es / 6);
//  const int step = 5 + (ewch * 2);

    // set up the 5-alpha code book
    acodes[ERR5][2 + es] = (ccr8)(((4 - es) * aline[ERR5][ASTRT] + (1 + es) * aline[ERR5][ASTOP]) / 5) & (es == 4 ? 0 : 255) | (es == 5 ? 255 : 0);
    // set up the 7-alpha code book
    acodes[ERR7][2 + es] = (ccr8)(((6 - es) * aline[ERR7][ASTRT] + (1 + es) * aline[ERR7][ASTOP]) / 7);
  }

  // make it visible
  tile_static_memory_fence(barrier);

  /* fit the data to both code books
   *
   *  erros[ERR5] = FitCodes(barrier, thread, rgba, mask, acodes[ERR5], matched[ERR5]);
   *  erros[ERR7] = FitCodes(barrier, thread, rgba, mask, acodes[ERR7], matched[ERR7]);
   */
  threaded_for(ef, ERRS) {
    error[ef] = FitCodes(barrier, thread, rgba, mask, acodes[ef], matched[ef]);
  }

  // AMP: final reduction (make the choice)
  {
    // save the block with least erros
    if (error[ERR5] <= error[ERR7])
      WriteAlphaBlock5(barrier, thread, aline[ERR5], matched[ERR5], block, yArr);
    else
      WriteAlphaBlock7(barrier, thread, aline[ERR7], matched[ERR7], block, yArr);
  }
}

#undef	ERR5
#undef	ERR7
#undef	ERRS
#undef	CNUMS
#undef	CBITS
#undef	CMASK

#if	defined(USE_COMPUTE)
  tile_static unsigned int frags[4];
#endif

/* C++ AMP version */
void CompressAlphaDxt3(tile_barrier barrier, const int thread, pixel16 rgba, int mask, code64 block,
		         IndexBlockLUT yArr) amp_restricted
{
#if	!defined(USE_COMPUTE)
  tile_static unsigned int frags[4];
#endif

  // quantise and pack the alpha values quad-wise (AMP: prefer vectorization over parallelism)
  threaded_for(qs, 4) {
    const int opos = qs << 2;
    const int hilo = qs >> 1;

    // quantise down to 4 bits
    float4 alpha = float4(
      (float)rgba[opos + 0][0],
      (float)rgba[opos + 1][0],
      (float)rgba[opos + 2][0],
      (float)rgba[opos + 3][0]
    );

    int4 quant = FloatToInt(alpha * 15.0f / 255.0f, 15);

    // set alpha to zero where masked
    quant.x = ((mask >> (opos + 0)) & 1) ? quant.x : 0;
    quant.y = ((mask >> (opos + 1)) & 1) ? quant.y : 0;
    quant.z = ((mask >> (opos + 2)) & 1) ? quant.z : 0;
    quant.w = ((mask >> (opos + 3)) & 1) ? quant.w : 0;

    // pack into the short
    frags[qs] = (
      (quant.x <<  0) |
      (quant.y <<  4) |
      (quant.z <<  8) |
      (quant.w << 12)
    ) << (hilo * 16);
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // AMP: final reduction (make the merge)
  threaded_cse(0) {
    // save the block
    block[0] = frags[1] | frags[0];
    block[1] = frags[3] | frags[2];
  }
}
#endif

} // namespace squish