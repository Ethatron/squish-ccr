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

#ifndef SQUISH_H
#define SQUISH_H

#include "config.h"

#if	defined(USE_COMPUTE) || defined(USE_AMP)
#include "singlecolourlookup_ccr.inl"
#include "degeneracy_ccr.inl"
#endif

//! All squish API functions live in this namespace.
namespace squish {

// -----------------------------------------------------------------------------

// Skip enum-type in HLSL
#if	!defined(USE_COMPUTE)
enum
{
	//! Use DXT1/BC1/BC1 compression.
	kBtc1 = ( 1 << 0 ),
	//! Use DXT3/BC2/BC2 compression.
	kBtc2 = ( 1 << 1 ),
	//! Use DXT5/BC3/BC3 compression.
	kBtc3 = ( 1 << 2 ),
	//! Use ATI1/BC4 compression.
	kBtc4 = ( 1 << 3 ),
	//! Use ATI2/BC5 compression.
	kBtc5 = ( 1 << 4 ),
	//! Use BC6H compression.
	kBtc6 = ( 1 << 5 ),
	//! Use BC7 compression.
	kBtc7 = ( 1 << 6 ),

	//! Use a uniform metric for colour error.
	kColourMetricUniform = ( 1 << 8 ),
	//! Use a perceptual metric for colour error (the default).
	kColourMetricPerceptual = ( 1 << 9 ),
	//! Use a unit metric for colour error.
	kColourMetricUnit = ( 1 << 10 ),

	//! Weight the colour by alpha during cluster fit (disabled by default).
	kWeightColourByAlpha = ( 1 << 11 ),
	//! Don't code alpha, set alpha to 255 after weighting (disabled by default).
	kExcludeAlphaFromPalette = ( 1 << 12 ),

	//! Transform input values/points from sRGB to linear RGB (disabled by default).
	kSrgbIn = ( 1 << 13 ),
	//! Transform output points/values from linear RGB to sRGB (disabled by default).
	kSrgbOut = ( 1 << 14 ),

	//! Use a fast but low quality colour compressor.
	kColourRangeFit	= ( 1 << 15 ),
	//! Use a slow but high quality colour compressor (the default).
	kColourClusterFit = ( 1 << 16 ),
	//! Use a very slow but very high quality colour compressor.
	kColourIterativeClusterFit  = (  8 << 16 ),
	//! Specify the number of iterations explicitly. You can go until 15.
	kColourIterativeClusterFit1 = (  1 << 16 ),
	kColourIterativeClusterFit2 = (  2 << 16 ),
	kColourIterativeClusterFit4 = (  4 << 16 ),
	kColourIterativeClusterFit8 = (  8 << 16 ),
	kColourIterativeClusterFits = ( 15 << 16 ),

	//! Use to code a specific BC6/7 mode, coded as "1 + mode-number" (disabled by default).
	kVariableCodingMode1  = (  1 << 24 ),
	kVariableCodingMode2  = (  2 << 24 ),
	kVariableCodingMode3  = (  3 << 24 ),
	kVariableCodingMode4  = (  4 << 24 ),
	kVariableCodingMode5  = (  5 << 24 ),
	kVariableCodingMode6  = (  6 << 24 ),
	kVariableCodingMode7  = (  7 << 24 ),
	kVariableCodingMode8  = (  8 << 24 ),
	kVariableCodingMode9  = (  9 << 24 ),
	kVariableCodingMode10 = ( 10 << 24 ),
	kVariableCodingMode11 = ( 11 << 24 ),
	kVariableCodingMode12 = ( 12 << 24 ),
	kVariableCodingMode13 = ( 13 << 24 ),
	kVariableCodingMode14 = ( 14 << 24 ),
	kVariableCodingModes  = ( 15 << 24 ),
};
#endif

// -----------------------------------------------------------------------------

#if	!defined(USE_PRE)
//! Typedef a quantity that is a single unsigned byte.
typedef unsigned char u8;
//! Typedef a quantity that is a single unsigned short.
typedef unsigned short u16;

// -----------------------------------------------------------------------------

/*! @brief Compresses a 4x4 block of pixels.

	@param rgba	The rgba values of the 16 source pixels.
	@param block	Storage for the compressed DXT/BTC block.
	@param flags	Compression flags.

	The source pixels should be presented as a contiguous array of 16 rgba
	values, with each component as 1 byte each. In memory this should be:

		{ r1, g1, b1, a1, .... , r16, g16, b16, a16 }

	The flags parameter should specify either kBtc1, kBtc2 or kBtc3 compression,
	however, DXT1/BC1 will be used by default if none is specified. When using DXT1/BC1
	compression, 8 bytes of storage are required for the compressed DXT/BTC block.
	DXT3/BC2 and DXT5/BC3 compression require 16 bytes of storage per block.

	The flags parameter can also specify a preferred colour compressor and
	colour error metric to use when fitting the RGB components of the data.
	Possible colour compressors are: kColourClusterFit (the default),
	kColourRangeFit or kColourIterativeClusterFit. Possible colour error metrics
	are: kColourMetricPerceptual (the default) or kColourMetricUniform. If no
	flags are specified in any particular category then the default will be
	used. Unknown flags are ignored.

	When using kColourClusterFit, an additional flag can be specified to
	weight the colour of each pixel by its alpha value. For images that are
	rendered using alpha blending, this can significantly increase the
	perceived quality.
*/
void Compress( u8 const* rgba, void* block, int flags );

// -----------------------------------------------------------------------------

/*! @brief Compresses a 4x4 block of pixels.

	@param rgba	The rgba values of the 16 source pixels.
	@param mask	The valid pixel mask.
	@param block	Storage for the compressed DXT/BTC block.
	@param flags	Compression flags.

	The source pixels should be presented as a contiguous array of 16 rgba
	values, with each component as 1 byte each. In memory this should be:

		{ r1, g1, b1, a1, .... , r16, g16, b16, a16 }

	The mask parameter enables only certain pixels within the block. The lowest
	bit enables the first pixel and so on up to the 16th bit. Bits beyond the
	16th bit are ignored. Pixels that are not enabled are allowed to take
	arbitrary colours in the output block. An example of how this can be used
	is in the CompressImage function to disable pixels outside the bounds of
	the image when the width or height is not divisible by 4.

	The flags parameter should specify either kBtc1, kBtc2 or kBtc3 compression,
	however, DXT1/BC1 will be used by default if none is specified. When using DXT1/BC1
	compression, 8 bytes of storage are required for the compressed DXT/BTC block.
	DXT3/BC2 and DXT5/BC3 compression require 16 bytes of storage per block.

	The flags parameter can also specify a preferred colour compressor and
	colour error metric to use when fitting the RGB components of the data.
	Possible colour compressors are: kColourClusterFit (the default),
	kColourRangeFit or kColourIterativeClusterFit. Possible colour error metrics
	are: kColourMetricPerceptual (the default) or kColourMetricUniform. If no
	flags are specified in any particular category then the default will be
	used. Unknown flags are ignored.

	When using kColourClusterFit, an additional flag can be specified to
	weight the colour of each pixel by its alpha value. For images that are
	rendered using alpha blending, this can significantly increase the
	perceived quality.
*/
void CompressMasked( u8 const* rgba, int mask, void* block, int flags );

// -----------------------------------------------------------------------------

/*! @brief Decompresses a 4x4 block of pixels.

	@param rgba	Storage for the 16 decompressed pixels.
	@param block	The compressed DXT/BTC block.
	@param flags	Compression flags.

	The decompressed pixels will be written as a contiguous array of 16 rgba
	values, with each component as 1 byte each. In memory this is:

		{ r1, g1, b1, a1, .... , r16, g16, b16, a16 }

	The flags parameter should specify either kBtc1, kBtc2 or kBtc3 compression,
	however, DXT1/BC1 will be used by default if none is specified. All other flags
	are ignored.
*/
void Decompress( u8* rgba, void const* block, int flags );

// -----------------------------------------------------------------------------

/*! @brief Computes the amount of compressed storage required.

	@param width	The width of the image.
	@param height	The height of the image.
	@param flags	Compression flags.

	The flags parameter should specify either kBtc1, kBtc2 or kBtc3 compression,
	however, DXT1/BC1 will be used by default if none is specified. All other flags
	are ignored.

	Most DXT/BTC images will be a multiple of 4 in each dimension, but this
	function supports arbitrary size images by allowing the outer blocks to
	be only partially used.
*/
int GetStorageRequirements( int width, int height, int flags );

// -----------------------------------------------------------------------------

/*! @brief Compresses an image in memory.

	@param rgba	The pixels of the source.
	@param width	The width of the source image.
	@param height	The height of the source image.
	@param blocks	Storage for the compressed output.
	@param flags	Compression flags.

	The source pixels should be presented as a contiguous array of width*height
	rgba values, with each component as 1 byte each. In memory this should be:

		{ r1, g1, b1, a1, .... , rn, gn, bn, an } for n = width*height

	The flags parameter should specify either kBtc1, kBtc2 or kBtc3 compression,
	however, DXT1/BC1 will be used by default if none is specified. When using DXT1/BC1
	compression, 8 bytes of storage are required for each compressed DXT/BTC block.
	DXT3/BC2 and DXT5/BC3 compression require 16 bytes of storage per block.

	The flags parameter can also specify a preferred colour compressor and
	colour error metric to use when fitting the RGB components of the data.
	Possible colour compressors are: kColourClusterFit (the default),
	kColourRangeFit or kColourIterativeClusterFit. Possible colour error metrics
	are: kColourMetricPerceptual (the default) or kColourMetricUniform. If no
	flags are specified in any particular category then the default will be
	used. Unknown flags are ignored.

	When using kColourClusterFit, an additional flag can be specified to
	weight the colour of each pixel by its alpha value. For images that are
	rendered using alpha blending, this can significantly increase the
	perceived quality.

	Internally this function calls squish::Compress for each block. To see how
	much memory is required in the compressed image, use
	squish::GetStorageRequirements.
*/
void CompressImage( u8 const* rgba, int width, int height, void* blocks, int flags );

// -----------------------------------------------------------------------------

/*! @brief Decompresses an image in memory.

	@param rgba	Storage for the decompressed pixels.
	@param width	The width of the source image.
	@param height	The height of the source image.
	@param blocks	The compressed DXT/BTC blocks.
	@param flags	Compression flags.

	The decompressed pixels will be written as a contiguous array of width*height
	16 rgba values, with each component as 1 byte each. In memory this is:

		{ r1, g1, b1, a1, .... , rn, gn, bn, an } for n = width*height

	The flags parameter should specify either kBtc1, kBtc2 or kBtc3 compression,
	however, DXT1/BC1 will be used by default if none is specified. All other flags
	are ignored.

	Internally this function calls squish::Decompress for each block.
*/
void DecompressImage( u8* rgba, int width, int height, void const* blocks, int flags );
#endif

/* *****************************************************************************
 */
#if	defined(USE_AMP) || defined(USE_COMPUTE)
//! Typedef a quantity that is a single unsigned byte.
typedef unsigned int ccr8;

// Start- and end-point positions in arrays, for colors
#define	CSTRT	0
#define	CSTOP	1
#define	CVALS	2

// Start- and end-point positions in arrays, for alphas
#define	ASTRT	0
#define	ASTOP	1
#define	AVALS	2

/* -----------------------------------------------------------------------------
 * Provide types for the C++ and AMP version, these a only for passing function
 * arguments as references
 */
#if	!defined(USE_COMPUTE)
typedef const int (&pixel16)[16][DIM];
typedef float3 (&lineC2)[CVALS];
typedef int (&lineI2)[CVALS];
typedef int (&lineA2)[AVALS];
typedef float3 (&point16)[16];
typedef float (&weight16)[16];
typedef ccr8 (&index16)[16];
typedef ccr8 (&index16x2)[2][16];
typedef ccr8 (&index8)[8];
typedef unsigned int (&code64)[2];
typedef float3 &float3r;
/* -----------------------------------------------------------------------------
 * Provide types for the DirectCompute version, these a only for passing function
 * arguments as "references" (no such thing as pointers and refs in HLSL)
 */
#else
typedef const int pixel16[16][DIM];
typedef float3 lineC2[CVALS];
typedef int lineI2[CVALS];
typedef int lineA2[AVALS];
typedef float3 point16[16];
typedef float weight16[16];
typedef ccr8 index16[16];
typedef ccr8 index16x2[2][16];
typedef ccr8 index8[8];
typedef unsigned int code64[2];
typedef float3 float3r;
#endif

// -----------------------------------------------------------------------------

#define	SQUISH_METRIC_UNIFORM		0
#define	SQUISH_METRIC_PERCEPTUAL	1
#define	SQUISH_METRIC_UNIT		2

#define	SQUISH_FIT_RANGE		0
#define	SQUISH_FIT_CLUSTER		1
#define	SQUISH_FIT_CLUSTERITERATIVE	8

// -----------------------------------------------------------------------------

#if	defined(USE_AMP)
struct ColourSet_CCR;

#if	!defined(USE_COMPUTE)
typedef ColourSet_CCR &ColourSet_CCRr;
#else
typedef ColourSet_CCR  ColourSet_CCRr;
#endif

void CompressColorBtc (tile_barrier barrier, const int thread,
		       pixel16 rgba, ColourSet_CCRr colours, out code64 block,
		       int metric, bool trans, int fit,
		       IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted;
void CompressColorBtc1(tile_barrier barrier, const int thread,
		       pixel16 rgba, int mask, out code64 block,
		       int metric, bool trans, int fit,
		       IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted;
void CompressColorBtc2(tile_barrier barrier, const int thread,
		       pixel16 rgba, int mask, out code64 block,
		       int metric, bool trans, int fit,
		       IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted;
void CompressColorBtc3(tile_barrier barrier, const int thread,
		       pixel16 rgba, int mask, out code64 block,
		       int metric, bool trans, int fit,
		       IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted;
void CompressColorBtc4(tile_barrier barrier, const int thread,
		       pixel16 z, int mask, out code64 block,
		       int metric, bool trans, int fit,
		       IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted;
void CompressColorBtc5(tile_barrier barrier, const int thread,
		       pixel16 xy, int mask, out code64 block,
		       int metric, bool trans, int fit,
		       IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted;
void CompressMixedBtc (tile_barrier barrier, const int thread,
		       pixel16 rgba, ColourSet_CCRr colours, out code64 block,
		       int metric, bool trans, int fit,
		       IndexBlockLUT yArr, SingleColourLUT lArr) amp_restricted;

#endif
#endif

} // namespace squish

#endif // ndef SQUISH_H

