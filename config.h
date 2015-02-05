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

#ifndef SQUISH_CONFIG_H
#define SQUISH_CONFIG_H

/* define the algorithm to determine the best start/end-point
 * for a single colour
 * as Ryg already observed, it's possible to always force the highest
 * precision interpolant to be at index "1" for the euclidian error-
 * distance case [the sole exception being the 4u1, but the improvement
 * there is too mior to justify the loss in speed and the growth of the
 * lookup tables]
 * if you want to use non-euclidian error-distances you can swap your
 * own LUTs in and go back to the iterative fit functions
 *
 * - ColourSingleSnap, ColourSingleFit
 * - PaletteSingleSnap, PaletteSingleFit
 */
#define	ColourSingleMatch	ColourSingleSnap
#define	PaletteSingleMatch	PaletteSingleSnap
#define	HDRSingleMatch		HDRSingleSnap

/* use the power-method to estimate the principle component axis
 * should be more precise if not faster
 */
#define	FEATURE_POWERESTIMATE

/* totally blank out any colours in the alpha==0 case when colours
 * are weighted by alpha (which is the indicator too assure the third
 * channel actually really is to be used for alpha or at least contains
 * the importance-factor)
 *
 * we assume the blanked out value receives at least a quantized and
 * rounded alpha value that makes the pixel unperceivable in case it
 * won't be assigned 0 by the matcher (very optimistic but probably
 * better than the "damage" done to the fit-algorithms in case the
 * colour gets through)
 *   
 * Results (ARGB):
 *   
 *				  omit					     keep
 *   #cluster		R       G       B	J		R       G       B	J
 *   RMSE no m.		2.0099	1.7891	1.9644	1.9235		2.0174	1.7958	1.9753	1.9319
 *   RMSE metric	2.0846	1.6368	2.1797	1.9812		2.0949	1.6444	2.1911	1.9911
 *   SSIM no m.			0.0014593				    0.0014612
 *   SSIM metric		0.0012824				    0.0012850
 *   
 *				  omit					     keep
 *   #range		R       G       B	J		R       G       B	J
 *   RMSE no m.		3.4200	2.8293	3.1534	3.1435		3.2591	2.7265	3.0184	3.0092
 *   RMSE metric	4.0410	2.3455	3.6444	3.4211		4.0410	2.3455	3.6444	3.4211
 *   SSIM no m.			0.0021402				    0.0023996
 *   SSIM metric		0.0020360				    0.0023973
 */
#define	FEATURE_IGNORE_ALPHA0

/* project the colour-values onto the principle component (which is
 * anchored at the centroid) instead of using the just the principle
 * direction
 * problems with the principle direction: if black exist it will always
 * end up as one endpoint (dot-product being 0)
 */
#define	FEATURE_RANGEFIT_PROJECT
#define	FEATURE_NORMALFIT_PROJECT

/* guarantee that the results of all z-complementing normal-fit algorithms
 * stay inside the unit-sphere
 * end-points, if not selected as index, may still go outside the unit-sphere
 */
#define	FEATURE_NORMALFIT_UNITGUARANTEE

/* inline the index-fit error-check, makes it use the minimal amount
 * of instructions possible at the cost of more code (~2x)
 * search all heuristically chosen interpolants for the best one, works
 * only if all checks are inlined
 *
 * first hit improvements:    thorough search:
 *   ib=2: x2.64		ib=2: x8.27	mse: 1.135e-4 vs. 4.221e-5
 *   ib=3: x3.39		ib=3: x7.22	mse: 1.360e-3 vs. 4.924e-4
 *   ib=4: x3.45		ib=4: x47.5	mse: 2.768e-5 vs. 7.118e-6
 */
#define	FEATURE_INDEXFIT_INLINED
#define	FEATURE_INDEXFIT_THOROUGH	true

/* - sqrt() the weights in the colourset, affects all fits
 */
#undef	FEATURE_WEIGHTS_ROOTED		// SSIM: 0.0012824 off, 0.0013257 on - RMSE: 1.9235 off, 1.9540 on

/* - use the metric building the covariance-matrix, affects all fits
 * - sqrt() the metric in the colourset, affects all fits
 *   use linear ( r²*rmetric  +  g²*gmetric  +  b²*bmetric  +  a²*ametric)
 *   instead of ((r*rmetric)² + (g*gmetric)² + (b*bmetric)² + (a*ametric)²)
 * - sqr() the metric for least-squares, affects only cluster-fit
 */
#undef	FEATURE_METRIC_COVARIANCE	// SSIM: 0.0013466 off, 0.0013596 on
#define	FEATURE_METRIC_ROOTED		// SSIM: 0.0006231 off, 0.0006097 on
#define	FEATURE_METRIC_SQUARED		// SSIM: 0.0006231 off, 0.0006102 on

/* push start/end values away from the midpoint if the codebook contains
 * less unique entries than possible indices
 * to fill four indices at least one axis need to have an interval of 4/255th
 */
#undef	FEATURE_ELIMINATE_FLATBOOKS

/* brute force search for the shared bits with the lowest error
 *
 * 0) ensure shared bit of 1 for opaque merged-alpha cases only (*1)
 * 1) trial shared bits for all cases when block is transparent (*4,*4*2^6)
 * 2) trial shared bits for transparent and low-precision cases (*4,*4*2^6,*4^3*2^4)
 * 3) trial all shared bits
 * 4) check all start/stop up/down combinations (*256 tries), incomplete implementation!
 */
#define	SHAREDBITS_TRIAL_ALPHAONLYOPAQUE	0
#define	SHAREDBITS_TRIAL_ALPHAONLY		1
#define	SHAREDBITS_TRIAL_LOWPRC			2
#define	SHAREDBITS_TRIAL_ALL			3
#define	SHAREDBITS_TRIAL_PERMUTE		4

#define	FEATURE_SHAREDBITS_TRIALS		SHAREDBITS_TRIAL_LOWPRC
 
/* enable the set-builder to detect straight lines of values
 */
#define	FEATURE_TEST_LINES

/* .............................................................................
 */

#ifndef NDEBUG

// adjustments working in "Debug" or "Release" builds:
// throw the quantized rgba values back into the input.image
#undef	VERIFY_QUANTIZER
// throw the decoded rgba values back into the input-image
#undef	VERIFY_ENCODER

// adjustments working only in "Debug" builds:
// code only a specific mode-setting
#undef	DEBUG_SETTING
// print out lots of information about the search
#undef	DEBUG_DETAILS

// print out lots of information about the algorithm behaviour
#undef	TRACK_STATISTICS

#endif // NDEBUG

#if defined(TRACK_STATISTICS)
namespace squish {
  extern struct statistics {
    int num_counts[8][64][4][16];
    int num_channels[8][64][4][5];
    int btr_index[5][2];
    int btr_cluster[8][2];
    int win_partition[8][64];
    int win_rotation[8][4];
    int win_swap[8][4][2];
    int win_cluster[8][2];
    int win_mode[8];
//  int num_lines[4];
    int has_countsets[4];
    int has_noweightsets[8][4][2];
    int num_poweritrs[64];
    int alpha[6];
    float err_index[5][2];
  } gstat;
}
#endif

/* -----------------------------------------------------------------------------
 */

// Set to 1 when building squish to use Altivec instructions.
#ifndef SQUISH_USE_ALTIVEC
#define SQUISH_USE_ALTIVEC  0
#endif

// Set to 1 or 2 or 3 or 4 when building squish to use SSE or SSE2, SSE3 or SSE4 instructions.
#ifndef SQUISH_USE_SSE
#define SQUISH_USE_SSE	    0
#endif

// Set to 3 or 4 when building squish to use SSSE3 or SSE4A instructions.
#ifndef SQUISH_USE_XSSE
#define SQUISH_USE_XSSE	    0
#endif

// Internally et SQUISH_USE_SIMD when either Altivec or SSE is available.
#if SQUISH_USE_ALTIVEC && SQUISH_USE_SSE
#error "Cannot enable both Altivec and SSE!"
#endif
#if SQUISH_USE_ALTIVEC || SQUISH_USE_SSE
#define SQUISH_USE_SIMD 1
#if defined(__GNUC__)
#define a16		__attribute__ ((__aligned__ (16)))
#elif defined(_MSC_VER)
#define a16		__declspec(align(16))
#else
#define a16		
#endif
#else
#define SQUISH_USE_SIMD 0
#if defined(__GNUC__)
#define a16		__attribute__ ((__aligned__ (4)))
#elif defined(_MSC_VER)
#define a16		__declspec(align(4))
#else
#define a16		
#endif
#endif

#if defined(__GNUC__)
typedef long long	__int64;
typedef unsigned long long unsigned__int64;
#elif defined(_MSC_VER)
typedef unsigned __int64 unsigned__int64;
#else
typedef long long	__int64;
typedef unsigned long long unsigned__int64;
#endif

/* *****************************************************************************
 * Turn explicit vectorization of in case AMP or DirectCompute is requested
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#undef	SQUISH_USE_ALTIVEC
#undef	SQUISH_USE_SSE
#undef	SQUISH_USE_SIMD

#define SQUISH_USE_ALTIVEC	0
#define SQUISH_USE_SSE		0
#define SQUISH_USE_SIMD		0
#endif

/* -----------------------------------------------------------------------------
 * provide common-interface implementations for the debug-version of the AMP
 * implementation
 */
#if	defined(SQUISH_USE_AMP) && defined(USE_AMP_DEBUG)
#define	tile_static						// void attribute
#define	tile_barrier	int					// dummy var
#define	tile_static_memory_fence(x)				// void function

#define	amp_restricted
#define	ccr_restricted
#define	static_hlsl
#define	extern_hlsl	extern
#define	inherit_hlsl	public
#define	public_hlsl	public:
#define	protected_hlsl	protected:
#define	private_hlsl	private:

#define	threaded_max(a, b)	if (a < b) a = b
#define	threaded_min(a, b)	if (a > b) a = b
#define	threaded_set(a, b, c)	if (a == b) a = c;
#define	threaded_xch(a, b)	a = b;
#define	threaded_add(a, b)	a += b;
#define	threaded_lor(a, b)	a |= b;
#define	threaded_inc(a, b)	a += 1;
#define	threaded_cse(n)						// passthrough
#define	threaded_for(i, n)	for (int i = 0; i < (n); i++)	// loop for real
#define	wavefrnt_for(i, n)	for (int i = 0; i < (n); i++)	// loop for real

#define	out							// compute attribute
#define	inout							// compute attribute
/* -----------------------------------------------------------------------------
 * provide common-interface implementations for the AMP implementation
 */
#elif	defined(SQUISH_USE_AMP)
#include <amp.h>
#include <amp_math.h>
#include <amp_graphics.h>
#include <amp_short_vectors.h>

#define	ccr_restricted	restrict(cpu,amp)
#define	amp_restricted	restrict(amp)
#define	static_hlsl
#define	extern_hlsl	extern
#define	inherit_hlsl	public
#define	public_hlsl	public:
#define	protected_hlsl	protected:
#define	private_hlsl	private:

#define	threaded_max(a, b)	atomic_fetch_max(&a, b)
#define	threaded_min(a, b)	atomic_fetch_min(&a, b)
#define	threaded_set(a, b, c)	atomic_compare_exchange(&a, &b, c)
#define	threaded_xch(a, b)	atomic_exchange(&a, b)
#define	threaded_add(a, b)	atomic_fetch_add(&a, b)
#define	threaded_lor(a, b)	atomic_fetch_or(&a, b)
#define	threaded_inc(a)		atomic_fetch_inc(&a)
#define	threaded_cse(n)		if (thread == n)
#define	threaded_for(i, n)	const int i = thread; if (i < (n))
#define	wavefrnt_for(i, n)	const int i = thread;

#define	out							// compute attribute
#define	inout							// compute attribute

typedef	::Concurrency::graphics::int_3		int3;
typedef	::Concurrency::graphics::int_4		int4;
typedef	/*::Concurrency::graphics::uint_3*/ unsigned int uint3[3];  /* array-access necessary */
typedef	/*::Concurrency::graphics::uint_4*/ unsigned int uint4[4];  /* array-access necessary */
typedef	::Concurrency::graphics::float_3	float3;
typedef	::Concurrency::graphics::float_4	float4;

using namespace ::Concurrency;
/* -----------------------------------------------------------------------------
 * provide common-interface implementations for the DirectCompute implementation
 */
#elif	defined(SQUISH_USE_COMPUTE)
#define	tile_static	groupshared
#define	tile_barrier	int
#define	tile_static_memory_fence(x)	GroupMemoryBarrier()

#define	amp_restricted
#define	ccr_restricted
#define	static_hlsl	static
#define	extern_hlsl
#define	inherit_hlsl	/*public*/				// C++ inheritance
#define	public_hlsl	/*public:*/				// C++ section
#define	protected_hlsl	/*protected:*/				// C++ section
#define	private_hlsl	/*private:*/				// C++ section

#define	threaded_max(a, b)	InterlockedMax(a, b)
#define	threaded_min(a, b)	InterlockedMin(a, b)
#define	threaded_set(a, b, c)	{ int o; InterlockedCompareExchange(a, b, c, o); }
#define	threaded_xch(a, b)	{ int o; InterlockedExchange(a, b, o); }
#define	threaded_add(a, b)	InterlockedAdd(a, b)
#define	threaded_lor(a, b)	InterlockedOr(a, b)
#define	threaded_inc(a)		InterlockedAdd(a, 1)
#define	threaded_cse(n)		if (thread == n)
#define	threaded_for(i, n)	const int i = thread; if (i < (n))
#define	wavefrnt_for(i, n)	const int i = thread;

// map to intrinsics
#define	powf		pow
#define	sqrtf		sqrt
#define	atan2f		atan2
#define	floorf		floor
#define	fabsf		abs
#define	ceilf		ceil
#define	sinf		sin
#define	cosf		cos

// map to intrinsics
#define	maximum		max
#define	minimum		min
#define	minimax		clamp
#define	muladd		mad
#define	recip		rcp
#define	truncate	trunc

#define FLT_MAX         3.402823466e+38F        /* max value */
#define FLT_EPSILON     1.192092896e-07F        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
/* -----------------------------------------------------------------------------
 */
#else
#define	ccr_restricted
#define	amp_restricted

#define	static_hlsl
#define	extern_hlsl	extern
#define	inherit_hlsl	public
#define	public_hlsl	public:
#define	protected_hlsl	protected:
#define	private_hlsl	private:

#if	!defined(SQUISH_USE_CPP)
#pragma message( "You may need to select a configuration." )
#pragma message( "These enable mutual exclusive implementations:" )
#pragma message( " SQUISH_USE_CPP      - if you are happy with the regular C++ code and want to omit this message" )
#pragma message( " SQUISH_USE_AMP      - if you want the AMP code" )
#pragma message( " SQUISH_USE_AMP_DEBG - if you want the AMP code without using AMP" )
#pragma message( " SQUISH_USE_COMPUTE  - if you want the DirectCompute code" )
#pragma message( "These disable implementations enabled by default:" )
#pragma message( " SQUISH_USE_PRE      - if you don't want the regular C++ code" )
#endif
#endif

/* *****************************************************************************
*/
#if	defined(SQUISH_USE_AMP)
#ifndef	DIM
#pragma message( "Input-array dimensionality has been set to 4" )

// at least RGBA, which is 4 components, can be more
#define	DIM	4
#endif
#endif

#ifdef __GNUC__
#define assume
#define doinline
#define	passreg		__fastcall
#else
#define assume		__assume
#define doinline	__forceinline
#define	passreg		__fastcall
#endif

#endif // ndef SQUISH_CONFIG_H
