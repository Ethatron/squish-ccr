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

/* use frequencies of colour/palette-entries to calculate errors
 * this is more exact but a bit slower and uses more memory
 */
#define	FEATURE_EXACT_ERROR

/* use the power-method to estimate the principle component axis
 * should be more precise if not faster
 */
#undef	FEATURE_POWERESTIMATE

/* push start/end values away from the midpoint if the codebook contains
 * less unique entries than possible indices
 * to fill four indices at least one axis need to have an interval of 4/255th
 */
#undef	FEATURE_ELIMINATE_FLATBOOKS

/* brute force search for the shared bits with the lowest error
 *
 * 1) don't care about the start/end point skew
 * 2) if start is lower and goes down, make stop go up and vice versa (heuristic)
 * 3) check all start/stop up/down combinations (*256 tries), incomplete implementation!
 *
 * normally this is not worth it in the current state
 */
#undef	FEATURE_SHAREDBITS_TRIALS

/* .............................................................................
 */

#ifndef NDEBUG

// adjustments working in "Debug" or "Release" builds:
// throw the quantized rgba values back into the input.image
#define	VERIFY_QUANTIZER
// throw the decoded rgba values back into the input-image
#undef	VERIFY_ENCODER

// print out lots of information about the algorithm behaviour
#define	TRACK_STATISTICS

// adjustments working only in "Debug" builds:
// code only a specific mode-setting
#undef	DEBUG_SETTING
// print out lots of information about the search
#undef	DEBUG_DETAILS

#endif // NDEBUG

#if defined(TRACK_STATISTICS)
namespace squish {
  extern struct statistics {
    int num_counts[8][64][4][16];
    int win_partition[8][64];
    int win_rotation[8][4];
    int win_swap[8][2][2];
    int win_cluster[8][2];
    int win_mode[8];
    int has_countsets[4];
    int has_noweightsets[8][4][2];
  } gstat;
}
#endif

/* -----------------------------------------------------------------------------
 */

// Set to 1 when building squish to use Altivec instructions.
#ifndef SQUISH_USE_ALTIVEC
#define SQUISH_USE_ALTIVEC 0
#endif

// Set to 1 or 2 or 3 or 4 when building squish to use SSE or SSE2, SSE3 or SSE4 instructions.
#ifndef SQUISH_USE_SSE
#define SQUISH_USE_SSE 0
#endif

// Set to 3 or 4 when building squish to use SSSE3 or SSE4A instructions.
#ifndef SQUISH_USE_XSSE
#define SQUISH_USE_XSSE 0
#endif

// Internally et SQUISH_USE_SIMD when either Altivec or SSE is available.
#if SQUISH_USE_ALTIVEC && SQUISH_USE_SSE
#error "Cannot enable both Altivec and SSE!"
#endif
#if SQUISH_USE_ALTIVEC || SQUISH_USE_SSE
#define SQUISH_USE_SIMD 1
#ifdef __GNUC__
#define a16		__attribute__ ((__aligned__ (16)))
typedef long long	__int64;
#else
#define a16		__declspec(align(16))
#endif
#else
#define SQUISH_USE_SIMD 0
#ifdef __GNUC__
#define a16		__attribute__ ((__aligned__ (4)))
typedef unsigned long long unsigned__int64;
#else
#define a16		__declspec(align(4))
typedef unsigned __int64 unsigned__int64;
#endif
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
