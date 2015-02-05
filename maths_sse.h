/* -----------------------------------------------------------------------------

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

#ifndef SQUISH_MATH_SSE_H
#define SQUISH_MATH_SSE_H

#include <xmmintrin.h>
#if ( SQUISH_USE_SSE > 1 )
#include <emmintrin.h>
#endif
#if ( SQUISH_USE_SSE >= 3 )
#include <pmmintrin.h>
#endif
#if ( SQUISH_USE_SSE >= 4 )
#include <smmintrin.h>
#endif
#if ( SQUISH_USE_XSSE == 3 )
#include <tmmintrin.h>
#endif
#if ( SQUISH_USE_XSSE == 4 )
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>
#endif
#endif

#pragma warning(disable: 4127)

namespace squish {
namespace math {
  
	static doinline float rcp(float in) {
		__m128 s, e, d;
		float r;

		s = _mm_load_ss(&in);
		e = _mm_rcp_ss(s);
		d = _mm_sub_ss(_mm_set1_ps(1.0f), _mm_mul_ss(e, s));
		s = _mm_add_ss(_mm_mul_ss(d, e), e);

		_mm_store_ss(&r, s);
		return r;
	}

	// 5-6% of execution time as std::sqrt
	static doinline float sqrt(float in) {
		__m128 s;
		float r;

		s = _mm_load_ss(&in);
		s = _mm_sqrt_ss(s);

		_mm_store_ss(&r, s);
		return r;
	}

	// 3-5% of execution time as std::pow
	static doinline float cbrt(float in) {
		__m128 n, x, c, u, v;
		float r;
	//	float check = std::pow(in, 1.0f/3.0f);

	//	n = _mm_load_ss(&check);
		n = _mm_load_ss(&in);
		// initial guess, poor: sqrt( n ) + 1e-15
		x = _mm_add_ss(_mm_sqrt_ss( n ), _mm_set1_ps(1e-15f));
		// initial guess, cool hack: Y = (((Y >> 17) * 0xAAAB)) + (709921077L << 0);
		x = _mm_castsi128_ps(_mm_add_epi32(_mm_mul_epu32(_mm_srli_epi32(_mm_castps_si128( n ), 17), _mm_set1_epi32(0xAAAB)), _mm_set1_epi32(709921077L)));

		for (int i = 0; i < 1; i++) {
		  // X * X * X
		  c = _mm_mul_ss(x, _mm_mul_ss(x, x));

		  // ((c * t) + (n * f))
		  u = _mm_add_ss(_mm_mul_ss(c, _mm_set1_ps(2)), _mm_mul_ss(n, _mm_set1_ps(4)));
		  // ((c * f) + (n * t))
		  v = _mm_add_ss(_mm_mul_ss(c, _mm_set1_ps(4)), _mm_mul_ss(n, _mm_set1_ps(2)));

		  // 1 / u
		  __m128 estimate = _mm_rcp_ss(v);
		  __m128 diff = _mm_sub_ss(_mm_set1_ps(1.0f), _mm_mul_ss(estimate, v));
		  v = _mm_add_ss(_mm_mul_ss(diff, estimate), estimate);

		  // ((c * t) + (n * f)) * (x) / ((c * f) + (n * t))
		  x = _mm_mul_ss(x, _mm_mul_ss(u, v));
		}

		_mm_store_ss(&r, x);
		return r;
	}

}
}

#endif
