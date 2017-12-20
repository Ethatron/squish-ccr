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

#ifndef SQUISH_SIMD_H
#define SQUISH_SIMD_H

#include "maths.h"

namespace squish {

  // FloatTo...
  extern unsigned short uhLUTb[1 << 9];
  extern char           uhLUTs[1 << 9];
  extern unsigned short shLUTb[1 << 9];
  extern char           shLUTs[1 << 9];

  static inline u16 FloatToUHalf(f23 c) {
    unsigned int f = *((unsigned int *)&c); return uhLUTb[(f >> 23) & 0x01FF] + (u16)((f & 0x007FFFFF) >> uhLUTs[(f >> 23) & 0x01FF]); }
  static inline u16 FloatToSHalf(f23 c) {
    unsigned int f = *((unsigned int *)&c); return shLUTb[(f >> 23) & 0x01FF] + (u16)((f & 0x007FFFFF) >> shLUTs[(f >> 23) & 0x01FF]); }

  // ...ToFloat
  extern unsigned int   uhLUTo[1 << 5];
  extern unsigned int   uhLUTm[1 << 12];
  extern unsigned int   uhLUTe[1 << 5];
  extern unsigned int   shLUTo[1 << 6];
  extern unsigned int   shLUTm[1 << 11];
  extern unsigned int   shLUTe[1 << 6];

  static inline f23 UHalfToFloat(u16 h) {
    unsigned int c = uhLUTm[uhLUTo[h >> 11] + (h & 0x07FF)] + uhLUTe[h >> 11]; return *((float *)&c); }
  static inline f23 SHalfToFloat(u16 h) {
    unsigned int c = shLUTm[shLUTo[h >> 10] + (h & 0x03FF)] + shLUTe[h >> 10]; return *((float *)&c); }

};

#if	SQUISH_USE_ALTIVEC
#include "simd_ve.h"
#elif	SQUISH_USE_SSE
#include "simd_sse.h"
#else
#include "simd_float.h"
#endif

namespace squish {

#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#if	defined(USE_AMP_DEBUG)
typedef	Vec4	float4;
typedef	Col4	int4;
#endif

#if	!defined(SQUISH_USE_COMPUTE)
namespace Concurrency {
namespace vector_math {
#endif

#if	!defined(SQUISH_USE_COMPUTE)
  float4 minimum( float4 left, float4 right ) amp_restricted
  {
    return float4(
      left.x < right.x ? left.x : right.x,
      left.y < right.y ? left.y : right.y,
      left.z < right.z ? left.z : right.z,
      left.w < right.w ? left.w : right.w
    );
  }

  int4 minimum( int4 left, int4 right ) amp_restricted
  {
    return int4(
      left.r < right.r ? left.r : right.r,
      left.g < right.g ? left.g : right.g,
      left.b < right.b ? left.b : right.b,
      left.a < right.a ? left.a : right.a
    );
  }

  float4 maximum( float4 left, float4 right ) amp_restricted
  {
    return float4(
      left.x > right.x ? left.x : right.x,
      left.y > right.y ? left.y : right.y,
      left.z > right.z ? left.z : right.z,
      left.w > right.w ? left.w : right.w
    );
  }

  int4 maximum( int4 left, int4 right ) amp_restricted
  {
    return int4(
      left.r > right.r ? left.r : right.r,
      left.g > right.g ? left.g : right.g,
      left.b > right.b ? left.b : right.b,
      left.a > right.a ? left.a : right.a
    );
  }

  float4 minimax( float4 center, float4 left, float4 right ) amp_restricted
  {
    return minimum(maximum(center, left), right);
  }

  int4 minimax( int4 center, int4 left, int4 right ) amp_restricted
  {
    return minimum(maximum(center, left), right);
  }

  float4 saturate( float4 center ) amp_restricted
  {
    return minimax(center, 0.0f, 1.0f);
  }

  float4 recip( float4 center ) amp_restricted
  {
    return float4(
      1.0f / center.x,
      1.0f / center.y,
      1.0f / center.z,
      1.0f / center.w
    );
  }

  float4 muladd( float4 left, float4 right, float4 offset ) amp_restricted
  {
    return float4(
      left.x * right.x + offset.x,
      left.y * right.y + offset.y,
      left.z * right.z + offset.z,
      left.w * right.w + offset.w
    );
  }

  float4 truncate( float4 center ) amp_restricted
  {
    return float4(
      center.x > 0.0f ? floorf( center.x ) : ceilf( center.x ),
      center.y > 0.0f ? floorf( center.y ) : ceilf( center.y ),
      center.z > 0.0f ? floorf( center.z ) : ceilf( center.z ),
      center.w > 0.0f ? floorf( center.w ) : ceilf( center.w )
    );
  }

  float4 round( float4 center ) amp_restricted
  {
    return float4(
      floorf( center.x + 0.5f ),
      floorf( center.y + 0.5f ),
      floorf( center.z + 0.5f ),
      floorf( center.w + 0.5f )
    );
  }
#endif

  float4 submul( float4 left, float4 right, float4 offset ) amp_restricted
  {
    return muladd(-left, right, offset);
  }

  bool less( float4 left, float4 right ) amp_restricted
  {
    return left.x < right.x
      || left.y < right.y
      || left.z < right.z
      || left.w < right.w;
  }

  bool less( float left, float right ) amp_restricted
  {
    return left < right;
  }

#if	!defined(SQUISH_USE_COMPUTE)
}
}
#endif
#endif

}


#endif // ndef SQUISH_SIMD_H
