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

#ifndef SQUISH_SIMD_FLOAT_H
#define SQUISH_SIMD_FLOAT_H

#if	!defined(USE_COMPUTE)
#include <cmath>
#include <algorithm>
#endif

namespace squish {

#if	!defined(USE_COMPUTE)
#define VEC4_CONST( X ) Vec4( X )

class Vec4
{
public:
  typedef Vec4 const& Arg;

  Vec4() {}

  explicit Vec4( float _s )
    : x( _s ),
    y( _s ),
    z( _s ),
    w( _s )
  {
  }

  Vec4( float _x, float _y, float _z, float _w )
    : x( _x ),
    y( _y ),
    z( _z ),
    w( _w )
  {
  }

  Vec4( Vec3 _v, float _w )
    : x( _v.x ),
    y( _v.y ),
    z( _v.z ),
    w( _w )
  {
  }

  Vec3 GetVec3() const
  {
    return Vec3( x, y, z );
  }

  Vec4 SplatX() const { return Vec4( x ); }
  Vec4 SplatY() const { return Vec4( y ); }
  Vec4 SplatZ() const { return Vec4( z ); }
  Vec4 SplatW() const { return Vec4( w ); }

  operator Vec3() const
  {
    Vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
  }

  Vec4 operator-() const
  {
    Vec4 v;
    v.x = -x;
    v.y = -y;
    v.z = -z;
    v.w = -w;
    return v;
  }

  Vec4& operator=(const float &f)
  {
    x = f;
    y = f;
    z = f;
    w = f;
    return *this;
  }

  Vec4& operator+=( Arg v )
  {
    x += v.x;
    y += v.y;
    z += v.z;
    w += v.w;
    return *this;
  }

  Vec4& operator-=( Arg v )
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    w -= v.w;
    return *this;
  }

  Vec4& operator*=( Arg v )
  {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    w *= v.w;
    return *this;
  }

  friend Vec4 operator+( Vec4::Arg left, Vec4::Arg right  )
  {
    Vec4 copy( left );
    return copy += right;
  }

  friend Vec4 operator-( Vec4::Arg left, Vec4::Arg right  )
  {
    Vec4 copy( left );
    return copy -= right;
  }

  friend Vec4 operator*( Vec4::Arg left, Vec4::Arg right  )
  {
    Vec4 copy( left );
    return copy *= right;
  }

  friend Vec4 operator*( Vec4::Arg left, float right  )
  {
    Vec4 copy( left );
    copy.x *= right;
    copy.y *= right;
    copy.z *= right;
    copy.w *= right;
    return copy;
  }

  //! Returns a*b + c
  friend Vec4 MultiplyAdd( Vec4::Arg a, Vec4::Arg b, Vec4::Arg c )
  {
    return a*b + c;
  }

  //! Returns -( a*b - c )
  friend Vec4 NegativeMultiplySubtract( Vec4::Arg a, Vec4::Arg b, Vec4::Arg c )
  {
    return c - a*b;
  }

  friend Vec4 Reciprocal( Vec4::Arg v )
  {
    return Vec4(
      1.0f/v.x,
      1.0f/v.y,
      1.0f/v.z,
      1.0f/v.w
      );
  }

  friend Vec4 Min( Vec4::Arg left, Vec4::Arg right )
  {
    return Vec4(
      std::min<float>( left.x, right.x ),
      std::min<float>( left.y, right.y ),
      std::min<float>( left.z, right.z ),
      std::min<float>( left.w, right.w )
      );
  }

  friend Vec4 Max( Vec4::Arg left, Vec4::Arg right )
  {
    return Vec4(
      std::max<float>( left.x, right.x ),
      std::max<float>( left.y, right.y ),
      std::max<float>( left.z, right.z ),
      std::max<float>( left.w, right.w )
      );
  }

  friend Vec4 Truncate( Vec4::Arg v )
  {
    return Vec4(
      v.x > 0.0f ? std::floor( v.x ) : std::ceil( v.x ),
      v.y > 0.0f ? std::floor( v.y ) : std::ceil( v.y ),
      v.z > 0.0f ? std::floor( v.z ) : std::ceil( v.z ),
      v.w > 0.0f ? std::floor( v.w ) : std::ceil( v.w )
      );
  }

  friend bool CompareAnyLessThan( Vec4::Arg left, Vec4::Arg right )
  {
    return left.x < right.x
      || left.y < right.y
      || left.z < right.z
      || left.w < right.w;
  }

#if	!defined(USE_AMP)
private:
#endif
  float x;
  float y;
  float z;
  float w;
};

class Col4
{
public:
  typedef Col4 const& Arg;

  Col4() {}

  explicit Col4( int _s )
    : r( _s ),
    g( _s ),
    b( _s ),
    a( _s )
  {
  }

  Col4( int _r, int _g, int _b, int _a )
    : r( _r ),
    g( _g ),
    b( _b ),
    a( _a )
  {
  }

  Col4( Col3 _v, int _w )
    : r( _v.r ),
    g( _v.g ),
    b( _v.b ),
    a( _w )
  {
  }

  Col3 GetCol3() const
  {
    return Col3( r, g, b );
  }

  Col4 SplatX() const { return Col4( r ); }
  Col4 SplatY() const { return Col4( g ); }
  Col4 SplatZ() const { return Col4( b ); }
  Col4 SplatW() const { return Col4( a ); }

  operator Col3() const
  {
    Col3 v;
    v.r = r;
    v.g = g;
    v.b = b;
    return v;
  }

  Col4 operator-() const
  {
    Col4 v;
    v.r = -r;
    v.g = -g;
    v.b = -b;
    v.a = -a;
    return v;
  }

  Col4& operator=(const int &f)
  {
    r = f;
    g = f;
    b = f;
    a = f;
    return *this;
  }

  Col4& operator+=( Arg v )
  {
    r += v.r;
    g += v.g;
    b += v.b;
    a += v.a;
    return *this;
  }

  Col4& operator-=( Arg v )
  {
    r -= v.r;
    g -= v.g;
    b -= v.b;
    a -= v.a;
    return *this;
  }

  Col4& operator*=( Arg v )
  {
    r *= v.r;
    g *= v.g;
    b *= v.b;
    a *= v.a;
    return *this;
  }

  friend Col4 operator+( Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy += right;
  }

  friend Col4 operator-( Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy -= right;
  }

  friend Col4 operator*( Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy *= right;
  }

  friend Col4 operator*( Col4::Arg left, int right  )
  {
    Col4 copy( left );
    copy.r *= right;
    copy.g *= right;
    copy.b *= right;
    copy.a *= right;
    return copy;
  }

  //! Returns a*b + c
  friend Col4 MultiplyAdd( Col4::Arg a, Col4::Arg b, Col4::Arg c )
  {
    return a*b + c;
  }

  //! Returns -( a*b - c )
  friend Col4 NegativeMultiplySubtract( Col4::Arg a, Col4::Arg b, Col4::Arg c )
  {
    return c - a*b;
  }

  friend Col4 Min( Col4::Arg left, Col4::Arg right )
  {
    return Col4(
      std::min<int>( left.r, right.r ),
      std::min<int>( left.g, right.g ),
      std::min<int>( left.b, right.b ),
      std::min<int>( left.a, right.a )
      );
  }

  friend Col4 Max( Col4::Arg left, Col4::Arg right )
  {
    return Col4(
      std::max<int>( left.r, right.r ),
      std::max<int>( left.g, right.g ),
      std::max<int>( left.b, right.b ),
      std::max<int>( left.a, right.a )
      );
  }

  friend bool CompareAnyLessThan( Col4::Arg left, Col4::Arg right )
  {
    return left.r < right.r
      || left.g < right.g
      || left.b < right.b
      || left.a < right.a;
  }

#if	!defined(USE_AMP)
private:
#endif
  int r;
  int g;
  int b;
  int a;
};
#endif

} // namespace squish

#endif // ndef SQUISH_SIMD_FLOAT_H

