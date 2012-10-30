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

#pragma warning(disable: 4127)
#pragma warning(disable: 4293)	// line 459
#pragma warning(disable: 4505)	// line 954

#if	!defined(SQUISH_USE_COMPUTE)
#include <cmath>
#include <algorithm>
#endif

namespace squish {

class Col3
{
public:
  typedef Col3 const& Arg;
  typedef Col3 & aArg;

  Col3()
  {
  }

  explicit Col3( int _s )
  {
    x = _s;
    y = _s;
    z = _s;
  }

  Col3( int _x, int _y, int _z )
  {
    x = _x;
    y = _y;
    z = _z;
  }

  int X() const { return x; }
  int Y() const { return y; }
  int Z() const { return z; }

  Col3 operator-() const
  {
    return Col3( -x, -y, -z );
  }

  Col3& operator+=(Arg v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  Col3& operator-=(Arg v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  Col3& operator*=(Arg v)
  {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
  }

  Col3& operator*=( int s )
  {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }

  Col3& operator/=(Arg v)
  {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return *this;
  }

  Col3& operator/=( int s )
  {
    int t = s;
    x /= t;
    y /= t;
    z /= t;
    return *this;
  }

  friend Col3 operator+(Arg left, Arg right)
  {
    Col3 copy( left );
    return copy += right;
  }

  friend Col3 operator-(Arg left, Arg right)
  {
    Col3 copy( left );
    return copy -= right;
  }

  friend Col3 operator*(Arg left, Arg right)
  {
    Col3 copy( left );
    return copy *= right;
  }

  friend Col3 operator*(Arg left, int right )
  {
    Col3 copy( left );
    return copy *= right;
  }

  friend Col3 operator*( int left, Arg right)
  {
    Col3 copy( right );
    return copy *= left;
  }

  friend Col3 operator/(Arg left, Arg right)
  {
    Col3 copy( left );
    return copy /= right;
  }

  friend Col3 operator/(Arg left, int right )
  {
    Col3 copy( left );
    return copy /= right;
  }
  
  friend Col3 Dot(Arg left, Arg right)
  {
    return Col3(left.x*right.x + left.y*right.y + left.z*right.z);
  }

  friend void Dot(Arg left, Arg right, int *r )
  {
    *r = left.x*right.x + left.y*right.y + left.z*right.z;
  }

  friend Col3 Min(Arg left, Arg right)
  {
    return Col3(
      std::min<int>( left.x, right.x ),
      std::min<int>( left.y, right.y ),
      std::min<int>( left.z, right.z )
      );
  }

  friend Col3 Max(Arg left, Arg right)
  {
    return Col3(
      std::max<int>( left.x, right.x ),
      std::max<int>( left.y, right.y ),
      std::max<int>( left.z, right.z )
      );
  }

  friend class Col4;
  friend class Vec3;

#if	!defined(SQUISH_USE_AMP)
private:
#endif
  union { int x; int r; };
  union { int y; int g; };
  union { int z; int b; };
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

  explicit Col4( float _s )
    : r( (int)_s ),
      g( (int)_s ),
      b( (int)_s ),
      a( (int)_s )
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

  explicit Col4( unsigned int _s )
    : r( _s ),
      g( _s ),
      b( _s ),
      a( _s )
  {
  }

  explicit Col4( const unsigned int (&_rgba)[4] )
    : r( _rgba[0] ),
      g( _rgba[1] ),
      b( _rgba[2] ),
      a( _rgba[3] )
  {
  }

  explicit Col4( u8 const *_source )
    : r( _source[0] ),
      g( _source[1] ),
      b( _source[2] ),
      a( _source[3] )
  {
  }

  Col3 GetCol3() const
  {
    return Col3( r, g, b );
  }

  int GetLong() const
  {
    return r;
  }

  Col4 SetLong( int v ) const
  {
    return Col4 ( v, 0, 0, 0 );
  }

  int R() const { return r; }
  int G() const { return g; }
  int B() const { return b; }
  int A() const { return a; }

  Col4 SplatR() const { return Col4( r ); }
  Col4 SplatG() const { return Col4( g ); }
  Col4 SplatB() const { return Col4( b ); }
  Col4 SplatA() const { return Col4( a ); }

  template<const int inv>
  void SetRGBA( int _r, int _g, int _b, int _a ) {
    r = (inv ? inv - _r : _r);
    g = (inv ? inv - _g : _g);
    b = (inv ? inv - _b : _b);
    a = (inv ? inv - _a : _a);
  }

  template<const int inv>
  void SetRGBApow2( int _r, int _g, int _b, int _a ) {
    r = 1 << (inv ? inv - _r : _r);
    g = 1 << (inv ? inv - _g : _g);
    b = 1 << (inv ? inv - _b : _b);
    a = 1 << (inv ? inv - _a : _a);
  }

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

  Col4& operator&=(Arg v)
  {
    r &= v.r;
    g &= v.g;
    b &= v.b;
    a &= v.a;

    return *this;
  }

  Col4& operator^=(Arg v)
  {
    r ^= v.r;
    g ^= v.g;
    b ^= v.b;
    a ^= v.a;

    return *this;
  }

  Col4& operator|=(Arg v)
  {
    r |= v.r;
    g |= v.g;
    b |= v.b;
    a |= v.a;

    return *this;
  }

  Col4& operator>>=( const int n )
  {
    r >>= n;
    g >>= n;
    b >>= n;
    a >>= n;

    return *this;
  }

  Col4& operator<<=( const int n )
  {
    r <<= n;
    g <<= n;
    b <<= n;
    a <<= n;

    return *this;
  }

  Col4& operator+=(Arg v)
  {
    r += v.r;
    g += v.g;
    b += v.b;
    a += v.a;

    return *this;
  }

  Col4& operator-=(Arg v)
  {
    r -= v.r;
    g -= v.g;
    b -= v.b;
    a -= v.a;

    return *this;
  }

  Col4& operator*=(Arg v)
  {
    r *= v.r;
    g *= v.g;
    b *= v.b;
    a *= v.a;

    return *this;
  }

  friend Col4 operator&( Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy &= right;
  }

  friend Col4 operator^( Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy ^= right;
  }

  friend Col4 operator|( Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy |= right;
  }

  friend Col4 operator>>( Col4::Arg left, int n  )
  {
    Col4 copy( left );
    return copy >>= n;
  }

  friend Col4 operator<<( Col4::Arg left, int n  )
  {
    Col4 copy( left );
    return copy <<= n;
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

  friend Col4 ShiftLeft( Col4::Arg a, int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    if (n >= 64) {
      bb[1] = (ba[0] << (n - 64));
      bb[0] = 0;
    }
    else {
      bb[1] = (ba[1] << (n)) | (ba[0] >> (64 - n));
      bb[0] = (ba[0] << (n));
    }

    return b;
  }

  template<const int n>
  friend Col4 ShiftLeft( Col4::Arg a )
  {
    return ShiftLeft( a, n );
  }

  friend Col4 ShiftRight( Col4::Arg a, int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    if (n >= 64) {
      bb[0] = (ba[1] >> (n - 64));
      bb[1] = 0;
    }
    else {
      bb[0] = (ba[0] >> (n)) | (ba[1] << (64 - n));
      bb[1] = (ba[1] >> (n));
    }

    return b;
  }

  template<const int n>
  friend Col4 ShiftRight( Col4::Arg a )
  {
    return ShiftRight( a, n );
  }

  template<const int n>
  friend Col4 ShiftRightHalf( Col4::Arg a )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] >> (n));
    bb[1] = (ba[1] >> (n));

    return b;
  }

  friend Col4 ShiftRightHalf( Col4::Arg a, int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] >> (n));
    bb[1] = (ba[1] >> (n));

    return b;
  }

  friend Col4 ShiftRightHalf( Col4::Arg a, Col4::Arg b )
  {
    Col4 c;

    unsigned__int64 *bc = (unsigned__int64 *)&c.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bc[0] = (ba[0] >> (b.r));
    bc[1] = (ba[1] >> (b.r));

    return c;
  }

  template<const int n>
  friend Col4 ShiftLeftHalf( Col4::Arg a )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] << (n));
    bb[1] = (ba[1] << (n));

    return b;
  }

  friend Col4 ShiftLeftHalf( Col4::Arg a, const int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] << (n));
    bb[1] = (ba[1] << (n));

    return b;
  }

  template<const int r, const int g, const int b, const int a>
  friend Col4 ShiftLeftLo( Col4::Arg v )
  {
    Col4 r;

    r.r <<= v.r;
    r.g <<= v.g;
    r.b <<= v.b;
    r.a <<= v.a;

    return r;
  }

  template<const int n, const int p>
  friend Col4 MaskBits( Col4::Arg a )
  {
    if ((p + n) <= 0)
      return Col4(0);
    if ((p + n) >= 64)
      return a;

    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] & (~(0xFFFFFFFFFFFFFFFFLL << (p + n))));
    bb[1] = (ba[1] &    0                                );

    return b;
  }

  friend Col4 MaskBits( Col4::Arg a, const int n, const int p )
  {
    if ((p + n) >= 64)
      return a;

    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] & (~(0xFFFFFFFFFFFFFFFFLL << (p + n))));
    bb[1] = (ba[1] &    0                                );

    return b;
  }

  template<const int n, const int p>
  friend Col4 CopyBits( Col4::Arg left, Col4::Arg right)
  {
    if (!n)
      return left;
    if (!p)
      return MaskBits<n, 0>(right);
    if ((p + n) >= 64)
      return (left) | ShiftLeftHalf<p>(right);

    return MaskBits<p, 0>(left) | MaskBits<n, p>(ShiftLeftHalf<p>(right));
  }

  friend Col4 CopyBits( Col4::Arg left, Col4 right, const int n, const int p )
  {
    return MaskBits(left, p, 0) | MaskBits(ShiftLeftHalf(right, p), n, p);
  }

  template<const int n, const int p>
  friend Col4 ExtrBits( Col4::Arg a )
  {
    if (!n)
      return Col4(0);
    if (!p)
      return MaskBits<n, 0>(a);
    if ((n + p) >= 64)
      return ShiftRightHalf<p>(a);

    return MaskBits<n, 0>(ShiftRightHalf<p>(a));
  }

  friend Col4 ExtrBits( Col4::Arg a, const int n, const int p )
  {
    return MaskBits(ShiftRightHalf(a, p), n, 0);
  }

  template<const int n, const int p>
  friend void ExtrBits( Col4::Arg left, Col4 &right )
  {
    right  = ExtrBits<n, p>( left );
  }

  template<const int n, const int p>
  friend void ConcBits( Col4::Arg left, Col4 &right )
  {
    right  = ShiftLeft<32>( right );
    if (n > 0)
      right |= ExtrBits<n, p>( left );
  }

  template<const int n, const int p>
  friend void ReplBits( Col4::Arg left, Col4 &right )
  {
    if (!n)
      return;
    if ((n < 0)) {
      right = ExtrBits<-n, p>( left );
      right = Col4(right.r, right.r, right.r, right.a);
    }
    else {
      right = ExtrBits< n, p>( left );
      right = Col4(right.r, right.r, right.r, right.r);
    }
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

  template<const int f, const int t>
  friend Col4 Shuffle( Arg a )
  {
    Col4 b = a;

    int *bb = (int *)&b.r;
    int *ba = (int *)&a.r;

    bb[t] = ba[f];

    return b;
  }

  template<const int f, const int t>
  friend Col4 Exchange( Arg a )
  {
    Col4 b = a;

    int *bb = (int *)&b.r;
    int *ba = (int *)&a.r;

    std::swap(bb[t], ba[f]);

    return b;
  }

  friend Col4 HorizontalAdd( Arg a )
  {
    return Col4( a.r + a.g + a.b + a.a );
  }

  friend Col4 HorizontalAdd( Arg a, Arg b )
  {
    return HorizontalAdd( a ) + HorizontalAdd( b );
  }

  friend Col4 HorizontalAddTiny( Arg a )
  {
    return HorizontalAdd( a );
  }

  friend Col4 HorizontalAddTiny( Arg a, Arg b )
  {
    return HorizontalAdd( a, b );
  }

  friend Col4 Dot(Arg left, Arg right)
  {
    return HorizontalAdd( left*right );
  }

  friend Col4 DotTiny(Arg left, Arg right)
  {
    return HorizontalAddTiny( left*right );
  }
  
  friend void Dot(Arg left, Arg right, int *r )
  {
    Col4 res = Dot( left, right );

    *r = res.R();
  }

  friend Col4 Min( Col4::Arg left, Col4::Arg right)
  {
    return Col4(
      std::min<int>( left.r, right.r ),
      std::min<int>( left.g, right.g ),
      std::min<int>( left.b, right.b ),
      std::min<int>( left.a, right.a )
    );
  }

  friend Col4 Max( Col4::Arg left, Col4::Arg right)
  {
    return Col4(
      std::max<int>( left.r, right.r ),
      std::max<int>( left.g, right.g ),
      std::max<int>( left.b, right.b ),
      std::max<int>( left.a, right.a )
    );
  }

  friend bool CompareAnyLessThan( Col4::Arg left, Col4::Arg right)
  {
    return
      left.r < right.r ||
      left.g < right.g ||
      left.b < right.b ||
      left.a < right.a;
  }

  friend bool CompareAllEqualTo( Col4::Arg left, Col4::Arg right)
  {
    return
      left.r == right.r &&
      left.g == right.g &&
      left.b == right.b &&
      left.a == right.a;
  }

  friend Col4 IsNotZero( Col4::Arg v )
  {
    return Col4(
      v.r > 0 ? ~0 : 0,
      v.g > 0 ? ~0 : 0,
      v.b > 0 ? ~0 : 0,
      v.a > 0 ? ~0 : 0
    );
  }

  friend Col4 IsOne( Col4::Arg v )
  {
    return Col4(
      v.r == 0xFF ? ~0 : 0,
      v.g == 0xFF ? ~0 : 0,
      v.b == 0xFF ? ~0 : 0,
      v.a == 0xFF ? ~0 : 0
    );
  }

  friend Col4 TransferA( Col4::Arg left, Col4::Arg right)
  {
    return Col4( left.r, left.g, left.b, right.a );
  }

  friend Col4 KillA( Col4::Arg left)
  {
    return Col4( left.r, left.g, left.b, 0xFF );
  }

  friend void PackBytes( Col4::Arg v, int &loc )
  {
    loc  = v.a; loc <<= 8;
    loc |= v.b; loc <<= 8;
    loc |= v.g; loc <<= 8;
    loc |= v.r; loc <<= 0;
  }

  // clamp the output to [0, 1]
  Col4 Clamp() const {
    Col4 const one (0xFF);
    Col4 const zero(0x00);

    return Min(one, Max(zero, *this));
  }

  friend void LoadAligned( Col4 &a, Col4 &b, Col4::Arg c )
  {
    a.r = c.r;
    a.g = c.g;

    b.b = c.b;
    b.a = c.a;
  }

  friend void LoadAligned( Col4 &a, void const *source )
  {
    a.r = ((int *)source)[0];
    a.g = ((int *)source)[1];

    a.b = ((int *)source)[2];
    a.a = ((int *)source)[3];
  }

  friend void LoadAligned( Col4 &a, Col4 &b, void const *source )
  {
    a.r = ((int *)source)[0];
    a.g = ((int *)source)[1];

    b.b = ((int *)source)[2];
    b.a = ((int *)source)[3];
  }

  friend void LoadUnaligned( Col4 &a, Col4 &b, void const *source )
  {
    a.r = ((int *)source)[0];
    a.g = ((int *)source)[1];

    b.b = ((int *)source)[2];
    b.a = ((int *)source)[3];
  }

  friend void StoreAligned( Col4::Arg a, Col4::Arg b, Col4 &c )
  {
    c.r = a.r;
    c.g = a.g;

    c.b = b.b;
    c.a = b.a;
  }

  friend void StoreAligned( Col4::Arg a, void *destination )
  {
    ((int *)destination)[0] = a.r;
    ((int *)destination)[1] = a.g;

    ((int *)destination)[2] = a.b;
    ((int *)destination)[3] = a.a;
  }

  friend void StoreAligned( Col4::Arg a, Col4::Arg b, void *destination )
  {
    ((int *)destination)[0] = a.r;
    ((int *)destination)[1] = a.g;

    ((int *)destination)[2] = b.b;
    ((int *)destination)[3] = b.a;
  }

  friend void StoreUnaligned( Col4::Arg a, Col4::Arg b, void *destination )
  {
    ((int *)destination)[0] = a.r;
    ((int *)destination)[1] = a.g;

    ((int *)destination)[2] = b.b;
    ((int *)destination)[3] = b.a;
  }

#if	!defined(SQUISH_USE_AMP)
private:
#endif
  int r;
  int g;
  int b;
  int a;
};

// scalar types
typedef	float Scr3;
typedef	float Scr4;

static float Reciprocal(float v) {
  return math::rcp(v);
}

static float ReciprocalSqrt(float v) {
  return math::rcp(math::sqrt(v));
}

static float Min(float a, float b) {
  return std::min(a, b);
}

static float Abs(float v) {
  return std::abs(v);
}

#if	!defined(SQUISH_USE_COMPUTE)
#define VEC4_CONST( X ) Vec4( X )

class Vec3
{
public:
  typedef Vec3 const& Arg;
  typedef Vec3 & aArg;

  Vec3()
  {
  }

  explicit Vec3( float _s )
  {
    x = _s;
    y = _s;
    z = _s;
  }

  Vec3( float _x, float _y, float _z )
  {
    x = _x;
    y = _y;
    z = _z;
  }
  
  Vec3( Vec3 _x, Vec3 _y, Vec3 _z )
  {
    x = _x.x;
    y = _y.x;
    z = _z.x;
  }

  float X() const { return x; }
  float Y() const { return y; }
  float Z() const { return z; }

  Vec3 operator-() const
  {
    return Vec3( -x, -y, -z );
  }

  Vec3& operator+=(Arg v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  Vec3& operator+=(float v)
  {
    x += v;
    y += v;
    z += v;
    return *this;
  }

  Vec3& operator-=(Arg v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  Vec3& operator*=(Arg v)
  {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
  }

  Vec3& operator*=( float s )
  {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }

  Vec3& operator/=(Arg v)
  {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return *this;
  }

  Vec3& operator/=( float s )
  {
    float t = 1.0f/s;
    x *= t;
    y *= t;
    z *= t;
    return *this;
  }
  
  Vec3& operator/=( int s )
  {
    float t = 1.0f/s;
    x *= t;
    y *= t;
    z *= t;
    return *this;
  }
  
  friend int operator<( Vec3::Arg left, Vec3::Arg right  )
  {
    return CompareFirstLessThan(left, right);
  }
	
  friend int operator>( Vec3::Arg left, Vec3::Arg right  )
  {
    return CompareFirstGreaterThan(left, right);
  }

  friend Vec3 operator+(Arg left, Arg right)
  {
    Vec3 copy( left );
    return copy += right;
  }

  friend Vec3 operator+(Arg left, float right )
  {
    Vec3 copy( left );
    return copy += right;
  }

  friend Vec3 operator-(Arg left, Arg right)
  {
    Vec3 copy( left );
    return copy -= right;
  }

  friend Vec3 operator*(Arg left, Arg right)
  {
    Vec3 copy( left );
    return copy *= right;
  }

  friend Vec3 operator*(Arg left, float right )
  {
    Vec3 copy( left );
    return copy *= right;
  }

  friend Vec3 operator*( float left, Arg right)
  {
    Vec3 copy( right );
    return copy *= left;
  }

  friend Vec3 operator/(Arg left, Arg right)
  {
    Vec3 copy( left );
    return copy /= right;
  }

  friend Vec3 operator/(Arg left, float right )
  {
    Vec3 copy( left );
    return copy /= right;
  }
  
  friend Vec3 Reciprocal( Vec3::Arg v )
  {
    return Vec3(
      1.0f / v.x,
      1.0f / v.y,
      1.0f / v.z
    );
  }
  
  friend Vec3 ReciprocalSqrt( Vec3::Arg v )
  {
    return Vec3(
      math::rsqrt(v.x),
      math::rsqrt(v.y),
      math::rsqrt(v.z)
    );
  }

  friend Scr3 HorizontalMax( Arg a )
  {
    return Scr3(
      std::max<float>( std::max<float>( a.x, a.y ), a.z )
    );
  }

  friend Scr3 HorizontalAdd( Arg a )
  {
    return Scr3( a.x + a.y + a.z );
  }

  friend Scr3 HorizontalAdd( Arg a, Arg b )
  {
    return HorizontalAdd( a + b );
  }

  friend Vec3 Normalize(Arg left)
  {
    Vec3 sum = (left * left);
    float rsq = math::rsqrt(sum.x + sum.y + sum.z);

    return left * rsq;
  }

  friend Scr3 Dot(Arg left, Arg right)
  {
    return HorizontalAdd( left*right );
  }

  friend void Dot(Arg left, Arg right, float *r )
  {
    *r = HorizontalAdd( left*right );
  }
  
  friend Vec3 Abs( Vec3::Arg v )
  {
    return Vec3(
      abs( v.x ),
      abs( v.y ),
      abs( v.z )
    );
  }
	
  friend Vec3 Min(Arg left, Arg right)
  {
    return Vec3(
      std::min<float>( left.x, right.x ),
      std::min<float>( left.y, right.y ),
      std::min<float>( left.z, right.z )
    );
  }

  friend Vec3 Max(Arg left, Arg right)
  {
    return Vec3(
      std::max<float>( left.x, right.x ),
      std::max<float>( left.y, right.y ),
      std::max<float>( left.z, right.z )
    );
  }

  // clamp the output to [0, 1]
  Vec3 Clamp() const {
    Vec3 const one (1.0f);
    Vec3 const zero(0.0f);

    return Min(one, Max(zero, *this));
  }

  friend Vec3 Truncate(Arg v)
  {
    return Vec3(
      v.x > 0.0f ? std::floor( v.x ) : std::ceil( v.x ),
      v.y > 0.0f ? std::floor( v.y ) : std::ceil( v.y ),
      v.z > 0.0f ? std::floor( v.z ) : std::ceil( v.z )
    );
  }
  
  friend Vec3 AbsoluteDifference( Vec3::Arg left, Vec3::Arg right )
  {
    return Abs( left - right );
  }

  friend Scr3 SummedAbsoluteDifference( Vec3::Arg left, Vec3::Arg right )
  {
    return HorizontalAdd( AbsoluteDifference( left, right ) );
  }

  friend Scr3 MaximumAbsoluteDifference( Vec3::Arg left, Vec3::Arg right )
  {
    return HorizontalMax( AbsoluteDifference( left, right ) );
  }

  friend int CompareEqualTo(Vec3::Arg left, Vec3::Arg right)
  {
    return
      (left.x == right.x ? 0x1 : 0x0) |
      (left.y == right.y ? 0x2 : 0x0) |
      (left.z == right.z ? 0x4 : 0x0);
  }

  friend bool CompareAnyLessThan(Vec3::Arg left, Vec3::Arg right)
  {
    return
      left.x < right.x ||
      left.y < right.y ||
      left.z < right.z;
  }
  
  friend bool CompareAnyGreaterThan(Vec3::Arg left, Vec3::Arg right)
  {
    return
      left.x > right.x ||
      left.y > right.y ||
      left.z > right.z;
  }

  friend int CompareFirstLessThan( Vec3::Arg left, Vec3::Arg right)
  {
    return left.x < right.x;
  }

  friend int CompareFirstGreaterThan( Vec3::Arg left, Vec3::Arg right)
  {
    return left.x > right.x;
  }

  friend class Col3;
  friend class Vec4;

#if	!defined(SQUISH_USE_AMP) && !defined(SQUISH_USE_COMPUTE)
private:
#endif
  union { float x; float r; };
  union { float y; float g; };
  union { float z; float b; };
};

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

  explicit Vec4( int _s )
    : x( (float) _s ),
      y( (float) _s ),
      z( (float) _s ),
      w( (float) _s )
  {
  }

  Vec4( float _x, float _y, float _z, float _w )
    : x( _x ),
      y( _y ),
      z( _z ),
      w( _w )
  {
  }

  Vec4( float _x, float _y, float _z )
    : x( _x ),
      y( _y ),
      z( _z ),
      w( 0.0f )
  {
  }

  Vec4( float _x, float _y )
    : x( _x ),
      y( _y ),
      z( 0.0f ),
      w( 0.0f )
  {
  }

  Vec4( Vec4 _x, Vec4 _y, Vec4 _z, Vec4 _w )
    : x( _x.x ),
      y( _y.x ),
      z( _z.x ),
      w( _w.x )
  {
  }

  Vec4( Vec4 _x, Vec4 _y, Vec4 _z )
    : x( _x.x ),
      y( _y.x ),
      z( _z.x ),
      w( 0.0f )
  {
  }

  Vec4( Vec4 _x, Vec4 _y )
    : x( _x.x ),
      y( _y.x ),
      z( 0.0f ),
      w( 0.0f )
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

  float X() const { return x; }
  float Y() const { return y; }
  float Z() const { return z; }
  float W() const { return w; }

  float &GetX() { return x; }
  float &GetY() { return y; }
  float &GetZ() { return z; }
  float &GetW() { return w; }
  // let the compiler figure this one out, probably spills to memory
  float &GetO(int o) { return ((float *)this)[o]; }
  
  Vec4 Swap  () const { return Vec4( z, w, x, y ); }
  Vec4 SplatX() const { return Vec4( x ); }
  Vec4 SplatY() const { return Vec4( y ); }
  Vec4 SplatZ() const { return Vec4( z ); }
  Vec4 SplatW() const { return Vec4( w ); }

  template<const int inv>
  void SetXYZW( int _x, int _y, int _z, int _w ) {
    x = (float)(inv ? inv - _x : _x);
    y = (float)(inv ? inv - _y : _y);
    z = (float)(inv ? inv - _z : _z);
    w = (float)(inv ? inv - _w : _w);
  }

  template<const int inv>
  void SetXYZWpow2( int _x, int _y, int _z, int _w ) {
    x = (float)(1 << (inv ? inv - _x : _x));
    y = (float)(1 << (inv ? inv - _y : _y));
    z = (float)(1 << (inv ? inv - _z : _z));
    w = (float)(1 << (inv ? inv - _w : _w));
  }

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

  Vec4& operator+=(Arg v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
    w += v.w;

    return *this;
  }

  Vec4& operator-=(Arg v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    w -= v.w;

    return *this;
  }

  Vec4& operator*=(Arg v)
  {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    w *= v.w;

    return *this;
  }
  
  Vec4& operator*=(float v)
  {
    x *= v;
    y *= v;
    z *= v;
    w *= v;

    return *this;
  }

  Vec4& operator/=(float v)
  {
    *this *= Reciprocal( Vec4( v ) );
    return *this;
  }
  
  Vec4& operator/=(int v)
  {
    *this *= Reciprocal( Vec4( v ) );
    return *this;
  }
  
  friend int operator<(Vec4::Arg left, Vec4::Arg right  )
  {
    return CompareFirstLessThan(left, right);
  }
	
  friend int operator>(Vec4::Arg left, Vec4::Arg right  )
  {
    return CompareFirstGreaterThan(left, right);
  }

  friend Vec4 operator&(Vec4::Arg left, Vec4::Arg right  )
  {
    Vec4 copy( left );

    *((int *)&copy.x) &= *((int *)&right.x);
    *((int *)&copy.y) &= *((int *)&right.y);
    *((int *)&copy.z) &= *((int *)&right.z);
    *((int *)&copy.w) &= *((int *)&right.w);

    return copy;
  }

  friend Vec4 operator+(Vec4::Arg left, Vec4::Arg right  )
  {
    Vec4 copy( left );
    return copy += right;
  }

  friend Vec4 operator-(Vec4::Arg left, Vec4::Arg right  )
  {
    Vec4 copy( left );
    return copy -= right;
  }

  friend Vec4 operator*(Vec4::Arg left, Vec4::Arg right  )
  {
    Vec4 copy( left );
    return copy *= right;
  }

  friend Vec4 operator*(Vec4::Arg left, float right  )
  {
    Vec4 copy( left );

    copy.x *= right;
    copy.y *= right;
    copy.z *= right;
    copy.w *= right;

    return copy;
  }

  friend Vec4 operator*( float left, Vec4::Arg right  )
  {
    Vec4 copy( right );
    return copy * left;
  }

  friend Vec4 operator/(Vec4::Arg left, float right  )
  {
    Vec4 copy( left );
    copy /= right;
    return copy;
  }

  friend Vec4 operator*(Vec4::Arg left, int right  )
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

  template<const int f, const int t>
  friend Vec4 Shuffle( Arg a )
  {
    Vec4 b = a;

    float *bb = (float *)&b.x;
    float *ba = (float *)&a.x;

    bb[t] = ba[f];

    return b;
  }

  template<const int f, const int t>
  friend Vec4 Exchange( Arg a )
  {
    Vec4 b = a;

    float *bb = (float *)&b.x;
    float *ba = (float *)&a.x;

    std::swap(bb[t], ba[f]);

    return b;
  }

  friend Vec4 Reciprocal(Vec4::Arg v)
  {
    return Vec4(
      1.0f / v.x,
      1.0f / v.y,
      1.0f / v.z,
      1.0f / v.w
    );
  }
  
  friend Vec4 ReciprocalSqrt(Vec4::Arg v)
  {
    return Vec4(
      math::rsqrt(v.x),
      math::rsqrt(v.y),
      math::rsqrt(v.z),
      math::rsqrt(v.w)
    );
  }

  friend Scr4 HorizontalMin( Arg a )
  {
    return Scr4(
      std::min<float>( std::min<float>( a.x, a.y ), std::min<float>( a.z, a.w ) )
    );
  }
  
  friend Scr4 HorizontalMax( Arg a )
  {
    return Scr4(
      std::max<float>( std::max<float>( a.x, a.y ), std::max<float>( a.z, a.w ) )
    );
  }

  friend Scr4 HorizontalAdd( Arg a )
  {
    return Scr4(a.x + a.y + a.z + a.w);
  }

  friend Scr4 HorizontalAdd( Arg a, Arg b )
  {
    return HorizontalAdd(a + b);
  }
  
  friend Scr4 Length(Arg left)
  {
    Vec4 sum = (left * left);
    return math::sqrt(sum.x + sum.y + sum.z + sum.w);
  }
  
  friend Scr4 ReciprocalLength(Arg left)
  {
    Vec4 sum = (left * left);
    return math::rsqrt(sum.x + sum.y + sum.z + sum.w);
  }
  
  friend Vec4 Normalize(Arg left)
  {
    return left * ReciprocalLength(left);
  }

  friend Scr4 Dot(Arg left, Arg right)
  {
    return HorizontalAdd(left * right);
  }

  friend void Dot(Arg left, Arg right, float *r )
  {
    *r = HorizontalAdd(left * right);
  }
  
  friend Vec4 Abs( Vec4::Arg v )
  {
    return Vec4(
      abs( v.x ),
      abs( v.y ),
      abs( v.z ),
      abs( v.w )
    );
  }
	
  friend Vec4 Min(Vec4::Arg left, Vec4::Arg right)
  {
    return Vec4(
      std::min<float>( left.x, right.x ),
      std::min<float>( left.y, right.y ),
      std::min<float>( left.z, right.z ),
      std::min<float>( left.w, right.w )
    );
  }

  friend Vec4 Max(Vec4::Arg left, Vec4::Arg right)
  {
    return Vec4(
      std::max<float>( left.x, right.x ),
      std::max<float>( left.y, right.y ),
      std::max<float>( left.z, right.z ),
      std::max<float>( left.w, right.w )
    );
  }

  // clamp the output to [0, 1]
  Vec4 Clamp() const {
    Vec4 const one (1.0f);
    Vec4 const zero(0.0f);

    return Min(one, Max(zero, *this));
  }

  template<const bool round>
  friend Col4 FloatToInt(Vec4::Arg v)
  {
    return Col4(
      (int)(v.x > 0.0f ? std::floor( v.x + (round ? 0.5f : 0.0f) ) : std::ceil( v.x - (round ? 0.5f : 0.0f) )),
      (int)(v.y > 0.0f ? std::floor( v.y + (round ? 0.5f : 0.0f) ) : std::ceil( v.y - (round ? 0.5f : 0.0f) )),
      (int)(v.z > 0.0f ? std::floor( v.z + (round ? 0.5f : 0.0f) ) : std::ceil( v.z - (round ? 0.5f : 0.0f) )),
      (int)(v.w > 0.0f ? std::floor( v.w + (round ? 0.5f : 0.0f) ) : std::ceil( v.w - (round ? 0.5f : 0.0f) ))
    );
  }

  friend Vec4 Truncate(Vec4::Arg v)
  {
    return Vec4(
      v.x > 0.0f ? std::floor( v.x ) : std::ceil( v.x ),
      v.y > 0.0f ? std::floor( v.y ) : std::ceil( v.y ),
      v.z > 0.0f ? std::floor( v.z ) : std::ceil( v.z ),
      v.w > 0.0f ? std::floor( v.w ) : std::ceil( v.w )
    );
  }
  
  friend Vec4 AbsoluteDifference( Vec4::Arg left, Vec4::Arg right )
  {
    return Abs( left - right );
  }

  friend Scr4 SummedAbsoluteDifference( Vec4::Arg left, Vec4::Arg right )
  {
    return HorizontalAdd( AbsoluteDifference( left, right ) );
  }

  friend Scr4 MaximumAbsoluteDifference( Vec4::Arg left, Vec4::Arg right )
  {
    return HorizontalMax( AbsoluteDifference( left, right ) );
  }

  friend int CompareEqualTo(Vec4::Arg left, Vec4::Arg right)
  {
    return
      (left.x == right.x ? 0x1 : 0x0) |
      (left.y == right.y ? 0x2 : 0x0) |
      (left.z == right.z ? 0x4 : 0x0) |
      (left.w == right.w ? 0x8 : 0x0);
  }

  friend bool CompareAnyLessThan(Vec4::Arg left, Vec4::Arg right)
  {
    return
      left.x < right.x ||
      left.y < right.y ||
      left.z < right.z ||
      left.w < right.w;
  }
  
  friend bool CompareAnyGreaterThan(Vec4::Arg left, Vec4::Arg right)
  {
    return
      left.x > right.x ||
      left.y > right.y ||
      left.z > right.z ||
      left.w > right.w;
  }

  friend int CompareFirstLessThan(Vec4::Arg left, Vec4::Arg right)
  {
    return left.x < right.x;
  }

  friend int CompareFirstGreaterThan(Vec4::Arg left, Vec4::Arg right)
  {
    return left.x > right.x;
  }

  Vec4 IsOne( ) const
  {
    Vec4 m;

    *((int *)&m.x) = (x == 1.0f ? ~0 : 0);
    *((int *)&m.y) = (y == 1.0f ? ~0 : 0);
    *((int *)&m.z) = (z == 1.0f ? ~0 : 0);
    *((int *)&m.w) = (w == 1.0f ? ~0 : 0);

    return m;
  }

  Vec4 IsNotOne( ) const
  {
    Vec4 m;

    *((int *)&m.x) = (x != 1.0f ? ~0 : 0);
    *((int *)&m.y) = (y != 1.0f ? ~0 : 0);
    *((int *)&m.z) = (z != 1.0f ? ~0 : 0);
    *((int *)&m.w) = (w != 1.0f ? ~0 : 0);

    return m;
  }

  friend Vec4 TransferZW(Vec4::Arg left, Vec4::Arg right)
  {
    return Vec4( left.x, left.y, right.z, right.w );
  }

  friend Vec4 TransferW(Vec4::Arg left, Vec4::Arg right)
  {
    return Vec4( left.x, left.y, left.z, right.w );
  }

  friend Vec4 KillW(Vec4::Arg left)
  {
    return Vec4( left.x, left.y, left.z, 0.0f );
  }

  friend Vec4 OnlyW(Vec4::Arg left)
  {
    return Vec4( 0.0f, 0.0f, 0.0f, left.w );
  }

  void SwapXYZW( Vec4 &with )
  {
    std::swap(x, with.x);
    std::swap(y, with.y);
    std::swap(z, with.z);
    std::swap(w, with.w);
  }

  void SwapXYZ ( Vec4 &with )
  {
    std::swap(x, with.x);
    std::swap(y, with.y);
    std::swap(z, with.z);
  }

  void SwapW   ( Vec4 &with )
  {
    std::swap(w, with.w);
  }

#if	!defined(SQUISH_USE_AMP)
private:
#endif
  float x;
  float y;
  float z;
  float w;
};

#if	!defined(SQUISH_USE_PRE)
inline Scr3 LengthSquared( Vec3::Arg v )
{
  return Dot( v, v );
}

inline void LengthSquared( Vec3::Arg v , float *r )
{
  Dot( v, v, r );
}

inline Scr3 LengthSquared(Vec4::Arg v)
{
  return Dot( v, v );
}

inline void LengthSquared( Vec4::Arg v , float *r )
{
  Dot( v, v, r );
}
#endif
#endif // ndef SQUISH_USE_COMPUTE

} // namespace squish

#endif // ndef SQUISH_SIMD_FLOAT_H
