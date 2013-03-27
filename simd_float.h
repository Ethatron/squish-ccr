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
#pragma warning(disable: 4172)	// line 1130

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

  explicit Col3(int _s)
  {
    x = _s;
    y = _s;
    z = _s;
  }

  Col3(int _x, int _y, int _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }
  
  Col3(int _x, int _y)
  {
    x = _x;
    y = _y;
    z = 0;
  }

  int R() const { return x; }
  int G() const { return y; }
  int B() const { return z; }

  Col3 operator-() const
  {
    return Col3( -x, -y, -z );
  }

  Col3& operator>>=( const int n )
  {
    x >>= n;
    y >>= n;
    z >>= n;
    return *this;
  }

  Col3& operator<<=( const int n )
  {
    x <<= n;
    y <<= n;
    z <<= n;
    return *this;
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

  friend Col3 operator>>( Col3::Arg left, int n  )
  {
    Col3 copy( left );
    return copy >>= n;
  }

  friend Col3 operator<<( Col3::Arg left, int n  )
  {
    Col3 copy( left );
    return copy <<= n;
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

  friend void PackBytes(Col3::Arg v, int &loc )
  {
    loc  = v.b; loc <<= 8;
    loc += v.g; loc <<= 8;
    loc += v.r; loc <<= 0;
  }

  // clamp the output to [0, 1]
  Col3 Clamp() const {
    Col3 const one (0xFF);
    Col3 const zero(0x00);

    return Min(one, Max(zero, *this));
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

  explicit Col4(int _s)
    : r(_s),
      g(_s),
      b(_s),
      a(_s)
  {
  }

  explicit Col4(float _s)
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

  Col4( int _r, int _g, int _b )
    : r( _r ),
      g( _g ),
      b( _b ),
      a( 0 )
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
    : r(_s),
      g(_s),
      b(_s),
      a(_s)
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
  
  int GetM8() const
  {
    return
      (r & 0x000000FF ? 0x0001 : 0x0000) +
      (r & 0x0000FF00 ? 0x0002 : 0x0000) +
      (r & 0x00FF0000 ? 0x0004 : 0x0000) +
      (r & 0xFF000000 ? 0x0008 : 0x0000) +
      (g & 0x000000FF ? 0x0010 : 0x0000) +
      (g & 0x0000FF00 ? 0x0020 : 0x0000) +
      (g & 0x00FF0000 ? 0x0040 : 0x0000) +
      (g & 0xFF000000 ? 0x0080 : 0x0000) +
      (b & 0x000000FF ? 0x0100 : 0x0000) +
      (b & 0x0000FF00 ? 0x0200 : 0x0000) +
      (b & 0x00FF0000 ? 0x0400 : 0x0000) +
      (b & 0xFF000000 ? 0x0800 : 0x0000) +
      (a & 0x000000FF ? 0x1000 : 0x0000) +
      (a & 0x0000FF00 ? 0x2000 : 0x0000) +
      (a & 0x00FF0000 ? 0x4000 : 0x0000) +
      (a & 0xFF000000 ? 0x8000 : 0x0000);
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
  void SetRGBA(int _r, int _g, int _b, int _a) {
    r = (inv ? inv - _r : _r);
    g = (inv ? inv - _g : _g);
    b = (inv ? inv - _b : _b);
    a = (inv ? inv - _a : _a);
  }

  template<const int inv>
  void SetRGBApow2(int _r, int _g, int _b, int _a) {
    r = 1 << (inv ? inv - _r : _r);
    g = 1 << (inv ? inv - _g : _g);
    b = 1 << (inv ? inv - _b : _b);
    a = 1 << (inv ? inv - _a : _a);
  }
  
  template<const int inv>
  void SetRGBApow2(int _v) {
    r = 
    g = 
    b = 
    a = 1 << (inv ? inv - _v : _v);
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

  friend Col4 operator&(Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy &= right;
  }

  friend Col4 operator^(Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy ^= right;
  }

  friend Col4 operator|(Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy |= right;
  }

  friend Col4 operator>>(Col4::Arg left, int n  )
  {
    Col4 copy( left );
    return copy >>= n;
  }

  friend Col4 operator<<(Col4::Arg left, int n  )
  {
    Col4 copy( left );
    return copy <<= n;
  }

  friend Col4 operator+(Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy += right;
  }

  friend Col4 operator-(Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy -= right;
  }

  friend Col4 operator*(Col4::Arg left, Col4::Arg right  )
  {
    Col4 copy( left );
    return copy *= right;
  }

  friend Col4 operator*(Col4::Arg left, int right  )
  {
    Col4 copy( left );

    copy.r *= right;
    copy.g *= right;
    copy.b *= right;
    copy.a *= right;

    return copy;
  }

  friend Col4 ShiftLeft(Col4::Arg a, int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    if (n >= 64) {
      bb[1] = (ba[0] << (n - 64));
      bb[0] = 0;
    }
    else {
      bb[1] = (ba[1] << (n)) + (ba[0] >> (64 - n));
      bb[0] = (ba[0] << (n));
    }

    return b;
  }

  template<const int n>
  friend Col4 ShiftLeft(Col4::Arg a )
  {
    return ShiftLeft( a, n );
  }

  friend Col4 ShiftRight(Col4::Arg a, int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    if (n >= 64) {
      bb[0] = (ba[1] >> (n - 64));
      bb[1] = 0;
    }
    else {
      bb[0] = (ba[0] >> (n)) + (ba[1] << (64 - n));
      bb[1] = (ba[1] >> (n));
    }

    return b;
  }

  template<const int n>
  friend Col4 ShiftRight(Col4::Arg a )
  {
    return ShiftRight( a, n );
  }

  template<const int n>
  friend Col4 ShiftRightHalf(Col4::Arg a )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] >> (n));
    bb[1] = (ba[1] >> (n));

    return b;
  }

  friend Col4 ShiftRightHalf(Col4::Arg a, int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] >> (n));
    bb[1] = (ba[1] >> (n));

    return b;
  }

  friend Col4 ShiftRightHalf(Col4::Arg a, Col4::Arg b )
  {
    Col4 c;

    unsigned__int64 *bc = (unsigned__int64 *)&c.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bc[0] = (ba[0] >> (b.r));
    bc[1] = (ba[1] >> (b.r));

    return c;
  }

  template<const int n>
  friend Col4 ShiftLeftHalf(Col4::Arg a )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] << (n));
    bb[1] = (ba[1] << (n));

    return b;
  }

  friend Col4 ShiftLeftHalf(Col4::Arg a, const int n )
  {
    Col4 b;

    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] << (n));
    bb[1] = (ba[1] << (n));

    return b;
  }

  template<const int r, const int g, const int b, const int a>
  friend Col4 ShiftLeftLo(Col4::Arg v )
  {
    Col4 r;

    r.r <<= v.r;
    r.g <<= v.g;
    r.b <<= v.b;
    r.a <<= v.a;

    return r;
  }

  template<const int n, const int p>
  friend Col4 MaskBits(Col4::Arg a )
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

  friend Col4 MaskBits(Col4::Arg a, const int n, const int p )
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
  friend Col4 CopyBits(Col4::Arg left, Col4::Arg right)
  {
    if (!n)
      return left;
    if (!p)
      return MaskBits<n, 0>(right);
    if ((p + n) >= 64)
      return (left) + ShiftLeftHalf<p>(right);

    return MaskBits<p, 0>(left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
  }

  friend Col4 CopyBits(Col4::Arg left, Col4 right, const int n, const int p )
  {
    return MaskBits(left, p, 0) + MaskBits(ShiftLeftHalf(right, p), n, p);
  }

  template<const int n, const int p>
  friend Col4 KillBits(Col4::Arg a )
  {
    if ((p + n) <= 0)
      return Col4(0);
    if ((p + n) >= 64)
      return a;

    Col4 b;

    unsigned__int64 base1 =  (0xFFFFFFFFFFFFFFFFULL << (     (p + 0) & 63));
    unsigned__int64 base2 =  (0xFFFFFFFFFFFFFFFFULL >> (64 - (p + n) & 63));
    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] & (base1 ^ base2));
    bb[1] = (ba[1] &        0       );

    return b;
  }

  friend Col4 KillBits(Col4::Arg a, const int n, const int p )
  {
    if ((p + n) >= 64)
      return a;

    Col4 b;

    unsigned__int64 base1 =  (0xFFFFFFFFFFFFFFFFULL << (     (p + 0) & 63));
    unsigned__int64 base2 =  (0xFFFFFFFFFFFFFFFFULL >> (64 - (p + n) & 63));
    unsigned__int64 *bb = (unsigned__int64 *)&b.r;
    unsigned__int64 *ba = (unsigned__int64 *)&a.r;

    bb[0] = (ba[0] & (base1 ^ base2));
    bb[1] = (ba[1] &        0       );

    return b;
  }

  template<const int n, const int p>
  friend Col4 InjtBits(Col4::Arg left, Col4::Arg right)
  {
    if (!n)
      return left;
    if (!p)
      return MaskBits<n, 0>(right);
    if ((p + n) >= 64)
      return (left) + ShiftLeftHalf<p>(right);

    return KillBits<n, p>(left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
  }

  friend Col4 InjtBits(Col4::Arg left, Col4 right, const int n, const int p )
  {
    return KillBits(left, n, p) + MaskBits(ShiftLeftHalf(right, p), n, p);
  }

  template<const int n, const int p>
  friend Col4 ExtrBits(Col4::Arg a )
  {
    if (!n)
      return Col4(0);
    if (!p)
      return MaskBits<n, 0>(a);
    if ((n + p) >= 64)
      return ShiftRightHalf<p>(a);

    return MaskBits<n, 0>(ShiftRightHalf<p>(a));
  }

  friend Col4 ExtrBits(Col4::Arg a, const int n, const int p )
  {
    return MaskBits(ShiftRightHalf(a, p), n, 0);
  }

  template<const int n, const int p>
  friend void ExtrBits(Col4::Arg left, Col4 &right )
  {
    right  = ExtrBits<n, p>( left );
  }

  template<const int n, const int p>
  friend void ConcBits(Col4::Arg left, Col4 &right )
  {
    right  = ShiftLeft<32>( right );
    if (n > 0)
      right += ExtrBits<n, p>( left );
  }

  template<const int n, const int p>
  friend void ReplBits(Col4::Arg left, Col4 &right )
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
  friend Col4 MultiplyAdd(Col4::Arg a, Col4::Arg b, Col4::Arg c )
  {
    return a*b + c;
  }

  //! Returns -( a*b - c )
  friend Col4 NegativeMultiplySubtract(Col4::Arg a, Col4::Arg b, Col4::Arg c )
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

  friend Col4 Min(Col4::Arg left, Col4::Arg right)
  {
    return Col4(
      std::min<int>( left.r, right.r ),
      std::min<int>( left.g, right.g ),
      std::min<int>( left.b, right.b ),
      std::min<int>( left.a, right.a )
    );
  }

  friend Col4 Max(Col4::Arg left, Col4::Arg right)
  {
    return Col4(
      std::max<int>( left.r, right.r ),
      std::max<int>( left.g, right.g ),
      std::max<int>( left.b, right.b ),
      std::max<int>( left.a, right.a )
    );
  }

  friend bool CompareAnyLessThan(Col4::Arg left, Col4::Arg right)
  {
    return
      left.r < right.r ||
      left.g < right.g ||
      left.b < right.b ||
      left.a < right.a;
  }

  friend bool CompareAllEqualTo(Col4::Arg left, Col4::Arg right)
  {
    return
      left.r == right.r &&
      left.g == right.g &&
      left.b == right.b &&
      left.a == right.a;
  }
  
  friend Col4 IsNotZero(Col4::Arg v )
  {
    return Col4(
      v.r > 0 ? ~0 : 0,
      v.g > 0 ? ~0 : 0,
      v.b > 0 ? ~0 : 0,
      v.a > 0 ? ~0 : 0
    );
  }

  friend Col4 IsOne(Col4::Arg v )
  {
    return Col4(
      v.r == 0xFF ? ~0 : 0,
      v.g == 0xFF ? ~0 : 0,
      v.b == 0xFF ? ~0 : 0,
      v.a == 0xFF ? ~0 : 0
    );
  }
  
  friend Col4 IsZero(Col4::Arg v )
  {
    return Col4(
      v.r == 0x00 ? ~0 : 0,
      v.g == 0x00 ? ~0 : 0,
      v.b == 0x00 ? ~0 : 0,
      v.a == 0x00 ? ~0 : 0
    );
  }

  friend Col4 TransferA(Col4::Arg left, Col4::Arg right)
  {
    return Col4( left.r, left.g, left.b, right.a );
  }

  friend Col4 KillA(Col4::Arg left)
  {
    return Col4( left.r, left.g, left.b, 0xFF );
  }

  friend void PackBytes(Col4::Arg v, int &loc )
  {
    loc  = v.a; loc <<= 8;
    loc += v.b; loc <<= 8;
    loc += v.g; loc <<= 8;
    loc += v.r; loc <<= 0;
  }

  // clamp the output to [0, 1]
  Col4 Clamp() const {
    Col4 const one (0xFF);
    Col4 const zero(0x00);

    return Min(one, Max(zero, *this));
  }

  friend void LoadAligned(Col4 &a, Col4 &b, Col4::Arg c )
  {
    a.r = c.r;
    a.g = c.g;

    b.b = c.b;
    b.a = c.a;
  }

  friend void LoadAligned(Col4 &a, void const *source )
  {
    a.r = ((int *)source)[0];
    a.g = ((int *)source)[1];

    a.b = ((int *)source)[2];
    a.a = ((int *)source)[3];
  }

  friend void LoadAligned(Col4 &a, Col4 &b, void const *source )
  {
    a.r = ((int *)source)[0];
    a.g = ((int *)source)[1];

    b.b = ((int *)source)[2];
    b.a = ((int *)source)[3];
  }

  friend void LoadUnaligned(Col4 &a, Col4 &b, void const *source )
  {
    a.r = ((int *)source)[0];
    a.g = ((int *)source)[1];

    b.b = ((int *)source)[2];
    b.a = ((int *)source)[3];
  }

  friend void StoreAligned(Col4::Arg a, Col4::Arg b, Col4 &c )
  {
    c.r = a.r;
    c.g = a.g;

    c.b = b.b;
    c.a = b.a;
  }

  friend void StoreAligned(Col4::Arg a, void *destination )
  {
    ((int *)destination)[0] = a.r;
    ((int *)destination)[1] = a.g;

    ((int *)destination)[2] = a.b;
    ((int *)destination)[3] = a.a;
  }

  friend void StoreAligned(Col4::Arg a, Col4::Arg b, void *destination )
  {
    ((int *)destination)[0] = a.r;
    ((int *)destination)[1] = a.g;

    ((int *)destination)[2] = b.b;
    ((int *)destination)[3] = b.a;
  }

  friend void StoreUnaligned(Col4::Arg a, Col4::Arg b, void *destination )
  {
    ((int *)destination)[0] = a.r;
    ((int *)destination)[1] = a.g;

    ((int *)destination)[2] = b.b;
    ((int *)destination)[3] = b.a;
  }

  friend class Vec4;
  friend class Col8;

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
  return math::rsqrt(v);
}

static float Sqrt(float v) {
  return math::sqrt(v);
}

static float Min(float a, float b) {
  return std::min(a, b);
}

static float Max(float a, float b) {
  return std::max(a, b);
}

static float Abs(float v) {
  return std::abs(v);
}

class Col8
{
public:
  typedef Col8 const& Arg;

  Col8() {}

  Col8( Col8 const& arg ) {
    int i = 7; do { s[i] = arg.s[i]; } while (--i >= 0); }

  Col8& operator=( Col8 const& arg ) {
    int i = 7; do { s[i] = arg.s[i]; } while (--i >= 0);
    return *this;
  }

  explicit Col8(Col4 &v) {
    s[0] = s[1] = (short)v.r;
    s[2] = s[3] = (short)v.g;
    s[4] = s[5] = (short)v.b;
    s[6] = s[7] = (short)v.a;
  }

  explicit Col8(int v) {
    s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = (short)v; }
  explicit Col8(short v) {
    s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = v; }

  Col8( int a, int b, int c, int d, int e, int f, int g, int h ) {
    s[0] = (short)a; s[1] = (short)b; s[2] = (short)c; s[3] = (short)d;
    s[4] = (short)e; s[5] = (short)f; s[6] = (short)g; s[7] = (short)h; }
  Col8( u16 a, u16 b, u16 c, u16 d, u16 e, u16 f, u16 g, u16 h ) {
    s[0] = a; s[1] = b; s[2] = c; s[3] = d; s[4] = e; s[5] = f; s[6] = g; s[7] = h; }

  int Get0() const
  {
    return s[0];
  }

  const u16 &operator[]( int pos ) const
  {
    return s[pos];
  }

  Col8& operator*=( Arg v ) {
    int i = 7; do { s[i] *= v.s[i]; } while (--i >= 0);
    return *this;
  }

  friend Col8 operator>>( Col8::Arg left, int right ) {
    Col8 res;
    int i = 7; do { res.s[i] = left.s[i] >> right; } while (--i >= 0);
    return res;
  }

  friend Col8 operator<<( Col8::Arg left, int right ) {
    Col8 res;
    int i = 7; do { res.s[i] = left.s[i] << right; } while (--i >= 0);
    return res;
  }

  friend Col8 operator+( Col8::Arg left, Col8::Arg right ) {
    Col8 res;
    int i = 7; do { res.s[i] = left.s[i] + right.s[i]; } while (--i >= 0);
    return res;
  }

  friend Col8 operator-( Col8::Arg left, Col8::Arg right ) {
    Col8 res;
    int i = 7; do { res.s[i] = left.s[i] - right.s[i]; } while (--i >= 0);
    return res;
  }

  friend Col8 operator*( Col8::Arg left, Col8::Arg right ) {
    Col8 res;
    int i = 7; do { res.s[i] = left.s[i] * right.s[i]; } while (--i >= 0);
    return res;
  }

  friend Col8 operator*( Col8::Arg left, int right ) {
    Col8 res;
    int i = 7; do { res.s[i] = (short)(left.s[i] * right); } while (--i >= 0);
    return res;
  }

  friend Col8 HorizontalMin( Arg a )
  {
    Col8 res;

    res.s[0] = std::min(a.s[0], a.s[1]);
    res.s[2] = std::min(a.s[2], a.s[3]);
    res.s[4] = std::min(a.s[4], a.s[5]);
    res.s[6] = std::min(a.s[6], a.s[7]);

    res.s[0] = std::min(res.s[0], res.s[2]);
    res.s[4] = std::min(res.s[4], res.s[6]);

    res.s[0] = std::min(res.s[0], res.s[4]);
    res.s[1] = res.s[2] = res.s[3] = res.s[4] = res.s[5] = res.s[6] = res.s[7] = res.s[0];

    return res;
  }

  friend Col8 HorizontalMax( Arg a )
  {
    Col8 res;

    res.s[0] = std::max(a.s[0], a.s[1]);
    res.s[2] = std::max(a.s[2], a.s[3]);
    res.s[4] = std::max(a.s[4], a.s[5]);
    res.s[6] = std::max(a.s[6], a.s[7]);

    res.s[0] = std::max(res.s[0], res.s[2]);
    res.s[4] = std::max(res.s[4], res.s[6]);

    res.s[0] = std::max(res.s[0], res.s[4]);
    res.s[1] = res.s[2] = res.s[3] = res.s[4] = res.s[5] = res.s[6] = res.s[7] = res.s[0];

    return res;
  }

  friend Col4 Expand(Arg a, int ia) {
    return Col4(
      a.s[ia - 0],
      a.s[ia - 1],
      a.s[ia - 2],
      a.s[ia - 3]);
  }

  friend Col4 Repeat(Arg a, int ia) {
    return Col4(
      a.s[ia],
      a.s[ia],
      a.s[ia],
      a.s[ia]);
  }

  friend Col4 Interleave(Arg a, Arg b, int ia, int ib) {
    return Col4(
      a.s[ia],
      b.s[ib],
      a.s[ia],
      b.s[ib]);
  }

  friend Col4 Replicate(Arg a, Arg b, int ia, int ib) {
    return Col4(
      a.s[ia],
      a.s[ia],
      b.s[ib],
      b.s[ib]);
  }

#if	!defined(SQUISH_USE_AMP)
private:
#endif
  short s[8];
};

#if	!defined(SQUISH_USE_COMPUTE)
#define VEC4_CONST( X ) Vec4( X )

class Vec3
{
public:
  typedef Vec3 const& Arg;
  typedef Vec3 & aArg;

  Vec3() {}

  explicit Vec3(float _s) { x = _s; y = _s; z = _s; }

  Vec3(const float *_x, const float *_y, const float *_z) { x = *_x;  y = *_y;  z = *_z; }
  Vec3(float  _x, float  _y, float  _z) { x = _x;   y = _y;   z = _z; }
  Vec3(Vec3   _x, Vec3   _y, Vec3   _z) { x = _x.x; y = _y.x; z = _z.x; }
  
  Vec3(Col3 &c) { x = (float)c.r; y = (float)c.g; z = (float)c.b; }

  void StoreX(float *_x) const { *_x = x; }
  void StoreY(float *_y) const { *_y = y; }
  void StoreZ(float *_z) const { *_z = z; }

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

  friend int operator<(Vec3::Arg left, Vec3::Arg right  )
  {
    return CompareFirstLessThan(left, right);
  }

  friend int operator>(Vec3::Arg left, Vec3::Arg right  )
  {
    return CompareFirstGreaterThan(left, right);
  }

  friend int operator==(Vec3::Arg left, Vec3::Arg right  )
  {
    return CompareFirstEqualTo(left, right);
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

  friend Vec3 Reciprocal(Vec3::Arg v )
  {
    return Vec3(
      1.0f / v.x,
      1.0f / v.y,
      1.0f / v.z
    );
  }

  friend Vec3 ReciprocalSqrt(Vec3::Arg v )
  {
    return Vec3(
      math::rsqrt(v.x),
      math::rsqrt(v.y),
      math::rsqrt(v.z)
    );
  }

  friend Vec3 Sqrt(Vec3::Arg v )
  {
    return Vec3(
      math::sqrt(v.x),
      math::sqrt(v.y),
      math::sqrt(v.z)
    );
  }

  friend Vec3 Select( Arg a, Arg b, Scr3 c )
  {
    if (b.x == c)
      return Vec3(a.x);
    if (b.y == c)
      return Vec3(a.y);
//  if (b.z == c)
      return Vec3(a.z);
  }

  template<const int n>
  friend Vec3 RotateLeft( Arg a )
  {
    return Vec3(
      n == 1 ? a.y : (n == 2 ? a.z : a.x),
      n == 1 ? a.z : (n == 2 ? a.x : a.y),
      n == 1 ? a.x : (n == 2 ? a.y : a.z)
    );
  }

  friend Scr3 HorizontalMin( Arg a )
  {
    return Scr3(
      std::min<float>(std::min<float>(a.x, a.y), a.z)
    );
  }

  friend Scr3 HorizontalMax( Arg a )
  {
    return Scr3(
      std::max<float>(std::max<float>(a.x, a.y), a.z)
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

  friend Vec3 Abs(Vec3::Arg v )
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

  template<const bool round>
  friend Col3 FloatToInt(Vec3::Arg v)
  {
    return Col3(
      (int)(v.x > 0.0f ? std::floor( v.x + (round ? 0.5f : 0.0f) ) : std::ceil( v.x - (round ? 0.5f : 0.0f) )),
      (int)(v.y > 0.0f ? std::floor( v.y + (round ? 0.5f : 0.0f) ) : std::ceil( v.y - (round ? 0.5f : 0.0f) )),
      (int)(v.z > 0.0f ? std::floor( v.z + (round ? 0.5f : 0.0f) ) : std::ceil( v.z - (round ? 0.5f : 0.0f) ))
    );
  }

  friend Vec3 Truncate(Arg v)
  {
    return Vec3(
      v.x > 0.0f ? std::floor( v.x ) : std::ceil( v.x ),
      v.y > 0.0f ? std::floor( v.y ) : std::ceil( v.y ),
      v.z > 0.0f ? std::floor( v.z ) : std::ceil( v.z )
    );
  }

  friend Vec3 AbsoluteDifference(Vec3::Arg left, Vec3::Arg right )
  {
    return Abs( left - right );
  }

  friend Scr3 SummedAbsoluteDifference(Vec3::Arg left, Vec3::Arg right )
  {
    return HorizontalAdd( AbsoluteDifference( left, right ) );
  }

  friend Scr3 MaximumAbsoluteDifference(Vec3::Arg left, Vec3::Arg right )
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

  friend bool CompareAllEqualTo(Vec3::Arg left, Vec3::Arg right)
  {
    return
      left.x == right.x &&
      left.y == right.y &&
      left.z == right.z;
  }

  friend int CompareFirstLessThan(Vec3::Arg left, Vec3::Arg right)
  {
    return left.x < right.x;
  }

  friend int CompareFirstGreaterThan(Vec3::Arg left, Vec3::Arg right)
  {
    return left.x > right.x;
  }

  friend int CompareFirstEqualTo(Vec3::Arg left, Vec3::Arg right)
  {
    return left.x == right.x;
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

  explicit Vec4(float _s) : x(_s), y(_s), z(_s), w(_s) {}
  explicit Vec4(int   _s) : x( (float) _s ), y( (float) _s ), z( (float) _s ), w( (float) _s ) {}

  Vec4(const float *_x, const float *_y, const float *_z, const float *_w) : x(*_x), y(*_y), z(*_z), w(*_w) {}
  Vec4(const float *_x, const float *_y, const float *_z                 ) : x(*_x), y(*_y), z(*_z), w(0.f) {}
  Vec4(const float *_x, const float *_y                                  ) : x(*_x), y(*_y), z(0.f), w(0.f) {}
  Vec4(const float *_a                                                   ) : x(*_a), y(*_a), z(*_a), w(*_a) {}

  Vec4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
  Vec4(float _x, float _y, float _z          ) : x(_x), y(_y), z(_z), w(0.0f) {}
  Vec4(float _x, float _y                    ) : x(_x), y(_y), z(0.0f), w(0.0f) {}

  Vec4(Vec4 _x, Vec4 _y, Vec4 _z, Vec4 _w) : x(_x.x), y(_y.x), z(_z.x), w(_w.x) {}
  Vec4(Vec4 _x, Vec4 _y, Vec4 _z         ) : x(_x.x), y(_y.x), z(_z.x), w(0.0f) {}
  Vec4(Vec4 _x, Vec4 _y                  ) : x(_x.x), y(_y.x), z(0.0f), w(0.0f) {}

  Vec4(Vec3 _v, float _w) : x(_v.x), y(_v.y), z(_v.z), w(_w) {}

  Vec4(Col4 _c) : x((float)_c.r), y((float)_c.g), z((float)_c.b), w((float)_c.a) {}

  Vec3 GetVec3() const
  {
    return Vec3( x, y, z );
  }

  void StoreX(float *_x) const { *_x = x; }
  void StoreY(float *_y) const { *_y = y; }
  void StoreZ(float *_z) const { *_z = z; }
  void StoreW(float *_w) const { *_w = w; }

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
  void SetXYZW(int _x, int _y, int _z, int _w) {
    x = (float)(inv ? inv - _x : _x);
    y = (float)(inv ? inv - _y : _y);
    z = (float)(inv ? inv - _z : _z);
    w = (float)(inv ? inv - _w : _w);
  }

  template<const int inv>
  void SetXYZWpow2(int _x, int _y, int _z, int _w) {
    x = (float)(1 << (inv ? inv - _x : _x));
    y = (float)(1 << (inv ? inv - _y : _y));
    z = (float)(1 << (inv ? inv - _z : _z));
    w = (float)(1 << (inv ? inv - _w : _w));
  }

  template<const int p>
  void Set(const float val)
  {
    /**/ if (p == 0) x = val;
    else if (p == 1) y = val;
    else if (p == 2) z = val;
    else if (p == 3) w = val;
  }

  operator Vec3() const
  {
    Vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
  }
  
  operator Scr4() const
  {
    return x;
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

  Vec4& operator/=(Vec4 v)
  {
    *this *= Reciprocal( v );
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

  friend Vec4 Sqrt( Vec4::Arg v )
  {
    return Vec4(
      math::sqrt(v.x),
      math::sqrt(v.y),
      math::sqrt(v.z),
      math::sqrt(v.w)
    );
  }

  friend Vec4 Select( Arg a, Arg b, Scr4 c )
  {
    if (b.x == c)
      return Vec4(a.x);
    if (b.y == c)
      return Vec4(a.y);
    if (b.z == c)
      return Vec4(a.z);
//  if (b.w == c)
      return Vec4(a.w);
  }

  template<const int n>
  friend Vec4 RotateLeft(Arg a)
  {
    return Vec4(
      n == 1 ? a.y : (n == 2 ? a.z : (n == 3 ? a.w : a.x)),
      n == 1 ? a.z : (n == 2 ? a.w : (n == 3 ? a.x : a.y)),
      n == 1 ? a.w : (n == 2 ? a.x : (n == 3 ? a.y : a.z)),
      n == 1 ? a.x : (n == 2 ? a.y : (n == 3 ? a.z : a.w))
    );
  }
  
  friend Vec4 Threshold(Arg a, Arg b) {
    return Vec4(
      a.x >= b.x ? 1.0f : 0.0f,
      a.y >= b.y ? 1.0f : 0.0f,
      a.z >= b.z ? 1.0f : 0.0f,
      a.w >= b.w ? 1.0f : 0.0f
    );
  }

  friend Scr4 HorizontalMin(Arg a)
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
  
  friend Col4 CompareAllEqualTo_M8(Vec4::Arg left, Vec4::Arg right)
  {
    unsigned char *lc = (unsigned char *)&left.x;
    unsigned char *rc = (unsigned char *)&right.x;
    return Col4(
      (lc[ 0] == rc[ 0] ? 0x000000FF : 0x00000000) +
      (lc[ 1] == rc[ 1] ? 0x0000FF00 : 0x00000000) +
      (lc[ 2] == rc[ 2] ? 0x00FF0000 : 0x00000000) +
      (lc[ 3] == rc[ 3] ? 0xFF000000 : 0x00000000),
      (lc[ 4] == rc[ 4] ? 0x000000FF : 0x00000000) +
      (lc[ 5] == rc[ 5] ? 0x0000FF00 : 0x00000000) +
      (lc[ 6] == rc[ 6] ? 0x00FF0000 : 0x00000000) +
      (lc[ 7] == rc[ 7] ? 0xFF000000 : 0x00000000),
      (lc[ 8] == rc[ 8] ? 0x000000FF : 0x00000000) +
      (lc[ 9] == rc[ 9] ? 0x0000FF00 : 0x00000000) +
      (lc[10] == rc[10] ? 0x00FF0000 : 0x00000000) +
      (lc[11] == rc[11] ? 0xFF000000 : 0x00000000),
      (lc[12] == rc[12] ? 0x000000FF : 0x00000000) +
      (lc[13] == rc[13] ? 0x0000FF00 : 0x00000000) +
      (lc[14] == rc[14] ? 0x00FF0000 : 0x00000000) +
      (lc[15] == rc[15] ? 0xFF000000 : 0x00000000)
    );
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
  
  Vec4 IsZero( ) const
  {
    Vec4 m;

    *((int *)&m.x) = (x == 0.0f ? ~0 : 0);
    *((int *)&m.y) = (y == 0.0f ? ~0 : 0);
    *((int *)&m.z) = (z == 0.0f ? ~0 : 0);
    *((int *)&m.w) = (w == 0.0f ? ~0 : 0);

    return m;
  }

  Vec4 IsNotZero( ) const
  {
    Vec4 m;

    *((int *)&m.x) = (x != 0.0f ? ~0 : 0);
    *((int *)&m.y) = (y != 0.0f ? ~0 : 0);
    *((int *)&m.z) = (z != 0.0f ? ~0 : 0);
    *((int *)&m.w) = (w != 0.0f ? ~0 : 0);

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

  friend void LoadAligned(Vec4 &a, Vec4 &b, Vec4::Arg c )
  {
    a.x = c.x;
    a.y = c.y;

    b.z = c.z;
    b.w = c.w;
  }

  friend void LoadAligned(Vec4 &a, void const *source )
  {
    a.x = ((float *)source)[0];
    a.y = ((float *)source)[1];

    a.z = ((float *)source)[2];
    a.w = ((float *)source)[3];
  }

  friend void LoadAligned(Vec4 &a, Vec4 &b, void const *source )
  {
    a.x = ((float *)source)[0];
    a.y = ((float *)source)[1];

    b.z = ((float *)source)[2];
    b.w = ((float *)source)[3];
  }

  friend void LoadUnaligned(Vec4 &a, void const *source )
  {
    a.x = ((float *)source)[0];
    a.y = ((float *)source)[1];

    a.z = ((float *)source)[2];
    a.w = ((float *)source)[3];
  }

  friend void LoadUnaligned(Vec4 &a, Vec4 &b, void const *source )
  {
    a.x = ((float *)source)[0];
    a.y = ((float *)source)[1];

    b.z = ((float *)source)[2];
    b.w = ((float *)source)[3];
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
inline Scr3 LengthSquared(Vec3::Arg v )
{
  return Dot(v, v);
}

inline void LengthSquared(Vec3::Arg v, float *r)
{
  Dot(v, v, r);
}

inline Scr4 LengthSquared(Vec4::Arg v)
{
  return Dot(v, v);
}

inline void LengthSquared(Vec4::Arg v, float *r)
{
  Dot(v, v, r);
}

inline int CompareFirstLessThan(Scr4 left, Scr4 right)
{
  return left < right;
}

inline Scr4 Threshold(Scr4 a, Scr4 b)
{
  return a >= b ? 1.0f : 0.0f;
}

#endif
#endif // ndef SQUISH_USE_COMPUTE

} // namespace squish

#endif // ndef SQUISH_SIMD_FLOAT_H
