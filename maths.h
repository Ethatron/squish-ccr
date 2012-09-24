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

#ifndef SQUISH_MATHS_H
#define SQUISH_MATHS_H

#if	!defined(USE_COMPUTE)
#include <cmath>
#include <algorithm>
#include <float.h>
#endif

#include "squish.h"

namespace squish {

#if	!defined(USE_COMPUTE)
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

  float X() const { return x; }
  float Y() const { return y; }
  float Z() const { return z; }

  Vec3 operator-() const
  {
    return Vec3( -x, -y, -z );
  }

  Vec3& operator+=( Arg v )
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  Vec3& operator+=( float v )
  {
    x += v;
    y += v;
    z += v;
    return *this;
  }

  Vec3& operator-=( Arg v )
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  Vec3& operator*=( Arg v )
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

  Vec3& operator/=( Arg v )
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

  friend Vec3 operator+( Arg left, Arg right )
  {
    Vec3 copy( left );
    return copy += right;
  }

  friend Vec3 operator+( Arg left, float right )
  {
    Vec3 copy( left );
    return copy += right;
  }

  friend Vec3 operator-( Arg left, Arg right )
  {
    Vec3 copy( left );
    return copy -= right;
  }

  friend Vec3 operator*( Arg left, Arg right )
  {
    Vec3 copy( left );
    return copy *= right;
  }

  friend Vec3 operator*( Arg left, float right )
  {
    Vec3 copy( left );
    return copy *= right;
  }

  friend Vec3 operator*( float left, Arg right )
  {
    Vec3 copy( right );
    return copy *= left;
  }

  friend Vec3 operator/( Arg left, Arg right )
  {
    Vec3 copy( left );
    return copy /= right;
  }

  friend Vec3 operator/( Arg left, float right )
  {
    Vec3 copy( left );
    return copy /= right;
  }

  friend float Dot( Arg left, Arg right )
  {
    return left.x*right.x + left.y*right.y + left.z*right.z;
  }

  friend Vec3 Min( Arg left, Arg right )
  {
    return Vec3(
      std::min<float>( left.x, right.x ),
      std::min<float>( left.y, right.y ),
      std::min<float>( left.z, right.z )
    );
  }

  friend Vec3 Max( Arg left, Arg right )
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

  friend Vec3 Truncate( Arg v )
  {
    return Vec3(
      v.x > 0.0f ? std::floor( v.x ) : std::ceil( v.x ),
      v.y > 0.0f ? std::floor( v.y ) : std::ceil( v.y ),
      v.z > 0.0f ? std::floor( v.z ) : std::ceil( v.z )
    );
  }

  friend class Col3;
  friend class Vec4;

#if	!defined(USE_AMP) && !defined(USE_COMPUTE)
private:
#endif
  union { float x; float r; };
  union { float y; float g; };
  union { float z; float b; };
};

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

  Col3( Vec3 v )
  {
    x = (int)v.x;
    y = (int)v.y;
    z = (int)v.z;
  }

  int X() const { return x; }
  int Y() const { return y; }
  int Z() const { return z; }

  Col3 operator-() const
  {
    return Col3( -x, -y, -z );
  }

  Col3& operator+=( Arg v )
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  Col3& operator-=( Arg v )
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  Col3& operator*=( Arg v )
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

  Col3& operator/=( Arg v )
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

  friend Col3 operator+( Arg left, Arg right )
  {
    Col3 copy( left );
    return copy += right;
  }

  friend Col3 operator-( Arg left, Arg right )
  {
    Col3 copy( left );
    return copy -= right;
  }

  friend Col3 operator*( Arg left, Arg right )
  {
    Col3 copy( left );
    return copy *= right;
  }

  friend Col3 operator*( Arg left, int right )
  {
    Col3 copy( left );
    return copy *= right;
  }

  friend Col3 operator*( int left, Arg right )
  {
    Col3 copy( right );
    return copy *= left;
  }

  friend Col3 operator/( Arg left, Arg right )
  {
    Col3 copy( left );
    return copy /= right;
  }

  friend Col3 operator/( Arg left, int right )
  {
    Col3 copy( left );
    return copy /= right;
  }

  friend int Dot( Arg left, Arg right )
  {
    return left.x*right.x + left.y*right.y + left.z*right.z;
  }

  friend Col3 Min( Arg left, Arg right )
  {
    return Col3(
      std::min<int>( left.x, right.x ),
      std::min<int>( left.y, right.y ),
      std::min<int>( left.z, right.z )
      );
  }

  friend Col3 Max( Arg left, Arg right )
  {
    return Col3(
      std::max<int>( left.x, right.x ),
      std::max<int>( left.y, right.y ),
      std::max<int>( left.z, right.z )
      );
  }

  friend class Col4;
  friend class Vec3;

#if	!defined(USE_AMP)
private:
#endif
  union { int x; int r; };
  union { int y; int g; };
  union { int z; int b; };
};
#endif

#if	!defined(USE_PRE)
inline float LengthSquared( Vec3::Arg v )
{
  return Dot( v, v );
}
#endif

#if	defined(USE_AMP) || defined(USE_COMPUTE)
#if	defined(USE_AMP_DEBUG)
typedef	Vec3	float3;
typedef	Col3	int3;
#endif

#if	!defined(USE_COMPUTE)
namespace Concurrency {
namespace vector_math {
#endif

#if	!defined(USE_COMPUTE)
  float3 minimum( float3 left, float3 right ) amp_restricted
  {
    return float3(
      left.x < right.x ? left.x : right.x,
      left.y < right.y ? left.y : right.y,
      left.z < right.z ? left.z : right.z
    );
  }

  int3 minimum( int3 left, int3 right ) amp_restricted
  {
    return int3(
      left.r < right.r ? left.r : right.r,
      left.g < right.g ? left.g : right.g,
      left.b < right.b ? left.b : right.b
    );
  }

  float3 maximum( float3 left, float3 right ) amp_restricted
  {
    return float3(
      left.x > right.x ? left.x : right.x,
      left.y > right.y ? left.y : right.y,
      left.z > right.z ? left.z : right.z
    );
  }

  int3 maximum( int3 left, int3 right ) amp_restricted
  {
    return int3(
      left.r > right.r ? left.r : right.r,
      left.g > right.g ? left.g : right.g,
      left.b > right.b ? left.b : right.b
    );
  }

  float3 minimax( float3 center, float3 left, float3 right ) amp_restricted
  {
    return minimum(maximum(center, left), right);
  }

  int3 minimax( int3 center, int3 left, int3 right ) amp_restricted
  {
    return minimum(maximum(center, left), right);
  }

  float3 saturate( float3 center ) amp_restricted
  {
    return minimax(center, 0.0f, 1.0f);
  }

  float3 muladd( float3 left, float3 right, float3 offset ) amp_restricted
  {
    return float3(
      left.x * right.x + offset.x,
      left.y * right.y + offset.y,
      left.z * right.z + offset.z
    );
  }

  float3 truncate( float3 center ) amp_restricted
  {
    return float3(
      center.x > 0.0f ? floorf( center.x ) : ceilf( center.x ),
      center.y > 0.0f ? floorf( center.y ) : ceilf( center.y ),
      center.z > 0.0f ? floorf( center.z ) : ceilf( center.z )
    );
  }

  float3 round( float3 center ) amp_restricted
  {
    return float3(
      floorf( center.x + 0.5f ),
      floorf( center.y + 0.5f ),
      floorf( center.z + 0.5f )
    );
  }

  float dot( float3 left, float3 right ) amp_restricted
  {
    return left.x*right.x + left.y*right.y + left.z*right.z;
  }
#endif

  float lengthsquared( float3 center ) amp_restricted
  {
    return dot(center, center);
  }

#if	!defined(USE_COMPUTE)
}
}
#endif
#endif

#if	!defined(USE_COMPUTE)
class Sym2x2
{
public:
	Sym2x2() ccr_restricted
	{
	}

	Sym2x2( float s ) ccr_restricted
	{
		for( int i = 0; i < 3; ++i )
			m_x[i] = s;
	}

	float operator[]( int index ) const
	{
		return m_x[index];
	}

	float& operator[]( int index ) ccr_restricted
	{
		return m_x[index];
	}

private:
	float m_x[3];
};

class Sym3x3
{
public:
	Sym3x3() ccr_restricted
	{
	}

	Sym3x3( float s ) ccr_restricted
	{
		for( int i = 0; i < 6; ++i )
			m_x[i] = s;
	}

	float operator[]( int index ) const
	{
		return m_x[index];
	}

	float& operator[]( int index ) ccr_restricted
	{
		return m_x[index];
	}

private:
	float m_x[6];
};

class Sym4x4
{
public:
	Sym4x4() ccr_restricted
	{
	}

	Sym4x4( float s ) ccr_restricted
	{
		for( int i = 0; i < 10; ++i )
			m_x[i] = s;
	}

	float operator[]( int index ) const
	{
		return m_x[index];
	}

	float& operator[]( int index ) ccr_restricted
	{
		return m_x[index];
	}

private:
	float m_x[10];
};
#endif

} // namespace squish

#include "simd.h"

/* -------------------------------------------------------------------------- */

namespace squish {

#if	!defined(USE_PRE)
Sym3x3 ComputeWeightedCovariance3(int n, Vec3 const* points, float const* weights);
Sym2x2 ComputeWeightedCovariance2(int n, Vec4 const* points, float const* weights);
Sym3x3 ComputeWeightedCovariance3(int n, Vec4 const* points, float const* weights);
Sym4x4 ComputeWeightedCovariance4(int n, Vec4 const* points, float const* weights);
void   ComputePrincipleComponent(Sym3x3 const& smatrix, Vec3 &out);
void   ComputePrincipleComponent(Sym2x2 const& smatrix, Vec4 &out);
void   ComputePrincipleComponent(Sym3x3 const& smatrix, Vec4 &out);
void   ComputePrincipleComponent(Sym4x4 const& smatrix, Vec4 &out);
void   EstimatePrincipleComponent(Sym3x3 const& smatrix, Vec3 &out);
void   EstimatePrincipleComponent(Sym2x2 const& smatrix, Vec4 &out);
void   EstimatePrincipleComponent(Sym3x3 const& smatrix, Vec4 &out);
void   EstimatePrincipleComponent(Sym4x4 const& smatrix, Vec4 &out);

const float *ComputeGammaLUT(bool sRGB);
#endif

#if	defined(USE_AMP) || defined(USE_COMPUTE)
#if	!defined(USE_COMPUTE)
typedef Sym3x3 &Sym3x3r;
#else
typedef float4 Sym3x3[2];
typedef Sym3x3 Sym3x3r;
#endif

Sym3x3 ComputeWeightedCovariance3(tile_barrier barrier, const int thread, int n, point16 points, weight16 weights) amp_restricted;
float3 ComputePrincipleComponent(tile_barrier barrier, const int thread, Sym3x3r smatrix) amp_restricted;
#endif

/* -------------------------------------------------------------------------- */

#undef	OLD_QUANTIZER

template<const int rb, const int gb, const int bb>
class cQuantizer3 {

public:
  static const int rm = (1 << rb) - 1;
  static const int gm = (1 << gb) - 1;
  static const int bm = (1 << bb) - 1;

#ifdef	OLD_QUANTIZER
  Vec3 gridrcp;
  Vec3 grid;
  Vec3 gridoff;

  cQuantizer3() :
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm ),
    grid   ( 1.0f * rm, 1.0f * gm, 1.0f * bm ),
    gridoff( 0.5f     , 0.5f     , 0.5f      ) {
  }

  doinline Vec3 SnapToLattice(Vec3 const &val) {
    return Truncate(grid * val.Clamp() + gridoff) * gridrcp;
  }

  doinline Vec3 SnapToLatticeClamped(Vec3 const &val) {
    return Truncate(grid * val + gridoff) * gridrcp;
  }
#else
  Vec3 gridrcp;
  Vec3 grid;

  cQuantizer3() :
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm ),
    grid   ( 1.0f + rm, 1.0f + gm, 1.0f + bm ) {
  }

  doinline Vec3 SnapToLattice(Vec3 const &val) {
    return (Truncate(grid * val) * gridrcp).Clamp();
  }

  doinline Vec3 SnapToLatticeClamped(Vec3 const &val) {
    return Min(Truncate(grid * val) * gridrcp, Vec3(1.0f));
  }
#endif
};

template<const int rb, const int gb, const int bb, const int ab>
class cQuantizer4 {

public:
  static const int rm = (1 << rb) - 1;
  static const int gm = (1 << gb) - 1;
  static const int bm = (1 << bb) - 1;
  static const int am = (1 << ab) - 1;

#ifdef	OLD_QUANTIZER
  Vec4 gridrcp;
  Vec4 grid;
  Vec4 gridoff;

  cQuantizer4() :
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm, 0.0f ),
    grid   ( 1.0f * rm, 1.0f * gm, 1.0f * bm, 0.0f ),
    gridoff( 0.5f     , 0.5f     , 0.5f     , 0.0f ) {
  }

  doinline Vec4 SnapToLattice(Vec4 const &val) {
    return Truncate(MultiplyAdd(grid, val.Clamp(), gridoff)) * gridrcp;
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) {
    return Truncate(MultiplyAdd(grid, val, gridoff)) * gridrcp;
  }
#else
  Vec4 gridrcp;
  Vec4 grid;

  cQuantizer4() :
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm, 0.0f ),
    grid   ( 1.0f + rm, 1.0f + gm, 1.0f + bm, 0.0f ) {
  }

  doinline Vec4 SnapToLattice(Vec4 const &val) {
    return (Truncate(grid * val) * gridrcp).Clamp();
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) {
    return Min(Truncate(grid * val) * gridrcp, Vec4(1.0f));
  }
#endif
};

class vQuantizer {

public:

#ifdef	OLD_QUANTIZER
  Vec4 gridrcp;
  Vec4 grid;
  Vec4 gridoff;

  vQuantizer(const int rb, const int gb, const int bb, const int ab)/* :
    rm(1 << rb), gm(1 << gb), bm(1 << bb), am(1 << ab),
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm, ab ? 1.0f / am : 1.0f),
    grid   ( 1.0f * rm, 1.0f * gm, 1.0f * bm, ab ? 1.0f * am : 1.0f),
    gridoff( 0.5f     , 0.5f     , 0.5f     , ab ? 0.5f      : 0.0f)*/ {
/*
    grid.SetXYZWpow2<0>(rb, gb, bb, ab);
    Vec4 mask = grid.IsNotOne();

    grid   -= mask & Vec4(1.0f);
    gridrcp = Reciprocal(grid);
    gridoff = mask & Vec4(0.5f);
*/
    grid.SetXYZWpow2<0>(rb, gb, bb, ab ? ab : 8);

    grid   -= Vec4(1.0f);
    gridrcp = Reciprocal(grid);
    gridoff = Vec4(0.5f);
  }

  doinline Vec4 SnapToLattice(Vec4 const &val) {
    return Truncate(MultiplyAdd(grid, val.Clamp(), gridoff)) * gridrcp;
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) {
    return Truncate(MultiplyAdd(grid, val, gridoff)) * gridrcp;
  }

  doinline Col4 QuantizeToInt(Vec4 const &val) {
    // rcomplete[?] = FloatToInt(grid * colour[?].Clamp());
    return FloatToInt<false>(MultiplyAdd(grid, val.Clamp(), gridoff));
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val) {
    // rcomplete[?] = FloatToInt(grid * colour[?]);
    return FloatToInt<false>(MultiplyAdd(grid, val, gridoff));
  }
#else
  Vec4 gridrcp;
  Vec4 grid;
  Col4 gridinv;

#ifdef FEATURE_SHAREDBITS_TRIALS
  Vec4 gridhlf;		// grid * 0.5f
  Vec4 griddbl;		// gridrcp * 2.0f
#endif

  vQuantizer(const int rb, const int gb, const int bb, const int ab, const int sb = 0)/* :
    rm(1 << rb), gm(1 << gb), bm(1 << bb), am(1 << ab),
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm, ab ? 1.0f / am : 1.0f ),
    grid   ( 1.0f + rm, 1.0f + gm, 1.0f + bm, ab ? 1.0f + am : 1.0f )*/ {

    gridinv.SetRGBApow2<0>(rb, gb, bb, ab);
    grid.SetXYZWpow2<0>(rb, gb, bb, ab);

#if 0
    // set inactive channels to 1.0 (Truncate will zero the channel)
    Vec4 mask = grid.IsNotOne();
    Vec4 mone = mask & Vec4(1.0f);
         mask = mask % Vec4(1.0f);

    // quantization is (x*(1<<b))/(1<<b-1)
    // except if a value has zero bits, to prevent
    // the singularity we set rcp to 1 as well
    grid -= mone;
    gridrcp = Reciprocal(grid);
    grid += mone;
    grid -= gridrcp;
    grid  = Max(grid, mask);
#else
    // set inactive channels to 255.0 (Truncate will preserve the channel)
    Vec4 tff = Vec4(255.0f) & grid.IsOne();
    Vec4 one = Vec4(  1.0f);

    // quantization is (x*(1<<b))/(1<<b-1)
    // except if a value has zero bits, to prevent
    // the singularity we set rcp to 1 as well
    grid -= one;
    grid += tff;
    gridrcp = Reciprocal(grid);
    grid += one;
    grid -= gridrcp;
#endif

#ifdef FEATURE_SHAREDBITS_TRIALS
    // get a bit more off when truncating
    gridhlf = grid;
    griddbl = gridrcp;
    if (sb) {
      gridhlf *= Vec4(0.5f);
      griddbl *= Vec4(2.0f);
    }
#else
    // silence the compiler
    bool hb = !!sb; hb = false;
#endif
  }

  doinline Vec4 SnapToLattice(Vec4 const &val) {
    return (Truncate(val * grid) * gridrcp).Clamp();
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) {
    return Min(Truncate(val * grid) * gridrcp, Vec4(1.0f));
  }

  doinline Col4 QuantizeToInt(Vec4 const &val) {
    // [0,255]
    Vec4 Qf = SnapToLattice(val);

    // [0,1<<b-1]
    return (FloatToInt<true>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val) {
    // [0,255]
    Vec4 Qf = SnapToLatticeClamped(val);

    // [0,1<<b-1]
    return (FloatToInt<true>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }

#ifndef FEATURE_SHAREDBITS_TRIALS
#define gridhlf	grid
#define griddbl	gridrcp
#define bitrial	bittest && 0	// silence compiler
#else
#define bitrial	1
#endif

  // if "bitset & bittest" is a compiler constant it becomes 0 penalty
  // no-op vs. the regular version (no additional instructions)

  doinline Vec4 SnapToLattice(Vec4 const &val, int bitset, int bittest) {
    // truncate the last valid bit (half multiplier)
    Vec4 Qf = Truncate(val * gridhlf);

    // add it in if wanted (half 1.0f)
    if (bitrial && (bitset & bittest))
      Qf += Vec4(0.5f);

    // get the bit back in (double multiplier)
    return (Qf * griddbl).Clamp();
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val, int bitset, int bittest) {
    // truncate the last valid bit (half multiplier)
    Vec4 Qf = Truncate(val * gridhlf);

    // add it in if wanted (half 1.0f)
    if (bitrial && (bitset & bittest))
       Qf += Vec4(0.5f);

    // get the bit back in (double multiplier)
    return Min(Qf * griddbl, Vec4(1.0f));
  }

  doinline Col4 QuantizeToInt(Vec4 const &val, int bitset, int bittest) {
    // [0,255]
    Vec4 Qf = SnapToLattice(val, bitset, bittest) * Vec4(255.0f);

    // [0,1<<b-1]
    return (FloatToInt<true>(Qf) * gridinv) >> 8;
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val, int bitset, int bittest) {
    // [0,255]
    Vec4 Qf = SnapToLatticeClamped(val, bitset, bittest) * Vec4(255.0f);

    // [0,1<<b-1]
    return (FloatToInt<true>(Qf) * gridinv) >> 8;
  }

#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= 3)
  doinline Vec4 SnapToLattice(Vec4 const &val, int bitset, int bittest, int oppose) {
    // truncate the last valid bit (half multiplier)
    Vec4 Qf = Truncate(val * gridhlf);

    // add it in if wanted (half 1.0f)
    if (bitrial && (bitset & bittest)) {
      Qf += Vec4(0.5f);
      // but go down by two as well (half 2.0f) per component
      if (oppose)
	Qf -= Vec4(oppose & 0x01 ? 1.0f : 0.0f, oppose & 0x04 ? 1.0f : 0.0f, oppose & 0x10 ? 1.0f : 0.0f, oppose & 0x40 ? 1.0f : 0.0f);
    }
    // don't add it in if wanted (half 0.0f)
    else {
      // but go up by two as well (half 2.0f) per component
      if (oppose)
	Qf += Vec4(oppose & 0x01 ? 1.0f : 0.0f, oppose & 0x04 ? 1.0f : 0.0f, oppose & 0x10 ? 1.0f : 0.0f, oppose & 0x40 ? 1.0f : 0.0f);
    }

    // get the bit back in (double multiplier)
    return (Qf * griddbl).Clamp();
  }
#endif
#endif
};

} // namespace squish

#endif // ndef SQUISH_MATHS_H
