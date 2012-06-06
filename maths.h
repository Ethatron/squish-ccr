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

  friend Vec3 Truncate( Arg v )
  {
    return Vec3(
      v.x > 0.0f ? std::floor( v.x ) : std::ceil( v.x ),
      v.y > 0.0f ? std::floor( v.y ) : std::ceil( v.y ),
      v.z > 0.0f ? std::floor( v.z ) : std::ceil( v.z )
      );
  }

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

  float3 mad( float3 left, float3 right, float3 offset ) amp_restricted
  {
    return float3(
      left.x * right.x + offset.x,
      left.y * right.y + offset.y,
      left.z * right.z + offset.z
    );
  }

  float dot( float3 left, float3 right ) amp_restricted
  {
    return left.x*right.x + left.y*right.y + left.z*right.z;
  }

  float lengthsquared( float3 center ) amp_restricted
  {
    return dot(center, center);
  }

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
#endif

  float3 truncate( float3 center ) amp_restricted
  {
    return float3(
      center.x > 0.0f ? floorf( center.x ) : ceilf( center.x ),
      center.y > 0.0f ? floorf( center.y ) : ceilf( center.y ),
      center.z > 0.0f ? floorf( center.z ) : ceilf( center.z )
    );
  }

#if	!defined(USE_COMPUTE)
}
}
#endif
#endif

#if	!defined(USE_COMPUTE)
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
#endif

#if	!defined(USE_PRE)
Sym3x3 ComputeWeightedCovariance( int n, Vec3 const* points, float const* weights );
Vec3   ComputePrincipleComponent( Sym3x3 const& smatrix );
#endif

#if	defined(USE_AMP) || defined(USE_COMPUTE)
#if	!defined(USE_COMPUTE)
typedef Sym3x3 &Sym3x3r;
#else
typedef float4 Sym3x3[2];
typedef Sym3x3 Sym3x3r;
#endif

Sym3x3 ComputeWeightedCovariance(tile_barrier barrier, const int thread, int n, point16 points, weight16 weights) amp_restricted;
float3 ComputePrincipleComponent(tile_barrier barrier, const int thread, Sym3x3r smatrix) amp_restricted;
#endif

} // namespace squish

#endif // ndef SQUISH_MATHS_H
