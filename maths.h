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

#if	!defined(SQUISH_USE_COMPUTE)
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <float.h>
#endif

#include "squish.h"

#if	SQUISH_USE_ALTIVEC
//nclude "maths_ve.h"
#elif	SQUISH_USE_SSE
#include "maths_sse.h"
#else
#include "maths_std.h"
#endif

namespace squish {

#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#if	defined(USE_AMP_DEBUG)
typedef	Vec3	float3;
typedef	Col3	int3;
#endif

#if	!defined(SQUISH_USE_COMPUTE)
namespace Concurrency {
namespace vector_math {
#endif

#if	!defined(SQUISH_USE_COMPUTE)
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

#if	!defined(SQUISH_USE_COMPUTE)
}
}
#endif
#endif

#if	!defined(SQUISH_USE_COMPUTE)
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
#pragma warning(disable: 4100)

#if	!defined(SQUISH_USE_PRE)
void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec3 &centroid, int n, Vec3 const* points, Vec3 const &metric);
void ComputeWeightedCovariance2(Sym2x2 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric);
void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric);
void ComputeWeightedCovariance4(Sym4x4 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric);
void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec3 &centroid, int n, Vec3 const* points, Vec3 const &metric, float const* weights);
void ComputeWeightedCovariance2(Sym2x2 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, float const* weights);
void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, float const* weights);
void ComputeWeightedCovariance4(Sym4x4 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, float const* weights);
void  ComputePrincipleComponent(Sym3x3 const& smatrix, Vec3 &out);
void  ComputePrincipleComponent(Sym2x2 const& smatrix, Vec4 &out);
void  ComputePrincipleComponent(Sym3x3 const& smatrix, Vec4 &out);
void  ComputePrincipleComponent(Sym4x4 const& smatrix, Vec4 &out);
void EstimatePrincipleComponent(Sym3x3 const& smatrix, Vec3 &out);
void EstimatePrincipleComponent(Sym2x2 const& smatrix, Vec4 &out);
void EstimatePrincipleComponent(Sym3x3 const& smatrix, Vec4 &out);
void EstimatePrincipleComponent(Sym4x4 const& smatrix, Vec4 &out);

#ifdef FEATURE_POWERESTIMATE
#define GetPrincipleComponent(covariance, m_principle)	EstimatePrincipleComponent(covariance, m_principle)
#else
#define GetPrincipleComponent(covariance, m_principle)	 ComputePrincipleComponent(covariance, m_principle)
#endif

const float *ComputeGammaLUT(bool sRGB);
#endif

#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#if	!defined(SQUISH_USE_COMPUTE)
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

extern const a16 float qLUT_1[2];
extern const a16 float qLUT_2[4];
extern const a16 float qLUT_3[8];
extern const a16 float qLUT_4[16];
extern const a16 float qLUT_5[32];
extern const a16 float qLUT_6[64];
extern const a16 float qLUT_7[128];
extern const a16 float qLUT_8[256];

static const a16 float *qLUT_s[9] = {
  NULL,
  qLUT_1,
  qLUT_2,
  qLUT_3,
  qLUT_4,
  qLUT_5,
  qLUT_6,
  qLUT_7,
  qLUT_8,
};

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
  Col3 gridi;

  cQuantizer3() :
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm ),
    grid   ( 1.0f * rm, 1.0f * gm, 1.0f * bm ),
    gridi  (        rm,        gm,        bm ),
    gridoff( 0.5f     , 0.5f     , 0.5f      ) {
  }

  doinline Vec3 SnapToLattice(Vec3 const &val) {
    return Truncate(grid * val.Clamp() + gridoff) * gridrcp;
  }

  doinline Vec3 SnapToLatticeClamped(Vec3 const &val) {
    return Truncate(grid * val + gridoff) * gridrcp;
  }

  doinline Col3 QuantizeToInt(Vec3 const &val) {
    // [0,255]
    Vec3 Qf = SnapToLattice(val);

    // [0,1<<b-1]
    return FloatToInt<true>(grid * Qf);
  }

  doinline Col3 QuantizeToIntClamped(Vec3 const &val) {
    Vec3 Qf = SnapToLatticeClamped(val);

    // [0,1<<b-1]
    return FloatToInt<true>(grid * Qf);
  }

  doinline Col3 LatticeToInt(Vec3 const &val) {
    // [0,1<<b-1]
    return FloatToInt<true>(grid * val.Clamp());
  }

  doinline Col3 LatticeToIntClamped(Vec3 const &val) {
    // [0,1<<b-1]
    return FloatToInt<true>(grid * val);
  }
#else
  // remainders
  static const int rr = (1 << (8 - rb));
  static const int gr = (1 << (8 - gb));
  static const int br = (1 << (8 - bb));
  
  /* half remainders
  static const float rh = 0.5f * (1 << (8 - rb)) * (rm / 255.0f);
  static const float gh = 0.5f * (1 << (8 - gb)) * (gm / 255.0f);
  static const float bh = 0.5f * (1 << (8 - bb)) * (bm / 255.0f); */

  Vec3 grid;
  Vec3 gridgap;
  Col3 gridinv;

  cQuantizer3() :
    grid   ( 1.0f * rm, 1.0f * gm, 1.0f * bm ),
    gridinv(   1 << rb,   1 << gb,   1 << bb ),
    gridgap( (0.5f * rr * rm) / 255.0f,
	     (0.5f * gr * gm) / 255.0f,
	     (0.5f * br * bm) / 255.0f) {
  }

  doinline Vec3 SnapToLattice(Vec3 const &val) {
    Col3 p = FloatToInt<false>((grid * val.Clamp()) + gridgap);
#ifndef NDEBUG
    int lu; PackBytes(p, lu);
#endif

    // exact nearest least-error quantization values
    Vec3 r = Vec3(qLUT_s[rb][p.R() & rm]);
    Vec3 g = Vec3(qLUT_s[gb][p.G() & gm]);
    Vec3 b = Vec3(qLUT_s[bb][p.B() & bm]);
    
    assert(((int)(r.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((int)(g.X() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((int)(b.X() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return Vec3(r, g, b);
  }

  doinline Vec3 SnapToLatticeClamped(Vec3 const &val) {
    Col3 p = FloatToInt<false>((grid * val) + gridgap);
#ifndef NDEBUG
    int lu; PackBytes(p, lu);
#endif
    
    // exact nearest least-error quantization values
    Vec3 r = Vec3(qLUT_s[rb][p.R() & rm]);
    Vec3 g = Vec3(qLUT_s[gb][p.G() & gm]);
    Vec3 b = Vec3(qLUT_s[bb][p.B() & bm]);
    
    assert(((int)(r.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((int)(g.X() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((int)(b.X() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return Vec3(r, g, b);
  }

  doinline Col3 QuantizeToInt(Vec3 const &val) {
    // [0,255]
    Vec3 Qf = SnapToLattice(val);

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
  }

  doinline Col3 QuantizeToIntClamped(Vec3 const &val) {
    // [0,255]
    Vec3 Qf = SnapToLatticeClamped(val);

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
  }

  doinline Col3 LatticeToInt(Vec3 const &val) {
    // [0,255]
    Vec3 Qf = val.Clamp();

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
  }

  doinline Col3 LatticeToIntClamped(Vec3 const &val) {
    // [0,255]
    Vec3 Qf = val;

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
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
  // remainders
  static const int rr = (1 << (8 - rb));
  static const int gr = (1 << (8 - gb));
  static const int br = (1 << (8 - bb));
  
  /* half remainders
  static const float rh = 0.5f * (1 << (8 - rb)) * (rm / 255.0f);
  static const float gh = 0.5f * (1 << (8 - gb)) * (gm / 255.0f);
  static const float bh = 0.5f * (1 << (8 - bb)) * (bm / 255.0f); */

  Vec4 grid;
  Vec4 gridgap;
  Col4 gridinv;

  cQuantizer4() :
    grid   ( 1.0f * rm, 1.0f * gm, 1.0f * bm ),
    gridinv(   1 << rb,   1 << gb,   1 << bb ),
    gridgap( (0.5f * rr * rm) / 255.0f,
	     (0.5f * gr * gm) / 255.0f,
	     (0.5f * br * bm) / 255.0f) {
  }

  doinline Vec4 SnapToLattice(Vec4 const &val) {
    Col4 p = FloatToInt<false>((grid * val.Clamp()) + gridgap);
#ifndef NDEBUG
    int lu; PackBytes(p, lu);
#endif
    
    // exact nearest least-error quantization values
    Vec4 r = Vec4(qLUT_s[rb][p.R() & rm]);
    Vec4 g = Vec4(qLUT_s[gb][p.G() & gm]);
    Vec4 b = Vec4(qLUT_s[bb][p.B() & bm]);
    
    assert(((int)(r.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((int)(g.X() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((int)(b.X() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return Vec4(r, g, b);
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) {
    Col4 p = FloatToInt<false>((grid * val) + gridgap);
#ifndef NDEBUG
    int lu; PackBytes(p, lu);
#endif
    
    // exact nearest least-error quantization values
    Vec4 r = Vec4(qLUT_s[rb][p.R() & rm]);
    Vec4 g = Vec4(qLUT_s[gb][p.G() & gm]);
    Vec4 b = Vec4(qLUT_s[bb][p.B() & bm]);
    
    assert(((int)(r.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((int)(g.X() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((int)(b.X() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return Vec4(r, g, b);
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

  doinline Col4 LatticeToInt(Vec4 const &val) {
    // [0,255]
    Vec4 Qf = val.Clamp();

    // [0,1<<b-1]
    return (FloatToInt<true>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }

  doinline Col4 LatticeToIntClamped(Vec4 const &val) {
    // [0,255]
    Vec4 Qf = val;

    // [0,1<<b-1]
    return (FloatToInt<true>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }
#endif
};

#undef	OLD_QUANTIZER

class vQuantizer {

public:

#ifdef	OLD_QUANTIZER
  Vec4 gridrcp;
  Vec4 grid;
  Vec4 gridoff;

  vQuantizer(const int rb, const int gb, const int bb, const int ab, const int sb = 0)/* :
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

  doinline Vec4 SnapToLattice(Vec4 const &val, int bitset = 0, int bittest = 0) {
    return Truncate(MultiplyAdd(grid, val.Clamp(), gridoff)) * gridrcp;
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val, int bitset = 0, int bittest = 0) {
    return Truncate(MultiplyAdd(grid, val, gridoff)) * gridrcp;
  }

  doinline Col4 QuantizeToInt(Vec4 const &val, int bitset = 0, int bittest = 0) {
    // rcomplete[?] = FloatToInt(grid * colour[?].Clamp());
    return FloatToInt<false>(MultiplyAdd(grid, val.Clamp(), gridoff));
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val, int bitset = 0, int bittest = 0) {
    // rcomplete[?] = FloatToInt(grid * colour[?]);
    return FloatToInt<false>(MultiplyAdd(grid, val, gridoff));
  }
#else
  /* remainders
  static const int rr = (1 << (8 - rb));
  static const int gr = (1 << (8 - gb));
  static const int br = (1 << (8 - bb));
  static const int ar = (1 << (8 - ab));
  
  // half remainders
  static const float rh = 0.5f * (1 << (8 - rb)) * (rm / 255.0f);
  static const float gh = 0.5f * (1 << (8 - gb)) * (gm / 255.0f);
  static const float bh = 0.5f * (1 << (8 - bb)) * (bm / 255.0f);
  static const float ah = 0.5f * (1 << (8 - ab)) * (am / 255.0f); */

  const a16 float *qLUT_rgba[4];

  Vec4 gridrcp;
  Vec4 grid;
  Vec4 gridgap;
  Col4 gridinv;

#ifdef FEATURE_SHAREDBITS_TRIALS
  Vec4 gridhlf;		// grid * 0.5f
  Vec4 griddbl;		// gridrcp * 2.0f
#endif
  
  vQuantizer(const int rb, const int gb, const int bb, const int ab, const int sb = 0)/* :
    rm(1 << rb), gm(1 << gb), bm(1 << bb), am(1 << ab),
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm, ab ? 1.0f / am : 1.0f ),
    grid   ( 1.0f + rm, 1.0f + gm, 1.0f + bm, ab ? 1.0f + am : 1.0f )*/ {
    qLUT_rgba[0] = qLUT_s[rb];
    qLUT_rgba[1] = qLUT_s[gb];
    qLUT_rgba[2] = qLUT_s[bb];
    qLUT_rgba[3] = qLUT_s[ab];
      
    gridgap.SetXYZWpow2<8>(rb, gb, bb, ab);
    gridinv.SetRGBApow2<0>(rb, gb, bb, ab);
    grid.SetXYZWpow2<0>(rb, gb, bb, ab);
    
    // set inactive channels to 255.0 (Truncate will preserve the channel)
    Vec4 tff = Vec4(255.0f) & grid.IsOne();
    Vec4 one = Vec4(  1.0f);

    // if a value has zero bits, to prevent
    // the singularity we set rcp to 1 as well
    grid -= one;
    grid += tff;
    gridgap *= grid;
    gridrcp = Reciprocal(grid);
    gridgap *= Vec4(0.5f / 255.0f);
    grid += one;

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
    Col4 p = FloatToInt<false>((grid * val.Clamp()) + gridgap);
    
    // exact nearest least-error quantization values
    Vec4 r = Vec4(qLUT_rgba[0][p.R() & 0xFF]);
    Vec4 g = Vec4(qLUT_rgba[1][p.G() & 0xFF]);
    Vec4 b = Vec4(qLUT_rgba[2][p.B() & 0xFF]);
    Vec4 a = Vec4(qLUT_rgba[3][p.A() & 0xFF]);

    return Vec4(r, g, b, a);
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) {
    Col4 p = FloatToInt<false>((grid * val) + gridgap);
    
    // exact nearest least-error quantization values
    Vec4 r = Vec4(qLUT_rgba[0][p.R() & 0xFF]);
    Vec4 g = Vec4(qLUT_rgba[1][p.G() & 0xFF]);
    Vec4 b = Vec4(qLUT_rgba[2][p.B() & 0xFF]);
    Vec4 a = Vec4(qLUT_rgba[3][p.A() & 0xFF]);

    return Vec4(r, g, b, a);
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

//extern cQuantizer3<5,6,5  > q3_565;
//extern cQuantizer4<5,6,5,0> q4_565;

} // namespace squish

#endif // ndef SQUISH_MATHS_H
