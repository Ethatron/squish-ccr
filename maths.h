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
	
	const float* operator()( int index ) const
	{
		return &m_x[index];
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
	
	const float* operator()( int index ) const
	{
		return &m_x[index];
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
	
	const float* operator()( int index ) const
	{
		return &m_x[index];
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
void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec3 &centroid, int n, Vec3 const* points, Vec3 const &metric, Scr3 const* weights);
void ComputeWeightedCovariance2(Sym2x2 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, Vec4 const* weights);
void ComputeWeightedCovariance3(Sym3x3 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, Vec4 const* weights);
void ComputeWeightedCovariance4(Sym4x4 &covariance, Vec4 &centroid, int n, Vec4 const* points, Vec4 const &metric, Vec4 const* weights);
void  ComputePrincipleComponent(Sym3x3 const& smatrix, Vec3 &out);
void  ComputePrincipleComponent(Sym2x2 const& smatrix, Vec4 &out);
void  ComputePrincipleComponent(Sym3x3 const& smatrix, Vec4 &out);
void  ComputePrincipleComponent(Sym4x4 const& smatrix, Vec4 &out);
void EstimatePrincipleComponent(Sym3x3 const& smatrix, Vec3 &out);
void EstimatePrincipleComponent(Sym2x2 const& smatrix, Vec4 &out);
void EstimatePrincipleComponent(Sym3x3 const& smatrix, Vec4 &out);
void EstimatePrincipleComponent(Sym4x4 const& smatrix, Vec4 &out);
void GetPrincipleProjection(Vec3 &enter, Vec3 &leave, Vec3 const &principle, Vec3 const &centroid, int n, Vec3 const* points);
void GetPrincipleProjection(Vec4 &enter, Vec4 &leave, Vec4 const &principle, Vec4 const &centroid, int n, Vec4 const* points);

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

/* ################################################################################# */

#undef	OLD_QUANTIZER

extern const a16 float qLUT_1all[  2], qLUT_1clr[  2], qLUT_1set[  2];
extern const a16 float qLUT_2all[  4], qLUT_2clr[  4], qLUT_2set[  4];
extern const a16 float qLUT_3all[  8], qLUT_3clr[  8], qLUT_3set[  8];
extern const a16 float qLUT_4all[ 16], qLUT_4clr[ 16], qLUT_4set[ 16];
extern const a16 float qLUT_5all[ 32], qLUT_5clr[ 32], qLUT_5set[ 32];
extern const a16 float qLUT_6all[ 64], qLUT_6clr[ 64], qLUT_6set[ 64];
extern const a16 float qLUT_7all[128], qLUT_7clr[128], qLUT_7set[128];
extern const a16 float qLUT_8all[256], qLUT_8clr[256], qLUT_8set[256];

static const a16 float *qLUT_a[9] = {
  qLUT_8all,	qLUT_1all, qLUT_2all, qLUT_3all, qLUT_4all, qLUT_5all, qLUT_6all, qLUT_7all, qLUT_8all };
static const a16 float *qLUT_c[9] = {
  qLUT_8all,	qLUT_1clr, qLUT_2clr, qLUT_3clr, qLUT_4clr, qLUT_5clr, qLUT_6clr, qLUT_7clr, qLUT_8clr };
static const a16 float *qLUT_s[9] = {
  qLUT_8all,	qLUT_1set, qLUT_2set, qLUT_3set, qLUT_4set, qLUT_5set, qLUT_6set, qLUT_7set, qLUT_8set };

/* ********************************************************************************* */

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
  
  /* ------------------------------------------------------------------------------- */
  doinline Vec3 SnapToLattice(Vec3 const &val) const {
    return Truncate(grid * val.Clamp() + gridoff) * gridrcp;
  }

  doinline Vec3 SnapToLatticeClamped(Vec3 const &val) const {
    return Truncate(grid * val + gridoff) * gridrcp;
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col3 QuantizeToInt(Vec3 const &val) const {
    // [0,255]
    Vec3 Qf = SnapToLattice(val);

    // [0,1<<b-1]
    return FloatToInt<true>(grid * Qf);
  }

  doinline Col3 QuantizeToIntClamped(Vec3 const &val) const {
    Vec3 Qf = SnapToLatticeClamped(val);

    // [0,1<<b-1]
    return FloatToInt<true>(grid * Qf);
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col3 LatticeToInt(Vec3 const &val) const {
    // [0,1<<b-1]
    return FloatToInt<true>(grid * val.Clamp());
  }

  doinline Col3 LatticeToIntClamped(Vec3 const &val) const {
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
  
  doinline Vec3 LookUpLattice(int r, int g, int b) const {
    // exact nearest least-error quantization values
    return Vec3(
      &qLUT_a[rb][r],
      &qLUT_a[gb][g],
      &qLUT_a[bb][b]);
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Vec3 SnapToLattice(Vec3 const &val) const {
    Col3 p = FloatToInt<false>((grid * val.Clamp()) + gridgap);
#ifndef NDEBUG
    unsigned int lu; PackBytes(p, lu);
#endif
    int r = p.R(), g = p.G(), b = p.B();
    assert((r >= 0) && (r <= rm) &&
	   (g >= 0) && (g <= gm) &&
	   (b >= 0) && (b <= bm));

    // exact nearest least-error quantization values
    Vec3 rgb = LookUpLattice(r, g, b);

    assert(((unsigned int)(rgb.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((unsigned int)(rgb.Y() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((unsigned int)(rgb.Z() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return rgb;
  }

  doinline Vec3 SnapToLatticeClamped(Vec3 const &val) const {
    Col3 p = FloatToInt<false>((grid * val) + gridgap);
#ifndef NDEBUG
    unsigned int lu; PackBytes(p, lu);
#endif
    int r = p.R(), g = p.G(), b = p.B();
    assert((r >= 0) && (r <= rm) &&
	   (g >= 0) && (g <= gm) &&
	   (b >= 0) && (b <= bm));

    // exact nearest least-error quantization values
    Vec3 rgb = LookUpLattice(r, g, b);

    assert(((unsigned int)(rgb.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((unsigned int)(rgb.Y() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((unsigned int)(rgb.Z() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return rgb;
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col3 QuantizeToInt(Vec3 const &val) const {
    // [0,255]
    Vec3 Qf = SnapToLattice(val);

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
  }

  doinline Col3 QuantizeToIntClamped(Vec3 const &val) const {
    // [0,255]
    Vec3 Qf = SnapToLatticeClamped(val);

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col3 LatticeToInt(Vec3 const &val) const {
    // [0,255]
    Vec3 Qf = val.Clamp();

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
  }

  doinline Col3 LatticeToIntClamped(Vec3 const &val) const {
    // [0,255]
    Vec3 Qf = val;

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec3(255.0f)) * gridinv) >> 8;
  }
#endif
};

/* ********************************************************************************* */

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

  doinline Vec4 SnapToLattice(Vec4 const &val) const {
    return Truncate(MultiplyAdd(grid, val.Clamp(), gridoff)) * gridrcp;
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) const {
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
  
  doinline Vec4 LookUpLattice(int r, int g, int b) const {
    // exact nearest least-error quantization values
    return Vec4(
      &qLUT_a[rb][r],
      &qLUT_a[gb][g],
      &qLUT_a[bb][b]);
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Vec4 SnapToLattice(Vec4 const &val) const {
    Col4 p = FloatToInt<false>((grid * val.Clamp()) + gridgap);
#ifndef NDEBUG
    unsigned int lu; PackBytes(p, lu);
#endif
    int r = p.R(), g = p.G(), b = p.B();
    assert((r >= 0) && (r <= rm) &&
	   (g >= 0) && (g <= gm) &&
	   (b >= 0) && (b <= bm));

    // exact nearest least-error quantization values
    Vec4 rgb = LookUpLattice(r, g, b);

    assert(((unsigned int)(rgb.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((unsigned int)(rgb.Y() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((unsigned int)(rgb.Z() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return rgb;
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) const {
    Col4 p = FloatToInt<false>((grid * val) + gridgap);
#ifndef NDEBUG
    unsigned int lu; PackBytes(p, lu);
#endif
    int r = p.R(), g = p.G(), b = p.B();
    assert((r >= 0) && (r <= rm) &&
	   (g >= 0) && (g <= gm) &&
	   (b >= 0) && (b <= bm));

    // exact nearest least-error quantization values
    Vec4 rgb = LookUpLattice(r, g, b);

    assert(((unsigned int)(rgb.X() * 255.0f) >> (8 - rb)) == ((lu >>  0) & rm));
    assert(((unsigned int)(rgb.Y() * 255.0f) >> (8 - gb)) == ((lu >>  8) & gm));
    assert(((unsigned int)(rgb.Z() * 255.0f) >> (8 - bb)) == ((lu >> 16) & bm));

    return rgb;
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col4 QuantizeToInt(Vec4 const &val) const {
    // [0,255]
    Vec4 Qf = SnapToLattice(val);

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val) const {
    // [0,255]
    Vec4 Qf = SnapToLatticeClamped(val);

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col4 LatticeToInt(Vec4 const &val) const {
    // [0,255]
    Vec4 Qf = val.Clamp();

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }

  doinline Col4 LatticeToIntClamped(Vec4 const &val) const {
    // [0,255]
    Vec4 Qf = val;

    // [0,1<<b-1]
    return (FloatToInt<false>(Qf * Vec4(255.0f)) * gridinv) >> 8;
  }
#endif
};

/* ********************************************************************************* */

#undef	OLD_QUANTIZER
#define SMALLEST_MIDPOINT (0.5f / 255.0f)

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
  
  doinline Vec4 LookUpLattice(int r, int g, int b, int a) const {
    // exact nearest least-error quantization values
    return Vec4(
      r,
      g,
      b,
      a) * gridrcp;
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Vec4 SnapToLattice(Vec4 const &val, int bitset = 0, int bittest = 0) const {
    return Truncate(MultiplyAdd(grid, val.Clamp(), gridoff)) * gridrcp;
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val, int bitset = 0, int bittest = 0) const {
    return Truncate(MultiplyAdd(grid, val, gridoff)) * gridrcp;
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col4 QuantizeToInt(Vec4 const &val, int bitset = 0, int bittest = 0) const {
    // rcomplete[?] = FloatToInt(grid * colour[?].Clamp());
    return FloatToInt<false>(MultiplyAdd(grid, val.Clamp(), gridoff));
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val, int bitset = 0, int bittest = 0) const {
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

  const a16 float *qLUT_t[8];

  void ChangeShared(const int rb, const int gb, const int bb, const int ab, const int sb = -1) {
    qLUT_t[0] = qLUT_a[rb];
    qLUT_t[1] = qLUT_a[gb];
    qLUT_t[2] = qLUT_a[bb];
    qLUT_t[3] = qLUT_a[ab];

#ifdef FEATURE_SHAREDBITS_TRIALS
    // allow bailout if whole == -1 (bit-test always succeeds)
    if ((FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALL) || (~sb)) {
      qLUT_t[0] = qLUT_s[rb]; qLUT_t[4] = qLUT_c[rb]; 
      qLUT_t[1] = qLUT_s[gb]; qLUT_t[5] = qLUT_c[gb]; 
      qLUT_t[2] = qLUT_s[bb]; qLUT_t[6] = qLUT_c[bb]; 
      qLUT_t[3] = qLUT_s[ab]; qLUT_t[7] = qLUT_c[ab];
    }
#endif
  }
  
  Vec4 grid;
  Vec4 gridgap;
  Vec4 gridint;
//Col4 gridinv;
  
  vQuantizer(const int rb, const int gb, const int bb, const int ab, const int sb = -1)/* :
    rm(1 << rb), gm(1 << gb), bm(1 << bb), am(1 << ab),
    gridrcp( 1.0f / rm, 1.0f / gm, 1.0f / bm, ab ? 1.0f / am : 1.0f ),
    grid   ( 1.0f + rm, 1.0f + gm, 1.0f + bm, ab ? 1.0f + am : 1.0f )*/ {

    /* assign the LUTs */
    ChangeShared(rb, gb, bb, ab, sb);
    
//  gridinv.SetRGBApow2<0>(rb, gb, bb, ab);
    gridgap.SetXYZWpow2<8>(rb, gb, bb, ab);
    grid.SetXYZWpow2<0>(rb, gb, bb, ab);

    // set inactive channels to 255.0 (Truncate will preserve the channel)
    Vec4 tff = Vec4(255.0f) & grid.IsOne();
    Vec4 one = Vec4(  1.0f);
    
    // merge "(FloatToInt<false>(Qf * Vec4(255.0f)) * gridinv) >> 8"
    // into "FloatToInt<false>(Qf * gridinv * Vec4(255.0f) >> 8)"
    gridint  = grid;
    gridint *= Vec4(255.0f / 256.0f);
    gridint  = Max(gridint, tff);

    // if a value has zero bits, to prevent
    // the singularity we set rcp to 1 as well
    grid -= one;
    
    // results in "x * 0 + 255", if ab is "0"
    gridgap *= grid;
    gridgap *= Vec4(0.5f / 255.0f);
    gridgap  = Max(gridgap, tff);

    // silence the compiler
    bool hb = !!sb; hb = false;
  }
  
  doinline Vec4 LookUpLattice(int r, int g, int b, int a) const {
    // exact nearest least-error quantization values
    return Vec4(
      &qLUT_t[0][r],
      &qLUT_t[1][g],
      &qLUT_t[2][b],
      &qLUT_t[3][a]);
  }
  
  doinline Vec4 LookUpLattice(int r, int g, int b, int a, int o) const {
    // exact nearest least-error quantization values
    return Vec4(
      &qLUT_t[0 + o][r],
      &qLUT_t[1 + o][g],
      &qLUT_t[2 + o][b],
      &qLUT_t[3 + o][a]);
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Vec4 SnapToLattice(Vec4 const &val) const {
    Col4 p = FloatToInt<false>((grid * val.Clamp()) + gridgap);
    int r = p.R(), g = p.G(), b = p.B(), a = p.A();
    assert((r >= 0) && (r <= 0xFF) &&
	   (g >= 0) && (g <= 0xFF) &&
	   (b >= 0) && (b <= 0xFF) &&
	   (a >= 0) && (a <= 0xFF));

    // exact nearest least-error quantization values
    return LookUpLattice(r, g, b, a);
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val) const {
    Col4 p = FloatToInt<false>((grid * val) + gridgap);
    int r = p.R(), g = p.G(), b = p.B(), a = p.A();
    assert((r >= 0) && (r <= 0xFF) &&
	   (g >= 0) && (g <= 0xFF) &&
	   (b >= 0) && (b <= 0xFF) &&
	   (a >= 0) && (a <= 0xFF));

    // exact nearest least-error quantization values
    return LookUpLattice(r, g, b, a);
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col4 QuantizeToInt(Vec4 const &val) const {
    // [0,255]
    Vec4 Qf = SnapToLattice(val);

    // [0,1<<b-1]
    return FloatToInt<false>(Qf * gridint);
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val) const {
    // [0,255]
    Vec4 Qf = SnapToLatticeClamped(val);
    
    // [0,1<<b-1]
    return FloatToInt<false>(Qf * gridint);
  }
  
  /* -------------------------------------------------------------------------------
   * if "bitset & bittest" is a compiler constant it becomes 0 penalty
   * no-op vs. the regular version (no additional instructions)
   */
#if	!defined(FEATURE_SHAREDBITS_TRIALS)
#define bitrial	(bittest || bitset || 1)	// silence compiler
#else
#define bitrial	(bittest & bitset)
#endif

  doinline Vec4 SnapToLattice(Vec4 const &val, int bitset, int bittest) const {
    Col4 p = FloatToInt<false>((grid * val.Clamp()) + gridgap);
    int r = p.R(), g = p.G(), b = p.B(), a = p.A();
    assert((r >= 0) && (r <= 0xFF) &&
	   (g >= 0) && (g <= 0xFF) &&
	   (b >= 0) && (b <= 0xFF) &&
	   (a >= 0) && (a <= 0xFF));

    // exact nearest least-error quantization values
    return LookUpLattice(r, g, b, a, (bitrial ? 0 : 4));
  }

  doinline Vec4 SnapToLatticeClamped(Vec4 const &val, int bitset, int bittest) const {
    Col4 p = FloatToInt<false>((grid * val) + gridgap);
    int r = p.R(), g = p.G(), b = p.B(), a = p.A();
    assert((r >= 0) && (r <= 0xFF) &&
	   (g >= 0) && (g <= 0xFF) &&
	   (b >= 0) && (b <= 0xFF) &&
	   (a >= 0) && (a <= 0xFF));

    // exact nearest least-error quantization values
    return LookUpLattice(r, g, b, a, (bitrial ? 0 : 4));
  }
  
  /* ------------------------------------------------------------------------------- */
  doinline Col4 QuantizeToInt(Vec4 const &val, int bitset, int bittest) const {
    // [0,255]
    Vec4 Qf = SnapToLattice(val, bitset, bittest);
    
    // [0,1<<b-1]
    return FloatToInt<false>(Qf * gridint);
  }

  doinline Col4 QuantizeToIntClamped(Vec4 const &val, int bitset, int bittest) const {
    // [0,255]
    Vec4 Qf = SnapToLatticeClamped(val, bitset, bittest);
    
    // [0,1<<b-1]
    return FloatToInt<false>(Qf * gridint);
  }
  
  /* ------------------------------------------------------------------------------- */
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_PERMUTE)
  doinline Vec4 SnapToLattice(Vec4 const &val, int bitset, int bittest, int oppose) const {
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

#undef qLUT_t
#endif
};

//extern cQuantizer3<5,6,5  > q3_565;
//extern cQuantizer4<5,6,5,0> q4_565;

/* ********************************************************************************* */

#undef	OLD_QUANTIZER

class fQuantizer {

public:
  void ChangeField(const int tb, const int db) const {
  }

  int trunc;
  Col3 gridprc;
	Col3 gridrnd;
  Col3 griddlt;
  Col3 griddltp;
  Col3 griddltm;

  fQuantizer(const int tb, const int db) {
    /* assign the precision bits */
//  assert(tb > 0);
    trunc = tb;

    /* truncate: tb = 4, 1 << 4 = 0x10 - 1 = 0xF ~= 0xFFF0 */
		gridprc.SetRGBpow2<0>(tb);
		gridrnd  = gridprc;
		gridprc -= Col3(1);
		gridrnd >>= 1;
    gridprc = ~gridprc;

    if (db)
      /* TODO: joint load and unpack */
      griddlt.SetRGBpow2<0>((db >> 0) & 0xFF, (db >> 8) & 0xFF, (db >> 16) & 0xFF/*, (db >> 24) & 0xFF*/);
    else
      griddlt.SetRGBpow2<16>(tb);

    /* plus/minus delta range */
    griddltp = (griddlt >> 1) - Col3(1);
    griddltm = (~griddltp) + Col3(1);
  }

  doinline Col3 MaskLattice(Col3 &p) const {
    // exact nearest least-error quantization values
    // NOTE: rounding doesn't help ... probably because of the log-scale
    return (p & gridprc) + gridrnd;
  }

  doinline Col3 DeltaLattice(Col3 &p, Col3 &b) const {
    // bring values into signed word range first
    Col3 d = (p - b) >> trunc;

    // exact nearest least-error quantization values
    d = Max(Min(d, griddltp), griddltm);
    d = d << trunc;

    // mask to half-range
    return (d + b) & Col3(0x0000FFFF);
  }

	/* ------------------------------------------------------------------------------- */
#if 1
	doinline Col3 QuantizeToLattice(Vec3 const &val) const {
		Col3 p = FloatToSHalf<false>(val);
		int prec = (16 - trunc);

		// re-range X bits of precision to 0x7C00
		p <<= prec;
		p /= (0x7BFF + 1);

		// zero bottom bits
		p <<= trunc;
		// add half a unit (if not zero or INF)
		p += gridrnd & CompareAllEqualTo_M4(p, Col3(0x00000000));
		p |=           CompareAllEqualTo_M4(p, Col3(0x0000FFFF) & gridprc);

		return p;
	}

	doinline Vec3 UnquantizeFromLattice(Col3 const &val) const {
		Col3 p = val;

		// re-range 0xFFFF to X bits of precision
		// 0xFFFF * 32 >> 6 = 0x7BFF
		p *= 31 << (16 - 6);

		return SHalfToFloat(p);
	}

	doinline Col3 QuantizeToLattice(Vec3 const &val, Col3 const &bse) const {
		Col3 b =			 bse ;
		Col3 p = QuantizeToLattice(val);

		// exact nearest least-error quantization values
		p = DeltaLattice(p, b);

		return p;
	}

	doinline void QuantizeToLattice(Vec3 const (&val)[2], Col3 (&res)[2]) const {
		Col3 b = QuantizeToLattice(val[0]);
		Col3 p = QuantizeToLattice(val[1]);

		// exact nearest least-error quantization values
//	b = MaskLattice(b);
		p = DeltaLattice(p, b);

		res[0] = b;
		res[1] = p;
	}

#else
	doinline Col3 QuantizeToLattice(Vec3 const &val) const {
    Col3 p = FloatToUHalf<false>(val);
    assert((p.R() >= 0) && (p.R() <= 0xFFFF) &&
	   (p.G() >= 0) && (p.G() <= 0xFFFF) &&
	   (p.B() >= 0) && (p.B() <= 0xFFFF));

    // exact nearest least-error quantization values
    p = MaskLattice(p);

    return p;
  }

  doinline Col3 QuantizeToLattice(Vec3 const &val, Col3 const &bse) const {
    Col3 b =			 bse ;
    Col3 p = FloatToUHalf<false>(val);
    assert((p.R() >= 0) && (p.R() <= 0xFFFF) &&
	   (p.G() >= 0) && (p.G() <= 0xFFFF) &&
	   (p.B() >= 0) && (p.B() <= 0xFFFF));

    // exact nearest least-error quantization values
    p = DeltaLattice(p, b);

    return p;
  }

  doinline void QuantizeToLattice(Vec3 const (&val)[2], Col3 (&res)[2]) const {
    Col3 b = FloatToUHalf<false>(val[0]);
    Col3 p = FloatToUHalf<false>(val[1]);
    assert((p.R() >= 0) && (p.R() <= 0xFFFF) &&
	   (p.G() >= 0) && (p.G() <= 0xFFFF) &&
	   (p.B() >= 0) && (p.B() <= 0xFFFF));

    // exact nearest least-error quantization values
    b = MaskLattice(b);
    p = DeltaLattice(p, b);

    res[0] = b;
    res[1] = p;
  }

  /* ------------------------------------------------------------------------------- */
  doinline Vec3 SnapToLattice(Vec3 const &val) const {
    Col3 p = FloatToUHalf<false>(val);
    assert((p.R() >= 0) && (p.R() <= 0xFFFF) &&
	   (p.G() >= 0) && (p.G() <= 0xFFFF) &&
	   (p.B() >= 0) && (p.B() <= 0xFFFF));

    // exact nearest least-error quantization values
    p = MaskLattice(p);

    return UHalfToFloat(p);
  }

  doinline Vec3 SnapToLattice(Vec3 const &val, Vec3 const &bse) const {
    Col3 b = FloatToUHalf<false>(bse);
    Col3 p = FloatToUHalf<false>(val);
    assert((p.R() >= 0) && (p.R() <= 0xFFFF) &&
	   (p.G() >= 0) && (p.G() <= 0xFFFF) &&
	   (p.B() >= 0) && (p.B() <= 0xFFFF));

    // exact nearest least-error quantization values
    b = MaskLattice(b);
    p = DeltaLattice(p, b);

    return UHalfToFloat(p);
  }

  doinline void SnapToLattice(Vec3 const (&val)[2], Vec3 (&res)[2]) const {
    Col3 b = FloatToUHalf<false>(val[0]);
    Col3 p = FloatToUHalf<false>(val[1]);
    assert((p.R() >= 0) && (p.R() <= 0xFFFF) &&
	   (p.G() >= 0) && (p.G() <= 0xFFFF) &&
	   (p.B() >= 0) && (p.B() <= 0xFFFF));

    // exact nearest least-error quantization values
    b = MaskLattice(b);
    p = DeltaLattice(p, b);

    res[0] = UHalfToFloat(b);
    res[1] = UHalfToFloat(p);
  }

  doinline void SnapToLattice(Vec3 const (&val)[4], int base, Vec3 (&res)[4]) const {
    Col3 b = FloatToUHalf<false>(val[base]);
    assert((b.R() >= 0) && (b.R() <= 0xFFFF) &&
	   (b.G() >= 0) && (b.G() <= 0xFFFF) &&
	   (b.B() >= 0) && (b.B() <= 0xFFFF));

    // exact nearest least-error quantization values
    b = MaskLattice(b);

    for (int i; i < 4; i++) {
      if (i != base) {
	Col3 d = FloatToUHalf<false>(val[i]);
	assert((d.R() >= 0) && (d.R() <= 0xFFFF) &&
	       (d.G() >= 0) && (d.G() <= 0xFFFF) &&
	       (d.B() >= 0) && (d.B() <= 0xFFFF));

	// exact nearest least-error quantization values
	d = DeltaLattice(d, b);

	res[i] = UHalfToFloat(d);
      }
    }

    res[base] = UHalfToFloat(b);
	}
#endif

};

} // namespace squish

#endif // ndef SQUISH_MATHS_H
