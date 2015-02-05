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

#include "alpha.h"
#include "maths.h"
#include "simd.h"

#if	!defined(SQUISH_USE_COMPUTE)
#include <cmath>
#include <algorithm>
#endif

#include "inlineables.cpp"

// Codebook precision bits, up to 5 (+8)
#define CBLB	0	// low precision
#define CBHB	5	// high precision

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
template<typename dtyp>
static void CompressAlphaBtc2u(dtyp const* rgba, int mask, void* block, float scale)
{
  u8* bytes = reinterpret_cast< u8* >(block);

  // quantize and pack the alpha values pairwise
  for (int i = 0, imask = mask; i < 8; ++i, imask >>= 2) {
    // quantize down to 4 bits
    float alpha1 = (float)rgba[8 * i + 3] * scale;
    float alpha2 = (float)rgba[8 * i + 7] * scale;

    int quant1 = FloatToInt<true,false>(alpha1, 15);
    int quant2 = FloatToInt<true,false>(alpha2, 15);

    // set alpha to zero where masked
    if ((imask & 1) == 0)
      quant1 = 0;
    if ((imask & 2) == 0)
      quant2 = 0;

    // pack into the byte
    bytes[i] = (u8)((quant1 << 0) + (quant2 << 4));
  }
}

void CompressAlphaBtc2u(u8  const* rgba, int mask, void* block) {
  CompressAlphaBtc2u(rgba, mask, block, 15.0f / 255.0f); }
void CompressAlphaBtc2u(u16 const* rgba, int mask, void* block) {
  CompressAlphaBtc2u(rgba, mask, block, 15.0f / 65535.0f); }
void CompressAlphaBtc2u(f23 const* rgba, int mask, void* block) {
  CompressAlphaBtc2u(rgba, mask, block, 15.0f / 1.0f); }

template<typename dtyp>
static void DecompressAlphaBtc2u(dtyp* rgba, void const* block, dtyp scale)
{
  u8 const* bytes = reinterpret_cast< u8 const* >(block);

  // unpack the alpha values pairwise
  for (int i = 0; i < 8; ++i) {
    // quantize down to 4 bits
    u8 quant = bytes[i];

    // unpack the values
    u8 lo = quant & 0x0F;
    u8 hi = quant >> 4;

    // convert back up to bytes
    rgba[8 * i + 3] = (lo * 0x11) * scale;
    rgba[8 * i + 7] = (hi * 0x11) * scale;
  }
}

void DecompressAlphaBtc2u(u8 * rgba, void const* block) {
  DecompressAlphaBtc2u(rgba, block, (u8 )(255 / 255)); }
void DecompressAlphaBtc2u(u16* rgba, void const* block) {
  DecompressAlphaBtc2u(rgba, block, (u16)(65535 / 255)); }
void DecompressAlphaBtc2u(f23* rgba, void const* block) {
  DecompressAlphaBtc2u(rgba, block, (f23)(1.0f / 255.0f)); }

/* *****************************************************************************
 */
template<const int min, const int max, const int steps>
static void FixRange(int& minS, int& maxS)
{
  if (maxS - minS < steps)
    maxS = std::min<int>(minS + steps, max);
  if (maxS - minS < steps)
    minS = std::max<int>(min, maxS - steps);
}

/* -----------------------------------------------------------------------------
 */
template<const int mul, const int div, typename otyp, typename dtyp>
static Scr4 FitCodesS(dtyp const* rgba, int mask, Col8 const &codes, u8* indices)
{
  // fit each alpha value to the codebook
  int err = 0;
  for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
    // check this pixel is valid
    if ((imask & 1) == 0) {
      // use the first code
      indices[i] = 0;
      continue;
    }

    // find the least error and corresponding index
    // prefer lower indices over higher ones
    Col8 value = Col8(rgba[4 * i + 3]);
    
    // for all possible codebook-entries
    Col8 dst = codes - value; dst *= dst;
    
    // get best index, binary search
    int  index = 0;
    Col8 least = HorizontalMin(dst);
    int  match = CompareEqualTo(least, dst);

    if (!(match & 0xFF))
      index += 4, match >>= 8;
    if (!(match & 0x0F))
      index += 2, match >>= 4;
    if (!(match & 0x03))
      index += 1, match >>= 2;

    // save this index and accumulate the error
    indices[i] = (u8)index;
    err += least.Get0();
  }

  // return the total error
  return Scr4(err);
}

template<const int mul, const int div, typename otyp, typename dtyp>
static Scr4 FitCodesL(dtyp const* rgba, int mask, Col8 const &codes, u8* indices)
{
  const Vec4 codesl = LoVec4(codes, otyp(0));
  const Vec4 codesh = HiVec4(codes, otyp(0));
  const f23 rerange = f23(mul) / f23(div);

  // fit each alpha value to the codebook
  Scr4 err = Scr4(0.0f);
  for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
    // check this pixel is valid
    if ((imask & 1) == 0) {
      // use the first code
      indices[i] = 0;
      continue;
    }

    // find the least error and corresponding index
    // prefer lower indices over higher ones
    f23 value = rgba[4 * i + 3] * rerange;

    // find the closest code
    Vec4 val = Vec4(value);

    // for all possible codebook-entries
    Vec4 dsl = codesl - val; dsl *= dsl;
    Vec4 dsh = codesh - val; dsh *= dsh;
    Vec4 dst = Min(dsl, dsh);
    
    // get best index, binary search
    int  index = 0;
    Scr4 least = HorizontalMin(dst);
    int  match = 
      (CompareEqualTo(least, dsl) << 0) + 
      (CompareEqualTo(least, dsh) << 4);

    if (!(match & 0xF))
      index += 4, match >>= 4;
    if (!(match & 0x3))
      index += 2, match >>= 2;
    if (!(match & 0x1))
      index += 1, match >>= 1;

    // save this index and accumulate the error
    indices[i] = (u8)index;
    err += least;
  }

  // return the total error
  return Scr4(err);
}

template<const int prc, typename otyp>
static doinline Scr4 FitCodes(u8  const* rgba, int mask, Col8 const &codes, u8* indices) {
  return FitCodesS<1 << prc,1,otyp>(rgba, mask, codes, indices); }
template<const int prc, typename otyp>
static doinline Scr4 FitCodes(s8  const* rgba, int mask, Col8 const &codes, u8* indices) {
  return FitCodesS<1 << prc,1,otyp>(rgba, mask, codes, indices); }

template<const int prc, typename otyp>
static doinline Scr4 FitCodes(u16 const* rgba, int mask, Col8 const &codes, u8* indices) {
  return FitCodesL<1 << prc,257,otyp>(rgba, mask, codes, indices); }
template<const int prc, typename otyp>
static doinline Scr4 FitCodes(s16 const* rgba, int mask, Col8 const &codes, u8* indices) {
  return FitCodesL<1 << prc,258,otyp>(rgba, mask, codes, indices); }

template<const int prc, typename otyp, const int mul>
static doinline Scr4 FitCodes(f23 const* rgba, int mask, Col8 const &codes, u8* indices) {
  return FitCodesL<mul << prc,1,otyp>(rgba, mask, codes, indices); }

/* -----------------------------------------------------------------------------
 */
template<const int mul, const int div, typename otyp, typename dtyp>
static void GetErrorS(dtyp const* rgba, int mask,
		      Col8 const &codes5, Col8 const &codes7,
		      Scr4 &error5, Scr4 &error7,
		      float (&aaaa)[16])
{
  // initial values
  int err5 = 0;
  int err7 = 0;
  for (int v = 0; v < 16; v++) {
    int value = rgba[4 * v + 3];

    // find the closest code
    Col8 val = Col8(value);
      
    // create floating point vector
    aaaa[v] = float(value);

    // for all possible codebook-entries
    Col8 d5 = codes5 - val; d5 *= d5;
    Col8 d7 = codes7 - val; d7 *= d7;
      
    // accumulate the error
    err5 += HorizontalMin(d5).Get0();
    err7 += HorizontalMin(d7).Get0();
  }
  
  // return the total error
  error5 = Scr4(err5);
  error7 = Scr4(err7);
}

template<const int mul, const int div, typename otyp, typename dtyp>
static void GetErrorL(dtyp const* rgba, int mask,
		      Col8 const &codes5, Col8 const &codes7,
		      Scr4 &error5, Scr4 &error7,
		      float (&aaaa)[16])
{
  const Vec4 codes5l = LoVec4(codes5, otyp(0));
  const Vec4 codes5h = HiVec4(codes5, otyp(0));
  const Vec4 codes7l = LoVec4(codes7, otyp(0));
  const Vec4 codes7h = HiVec4(codes7, otyp(0));
  const f23 rerange = f23(mul) / f23(div);

  // initial values
  Scr4 err5 = Scr4(0.0f);
  Scr4 err7 = Scr4(0.0f);
  for (int v = 0; v < 16; v++) {
    f23 value = rgba[4 * v + 3] * rerange;

    // find the closest code
    Vec4 val = Vec4(value);
      
    // create floating point vector
    aaaa[v] = float(value);

    // for all possible codebook-entries
    Vec4 d5l = codes5l - val; d5l *= d5l;
    Vec4 d5h = codes5h - val; d5h *= d5h;
    Vec4 d7l = codes7l - val; d7l *= d7l;
    Vec4 d7h = codes7h - val; d7h *= d7h;
      
    // accumulate the error
    err5 += HorizontalMin(Min(d5l, d5h));
    err7 += HorizontalMin(Min(d7l, d7h));
  }
  
  // return the total error
  error5 = err5;
  error7 = err7;
}

template<const int prc, typename otyp>
static doinline void GetError(u8  const* rgba, int mask, Col8 const &codes5, Col8 const &codes7, Scr4 &error5, Scr4 &error7, float (&aaaa)[16]) {
  GetErrorS<1 << prc,1,otyp>(rgba, mask, codes5, codes7, error5, error7, aaaa); }
template<const int prc, typename otyp>
static doinline void GetError(s8  const* rgba, int mask, Col8 const &codes5, Col8 const &codes7, Scr4 &error5, Scr4 &error7, float (&aaaa)[16]) {
  GetErrorS<1 << prc,1,otyp>(rgba, mask, codes5, codes7, error5, error7, aaaa); }

template<const int prc, typename otyp>
static doinline void GetError(u16 const* rgba, int mask, Col8 const &codes5, Col8 const &codes7, Scr4 &error5, Scr4 &error7, float (&aaaa)[16]) {
  GetErrorL<1 << prc,257,otyp>(rgba, mask, codes5, codes7, error5, error7, aaaa); }
template<const int prc, typename otyp>
static doinline void GetError(s16 const* rgba, int mask, Col8 const &codes5, Col8 const &codes7, Scr4 &error5, Scr4 &error7, float (&aaaa)[16]) {
  GetErrorL<1 << prc,258,otyp>(rgba, mask, codes5, codes7, error5, error7, aaaa); }

template<const int prc, typename otyp, const int mul>
static doinline void GetError(f23 const* rgba, int mask, Col8 const &codes5, Col8 const &codes7, Scr4 &error5, Scr4 &error7, float (&aaaa)[16]) {
  GetErrorL<mul << prc,1,otyp>(rgba, mask, codes5, codes7, error5, error7, aaaa); }

/* -----------------------------------------------------------------------------
 */
#define FIT_THRESHOLD 1e-15f

template<const int min, const int max, const int prc, const int steps>
static Scr4 FitError(float const* aaaa, int &minS, int &maxS, Scr4 &errS) {
  const Vec4 mprc = Vec4(1.0f * (1 << prc));

  // for 5-step
  const Vec4 fouf = Vec4(4.0f * (1 << prc) / 5.0f);
  const Vec4 thrf = Vec4(3.0f * (1 << prc) / 5.0f);
  const Vec4 twof = Vec4(2.0f * (1 << prc) / 5.0f);
  const Vec4 onef = Vec4(1.0f * (1 << prc) / 5.0f);
  
  // for 7-step
  const Vec4 sixs = Vec4(6.0f * (1 << prc) / 7.0f);
  const Vec4 fivs = Vec4(5.0f * (1 << prc) / 7.0f);
  const Vec4 fous = Vec4(4.0f * (1 << prc) / 7.0f);
  const Vec4 thrs = Vec4(3.0f * (1 << prc) / 7.0f);
  const Vec4 twos = Vec4(2.0f * (1 << prc) / 7.0f);
  const Vec4 ones = Vec4(1.0f * (1 << prc) / 7.0f);

#if	(SQUISH_USE_SIMD > 0)
  int rng = (maxS - minS) >> 1;	// max 127

  Scr4 errC = errS;
  Vec4 s    = Vec4(minS);
  Vec4 e    = Vec4(maxS);
  Vec4 r    = Vec4(rng) * Vec4(1.0f, -1.0f, 0.0f, 0.0f);

  while (r > Vec4(0.0f)) {
    // s + os, s + os, s + os - r, s + os + r
    // e - oe - r, e - oe + r, e - oe, e - oe
    Vec4 cssmp = s + r       ; cssmp = Max(Min(cssmp, Vec4(max)), Vec4(min));
    Vec4 cmpee = e + r.Swap(); cmpee = Max(Min(cmpee, Vec4(max)), Vec4(min));
    
    // for all possible codebook-entries
    Vec4 cb0, cb1, cb2, cb3,
	 cb4, cb5, cb6, cb7;

    if (steps == 5) {
      cb0 = (cssmp * mprc)                 ; //0 = Truncate(cb0);
      cb1 = (cssmp * fouf) + (cmpee * onef); cb1 = Truncate(cb1);
      cb2 = (cssmp * thrf) + (cmpee * twof); cb2 = Truncate(cb2);
      cb3 = (cssmp * twof) + (cmpee * thrf); cb3 = Truncate(cb3);
      cb4 = (cssmp * onef) + (cmpee * fouf); cb4 = Truncate(cb4);
      cb5 =                  (cmpee * mprc); //5 = Truncate(cb5);
      cb6 = Vec4(min * (1 << prc));
      cb7 = Vec4(max * (1 << prc));
    }
    else if (steps == 7) {
      cb0 = (cssmp * mprc)                 ; //0 = Truncate(cb0);
      cb1 = (cssmp * sixs) + (cmpee * ones); cb1 = Truncate(cb1);
      cb2 = (cssmp * fivs) + (cmpee * twos); cb2 = Truncate(cb2);
      cb3 = (cssmp * fous) + (cmpee * thrs); cb3 = Truncate(cb3);
      cb4 = (cssmp * thrs) + (cmpee * fous); cb4 = Truncate(cb4);
      cb5 = (cssmp * twos) + (cmpee * fivs); cb5 = Truncate(cb5);
      cb6 = (cssmp * ones) + (cmpee * sixs); cb6 = Truncate(cb6);
      cb7 =                  (cmpee * mprc); //7 = Truncate(cb7);
    }

    Vec4 errors = Vec4(0.0f);
    for (int v = 0; v < 16; v++) {
      // find the closest code
      Vec4 val = Vec4(&aaaa[v]);
      Vec4 dists;

      Vec4 d0 = cb0 - val; d0 = d0 * d0;
      Vec4 d1 = cb1 - val; d1 = d1 * d1;
      Vec4 d2 = cb2 - val; d2 = d2 * d2;
      Vec4 d3 = cb3 - val; d3 = d3 * d3;
      Vec4 d4 = cb4 - val; d4 = d4 * d4;
      Vec4 d5 = cb5 - val; d5 = d5 * d5;
      Vec4 d6 = cb6 - val; d6 = d6 * d6;
      Vec4 d7 = cb7 - val; d7 = d7 * d7;

      // encourage OoO
      Vec4 da = Min(d0, d1);
      Vec4 db = Min(d2, d3);
      Vec4 dc = Min(d4, d5);
      Vec4 dd = Min(d6, d7);
      Vec4 de = Min(da, db);
      Vec4 df = Min(dc, dd);
      dists   = Min(de, df);

      // accumulate the error
      errors = errors + dists;
    }

    // check which range has smaller error
    int value = CompareLessThan(errors, errC);
      
#if defined(TRACK_STATISTICS)
			  gstat.alphacond[steps == 5 ? 0 : 1][0]++;
    /**/ if (value ==  0) gstat.alphacond[steps == 5 ? 0 : 1][1]++;
    else if (value & 0x1) gstat.alphacond[steps == 5 ? 0 : 1][2]++;
    else if (value & 0x2) gstat.alphacond[steps == 5 ? 0 : 1][3]++;
    else if (value & 0x4) gstat.alphacond[steps == 5 ? 0 : 1][4]++;
    else /* (value & 8)*/ gstat.alphacond[steps == 5 ? 0 : 1][5]++;
#endif

    /**/ if (value ==  0) { r = Truncate(r * Vec4(0.5f)); continue; }
    // e -= r;// range up
    else if (value & 0x8) errC = errors.SplatW(), e = cmpee.SplatW();
    // s += r;// range dn
    else if (value & 0x1) errC = errors.SplatX(), s = cssmp.SplatX();
    // e += r;// range up
    else if (value & 0x4) errC = errors.SplatZ(), e = cmpee.SplatZ();
    // s -= r;// range dn
    else /* (value & 2)*/ errC = errors.SplatY(), s = cssmp.SplatY();

    // lossless
    if (!(errC > Scr4(FIT_THRESHOLD)))
      break;
  }

#if defined(TRACK_STATISTICS)
  /* there is a clear skew towards unweighted clusterfit (all weights == 1)
    *
    * C == 3, numset ==
    *  [0]	0x00f796e0 {124800, 15616}
    */
  if (steps == 5) {
    gstat.alpha[0] += errS;
    gstat.alpha[2] += 1;
  }
  else if (steps == 7) {
    gstat.alpha[3] += errS;
    gstat.alpha[5] += 1;
  }
#endif

  // final match
  minS = FloatToInt<false>(Min(s, e)).GetLong();
  maxS = FloatToInt<false>(Max(s, e)).GetLong();
  errS = errC;

#if defined(TRACK_STATISTICS)
  /* there is a clear skew towards unweighted clusterfit (all weights == 1)
    *
    * C == 3, numset ==
    *  [0]	0x00f796e0 {124800, 15616}
    */
  if (steps == 5)
    gstat.alpha[1] += errS;
  else if (steps == 7)
    gstat.alpha[4] += errS;
#endif
#else
  // binary search, tangent-fitting
  int error0 = errS;
  int s = minS;
  int e = maxS;
  int r = (e - s) >> 1;	// max 127

  while (r != 0) {
    int ms = std::max(std::min(s - r, max), min);  // os - r
    int cs = std::max(std::min(s    , max), min);  // os
    int ps = std::max(std::min(s + r, max), min);  // os + r
    int me = std::max(std::min(e - r, max), min);  // oe + r
    int ce = std::max(std::min(e    , max), min);  // oe
    int pe = std::max(std::min(e + r, max), min);  // oe - r
    
    // construct code-books
    Col8 cb1, cb2, cb3, cb4;
    
    if (steps == 5) {
      Codebook6<min,max,prc>(codes5, Col8(cs), Col8(me));
      Codebook6<min,max,prc>(codes5, Col8(cs), Col8(pe));
      Codebook6<min,max,prc>(codes5, Col8(ms), Col8(ce));
      Codebook6<min,max,prc>(codes5, Col8(ps), Col8(ce));
    }
    else if (steps == 7) {
      Codebook8<min,max,prc>(codes5, Col8(cs), Col8(me));
      Codebook8<min,max,prc>(codes5, Col8(cs), Col8(pe));
      Codebook8<min,max,prc>(codes5, Col8(ms), Col8(ce));
      Codebook8<min,max,prc>(codes5, Col8(ps), Col8(ce));
    }

    int error1 = 0;
    int error2 = 0;
    int error3 = 0;
    int error4 = 0;
    for (int v = 0; v < 16; v++) {
      // find the closest code
      Col8 val = Col8(aaaa[v]);

      // for all possible codebook-entries
      Col8 d1 = cb1 - val; d1 *= d1;
      Col8 d2 = cb2 - val; d2 *= d2;
      Col8 d3 = cb3 - val; d3 *= d3;
      Col8 d4 = cb4 - val; d4 *= d4;
      
      // accumulate the error
      error1 += HorizontalMin(d1).Get0();
      error2 += HorizontalMin(d2).Get0();
      error3 += HorizontalMin(d3).Get0();
      error4 += HorizontalMin(d4).Get0();
    }

    int                merr = error0;
    if (merr > error1) merr = error1;
    if (merr > error4) merr = error4;
    if (merr > error2) merr = error2;
    if (merr > error3) merr = error3;

    /**/ if (merr == error0) r >>= 1;		// half range
    else if (merr == error1) e = me;//e -= r;	// range up
    else if (merr == error4) s = ps;//s += r;	// range dn
    else if (merr == error2) e = pe;//e += r;	// range up
    else if (merr == error3) s = ms;//s -= r;	// range dn

    // lossless
    error0 = merr;
    if (!error0)
      break;
  }

#if defined(TRACK_STATISTICS)
  /* way better!
    * [0]	17302780	int
    * [1]	3868483		int
    * [2]	69408		int
    */
  if (steps == 5) {
    gstat.alpha[0] += errS;
    gstat.alpha[2] += 1;
  }
  else if (steps == 7) {
    gstat.alpha[3] += errS;
    gstat.alpha[5] += 1;
  }
#endif

  // final match
  minS = std::min(s, e);
  maxS = std::max(s, e);
  errS = error0;

#if defined(TRACK_STATISTICS)
  if (steps == 5)
    gstat.alpha[1] += errS;
  else if (steps == 7)
    gstat.alpha[4] += errS;
#endif
#endif

  assert((minS >= min) && (minS <= max));
  assert((maxS >= min) && (maxS <= max));
  assert((minS <= maxS));

  return errS;
}

/* -----------------------------------------------------------------------------
 */
static void WriteAlphaBlock(int alpha0, int alpha1, u8 const* indices, void* block)
{
  u8* bytes = reinterpret_cast< u8* >(block);

  // write the first two bytes
  bytes[0] = (u8)alpha0;
  bytes[1] = (u8)alpha1;

  // pack the indices with 3 bits each
  u8* dest = bytes + 2;
  u8 const* src = indices;
  for (int i = 0; i < 2; ++i) {
    // pack 8 3-bit values
    int value = 0;
    for (int j = 0; j < 8; ++j) {
      int index = *src++;
      value += (index << (3 * j));
    }

    // store in 3 bytes
    for (int j = 0; j < 3; ++j) {
      int byte = (value >> (8 * j)) & 0xFF;
      *dest++ = (u8)byte;
    }

    // 77766655.54443332.22111000, FFFEEEDD.DCCCBBBA.AA999888
    // 22111000.54443332.77766655, AA999888.DCCCBBBA.FFFEEEDD
  }
}

static doinline void WriteAlphaBlock5(int alpha0, int alpha1, u8 const* indices, void* block)
{
  // check the relative values of the endpoints
  assert(alpha0 <= alpha1);
  
  if (alpha0 == alpha1) {
    assert((indices[ 0] % 6) < 2); assert((indices[ 1] % 6) < 2);
    assert((indices[ 2] % 6) < 2); assert((indices[ 3] % 6) < 2);
    assert((indices[ 4] % 6) < 2); assert((indices[ 5] % 6) < 2);
    assert((indices[ 6] % 6) < 2); assert((indices[ 7] % 6) < 2);
    assert((indices[ 8] % 6) < 2); assert((indices[ 9] % 6) < 2);
    assert((indices[10] % 6) < 2); assert((indices[11] % 6) < 2);
    assert((indices[12] % 6) < 2); assert((indices[13] % 6) < 2);
    assert((indices[14] % 6) < 2); assert((indices[15] % 6) < 2);
  }

  // write the block
  WriteAlphaBlock(alpha0, alpha1, indices, block);
}

static doinline void WriteAlphaBlock7(int alpha0, int alpha1, u8 const* indices, void* block)
{
  // check the relative values of the endpoints
  assert(alpha0 >= alpha1);
  
  if (alpha0 == alpha1) {
    assert(indices[ 0] < 2); assert(indices[ 1] < 2);
    assert(indices[ 2] < 2); assert(indices[ 3] < 2);
    assert(indices[ 4] < 2); assert(indices[ 5] < 2);
    assert(indices[ 6] < 2); assert(indices[ 7] < 2);
    assert(indices[ 8] < 2); assert(indices[ 9] < 2);
    assert(indices[10] < 2); assert(indices[11] < 2);
    assert(indices[12] < 2); assert(indices[13] < 2);
    assert(indices[14] < 2); assert(indices[15] < 2);
  }

  // write the block
  WriteAlphaBlock(alpha0, alpha1, indices, block);
}

/* -----------------------------------------------------------------------------
 */
template<const int min, const int max, const int prc, const int compress, typename otyp, typename dtyp>
static void CompressAlphaBtc3i(dtyp const* rgba, int mask, void* block, int flags)
{
  Col8 codes5, codes7;

  // get the range for 5-alpha and 7-alpha interpolation
  int min5 = max, max5 = min;
  int min7 = max, max7 = min;

  if (!((flags & kAlphaIterativeFit) && prc)) {
    for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
      // check this pixel is valid
      if ((imask & 1) == 0)
        continue;

      const int fudge = ((1 << compress) >> 1);

      // incorporate into the min/max
      {
        // create 8bit integer value, do round() with truncate() fudged
        int value = (rgba[4 * i + 3] + ((min < 0) && (rgba[4 * i + 3] < 0) ? -fudge : fudge)) >> compress;
        
        if (value < min7)
          min7 = value;
        if (value > max7)
          max7 = value;
        if ((value != min) && (value < min5))
          min5 = value;
        if ((value != max) && (value > max5))
          max5 = value;
   
      }
    }
  }
  else {
    for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
      // check this pixel is valid
      if ((imask & 1) == 0)
        continue;
    
      // calculate fudge value which guarantees "255.0f + fudge" truncate to 255
      const int fudge = ((1 << compress) - 1);

      // incorporate into the min/max
      {
        // create 8bit integer value, do floor() with truncate() fudged
        int value = (rgba[4 * i + 3] - ((min < 0) && (rgba[4 * i + 3] < 0) ? fudge : 0)) >> compress;
    
        if (value < min7)
          min7 = value;
        if ((value != min) && (value < min5))
          min5 = value;
      }
        
      {
        // create 8bit integer value, do ceil() with truncate() fudged
        int value = (rgba[4 * i + 3] + ((min < 0) && (rgba[4 * i + 3] < 0) ? 0 : fudge)) >> compress;
    
        if (value > max7)
          max7 = value;
        if ((value != max) && (value > max5))
          max5 = value;
      }
    }
  }

  // handle the case that no valid range was found
  if (min5 > max5) min5 = max5;
  if (min7 > max7) min7 = max7;

  // expand empty range in 5-code
  if (min5 == max5) min5 -= min5 > min, max5 += max5 < max;

  assert(min5 >= min); assert(max5 <= max);
  assert(min7 >= min); assert(max7 <= max);
  
  // fix the range to be the minimum in each case
//FixRange<min,max,5>(min5, max5);
//FixRange<min,max,5>(min7, max7);

  // construct code-book
  Codebook6<min,max,prc>(codes5, Col8(min5), Col8(max5));
  Codebook8<min,max,prc>(codes7, Col8(max7), Col8(min7));

  // do the iterative tangent search
  if (flags & kAlphaIterativeFit) {
    Scr4 err5, err7; float aaaa[16];

    GetError<prc,otyp>(rgba, mask, codes5, codes7, err5, err7, aaaa);

    // binary search, tangent-fitting
    Scr4 errM = Min(err5, err7);

    // !lossless, !lossless
    if (errM > Scr4(FIT_THRESHOLD)) {
      // if better, reconstruct code-book
      errM = FitError<min,max,prc,7>(aaaa, min7, max7, err7);
      Codebook8<min,max,prc>(codes7, Col8(max7), Col8(min7));

    if (errM > Scr4(FIT_THRESHOLD)) {
      // if better, reconstruct code-book
      errM = FitError<min,max,prc,5>(aaaa, min5, max5, err5);
      Codebook6<min,max,prc>(codes5, Col8(min5), Col8(max5));
    
    }}
  }

  // fit the data to both code books
  u8 indices5[16];
  u8 indices7[16];

  Scr4 err5 = FitCodes<prc,otyp>(rgba, mask, codes5, indices5);
  Scr4 err7 = FitCodes<prc,otyp>(rgba, mask, codes7, indices7);

  // save the block with least error (prefer 7)
  if (err7 > err5)
    WriteAlphaBlock5(min5, max5, indices5, block);
  else
    WriteAlphaBlock7(max7, min7, indices7, block);
}

template<const int min, const int max, const int prc, const int compress, typename otyp, typename dtyp>
static void CompressAlphaBtc3f(dtyp const* rgba, int mask, void* block, int flags)
{
  Col8 codes5, codes7;

  // get the range for 5-alpha and 7-alpha interpolation
  int min5 = max, max5 = min;
  int min7 = max, max7 = min;

  if (!((flags & kAlphaIterativeFit) && prc)) {
    for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
      // check this pixel is valid
      if ((imask & 1) == 0)
        continue;
    
      // incorporate into the min/max
      {
        // create 8bit integer value, do round() with truncate() fudged
        int value = int((rgba[4 * i + 3] * max) + ((min < 0) && (rgba[4 * i + 3] < 0) ? -0.5f : 0.5f));
        
        if (value < min7)
          min7 = value;
        if (value > max7)
          max7 = value;
        if ((value != min) && (value < min5))
          min5 = value;
        if ((value != max) && (value > max5))
          max5 = value;
      }
    }
  }
  else {
    for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
      // check this pixel is valid
      if ((imask & 1) == 0)
        continue;

      // calculate fudge value which guarantees "255.0f + fudge" truncate to 255
      const float fudge = 1.0f - 0.5f / (1 << prc);

      // incorporate into the min/max
      {
        // create 8bit integer value, do floor() with truncate() fudged
        int value = int((rgba[4 * i + 3] * max) - ((min < 0) && (rgba[4 * i + 3] < 0) ? fudge : 0));
        
        if (value < min7)
          min7 = value;
        if ((value != min) && (value < min5))
          min5 = value;
      }
    
      {
        // create 8bit integer value, do ceil() with truncate() fudged
        int value = int((rgba[4 * i + 3] * max) + ((min < 0) && (rgba[4 * i + 3] < 0) ? 0 : fudge));
        
        if (value > max7)
          max7 = value;
        if ((value != max) && (value > max5))
          max5 = value;
      }
    }
  }

  // handle the case that no valid range was found
  if (min5 > max5) min5 = max5;
  if (min7 > max7) min7 = max7;

  // expand empty range in 5-code
  if (min5 == max5) min5 -= min5 > min, max5 += max5 < max;
  
  assert(min5 >= min); assert(max5 <= max);
  assert(min7 >= min); assert(max7 <= max);
  
  // fix the range to be the minimum in each case
//FixRange<min,max,5>(min5, max5);
//FixRange<min,max,7>(min7, max7);

  // construct code-book
  Codebook6<min,max,prc>(codes5, Col8(min5), Col8(max5));
  Codebook8<min,max,prc>(codes7, Col8(max7), Col8(min7));

  // do the iterative tangent search
  if (flags & kAlphaIterativeFit) {
    Scr4 err5, err7; float aaaa[16];

    GetError<prc,otyp,max>(rgba, mask, codes5, codes7, err5, err7, aaaa);

    // binary search, tangent-fitting
    Scr4 errM = Min(err5, err7);

    // !lossless, !lossless
    if (errM > Scr4(FIT_THRESHOLD)) {
      // if better, reconstruct code-book
      errM = FitError<min,max,prc,7>(aaaa, min7, max7, err7);
      Codebook8<min,max,prc>(codes7, Col8(max7), Col8(min7));

    if (errM > Scr4(FIT_THRESHOLD)) {
      // if better, reconstruct code-book
      errM = FitError<min,max,prc,5>(aaaa, min5, max5, err5);
      Codebook6<min,max,prc>(codes5, Col8(min5), Col8(max5));
    
    }}
  }

  // fit the data to both code books
  u8 indices5[16];
  u8 indices7[16];

  Scr4 err5 = FitCodes<prc,otyp,max>(rgba, mask, codes5, indices5);
  Scr4 err7 = FitCodes<prc,otyp,max>(rgba, mask, codes7, indices7);

  // save the block with least error (prefer 7)
  if (err7 > err5)
    WriteAlphaBlock5(min5, max5, indices5, block);
  else
    WriteAlphaBlock7(max7, min7, indices7, block);
}

void CompressAlphaBtc3u(u8  const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<   0,255, CBLB, 0, unsigned>(rgba, mask, block, flags); }
void CompressAlphaBtc3s(s8  const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<-127,127, CBLB, 0,   signed>(rgba, mask, block, flags); }

void CompressAlphaBtc3u(u16 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<   0,255, CBLB, 8, unsigned>(rgba, mask, block, flags); }
void CompressAlphaBtc3s(s16 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<-127,127, CBLB, 8,   signed>(rgba, mask, block, flags); }

void CompressAlphaBtc3u(f23 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3f<   0,255, CBLB, 0, unsigned>(rgba, mask, block, flags); }
void CompressAlphaBtc3s(f23 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3f<-127,127, CBLB, 0,   signed>(rgba, mask, block, flags); }

void CompressDepthBtc4u(u8  const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<   0,255, CBLB, 0, unsigned>(rgba, mask, block, flags); }
void CompressDepthBtc4s(s8  const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<-127,127, CBLB, 0,   signed>(rgba, mask, block, flags); }

void CompressDepthBtc4u(u16 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<   0,255, CBHB, 8, unsigned>(rgba, mask, block, flags); }
void CompressDepthBtc4s(s16 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3i<-127,127, CBHB, 8,   signed>(rgba, mask, block, flags); }

void CompressDepthBtc4u(f23 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3f<   0,255, CBHB, 0, unsigned>(rgba, mask, block, flags); }
void CompressDepthBtc4s(f23 const* rgba, int mask, void* block, int flags) {
  CompressAlphaBtc3f<-127,127, CBHB, 0,   signed>(rgba, mask, block, flags); }

#undef	FIT_THRESHOLD

/* *****************************************************************************
 */
template<const int prc, typename ctyp, typename etyp>
static void ReadAlphaBlock(
  ctyp (&codes  )[ 8],
  u8   (&indices)[16],
  void const* block)
{
  // get the two alpha values
  etyp const* bytes = reinterpret_cast< etyp const* >(block);
  etyp alpha0 = bytes[0];
  etyp alpha1 = bytes[1];

  // compare the values to build the codebook
  codes[0] = (ctyp)alpha0 << prc;
  codes[1] = (ctyp)alpha1 << prc;

  // use 5-alpha or 7-alpha codebook
  Codebook6or8<prc>(codes, alpha0 <= alpha1);

  // decode the indices
  u8 const* src = (u8*)bytes + 2;
  u8* dest = indices;
  for (int i = 0; i < 2; ++i) {
    // grab 3 bytes
    int value = 0;
    for (int j = 0; j < 3; ++j) {
      int byte = *src++;
      value += (byte << 8 * j);
    }

    // unpack 8 3-bit values from it
    for (int j = 0; j < 8; ++j) {
      int index = (value >> 3 * j) & 0x7;
      *dest++ = (u8)index;
    }
  }
}

template<const int prc, typename dtyp, typename ctyp, typename etyp, const int scale, const int correction>
static void DecompressAlphaBtc3i(dtyp* rgba, void const* block)
{
  ctyp codes[8]; u8 indices[16];

  ReadAlphaBlock<prc,ctyp,etyp>(codes, indices, block);

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i)
    rgba[4 * i + 3] = dtyp((codes[indices[i]] * (scale)) >> (correction));
}

template<const int prc, typename dtyp, typename ctyp, typename etyp, const int scale>
static void DecompressAlphaBtc3f(dtyp* rgba, void const* block)
{
  ctyp codes[8]; u8 indices[16];

  ReadAlphaBlock<prc,ctyp,etyp>(codes, indices, block);

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i)
    rgba[4 * i + 3] = dtyp((codes[indices[i]]) * (1.0f / (scale)));
}

void DecompressAlphaBtc3u(u8* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBLB, u8 , u8 ,u8, 255 / 255, 0>(rgba, block); }
void DecompressAlphaBtc3s(s8* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBLB, s8 , s8 ,s8, 127 / 127, 0>(rgba, block); }

void DecompressAlphaBtc3u(u16* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBLB, u16, u8 ,u8, 65535 / 255, 0>(rgba, block); }
void DecompressAlphaBtc3s(s16* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBLB, s16, s8 ,s8, 0x1FFFFF / 127, 6>(rgba, block); }

void DecompressAlphaBtc3u(f23* rgba, void const* block, int flags) {
  DecompressAlphaBtc3f<CBLB, f23, u8 ,u8, 255>(rgba, block); }
void DecompressAlphaBtc3s(f23* rgba, void const* block, int flags) {
  DecompressAlphaBtc3f<CBLB, f23, s8 ,s8, 127>(rgba, block); }

void DecompressDepthBtc4u(u8* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBLB, u8 , u16,u8, 255 / 255, 0>(rgba, block); }
void DecompressDepthBtc4s(s8* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBLB, s8 , s16,s8, 127 / 127, 0>(rgba, block); }

void DecompressDepthBtc4u(u16* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBHB, u16, u16,u8, 65535 / 255, CBHB>(rgba, block); }
void DecompressDepthBtc4s(s16* rgba, void const* block, int flags) {
  DecompressAlphaBtc3i<CBHB, s16, s16,s8, 0x1FFFFF / 127, CBHB + 6>(rgba, block); }

void DecompressDepthBtc4u(f23* rgba, void const* block, int flags) {
  DecompressAlphaBtc3f<CBHB, f23, u16,u8, 255 * (1 << CBHB)>(rgba, block); }
void DecompressDepthBtc4s(f23* rgba, void const* block, int flags) {
  DecompressAlphaBtc3f<CBHB, f23, s16,s8, 127 * (1 << CBHB)>(rgba, block); }
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
/* C++ AMP version */
static void FixRange(out lineA2 aline, int steps) amp_restricted
{
  if ((aline[ASTOP] - aline[ASTRT]) < steps)
    aline[ASTOP] = ((aline[ASTRT] + steps) < 255 ? (aline[ASTRT] + steps) : 255);
  if ((aline[ASTOP] - aline[ASTRT]) < steps)
    aline[ASTRT] = (  0 > (aline[ASTOP] - steps) ? (aline[ASTOP] - steps) :   0);
}

#if	defined(SQUISH_USE_COMPUTE)
  tile_static int cerror;
#endif

/* C++ AMP version */
static int FitCodes(tile_barrier barrier, const int thread, pixel16 rgba, int mask, index8 codes, out index16 matched) amp_restricted
{
  // fit each alpha value to the codebook (AMP: prefer vectorization over parallelism)
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static int cerror;
#endif

  threaded_cse(0) {
    cerror = 0; }
  wavefrnt_for(fcscan, 16) {
    // use the first code by default
    matched[fcscan] = 0;

    // check this pixel is valid
    int bit = 1 << fcscan;
    if ((mask & bit) != 0) {
      // find the least error and corresponding index
      const int value = rgba[fcscan][0];
      int dists[8];

      dists[0] = value - (int)codes[0];
      dists[1] = value - (int)codes[1];
      dists[2] = value - (int)codes[2];
      dists[3] = value - (int)codes[3];

      dists[4] = value - (int)codes[4];
      dists[5] = value - (int)codes[5];
      dists[6] = value - (int)codes[6];
      dists[7] = value - (int)codes[7];

      // find the closest code (AMP: vectorized reduction, cset)
      int4 idx07, idx;

      idx07.x = 1 * (0 + (dists[1] < dists[0] ? 1 : 0));
      idx07.y = 1 * (2 + (dists[3] < dists[2] ? 1 : 0));
      idx07.z = 1 * (4 + (dists[5] < dists[4] ? 1 : 0));
      idx07.w = 1 * (6 + (dists[7] < dists[6] ? 1 : 0));

      idx.x = (dists[idx07.x] < dists[idx07.y] ? idx07.x : idx07.y);
      idx.y = (dists[idx07.z] < dists[idx07.w] ? idx07.z : idx07.w);

      idx.z = (dists[idx.x] < dists[idx.y] ? idx.x : idx.y);

      // save the index
      matched[fcscan] = (ccr8)idx.z;

      // accumulate the error
      threaded_add(cerror, dists[idx.z]);
    }
  }

  // return the total error
  return cerror;
}

/* C++ AMP version */
static void WriteAlphaBlock(tile_barrier barrier, const int thread, lineA2 alpha, index16 indices, out code64 block) amp_restricted
{
  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  threaded_cse(0) {
    // write the first two bytes
    block[0] =
      (alpha  [ASTRT] <<  0) +
      (alpha  [ASTOP] <<  8) +
      (indices[ 0] << 16) + (indices[ 1] << 19) +
      (indices[ 2] << 22) + (indices[ 3] << 25) +
      (indices[ 4] << 28) + (indices[ 5] << 31);

    // pack the indices with 3 bits each
    block[1] =
      (indices[ 5] >>  1) +
      (indices[ 6] <<  2) + (indices[ 7] <<  5) +
      (indices[ 8] <<  8) + (indices[ 9] << 11) +
      (indices[10] << 14) + (indices[11] << 17) +
      (indices[12] << 20) + (indices[13] << 23) +
      (indices[14] << 26) + (indices[15] << 29);
  }
}

/* C++ AMP version */
static void WriteAlphaBlock5(tile_barrier barrier, const int thread, lineA2 alpha, inout index16 indices, out code64 block,
			     IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 5, 4, 3, 2, 6, 7 };	// (a >  b)

  // degenerate case "start > stop" (AMP: local register, not enough for group-shared)
  int sorted[AVALS];

  // alpha[ASTRT] > alpha[ASTOP] := sorted[ASTRT] == alpha[ASTOP]
  sorted[ASTRT] = min(alpha[ASTRT], alpha[ASTOP]);
  sorted[ASTOP] = max(alpha[ASTRT], alpha[ASTOP]);

  // swap the indices
  wavefrnt_for(i, 16) {
#if	defined(SQUISH_USE_COMPUTE)
    indices[i] = (sorted[ASTRT] == alpha[ASTOP] ? yArr[IBL_ALPHA5][indices[i]].mapped : indices[i]);
#else
    indices[i] = (sorted[ASTRT] == alpha[ASTOP] ? yArr(IBL_ALPHA5, indices[i]).mapped : indices[i]);
#endif
  }

  // write the block
  WriteAlphaBlock(barrier, thread, sorted, indices, block);
}

/* C++ AMP version */
static void WriteAlphaBlock7(tile_barrier barrier, const int thread, lineA2 alpha, inout index16 indices, out code64 block,
			     IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 7, 6, 5, 4, 3, 2 };	// (a <  b)

  // degenerate case "start < stop" (AMP: local register, not enough for group-shared)
  int sorted[AVALS];

  // alpha[ASTRT] < alpha[ASTOP] := sorted[ASTRT] == alpha[ASTRT]
  sorted[ASTRT] = max(alpha[ASTRT], alpha[ASTOP]);
  sorted[ASTOP] = min(alpha[ASTRT], alpha[ASTOP]);

  // swap the indices
  wavefrnt_for(i, 16) {
#if	defined(SQUISH_USE_COMPUTE)
    indices[i] = (sorted[ASTRT] == alpha[ASTRT] ? yArr[IBL_ALPHA7][indices[i]].mapped : indices[i]);
#else
    indices[i] = (sorted[ASTRT] == alpha[ASTRT] ? yArr(IBL_ALPHA7, indices[i]).mapped : indices[i]);
#endif
  }

  // write the block
  WriteAlphaBlock(barrier, thread, sorted, indices, block);
}

#define	ERR5	0
#define	ERR7	1
#define	ERRS	2
#define CNUMS	4
#define CBITS	2	// number of bits
#define CMASK	3	// mask

#if	defined(SQUISH_USE_COMPUTE)
  tile_static ccr8 acodes [ERRS][ 8];
  tile_static ccr8 matched[ERRS][16];
  tile_static int aline[ERRS][AVALS];
//tile_static int error[ERRS];
  tile_static int error[4];
#endif

/* C++ AMP version */
void CompressAlphaBtc3u(tile_barrier barrier, const int thread, pixel16 rgba, int mask, out code64 block,
		       IndexBlockLUT yArr) amp_restricted
{
  // get the range for 5-alpha and 7-alpha interpolation
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static ccr8 acodes [ERRS][ 8];
  tile_static ccr8 matched[ERRS][16];
  tile_static int aline[ERRS][AVALS];
  tile_static int error[ERRS];
#endif

  // same for all
  const int ewch = thread & 1;
  const int cpos = thread >> 1;

  /* initialize:
   *
   *  aline[ERR7][ASTRT] = 255;
   *  aline[ERR7][ASTOP] = 0;
   *  aline[ERR5][ASTRT] = 255;
   *  aline[ERR5][ASTOP] = 0;
   */
  threaded_for(li, CNUMS) {
    aline[ewch][cpos] = (cpos ? 0 : 255);
  }

  wavefrnt_for(mmscan, 16) {
    // check this pixel is valid
    int bit = 1 << mmscan;
    if ((mask & bit) != 0) {
      // incorporate into the min/max
      const int value = rgba[mmscan][0];

      threaded_min(aline[ERR7][ASTRT], value);
      threaded_max(aline[ERR7][ASTOP], value);
      if ((value != 0) && (value != 255)) {
	threaded_min(aline[ERR5][ASTRT], value);
	threaded_max(aline[ERR5][ASTOP], value);
      }
    }
  }

  /* fix the range to be the minimum in each case
   *
   *  FixRange(aline[ERR5], 5);
   *  FixRange(aline[ERR7], 7);
   */
  threaded_for(ec, ERRS) {
    const int step = 5 + (ec * 2);

    /* handle the case that no valid range was found
     *
     *  threaded_min(aline[ERR5][ASTRT], aline[ERR5][ASTOP]);
     *  threaded_min(aline[ERR7][ASTRT], aline[ERR7][ASTOP]);
     */
    threaded_min(aline[ewch][ASTRT], aline[ewch][ASTOP]);

    FixRange(aline[ewch], step);
  }

  /* set up the 5-alpha code book
   * set up the 7-alpha code book
   *
   *  acodes[ERR5][0] = (ccr8)aline[ERR5][ASTRT];
   *  acodes[ERR5][1] = (ccr8)aline[ERR5][ASTOP];
   *
   *  acodes[ERR7][0] = (ccr8)aline[ERR7][ASTRT];
   *  acodes[ERR7][1] = (ccr8)aline[ERR7][ASTOP];
   */
  threaded_for(ei, CNUMS) {
    acodes[ewch][cpos] = (ccr8)aline[ewch][cpos];
  }

  /* set up the 5-alpha code book
   * set up the 7-alpha code book
   *
   *  acodes[ERR5][0] = (ccr8)aline[ERR5][ASTRT];
   *  acodes[ERR5][1] = (ccr8)aline[ERR5][ASTOP];
   *
   *  acodes[ERR7][0] = (ccr8)aline[ERR7][ASTRT];
   *  acodes[ERR7][1] = (ccr8)aline[ERR7][ASTOP];
   *
   *  acodes[ERR5][6] = 0;
   *  acodes[ERR5][7] = 255;
   */
  threaded_for(es, 6) {
//  const int ewch = (es / 6);
//  const int step = 5 + (ewch * 2);

    // set up the 5-alpha code book
    acodes[ERR5][2 + es] = (ccr8)(((4 - es) * aline[ERR5][ASTRT] + (1 + es) * aline[ERR5][ASTOP]) / 5) & (es == 4 ? 0 : 255) + (es == 5 ? 255 : 0);
    // set up the 7-alpha code book
    acodes[ERR7][2 + es] = (ccr8)(((6 - es) * aline[ERR7][ASTRT] + (1 + es) * aline[ERR7][ASTOP]) / 7);
  }

  // make it visible
  tile_static_memory_fence(barrier);

  /* fit the data to both code books
   *
   *  erros[ERR5] = FitCodes(barrier, thread, rgba, mask, acodes[ERR5], matched[ERR5]);
   *  erros[ERR7] = FitCodes(barrier, thread, rgba, mask, acodes[ERR7], matched[ERR7]);
   */
  threaded_for(ef, ERRS) {
    error[ewch] = FitCodes(barrier, thread, rgba, mask, acodes[ewch], matched[ewch]);
  }

  // AMP: final reduction (make the choice)
  {
    // save the block with least erros
    if (error[ERR5] <= error[ERR7])
      WriteAlphaBlock5(barrier, thread, aline[ERR5], matched[ERR5], block, yArr);
    else
      WriteAlphaBlock7(barrier, thread, aline[ERR7], matched[ERR7], block, yArr);
  }
}

#undef	ERR5
#undef	ERR7
#undef	ERRS
#undef	CNUMS
#undef	CBITS
#undef	CMASK

#if	defined(SQUISH_USE_COMPUTE)
  tile_static unsigned int frags[4];
#endif

/* C++ AMP version */
void CompressAlphaBtc2u(tile_barrier barrier, const int thread, pixel16 rgba, int mask, code64 block,
		         IndexBlockLUT yArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static unsigned int frags[4];
#endif

  // same for all
  const int opos = thread << 2;
  const int hilo = thread >> 1;

  // quantize and pack the alpha values quad-wise (AMP: prefer vectorization over parallelism)
  threaded_for(qs, 4) {
    // quantize down to 4 bits
    float4 alpha = float4(
      (float)rgba[opos + 0][0],
      (float)rgba[opos + 1][0],
      (float)rgba[opos + 2][0],
      (float)rgba[opos + 3][0]
    );

    int4 quant = QuantizeFloatToInt(alpha * 1.0f / 255.0f, 15);

    // set alpha to zero where masked
    quant.x = ((mask >> (opos + 0)) & 1) ? quant.x : 0;
    quant.y = ((mask >> (opos + 1)) & 1) ? quant.y : 0;
    quant.z = ((mask >> (opos + 2)) & 1) ? quant.z : 0;
    quant.w = ((mask >> (opos + 3)) & 1) ? quant.w : 0;

    // pack into the short
    frags[qs] = (
      (quant.x <<  0) |
      (quant.y <<  4) |
      (quant.z <<  8) |
      (quant.w << 12)
    ) << (hilo * 16);
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // AMP: final reduction (make the merge)
  threaded_cse(0) {
    // save the block
    block[0] = frags[1] + frags[0];
    block[1] = frags[3] + frags[2];
  }
}
#endif

} // namespace squish

#undef CBLB
#undef CBHB
