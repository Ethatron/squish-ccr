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

#include "paletteindexfit.h"
#include "paletteset.h"
#include "paletteblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
PaletteIndexFit::PaletteIndexFit(PaletteSet const* palette, int flags, int swap, int shared)
  : PaletteFit(palette, flags, swap, shared)
{
}

#ifdef	FEATURE_INDEXFIT_INLINED
void PaletteIndexFit::ErrorEndPoints(int set, Vec4 const &metric, vQuantizer &q, int sb, u8 (&closest)[16],
				     Vec4 const* values, Scr4 const* freq,
				     int ib, int idxs) {
  // snap floating-point-values to the integer-lattice
  Vec4 start = m_qstart = q.SnapToLatticeClamped(values[0], sb, 1 << SBSTART);
  Vec4 end   = m_qend   = q.SnapToLatticeClamped(values[1], sb, 1 << SBEND);

  // error-sum of the two values
  Scr4 serror = LengthSquared(metric * (values[0] - start)) * freq[0];
  Scr4 eerror = LengthSquared(metric * (values[1] - end  )) * freq[1];
      
  // return combined error
  m_qerror = Vec4(serror + eerror, serror, eerror);
  
  // final check if it really improved
  {
    m_berror = Scr4(m_qerror);

    // save the index
    closest[0] = 0; 
    closest[1] = (u8)idxs;

    // save the end-points
    m_start[set] = values[0]; 
    m_end  [set] = values[1];
  }
}

Scr4 PaletteIndexFit::ErrorInterpolants(Vec4 const &metric, vQuantizer &q, int sb,
					Vec4 const* values, Scr4 const* freq,
					int ib, int idxs, Vec4 &value0, Vec4 &value1, int closest0, int closest1) {
  // snap floating-point-values to the integer-lattice
  Vec4 start = q.SnapToLatticeClamped(value0   , sb, 1 << SBSTART);
  Vec4 end   = q.SnapToLatticeClamped(value1   , sb, 1 << SBEND);
    
  // created interpolated values
  Vec4 scode = weights_V4[ib][idxs - closest0] * start + weights_V4[ib][closest0] * end;
  Vec4 ecode = weights_V4[ib][idxs - closest1] * start + weights_V4[ib][closest1] * end;

  // error-sum of the two values
  Scr4 serror = LengthSquared(metric * (values[0] - scode)) * freq[0];
  Scr4 eerror = LengthSquared(metric * (values[1] - ecode)) * freq[1];
      
  // return combined error
  return (serror + eerror);
}

Scr4 PaletteIndexFit::ErrorInterpolantsS(Vec4 const &metric, vQuantizer &q, int sb,
					 Vec4 const* values, Scr4 const* freq,
					 int ib, int idxs, Vec4 &value0, int closest0) {
  // snap floating-point-values to the integer-lattice
  Vec4 start = q.SnapToLatticeClamped(value0   , sb, 1 << SBSTART);
  Vec4 end   = m_qend;
    
  // created interpolated values
  Vec4 scode = weights_V4[ib][idxs - closest0] * start + weights_V4[ib][closest0] * end;

  // error-sum of the two values
  Scr4 serror = LengthSquared(metric * (values[0] - scode)) * freq[0];
  Scr4 eerror = Scr4(m_qerror.SplatZ());
      
  // return combined error
  return (serror + eerror);
}

Scr4 PaletteIndexFit::ErrorInterpolantsE(Vec4 const &metric, vQuantizer &q, int sb,
					 Vec4 const* values, Scr4 const* freq,
					 int ib, int idxs, Vec4 &value1, int closest1) {
  // snap floating-point-values to the integer-lattice
  Vec4 start = m_qstart;
  Vec4 end   = q.SnapToLatticeClamped(value1   , sb, 1 << SBEND);
    
  // created interpolated values
  Vec4 ecode = weights_V4[ib][idxs - closest1] * start + weights_V4[ib][closest1] * end;

  // error-sum of the two values
  Scr4 serror = Scr4(m_qerror.SplatY());
  Scr4 eerror = LengthSquared(metric * (values[1] - ecode)) * freq[1];
      
  // return combined error
  return (serror + eerror);
}

void PaletteIndexFit::BetterInterpolants(int set, Vec4 const &metric, vQuantizer &q, int sb, u8 (&closest)[16],
					 Vec4 const* values, Scr4 const* freq,
				         int ib, int idxs, Vec4 &value0, Vec4 &value1, int closest0, int closest1) {
  Scr4 nerror = ErrorInterpolants(metric, q, sb, values, freq, ib, idxs, value0, value1, closest0, closest1);
  Scr4 berror = Min(nerror, m_berror);
  
  // final check if it really improved
  if (berror == nerror) {
    m_berror = nerror;

    // save the index
    closest[0] = (u8)closest0; 
    closest[1] = (u8)closest1;

    // save the end-points
    m_start[set] = value0; 
    m_end  [set] = value1;
  }
}
      
void PaletteIndexFit::BetterInterpolantsS(int set, Vec4 const &metric, vQuantizer &q, int sb, u8 (&closest)[16],
					 Vec4 const* values, Scr4 const* freq,
				         int ib, int idxs, Vec4 &value0, int closest0) {
  Scr4 nerror = ErrorInterpolantsS(metric, q, sb, values, freq, ib, idxs, value0, closest0);
  Scr4 berror = Min(nerror, m_berror);
  
  // final check if it really improved
  if (berror == nerror) {
    m_berror = nerror;

    // save the index
    closest[0] = (u8)closest0; 
    closest[1] = (1 << ib) - 1;

    // save the end-points
    m_start[set] = value0; 
    m_end  [set] = values[1];
  }
}

void PaletteIndexFit::BetterInterpolantsE(int set, Vec4 const &metric, vQuantizer &q, int sb, u8 (&closest)[16],
					 Vec4 const* values, Scr4 const* freq,
				         int ib, int idxs, Vec4 &value1, int closest1) {
  Scr4 nerror = ErrorInterpolantsE(metric, q, sb, values, freq, ib, idxs, value1, closest1);
  Scr4 berror = Min(nerror, m_berror);
  
  // final check if it really improved
  if (berror == nerror) {
    m_berror = nerror;

    // save the index
    closest[0] = 0; 
    closest[1] = (u8)closest1;

    // save the end-points
    m_start[set] = values[0]; 
    m_end  [set] = value1;
  }
}

#define ChooseInterpolants(_closest0, _value0, _closest1, _value1)	\
  BetterInterpolants(							\
    set, metric, q, sb, closest, values, freq, ib, idxs,		\
    _value0, _value1, _closest0, _closest1);				\
  if (!FEATURE_INDEXFIT_THOROUGH) break;

#define ChooseInterpolantsS(_closest0, _value0)				\
  BetterInterpolantsS(							\
    set, metric, q, sb, closest, values, freq, ib, idxs,		\
    _value0, _closest0);						\
  if (!FEATURE_INDEXFIT_THOROUGH) break;

#define ChooseInterpolantsE(_closest1, _value1)				\
  BetterInterpolantsE(							\
    set, metric, q, sb, closest, values, freq, ib, idxs,		\
    _value1, _closest1);						\
  if (!FEATURE_INDEXFIT_THOROUGH) break;
#else
#define ChooseInterpolants(_closest0, _value0, _closest1, _value1)	\
  closest0 = _closest0; value0 = _value0;				\
  closest1 = _closest1; value1 = _value1;				\
  break;

#define ChooseInterpolantsS(_closest0, _value0)				\
  closest0 = _closest0; value0 = _value0;				\
  break;

#define ChooseInterpolantsE(_closest1, _value1)				\
  closest1 = _closest1; value1 = _value1;				\
  break;
#endif

Scr4 PaletteIndexFit::StretchEndPoints(int set, Vec4 const &metric, vQuantizer &q, int sb, int ib, u8 (&closest)[16])
{
  // cache some values
  Vec4 const* values = m_palette->GetPoints(set);
  Scr4 const* freq = m_palette->GetWeights(set);
  
  // hold these in SSE-registers
  int gap, idxs = (1 << ib) - 1;
  Col4 igap, iidxs; iidxs.SetRGBApow2<0>(ib); iidxs = iidxs - Col4(1);

#ifdef FEATURE_INDEXFIT_INLINED
  // get the base error
  ErrorEndPoints(set, metric, q, sb, closest, values, freq, ib, idxs);
#endif

  // the idea is that quantization noise has less impact on the inner most
  // interpolated values, that is, the angle-deviation becomes smaller the
  // longer the line is
  Vec4 value0 = values[0];
  Vec4 value1 = values[1];
  Vec4 delta  = value1 - value0;

  int closest0 = 0;
  int closest1 = idxs;
      
  // try to find plain expansions into both directions first
  for (gap = 1, igap = Col4(1); gap < idxs; gap += 2, igap += Col4(2)) {
    // 3-1 >> 1 -> 1, 1/1
    // 7-1 >> 1 -> 3, 3/1
    // 7-3 >> 1 -> 2, 2/3
    // 7-5 >> 1 -> 1, 1/5
    // ...
    // value0 - delta * ((float)exp / gap);
    // value1 + delta * ((float)exp / gap);
    
    Col4 iexp = (iidxs - igap) >> 1;
    Vec4 fexp = Vec4(iexp);
    Vec4 fgap = Reciprocal(Vec4(igap));
    Vec4 stretch = delta * fexp * fgap;

    {
      Vec4 expA0 = value0 - stretch;
      Vec4 expA1 = value1 + stretch;

      Vec4 minA = Min(expA0, expA1);
      Vec4 maxA = Max(expA0, expA1);

      // expansion into both directions works
      if (!(CompareAnyGreaterThan(maxA, Vec4(1.0f + (0.25f / 255.0f))) ||
	    CompareAnyLessThan   (minA, Vec4(     -  0.25f / 255.0f)))) {
	int exp = (idxs - gap) >> 1;
	
	ChooseInterpolants(
	  closest0 + exp, expA0,
	  closest1 - exp, expA1
	);
      }
    }
  }

  /* nothing found */
  if (FEATURE_INDEXFIT_THOROUGH || (gap >= idxs)) {
    //  3>>1 -> 2/1, 1/2
    //  7>>1 -> 4/3, 3/4
    // 15>>1 -> 8/7, 7/8
    // ...
    // int eshort = (idxs + 0) >> 1;
    // int elong  = (idxs + 1) >> 1;
	  
    Col4 isht = (iidxs) >> 1;
    Col4 ilng = (iidxs - isht);
    
    do {
      Vec4 elong = Vec4(ilng);
      Vec4 eshort = Reciprocal(Vec4(isht));

      // ---------------------------------------------------------------------------------------
      Vec4 factorl = elong * eshort;
      Vec4 stretchl = delta * factorl;
      
      // expA0 = value0 - delta * ((float)elong / eshort);
      // expA1 = value1 + delta * ((float)elong / eshort);
      Vec4 expA0 = value0 - stretchl;
      Vec4 expA1 = value1 + stretchl;

      // expansion down works
      if (!CompareAnyGreaterThan(Abs(Vec4(0.5f) - expA0), Vec4(0.5f + (0.25f / 255.0f)))) {
	int eshort = (idxs + 0) >> 1; ChooseInterpolantsS(closest1 - eshort, expA0); }
      // expansion up works
      if (!CompareAnyGreaterThan(Abs(Vec4(0.5f) - expA1), Vec4(0.5f + (0.25f / 255.0f)))) {
	int eshort = (idxs + 0) >> 1; ChooseInterpolantsE(closest0 + eshort, expA1); }
      
      // ---------------------------------------------------------------------------------------
      Vec4 factors = Reciprocal(factorl);
      Vec4 stretchs = delta * factors;

      // expA0 = value0 - delta * ((float)eshort / elong);
      // expA1 = value1 + delta * ((float)eshort / elong);
      expA0 = value0 - stretchs;
      expA1 = value1 + stretchs;

      // expansion down works
      if (!CompareAnyGreaterThan(Abs(Vec4(0.5f) - expA0), Vec4(0.5f + (0.25f / 255.0f)))) {
	int elong  = (idxs + 1) >> 1; ChooseInterpolantsS(closest1 - elong, expA0); }
      // expansion up works
      if (!CompareAnyGreaterThan(Abs(Vec4(0.5f) - expA1), Vec4(0.5f + (0.25f / 255.0f)))) {
	int elong  = (idxs + 1) >> 1; ChooseInterpolantsE(closest0 + elong, expA1); }
    } while(0);
  }
  
#ifndef FEATURE_INDEXFIT_INLINED
  // snap floating-point-values to the integer-lattice
  Vec4 vstart = q.SnapToLatticeClamped(value0   , sb, 1 << SBSTART);
  Vec4 vend   = q.SnapToLatticeClamped(value1   , sb, 1 << SBEND);
  Vec4 estart = q.SnapToLatticeClamped(values[0], sb, 1 << SBSTART);
  Vec4 eend   = q.SnapToLatticeClamped(values[1], sb, 1 << SBEND);
    
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec4 vcode0 = weights_V4[ib][idxs - closest0] * vstart + weights_V4[ib][closest0] * vend;
  Vec4 vcode1 = weights_V4[ib][idxs - closest1] * vstart + weights_V4[ib][closest1] * vend;

  // error-sum of the two values
  Vec4 verror0 = LengthSquared(metric * (values[0] - vcode0)) * freq[0];
  Vec4 verror1 = LengthSquared(metric * (values[1] - vcode1)) * freq[1];
  Vec4 eerror0 = LengthSquared(metric * (values[0] - estart)) * freq[0];
  Vec4 eerror1 = LengthSquared(metric * (values[1] - eend  )) * freq[1];
      
  // final check if it really improved
  Vec4 verror = (verror0 + verror1);
  Vec4 eerror = (eerror0 + eerror1);
  Vec4 merror = Min(verror, eerror);
      
  // find the best code
  if (merror == eerror) {
    closest0 = (u8)   0;
    closest1 = (u8)idxs;

    value0 = values[0];
    value1 = values[1];
  }
      
  // save the index
  closest[0] = (u8)closest0; 
  closest[1] = (u8)closest1;

  // save the end-points
  m_start[set] = value0; 
  m_end  [set] = value1;

#if defined(TRACK_STATISTICS)
  gstat.btr_index[ib][merror == verror ? 1 : 0]++;
  gstat.err_index[ib][0] += merror.X();
  gstat.err_index[ib][1] += verror.X();
#endif

  // accumulate the error
  return merror;
#else

#if defined(TRACK_STATISTICS)
  gstat.btr_index[ib][m_qerror == m_berror ? 0 : 1]++;
  gstat.err_index[ib][0] += m_qerror.X();
  gstat.err_index[ib][1] += m_berror.X();
#endif

  // accumulate the error
  return m_berror;
#endif
}

#undef	ChooseInterpolants
#undef	ChooseInterpolantsS
#undef	ChooseInterpolantsE
#endif

} // namespace squish
