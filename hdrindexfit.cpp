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

#include "hdrindexfit.h"
#include "hdrset.h"
#include "hdrblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
HDRIndexFit::HDRIndexFit(HDRSet const* palette, int flags)
  : HDRFit(palette, flags)
{
}

#ifdef	FEATURE_INDEXFIT_INLINED
void HDRIndexFit::ErrorEndPoints(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
				     Vec3 const* values, Scr3 const* freq,
				     int ib, int idxs) {
  // snap floating-point-values to the integer-lattice
  Vec3 val[2];
  Col3 res[2];
  
  val[0] = values[0];
  val[1] = values[1];

  q.QuantizeToLattice(val, res);

  m_qstart = res[0];
  m_qend   = res[1];

  Vec3 start = q.UnquantizeFromLattice(m_qstart);
  Vec3 end   = q.UnquantizeFromLattice(m_qend);

  // error-sum of the two values
  Scr3 serror = LengthSquared(metric * (values[0] - start)) * freq[0];
  Scr3 eerror = LengthSquared(metric * (values[1] - end  )) * freq[1];
      
  // return combined error
  m_qerror = Vec3(serror + eerror, serror, eerror);
  
  // final check if it really improved
  {
    m_berror = Scr3(m_qerror);

    // save the index
    closest[0] = 0; 
    closest[1] = (u8)idxs;

    // save the end-points
    m_start[set] = values[0]; 
    m_end  [set] = values[1];
  }
}

Scr3 HDRIndexFit::ErrorInterpolants(Vec3 const &metric, fQuantizer &q,
					Vec3 const* values, Scr3 const* freq,
					int ib, int idxs, Vec3 &value0, Vec3 &value1, int closest0, int closest1) {
  // snap floating-point-values to the integer-lattice
  Vec3 val[2];
  Col3 res[2];
  
  val[0] = value0;
  val[1] = value1;

  q.QuantizeToLattice(val, res);

  Col3 start = res[0];
  Col3 end   = res[1];
    
  // created interpolated values
  Col3 scode = (Mul16x16u(weights_C4[ib][idxs - closest0].GetCol3(), start) + Mul16x16u(weights_C4[ib][closest0].GetCol3(), end) + Col3(32)) >> 6;
  Col3 ecode = (Mul16x16u(weights_C4[ib][idxs - closest1].GetCol3(), start) + Mul16x16u(weights_C4[ib][closest1].GetCol3(), end) + Col3(32)) >> 6;

  Vec3 sval = q.UnquantizeFromLattice(scode);
  Vec3 eval = q.UnquantizeFromLattice(ecode);

  // error-sum of the two values
  Scr3 serror = LengthSquared(metric * (values[0] - sval)) * freq[0];
  Scr3 eerror = LengthSquared(metric * (values[1] - eval)) * freq[1];
      
  // return combined error
  return (serror + eerror);
}

Scr3 HDRIndexFit::ErrorInterpolantsS(Vec3 const &metric, fQuantizer &q,
					 Vec3 const* values, Scr3 const* freq,
					 int ib, int idxs, Vec3 &value0, int closest0) {
  // snap floating-point-values to the integer-lattice
  Col3 start = q.QuantizeToLattice(value0, m_qend);
  Col3 end   = m_qend;
    
  // created interpolated values
  Col3 scode = (Mul16x16u(weights_C4[ib][idxs - closest0].GetCol3(), start) + Mul16x16u(weights_C4[ib][closest0].GetCol3(), end) + Col3(32)) >> 6;
  Vec3 sval  = q.UnquantizeFromLattice(scode);

  // error-sum of the two values
  Scr3 serror = LengthSquared(metric * (values[0] - sval)) * freq[0];
  Scr3 eerror = Scr3(m_qerror.SplatZ());
      
  // return combined error
  return (serror + eerror);
}

Scr3 HDRIndexFit::ErrorInterpolantsE(Vec3 const &metric, fQuantizer &q,
					 Vec3 const* values, Scr3 const* freq,
					 int ib, int idxs, Vec3 &value1, int closest1) {
  // snap floating-point-values to the integer-lattice
  Col3 start = m_qstart;
  Col3 end   = q.QuantizeToLattice(value1, m_qstart);
    
  // created interpolated values
  Col3 ecode = (Mul16x16u(weights_C4[ib][idxs - closest1].GetCol3(), start) + Mul16x16u(weights_C4[ib][closest1].GetCol3(), end) + Col3(32)) >> 6;
  Vec3 eval  = q.UnquantizeFromLattice(ecode);

  // error-sum of the two values
  Scr3 serror = Scr3(m_qerror.SplatY());
  Scr3 eerror = LengthSquared(metric * (values[1] - eval)) * freq[1];
      
  // return combined error
  return (serror + eerror);
}

void HDRIndexFit::BetterInterpolants(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
					 Vec3 const* values, Scr3 const* freq,
				         int ib, int idxs, Vec3 &value0, Vec3 &value1, int closest0, int closest1) {
  Scr3 nerror = ErrorInterpolants(metric, q, values, freq, ib, idxs, value0, value1, closest0, closest1);
  Scr3 berror = Min(nerror, m_berror);
  
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
      
void HDRIndexFit::BetterInterpolantsS(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
					 Vec3 const* values, Scr3 const* freq,
				         int ib, int idxs, Vec3 &value0, int closest0) {
  Scr3 nerror = ErrorInterpolantsS(metric, q, values, freq, ib, idxs, value0, closest0);
  Scr3 berror = Min(nerror, m_berror);
  
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

void HDRIndexFit::BetterInterpolantsE(int set, Vec3 const &metric, fQuantizer &q, u8 (&closest)[16],
					 Vec3 const* values, Scr3 const* freq,
				         int ib, int idxs, Vec3 &value1, int closest1) {
  Scr3 nerror = ErrorInterpolantsE(metric, q, values, freq, ib, idxs, value1, closest1);
  Scr3 berror = Min(nerror, m_berror);
  
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
    set, metric, q, closest, values, freq, ib, idxs,			\
    _value0, _value1, _closest0, _closest1);				\
  if (!FEATURE_INDEXFIT_THOROUGH) break;

#define ChooseInterpolantsS(_closest0, _value0)				\
  BetterInterpolantsS(							\
    set, metric, q, closest, values, freq, ib, idxs,			\
    _value0, _closest0);						\
  if (!FEATURE_INDEXFIT_THOROUGH) break;

#define ChooseInterpolantsE(_closest1, _value1)				\
  BetterInterpolantsE(							\
    set, metric, q, closest, values, freq, ib, idxs,			\
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

Scr3 HDRIndexFit::StretchEndPoints(int set, Vec3 const &metric, fQuantizer &q, int ib, u8 (&closest)[16])
{
  // cache some values
  Vec3 const* values = m_palette->GetPoints(set);
  Scr3 const* freq = m_palette->GetWeights(set);
  
  // hold these in SSE-registers
  int gap, idxs = (1 << ib) - 1;
  Col3 igap, iidxs; iidxs.SetRGBpow2<0>(ib); iidxs = iidxs - Col3(1);

#ifdef FEATURE_INDEXFIT_INLINED
  // get the base error
  ErrorEndPoints(set, metric, q, closest, values, freq, ib, idxs);
#endif

  // the idea is that quantization noise has less impact on the inner most
  // interpolated values, that is, the angle-deviation becomes smaller the
  // longer the line is
  Vec3 value0 = values[0];
  Vec3 value1 = values[1];
  Vec3 delta  = value1 - value0;

  int closest0 = 0;
  int closest1 = idxs;
      
  // try to find plain expansions into both directions first
  for (gap = 1, igap = Col3(1); gap < idxs; gap += 2, igap += Col3(2)) {
    // 3-1 >> 1 -> 1, 1/1
    // 7-1 >> 1 -> 3, 3/1
    // 7-3 >> 1 -> 2, 2/3
    // 7-5 >> 1 -> 1, 1/5
    // ...
    // value0 - delta * ((float)exp / gap);
    // value1 + delta * ((float)exp / gap);
    
    Col3 iexp = (iidxs - igap) >> 1;
    Vec3 fexp = Vec3(iexp);
    Vec3 fgap = Reciprocal(Vec3(igap));
    Vec3 stretch = delta * fexp * fgap;

    {
      Vec3 expA0 = value0 - stretch;
      Vec3 expA1 = value1 + stretch;

      Vec3 minA = Min(expA0, expA1);
      Vec3 maxA = Max(expA0, expA1);

      // expansion into both directions works
      if (!(CompareAnyGreaterThan(maxA, Vec3(1.0f + (0.25f / 255.0f))) ||
	    CompareAnyLessThan   (minA, Vec3(     -  0.25f / 255.0f)))) {
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
	  
    Col3 isht = (iidxs) >> 1;
    Col3 ilng = (iidxs - isht);
    
    do {
      Vec3 elong = Vec3(ilng);
      Vec3 eshort = Reciprocal(Vec3(isht));

      // ---------------------------------------------------------------------------------------
      Vec3 factorl = elong * eshort;
      Vec3 stretchl = delta * factorl;
      
      // expA0 = value0 - delta * ((float)elong / eshort);
      // expA1 = value1 + delta * ((float)elong / eshort);
      Vec3 expA0 = value0 - stretchl;
      Vec3 expA1 = value1 + stretchl;

      // expansion down works
      if (!CompareAnyGreaterThan(Abs(Vec3(0.5f) - expA0), Vec3(0.5f + (0.25f / 255.0f)))) {
	int eshort = (idxs + 0) >> 1; ChooseInterpolantsS(closest1 - eshort, expA0); }
      // expansion up works
      if (!CompareAnyGreaterThan(Abs(Vec3(0.5f) - expA1), Vec3(0.5f + (0.25f / 255.0f)))) {
	int eshort = (idxs + 0) >> 1; ChooseInterpolantsE(closest0 + eshort, expA1); }
      
      // ---------------------------------------------------------------------------------------
      Vec3 factors = Reciprocal(factorl);
      Vec3 stretchs = delta * factors;

      // expA0 = value0 - delta * ((float)eshort / elong);
      // expA1 = value1 + delta * ((float)eshort / elong);
      expA0 = value0 - stretchs;
      expA1 = value1 + stretchs;

      // expansion down works
      if (!CompareAnyGreaterThan(Abs(Vec3(0.5f) - expA0), Vec3(0.5f + (0.25f / 255.0f)))) {
	int elong  = (idxs + 1) >> 1; ChooseInterpolantsS(closest1 - elong, expA0); }
      // expansion up works
      if (!CompareAnyGreaterThan(Abs(Vec3(0.5f) - expA1), Vec3(0.5f + (0.25f / 255.0f)))) {
	int elong  = (idxs + 1) >> 1; ChooseInterpolantsE(closest0 + elong, expA1); }
    } while(0);
  }
  
#ifndef FEATURE_INDEXFIT_INLINED
  // snap floating-point-values to the integer-lattice
  Col3 vstart = q.QuantizeToLattice(value0);
  Col3 vend   = q.QuantizeToLattice(value1   , vstart);
  Col3 estart = q.QuantizeToLattice(values[0]);
  Col3 eend   = q.QuantizeToLattice(values[1], estart);

  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Col3 vcode0 = weights_C4[ib][idxs - closest0] * vstart + weights_C4[ib][closest0] * vend;
  Col3 vcode1 = weights_C4[ib][idxs - closest1] * vstart + weights_C4[ib][closest1] * vend;

  Vec3 vsval = q.UnquantizeFromLattice(vcode0);
  Vec3 veval = q.UnquantizeFromLattice(vcode1);
  Vec3 esval = q.UnquantizeFromLattice(estart);
  Vec3 eeval = q.UnquantizeFromLattice(eend);

  // error-sum of the two values
  Vec3 verror0 = LengthSquared(metric * (values[0] - vsval)) * freq[0];
  Vec3 verror1 = LengthSquared(metric * (values[1] - veval)) * freq[1];
  Vec3 eerror0 = LengthSquared(metric * (values[0] - esval)) * freq[0];
  Vec3 eerror1 = LengthSquared(metric * (values[1] - eeval)) * freq[1];
      
  // final check if it really improved
  Vec3 verror = (verror0 + verror1);
  Vec3 eerror = (eerror0 + eerror1);
  Vec3 merror = Min(verror, eerror);
      
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
