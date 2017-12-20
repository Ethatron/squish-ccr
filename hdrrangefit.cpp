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

#include "hdrrangefit.h"
#include "hdrset.h"
#include "hdrblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
HDRRangeFit::HDRRangeFit(HDRSet const* palette, int flags)
  : HDRSingleMatch(palette, flags)
  ,    HDRIndexFit(palette, flags)
  ,         HDRFit(palette, flags)
{
  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();

  assume((isets > 0) && (isets <= 2));

  for (int s = 0; s < isets; s++) {
    // cache some values
    bool const unweighted = m_palette->IsUnweighted(s);
    int const count = m_palette->GetCount(s);
    Vec3 const* values = m_palette->GetPoints(s);
    
    // the codebook endpoints
    Vec3 start(0.0f);
    Vec3 end(0.0f);

    // we don't do this for sparse sets
		if (count > 2) {
      Vec3 centroid;
      Vec3 principle;

      {
        Sym3x3 covariance;

        // get the covariance matrix
        if (unweighted)
	  ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric);
        else
	  ComputeWeightedCovariance3(covariance, centroid, count, values, m_metric, m_palette->GetWeights(s));

	// compute the principle component
	GetPrincipleComponent(covariance, principle);
      }

      // get the min and max range as the codebook endpoints
#ifdef	FEATURE_RANGEFIT_PROJECT
      Scr3 div = Reciprocal(Dot(principle, principle));
      Vec3 rec = Reciprocal(    principle            );
      Scr3 len, min, max;
      Vec3 chk;

      // compute the projection
      min = max = Dot(values[0] - centroid, principle);

      for (int i = 1; i < count; ++i) {
	len = Dot(values[i] - centroid, principle);
	min = Min(min, len);
	max = Max(max, len);
      }

      start = centroid + principle * min * div;
      end   = centroid + principle * max * div;

			// a) floating point values basically don't have any range-limits
			// b) unsigned floats can't be negative
			// c) typeless 16bit do have range-limits

      // take care that low magnitude overshoots can have very high
      // derivatives on the other axis (-0.0039 in R may be +1 in G,
      // thus we at least go to -0.0039/255 to get +1/255 -> 1/255²)

      if (1 /* IsTypeUH*/ ) {
	// intersect with negative axis-plane, clamp to 0.0
	chk = start;
	while (CompareAnyLessThan(chk, Vec3(-1.0f / (255.0f * 255.0f)))) {
	  Vec3 fct = chk * rec;
	  Vec3 hin = Select(fct, chk, HorizontalMin(chk));

	  start -= principle * hin;
	  chk = start;
	}

	// intersect negative undershoot with axis-plane(s), clamp to 0.0
	chk = end;
	while (CompareAnyLessThan(chk, Vec3(-1.0f / (255.0f * 255.0f)))) {
	  Vec3 fct = chk * rec;
	  Vec3 hin = Select(fct, chk, HorizontalMin(chk));

	  end -= principle * hin;
	  chk = end;
	}
      }

      if (0 /* IsTypeUH*/ ) {
	// intersect positive overshoot with axis-plane(s), clamp to 1.0
	chk = start - Vec3(1.0f);
	while (CompareAnyGreaterThan(chk, Vec3(1.0f / (255.0f * 255.0f)))) {
	  Vec3 fct = chk * rec;
	  Vec3 hax = Select(fct, chk, HorizontalMax(chk));

	  start -= principle * hax;
	  chk = start - Vec3(1.0f);
	}

	// intersect positive overshoot with axis-plane(s), clamp to 1.0
	chk = end - Vec3(1.0f);
	while (CompareAnyGreaterThan(chk, Vec3(1.0f / (255.0f * 255.0f)))) {
	  Vec3 fct = chk * rec;
	  Vec3 hax = Select(fct, chk, HorizontalMax(chk));

	  end -= principle * hax;
	  chk = end - Vec3(1.0f);
	}
      }

/*    assert(HorizontalMin(start).X() > -0.0001);
      assert(HorizontalMin(end  ).X() > -0.0001);
      assert(HorizontalMax(start).X() <  1.0001);
      assert(HorizontalMax(end  ).X() <  1.0001);  */
#else
      Scr4 min, max;

      // compute the range
      start = end = values[0];
      min = max = Dot(values[0], principle);

      for (int i = 1; i < count; ++i) {
	Scr4 val = Dot(values[i], principle);

	if (min > val) {
	  start = values[i];
	  min = val;
	}
	else if (max < val) {
	  end = values[i];
	  max = val;
	}
      }
#endif

      m_centroid[s] = centroid;
    }
    else if (count == 2) {
      start = values[1];
      end = values[0];
    }
    else if (count == 1) {
      start = 
      end = values[0];
    }

    // clamp the output to [0, 1]
    m_start[s] = Max(start, Vec3(0.0f));
    m_end  [s] = Max(end  , Vec3(0.0f));
  }

#ifdef	FEATURE_ELIMINATE_FLATBOOKS
  // backup original values, we don't now the precision yet
  for (int x = 0; x < (isets + asets); x++) {
    m_start_candidate[x] = m_start[x];
    m_end_candidate[x] = m_end[x];
  }
#endif
}

void HDRRangeFit::Compress(void* block, fQuantizer &q, int mode)
{
	int swaps = 0;
SwapSet1EndPointsAndRedoCalculation:

  int ib = GetIndexBits(mode);
  int tb = GetTruncationBits(mode);
  int db = GetDeltaBits(mode);

  q.ChangeField(tb, db);

  // match each point to the closest code
  Scr3 error = Scr3(0.0f);
  a16 u8 closest[4][16];

  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();

  assume((isets > 0) && (isets <= 3));

  // create a codebook
  Vec3 codes[1 << 4];

  /* TODO: part of the magic is to find a nicer range-reducer,
   * as we are capped by the mode's delta range, this is very
   * important for the performance, several approaches may be nice:
   * - centroid shrinking
   * - ...
   */

#if 0
  // NOTE: evil numbers:
  //
  // if we need to swap start/end of the first set, then
  // the whole calculation's result goes to hell, as the
  // anchor point for all deltas changes
  //
  // calculate both optimal range-shrinks upfront?
  Vec4 midpt = (m_start[0] + m_end[0]) * 0.5f;
  Col4 midhf = FloatToUHalf<false>(midpt);
  Col4 midqt = q.QuantizeToLattice(midpt);
  Col4 ctrqt = q.QuantizeToLattice(m_centroid[0]);

//int rng = ((1 << (db - 0)) - 0) << tb;
  int lwr = ((1 << (db - 1)) - 0) << tb;
  int upr = ((1 << (db - 1)) - 1) << tb;

  m_qstart[0] = q.QuantizeToLattice(m_start[0]);
  m_qend  [0] = q.QuantizeToLattice(m_end  [0]);

  Col4 range = MaxTiny(m_qstart[0], m_qend[0]) - MinTiny(m_qstart[0], m_qend[0]);
  Col4 maxr = HorizontalMaxTiny(range);

  if (maxr > Col4(upr)) {
    range = Div32x16(Mul16x16(range, Col4(upr)), maxr);
    midqt = (m_qstart[0] + m_qend[0]) >> 1;

    m_qstart[0] = midqt - (range >> 1);
    m_qend  [0] = m_qstart[0] + range;
  }

  m_qstart[0] = midqt - Col4(lwr);
  m_qend  [0] = midqt + Col4(upr);

  m_qstart[0] = ctrqt - Col4(lwr);
  m_qend  [0] = ctrqt + Col4(upr);

  m_qstart[0] = q.QuantizeToLattice(m_start[0]);
  m_qend  [0] = q.QuantizeToLattice(m_end  [0]);
  m_qstart[1] = q.QuantizeToLattice(m_start[1]);
  m_qend  [1] = q.QuantizeToLattice(m_end  [1]);

  Col4 range1 = (m_qstart[0] - m_qend  [0]);
  Col4 range2 = (m_qstart[0] - m_qstart[1]);
  Col4 range3 = (m_qstart[0] - m_qend  [0]);
#endif

  // if the range has been shrunk, we have to redo everything when
  // we do the post-swap of start/end
  if (db) {
#if 1
    // snap floating-point-values to the half's delta-lattice
    m_qstart[0] = q.QuantizeToLattice(m_start[0]);
    m_qend  [0] = q.QuantizeToLattice(m_end  [0]);

    do {
      Col3 delta = m_qend[0] - m_qstart[0];
      Col3 range = Abs(delta);
      Col3 qgrid = q.griddltp << tb;

      // iterative to allow different max-ranges per component
      int match = CompareGreaterThan(range, qgrid);
      if (!(match & 0x0FFF))
	break;

      Col3 qnorm = qgrid;
      Col3 rnorm = range;

      /**/ if (match & 0x000F)
	qnorm = qnorm.SplatR(),
	rnorm = rnorm.SplatR();
      else if (match & 0x00F0)
	qnorm = qnorm.SplatG(),
	rnorm = rnorm.SplatG();
      else // (match & 0x0F00)
	qnorm = qnorm.SplatB(),
	rnorm = rnorm.SplatB();

      // NOTE: rounding doesn't help ... probably because of the log-scale
      delta = Div32x16s(Mul16x16s(delta, qnorm), rnorm);
			delta &= q.gridprc;

#if 0
      Col3 m = (m_qstart[0] + m_qend[0]) >> 1;

      m_qstart[0] = m - (delta >> 1);
      m_qend  [0] = m_qstart[0] + delta;
#else
      Col3 s = HorizontalMaxTiny(m_qstart[0]);
      Col3 e = HorizontalMaxTiny(m_qend  [0]);

      if (s > e)
	m_qend[0] = m_qstart[0] + delta;
      else
	m_qstart[0] = m_qend[0] - delta;
#endif
    } while(1);
#else
    // snap floating-point-values to the half's delta-lattice
    m_qstart[0] = q.QuantizeToLattice(m_start[0]);
    m_qend  [0] = q.QuantizeToLattice(m_end  [0], m_qstart[0]);
#endif

#if 1
    if (isets > 1) {
      // snap floating-point-values to the half's delta-lattice
      m_qstart[1] = q.QuantizeToLattice(m_start[1]);
      m_qend  [1] = q.QuantizeToLattice(m_end  [1]);

      do {
	Col3 deltas = m_qstart[1] - m_qstart[0];
	Col3 deltae = m_qend  [1] - m_qstart[0];
	Col3 ranges = Abs(deltas);
	Col3 rangee = Abs(deltae);
	Col3 range = MaxTiny(ranges, rangee);
	Col3 qgrid = q.griddltp << tb;

	// iterative to allow different max-ranges per component
	int match = CompareGreaterThan(range, qgrid);
	if (!(match & 0x0FFF))
	  break;

	Col3 qnorm = qgrid;
	Col3 rnorm = range;

	/**/ if (match & 0x000F)
	  qnorm = qnorm.SplatR(),
	  rnorm = rnorm.SplatR();
	else if (match & 0x00F0)
	  qnorm = qnorm.SplatG(),
	  rnorm = rnorm.SplatG();
	else // (match & 0x0F00)
	  qnorm = qnorm.SplatB(),
	  rnorm = rnorm.SplatB();

	// NOTE: rounding doesn't help ... probably because of the log-scale
	deltas = Div32x16s(Mul16x16s(deltas, qnorm), rnorm);
	deltae = Div32x16s(Mul16x16s(deltae, qnorm), rnorm);
	deltas &= q.gridprc;
	deltae &= q.gridprc;

	m_qstart[1] = m_qstart[0] + deltas;
	m_qend  [1] = m_qstart[0] + deltae;
      } while(1);
    }
#else
    // loop over all sets and prepare initial delta-capped end-points
    for (int s = 1; s < isets; s++) {
      // snap floating-point-values to the half's delta-lattice
      m_qstart[s] = q.QuantizeToLattice(m_start[s], m_qstart[0]);
      m_qend  [s] = q.QuantizeToLattice(m_end  [s], m_qstart[0]);
    }
#endif
  }
  else {
    for (int s = 0; s < isets; s++) {
      // snap floating-point-values to the half's explicit-lattice
      m_qstart[s] = q.QuantizeToLattice(m_start[s]);
      m_qend  [s] = q.QuantizeToLattice(m_end  [s]);
    }
  }

  // loop over all sets
  for (int s = 0; s < isets; s++) {
    // cache some values
    int const count = m_palette->GetCount(s);
    Vec3 const* values = m_palette->GetPoints(s);
    Scr3 const* freq = m_palette->GetWeights(s);

    // in case of separate alpha the colors of the alpha-set have all been set to alpha
    Vec3 metric = m_metric;

#if 0
    // we do single entry fit for sparse sets
    if (count == 1) {
      // find the closest code
      Scr3 dist = ComputeEndPoints(s, metric, tb, db, ib, 0x7);

      // save the index (it's just a single one)
      closest[s][0] = GetIndex();
      
      // accumulate the error
      error += dist * freq[0];
    }
    // we do dual entry fit for sparse sets
    else if (count == 2) {
      // find the closest codes (it's just two)
      Scr3 dist = StretchEndPoints(s, metric, q, ib, closest[s]);
      
      // accumulate the error
      error += dist;
    }
    else
#endif
    {
      // create codebook as half (interpolation is in integer log-space)
      Col3 codeh[1 << 4];

      // then translate to float (for linear error measurement)
      int ccs = CodebookP(codeh, ib, m_qstart[s], m_qend[s]);
      for (int i = 0; i < ccs; ++i)
	codes[i] = metric * q.UnquantizeFromLattice(codeh[i]);

      for (int i = 0; i < count; ++i) {
	// find the closest code
	Scr3 dist = Scr3(FLT_MAX);
	Vec3 value = metric * values[i];
	int idx = 0;

	// this loop will always choose the lowest index in
	// case two of them are equal, this guarantees that
	// swapping of end-points quits leading 1s
	assume((ccs >= (1 << 2)) && (ccs <= (1 << 4)));
	int j = ccs - 1; do {
	  Scr3 d0 = LengthSquared(value - codes[j - 0]);
	  Scr3 d1 = LengthSquared(value - codes[j - 1]);
	  Scr3 d2 = LengthSquared(value - codes[j - 2]);
	  Scr3 d3 = LengthSquared(value - codes[j - 3]);

	  // encourage OoO
	  Scr3 da = Min(d0, d1);
	  Scr3 db = Min(d2, d3);
	  dist = Min(da, dist);
	  dist = Min(db, dist);

	  // will cause VS to make them all cmovs
	  if (d0 == dist) { idx = j; } j--;
	  if (d1 == dist) { idx = j; } j--;
	  if (d2 == dist) { idx = j; } j--;
	  if (d3 == dist) { idx = j; } j--;
	} while(j >= 0);

	// save the index
	closest[s][i] = (u8)idx;

	// accumulate the error
	error += dist * freq[i];
      }
    }

#if 0
    // kill early if this scheme looses
    if (!(error < m_besterror))
			return;
#endif
  }

  // remap the indices
  for (int s = 0; s < isets; s++)
    m_palette->RemapIndices(closest[s], m_indices, s);

#if 1
  // multiple sets, first set is always index 0
  if ((isets > 1) && (m_indices[0] & 0x4)) {
    m_start[0].SwapXYZ(m_end[0]);

		swaps++;
		if (swaps > 1)
			return;

		goto SwapSet1EndPointsAndRedoCalculation;
  }
#endif

	// kill late if this scheme looses (wait for the swap-correction)
	if (!(error < m_besterror))
		return;

  typedef u8 (&Ione)[1][16];

  typedef Col4 (&C4two)[2];
  typedef Col4 (&C4one)[1];

  // save the block
  int partition = m_palette->GetPartition();
  switch (mode) {
    case  0: WriteHDRBlock_m1(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  1: WriteHDRBlock_m2(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  2: WriteHDRBlock_m3(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  3: WriteHDRBlock_m4(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  4: WriteHDRBlock_m5(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  5: WriteHDRBlock_m6(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  6: WriteHDRBlock_m7(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  7: WriteHDRBlock_m8(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  8: WriteHDRBlock_m9(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case  9: WriteHDRBlock_mA(partition, (C4two)m_qstart, (C4two)m_qend, (Ione)m_indices, block); break;
    case 10: WriteHDRBlock_mB(           (C4one)m_qstart, (C4one)m_qend, (Ione)m_indices, block); break;
    case 11: WriteHDRBlock_mC(           (C4one)m_qstart, (C4one)m_qend, (Ione)m_indices, block); break;
    case 12: WriteHDRBlock_mD(           (C4one)m_qstart, (C4one)m_qend, (Ione)m_indices, block); break;
    case 13: WriteHDRBlock_mE(           (C4one)m_qstart, (C4one)m_qend, (Ione)m_indices, block); break;
  }

  // save the error
  m_besterror = error;
  m_best = true;
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
