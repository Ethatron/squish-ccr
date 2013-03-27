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

#include "colournormalfit.h"
#include "colourset.h"
#include "colourblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
ColourNormalFit::ColourNormalFit(ColourSet const* colours, int flags)
  : ColourFit(colours, flags)
{
  cQuantizer3<5,6,5> q = cQuantizer3<5,6,5>();

  // initialize the best error
  m_besterror = Scr3(FLT_MAX);

  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  float const* weights = m_colours->GetWeights();

#if 0
  u8 const* freq = m_colours->GetFrequencies();
  
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  // build the average normal vector
  Vec3 average = Vec3(0.0f);
  Vec3 norml[16];
  for (int i = 0; i < count; ++i) {
    norml[i] = Normalize(scale * (offset + values[i]));
    average += norml[i] * freq[i];
  }

  average = Normalize(average);
  average = average + Vec3(0.0f, 0.0f, 1.0f);
  average = Normalize(average);

  Vec3 rotX, rotY, rotZ;
  Vec3 invX, invY, invZ;
  Vec3 test, tEst;

  rotZ = average;
  rotY = average; rotY.GetZ() = 0.0f; rotY = Normalize(rotY);

  rotX.GetX() = rotZ.Y() * rotY.Z() - rotZ.Z() * rotY.Y();
  rotX.GetY() = rotZ.Z() * rotY.X() - rotZ.X() * rotY.Z();
  rotX.GetZ() = rotZ.X() * rotY.Y() - rotZ.Y() * rotY.X();
  
  rotY.GetX() = rotZ.Y() * rotX.Z() - rotZ.Z() * rotX.Y();
  rotY.GetY() = rotZ.Z() * rotX.X() - rotZ.X() * rotX.Z();
  rotY.GetZ() = rotZ.X() * rotX.Y() - rotZ.Y() * rotX.X();
  
  invX.GetX() = rotX.X(); invY.GetX() = rotX.Y(); invZ.GetX() = rotX.Z();
  invX.GetY() = rotY.X(); invY.GetY() = rotY.Y(); invZ.GetY() = rotY.Z();
  invX.GetZ() = rotZ.X(); invY.GetZ() = rotZ.Y(); invZ.GetZ() = rotZ.Z();
  
  test.GetX() = Dot(average, rotX).X();
  test.GetY() = Dot(average, rotY).Y();
  test.GetZ() = Dot(average, rotZ).Z();
  
  tEst.GetX() = Dot(test, invX).X();
  tEst.GetY() = Dot(test, invY).Y();
  tEst.GetZ() = Dot(test, invZ).Z();
  
  // transform vectors into tangent space
  Vec3 tspace[16];
  for (int i = 0; i < count; ++i) {
    tspace[i] = Vec3(0.0f);
    tspace[i].GetX() = Dot(norml[i], rotX).X() * 0.5f + 0.5f;
    tspace[i].GetY() = Dot(norml[i], rotY).Y() * 0.5f + 0.5f;
    tspace[i].GetZ() = Dot(norml[i], rotZ).Z() * 0.5f + 0.5f;
  }
  
  Sym3x3 covariance;
  Vec3 centroid;
  Vec3 principle;

  // get the covariance matrix
  if (m_colours->IsUnweighted())
    ComputeWeightedCovariance3(covariance, centroid, count, tspace, Vec3(1.0f));
  else
    ComputeWeightedCovariance3(covariance, centroid, count, tspace, Vec3(1.0f), weights);

  // compute the principle component
  GetPrincipleComponent(covariance, principle);

  // get the min and max normal as the codebook endpoints
  Vec3 start(0.0f);
  Vec3 end(0.0f);
  
//centroid = Vec3(2.0f) * (Vec3(-0.5f) + centroid);
//centroid = centroid * Vec3(1.0f, 1.0f, centroid.Z() * centroid.Z());
//centroid = (centroid * Vec3(0.5f)) + Vec3(0.5f);

//centroid = Vec3(0.5f, 0.5f, 0.00f);
//centroid = Vec3(0.5f, 0.5f, ((centroid.Z() - 0.5f) * 0.5f) + 0.5f);
//centroid = Vec3(0.5f, 0.5f, ((centroid.Z() * 2.0f) - 1.0f) + 0.5f);

  if (count > 0) {
#ifdef	FEATURE_NORMALFIT_PROJECT
    Scr3 div = Reciprocal(Dot(principle, principle));
    Vec3 rec = Reciprocal(    principle            );
    Scr3 len, min, max;
    Vec3 chk;

    // compute the projection
    min = max = Dot(tspace[0] - centroid, principle);

    for (int i = 1; i < count; ++i) {
      len = Dot(tspace[i] - centroid, principle);
      min = Min(min, len);
      max = Max(max, len);
    }

    start = centroid + principle * min * div;
    end   = centroid + principle * max * div;
    
//  start = (Complement(Vec3(2.0f) * (Vec3(-0.5f) + start)) * Vec3(0.5f)) + Vec3(0.5f);
//  end   = (Complement(Vec3(2.0f) * (Vec3(-0.5f) + end  )) * Vec3(0.5f)) + Vec3(0.5f);

    // intersect negative undershoot with axis-plane(s), clamp to 0.0
    chk = start;
    while (CompareAnyLessThan(chk, Vec3(-1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hin = Select(fct, chk, HorizontalMin(chk));

      start -= principle * hin;
      chk = start;
    }

    // intersect negative undershoot with axis-plane(s), clamp to 0.0
    chk = end;
    while (CompareAnyLessThan(chk, Vec3(-1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hin = Select(fct, chk, HorizontalMin(chk));

      end -= principle * hin;
      chk = end;
    }

    // intersect positive overshoot with axis-plane(s), clamp to 1.0
    chk = start - Vec3(1.0f);
    while (CompareAnyGreaterThan(chk, Vec3(1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hax = Select(fct, chk, HorizontalMax(chk));

      start -= principle * hax;
      chk = start - Vec3(1.0f);
    }

    // intersect positive overshoot with axis-plane(s), clamp to 1.0
    chk = end - Vec3(1.0f);
    while (CompareAnyGreaterThan(chk, Vec3(1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hax = Select(fct, chk, HorizontalMax(chk));

      end -= principle * hax;
      chk = end - Vec3(1.0f);
    }

/*  assert(HorizontalMin(start).X() > -0.0001);
    assert(HorizontalMin(end  ).X() > -0.0001);
    assert(HorizontalMax(start).X() <  1.0001);
    assert(HorizontalMax(end  ).X() <  1.0001);	 */
#else
    Scr3 min, max;

    // compute the normal
    start = end = tspace[0];
    min = max = Dot(tspace[0], principle);

    for (int i = 1; i < count; ++i) {
      Scr3 val = Dot(tspace[i], principle);

      if (val < min) {
	start = tspace[i];
	min = val;
      }
      else if (val > max) {
	end = tspace[i];
	max = val;
      }
    }
#endif
  }

  for (int i = 0; i < count; ++i) {
  for (int j = 0; j < count; ++j) {
  }
  }

  Vec3 istart;
  Vec3 iend;
  
  start = (start * 2.0f) - Vec3(1.0f);
  end   = (end   * 2.0f) - Vec3(1.0f);

  istart.GetX() = Dot(start, invX).X();
  istart.GetY() = Dot(start, invY).Y();
  istart.GetZ() = Dot(start, invZ).Z();
  
  iend.GetX() = Dot(end, invX).X();
  iend.GetY() = Dot(end, invY).Y();
  iend.GetZ() = Dot(end, invZ).Z();

//istart = Complement(istart);
//iend   = Complement(iend  );
//
//istart = Normalize(istart);
//iend   = Normalize(iend  );

  istart = (istart * 0.5f) + Vec3(0.5f);
  iend   = (iend   * 0.5f) + Vec3(0.5f);

//istart = Max(Vec3(0.0f), Min(Vec3(1.0f), istart));
//iend   = Max(Vec3(0.0f), Min(Vec3(1.0f), iend  ));

  // snap floating-point-values to the integer-lattice and save
  m_start = q.SnapToLattice(istart);
  m_end   = q.SnapToLattice(iend);
#else
  Sym3x3 covariance;
  Vec3 centroid;
  Vec3 principle;

  // get the covariance matrix
  if (m_colours->IsUnweighted())
    ComputeWeightedCovariance3(covariance, centroid, count, values, Vec3(1.0f));
  else
    ComputeWeightedCovariance3(covariance, centroid, count, values, Vec3(1.0f), weights);

  // compute the principle component
  GetPrincipleComponent(covariance, principle);

  // get the min and max normal as the codebook endpoints
  Vec3 start(0.0f);
  Vec3 end(0.0f);

  if (count > 0) {
#ifdef	FEATURE_NORMALFIT_PROJECT
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

    // intersect negative undershoot with axis-plane(s), clamp to 0.0
    chk = start;
    while (CompareAnyLessThan(chk, Vec3(-1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hin = Select(fct, chk, HorizontalMin(chk));

      start -= principle * hin;
      chk = start;
    }

    // intersect negative undershoot with axis-plane(s), clamp to 0.0
    chk = end;
    while (CompareAnyLessThan(chk, Vec3(-1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hin = Select(fct, chk, HorizontalMin(chk));

      end -= principle * hin;
      chk = end;
    }

    // intersect positive overshoot with axis-plane(s), clamp to 1.0
    chk = start - Vec3(1.0f);
    while (CompareAnyGreaterThan(chk, Vec3(1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hax = Select(fct, chk, HorizontalMax(chk));

      start -= principle * hax;
      chk = start - Vec3(1.0f);
    }

    // intersect positive overshoot with axis-plane(s), clamp to 1.0
    chk = end - Vec3(1.0f);
    while (CompareAnyGreaterThan(chk, Vec3(1.0f / 65536))) {
      Vec3 fct = chk * rec;
      Vec3 hax = Select(fct, chk, HorizontalMax(chk));

      end -= principle * hax;
      chk = end - Vec3(1.0f);
    }

/*  assert(HorizontalMin(start).X() > -0.0001);
    assert(HorizontalMin(end  ).X() > -0.0001);
    assert(HorizontalMax(start).X() <  1.0001);
    assert(HorizontalMax(end  ).X() <  1.0001);	 */
#else
    Scr3 min, max;

    // compute the normal
    start = end = values[0];
    min = max = Dot(values[0], principle);

    for (int i = 1; i < count; ++i) {
      Scr3 val = Dot(values[i], principle);

      if (val < min) {
	start = values[i];
	min = val;
      }
      else if (val > max) {
	end = values[i];
	max = val;
      }
    }
#endif
  }

  // snap floating-point-values to the integer-lattice and save
  m_start = q.SnapToLattice(start);
  m_end   = q.SnapToLattice(end  );
#endif
}

void ColourNormalFit::Compress3(void* block)
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  u8 const* freq = m_colours->GetFrequencies();

  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[3]; Codebook3(codes, m_start, m_end);

  codes[0] = Normalize(scale * (offset + codes[0]));
  codes[1] = Normalize(scale * (offset + codes[1]));
  codes[2] = Normalize(scale * (offset + codes[2]));

  // match each point to the closest code
  u8 closest[16];

  Scr3 error = Scr3(16.0f);
  for (int i = 0; i < count; ++i) {
    // find the closest code
    Scr3 dist;
    Vec3 value = Normalize(scale * (offset + values[i]));
    int idx = 0;

    {
      Scr3 d0 = Dot(value, codes[0]);
      Scr3 d1 = Dot(value, codes[1]);
      Scr3 d2 = Dot(value, codes[2]);

      // encourage OoO
      Scr3 da = Max(d0, d1);
      Scr3 db =    (d2    );
      dist    = Max(da, db);

      // will cause VS to make them all cmovs
      if (d2 == dist) { idx = 2; }
      if (d1 == dist) { idx = 1; }
      if (d0 == dist) { idx = 0; }
    }

    // save the index
    closest[i] = (u8)idx;

    // accumulate the error
    error -= dist * freq[i];
  }

  // save this scheme if it wins
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16];

    m_colours->RemapIndices(closest, indices);

    // save the block
    WriteColourBlock3(m_start, m_end, indices, block);
  }
}

void ColourNormalFit::Compress4(void* block)
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  u8 const* freq = m_colours->GetFrequencies();

  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[4]; Codebook4(codes, m_start, m_end);
  
  codes[0] = Normalize(scale * (offset + codes[0]));
  codes[1] = Normalize(scale * (offset + codes[1]));
  codes[2] = Normalize(scale * (offset + codes[2]));
  codes[3] = Normalize(scale * (offset + codes[3]));

  // match each point to the closest code
  u8 closest[16];

  Scr3 error = Scr3(16.0f);
  for (int i = 0; i < count; ++i) {
    // find the closest code
    Scr3 dist;
    Vec3 value = Normalize(scale * (offset + values[i]));
    int idx = 0;

    {
      Scr3 d0 = Dot(value, codes[0]);
      Scr3 d1 = Dot(value, codes[1]);
      Scr3 d2 = Dot(value, codes[2]);
      Scr3 d3 = Dot(value, codes[3]);

      // encourage OoO
      Scr3 da = Max(d0, d1);
      Scr3 db = Max(d2, d3);
      dist    = Max(da, db);

      // will cause VS to make them all cmovs
      if (d3 == dist) { idx = 3; }
      if (d2 == dist) { idx = 2; }
      if (d1 == dist) { idx = 1; }
      if (d0 == dist) { idx = 0; }
    }

    // save the index
    closest[i] = (u8)idx;

    // accumulate the error
    error -= dist * freq[i];
  }

  // save this scheme if it wins
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16];

    m_colours->RemapIndices(closest, indices);

    // save the block
    WriteColourBlock4(m_start, m_end, indices, block);
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish
