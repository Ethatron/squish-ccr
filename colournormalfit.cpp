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

  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* weights = m_colours->GetWeights();

#ifdef	FEATURE_NORMALFIT_PROJECT
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
  Vec3 start(127.5f, 127.5f, 255.0f);
  Vec3 end(127.5f, 127.5f, 255.0f);

  if (count > 0) {
#undef	FEATURE_NORMALFIT_PROJECT_NEAREST
#ifdef	FEATURE_NORMALFIT_PROJECT_NEAREST
    const Vec3 scale  = Vec3( 1.0f / 0.5f);
    const Vec3 offset = Vec3(-1.0f * 0.5f);
    const Vec3 scalei = Vec3( 1.0f * 0.5f);

    Vec3 centroidn = (scale * (offset + centroid));
    Vec3 rec = Reciprocal(principle);
    Scr3 min, max;
    Vec3 chk;

    // http://geomalgorithms.com/a07-_distance.html
    // compute the line parameters of the two closest points
    min = Scr3( FLT_MAX);
    max = Scr3(-FLT_MAX);

    for (int i = 0; i < count; ++i) {
      Vec3 valuenorm = Normalize(scale * (offset + values[i]));

      Vec3 u = principle;//L1.P1 - L1.P0;
      Vec3 v = valuenorm;//L2.P1 - L2.P0;
      Vec3 w = centroidn;//L1.P0 - L2.P0;
      Scr3 a = Dot(u, u);         // always >= 0
      Scr3 b = Dot(u, v);
      Scr3 c = Dot(v, v);         // always >= 0
      Scr3 d = Dot(u, w);
      Scr3 e = Dot(v, w);
      Scr3 D = a * c - b * b;     // always >= 0
      Scr3 sc, tc;

      // compute the line parameters of the two closest points
      if (D < Scr3(0.00001f)) {    // the lines are almost parallel
	sc = Scr3(0.0f);	   // use the largest denominator
	tc = (b > c 
	  ? d * Reciprocal(b) 
	  : e * Reciprocal(c)
	);    
      }
      else {
	D = Reciprocal(D);

	sc = (b * e - c * d) * D;
	tc = (a * e - b * d) * D;
      }

      // one dimension of the principle axis is 1
      // the maximum magnitude the principle axis
      // can move in the [-1,+1] cube is 1.41*2
      // without leaving the cube's boundaries
      sc = Min(sc, Scr3( 2.82842712474619f));
      sc = Max(sc, Scr3(-2.82842712474619f));

      min = Min(min, sc);
      max = Max(max, sc);
    }
    
    start = centroidn + principle * min;
    end   = centroidn + principle * max;

    start = Normalize(start);
    end   = Normalize(end  );
    
    start = (start * scalei) - offset;
    end   = (end   * scalei) - offset;
#else
    // compute the projection
    GetPrincipleProjection(start, end, principle, centroid, count, values);
#endif
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
  m_start_candidate = q.SnapToLattice(start);
  m_end_candidate   = q.SnapToLattice(end  );
}

void ColourNormalFit::kMeans3()
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);
  const Vec3 scalei = Vec3( 1.0f * 0.5f);
  
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();
  
  cQuantizer3<5,6,5> q = cQuantizer3<5,6,5>();
  Vec3 c_start = m_start, c_end = m_end;
  Vec3 l_start = m_start, l_end = m_end;
  Scr3 berror = Scr3(DEVIANCE_MAXSUM);
  
  int trie = 1 + (m_flags & kColourIterativeClusterFits) / kColourClusterFit;
  do {
    Vec3 means[3];

    means[0] = Vec3(0.0f);
    means[1] = Vec3(0.0f);
    means[2] = Vec3(0.0f);
    
    // create a codebook
    // resolve "metric * (value - code)" to "metric * value - metric * code"
    Vec3 codes[3]; Codebook3n(codes, c_start, c_end);

    Scr3 merror = Scr3(DEVIANCE_BASE);
    for (int i = 0; i < count; ++i) {
      int idx = 0;

      // find the closest code
      Vec3 value = Normalize(scale * (offset + values[i]));
      Scr3 dist; MinDeviance3<true>(dist, idx, value, codes);

      // accumulate the error
      AddDeviance(dist, merror, freq[i]);
      
      // accumulate the mean
      means[idx] += value * freq[i];
    }
    
    if (berror > merror) {
      berror = merror;

      m_start = c_start;
      m_end   = c_end;
    }

    means[0] = (Normalize(means[0]) * scalei) - offset;
    means[1] = (Normalize(means[1]) * scalei) - offset;
    
    l_start = c_start;
    l_end   = c_end;
    c_start = q.SnapToLattice(means[0]);
    c_end   = q.SnapToLattice(means[1]);

  } while(--trie && !(CompareAllEqualTo(c_start, l_start) && CompareAllEqualTo(c_end, l_end)));
}

void ColourNormalFit::kMeans4()
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);
  const Vec3 scalei = Vec3( 1.0f * 0.5f);
  
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();
  
  cQuantizer3<5,6,5> q = cQuantizer3<5,6,5>();
  Vec3 c_start = m_start, c_end = m_end;
  Vec3 l_start = m_start, l_end = m_end;
  Scr3 berror = Scr3(DEVIANCE_MAXSUM);
  
  int trie = 1 + (m_flags & kColourIterativeClusterFits) / kColourClusterFit;
  do {
    Vec3 means[4];

    means[0] = Vec3(0.0f);
    means[1] = Vec3(0.0f);
    means[2] = Vec3(0.0f);
    means[3] = Vec3(0.0f);
    
    // create a codebook
    Vec3 codes[4]; Codebook4n(codes, c_start, c_end);

    Scr3 merror = Scr3(DEVIANCE_BASE);
    for (int i = 0; i < count; ++i) {
      int idx = 0;

      // find the closest code
      Vec3 value = Normalize(scale * (offset + values[i]));
      Scr3 dist; MinDeviance4<true>(dist, idx, value, codes);
      
      // accumulate the error
      AddDeviance(dist, merror, freq[i]);
      
      // accumulate the mean
      means[idx] += value * freq[i];
    }
  
    if (berror > merror) {
      berror = merror;
      
      m_start = c_start;
      m_end   = c_end;
    }

    means[0] = (Normalize(means[0]) * scalei) - offset;
    means[1] = (Normalize(means[1]) * scalei) - offset;
    
    l_start = c_start;
    l_end   = c_end;
    c_start = q.SnapToLattice(means[0]);
    c_end   = q.SnapToLattice(means[1]);

  } while(--trie && !(CompareAllEqualTo(c_start, l_start) && CompareAllEqualTo(c_end, l_end)));
}

void ColourNormalFit::Permute3()
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);
  const Vec3 scalei = Vec3( 1.0f * 0.5f);
  
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();
  
  cQuantizer3<5,6,5> q = cQuantizer3<5,6,5>();
  Scr3 berror = Scr3(DEVIANCE_MAXSUM);
  
  Vec3 c_start = m_start;
  Vec3 c_end   = m_end;
  Scr3 l_start = LengthSquared(Normalize(scale * (offset + m_start)));
  Scr3 l_end   = LengthSquared(Normalize(scale * (offset + m_end)));
  Vec3 q_start = Reciprocal(q.grid + Vec3(1.0f));
  Vec3 q_end   = q_start;

  // adjust offset towards sphere-boundary
  if (!(l_start < Scr3(1.0f)))
    q_start = Vec3(0.0f) - q_start;
  if (!(l_end   < Scr3(1.0f)))
    q_end   = Vec3(0.0f) - q_end;
  
  int trie = 0x3F;
  do {
    // permute end-points +-1 towards sphere-boundary
    Vec3 p_start = q_start & Vec3(!(trie & 0x01), !(trie & 0x02), !(trie & 0x04));
    Vec3 p_end   = q_end   & Vec3(!(trie & 0x08), !(trie & 0x10), !(trie & 0x20));
    
    p_start = q.SnapToLattice(c_start + p_start);
    p_end   = q.SnapToLattice(c_end   + p_end);

    // create a codebook
    // resolve "metric * (value - code)" to "metric * value - metric * code"
    Vec3 codes[3]; Codebook3n(codes, p_start, p_end);

    Scr3 merror = Scr3(DEVIANCE_BASE);
    for (int i = 0; i < count; ++i) {
      // find the closest code
      Vec3 value = Normalize(scale * (offset + values[i]));
      Scr3 dist; MinDeviance3<false>(dist, i, value, codes);
      
      // accumulate the error
      AddDeviance(dist, merror, freq[i]);
    }
    
    if (berror > merror) {
      berror = merror;

      m_start = p_start;
      m_end   = p_end;
    }

  } while(--trie);
}

void ColourNormalFit::Permute4()
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);
  const Vec3 scalei = Vec3( 1.0f * 0.5f);
  
  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();
  
  cQuantizer3<5,6,5> q = cQuantizer3<5,6,5>();
  Scr3 berror = Scr3(DEVIANCE_MAXSUM);
  
  Vec3 c_start = m_start;
  Vec3 c_end   = m_end;
  Scr3 l_start = LengthSquared(Normalize(scale * (offset + m_start)));
  Scr3 l_end   = LengthSquared(Normalize(scale * (offset + m_end)));
  Vec3 q_start = Reciprocal(q.grid + Vec3(1.0f));
  Vec3 q_end   = q_start;

  // adjust offset towards sphere-boundary
  if (!(l_start < Scr3(1.0f)))
    q_start = Vec3(0.0f) - q_start;
  if (!(l_end   < Scr3(1.0f)))
    q_end   = Vec3(0.0f) - q_end;
  
  int trie = 0x3F;
  do {
    // permute end-points +-1 towards sphere-boundary
    Vec3 p_start = q_start & Vec3(!(trie & 0x01), !(trie & 0x02), !(trie & 0x04));
    Vec3 p_end   = q_end   & Vec3(!(trie & 0x08), !(trie & 0x10), !(trie & 0x20));
    
    p_start = q.SnapToLattice(c_start + p_start);
    p_end   = q.SnapToLattice(c_end   + p_end);

    // create a codebook
    Vec3 codes[4]; Codebook4n(codes, p_start, p_end);

    Scr3 merror = Scr3(DEVIANCE_BASE);
    for (int i = 0; i < count; ++i) {
      // find the closest code
      Vec3 value = Normalize(scale * (offset + values[i]));
      Scr3 dist; MinDeviance4<false>(dist, i, value, codes);
      
      // accumulate the error
      AddDeviance(dist, merror, freq[i]);
    }
    
    if (berror > merror) {
      berror = merror;

      m_start = p_start;
      m_end   = p_end;
    }

  } while(--trie);
}

void ColourNormalFit::Compress3(void* block)
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  // cache some values
  int const count = m_colours->GetCount();
  Vec3 const* values = m_colours->GetPoints();
  Scr3 const* freq = m_colours->GetWeights();
  
  // use a fitting algorithm
  m_start = m_start_candidate;
  m_end   = m_end_candidate;
  if ((m_flags & kColourIterativeClusterFits))
    kMeans3();
//if ((m_flags & kColourIterativeClusterFits) >= (kColourIterativeClusterFit))
//  Permute3();
  
  // create a codebook
  // resolve "metric * (value - code)" to "metric * value - metric * code"
  Vec3 codes[3]; Codebook3n(codes, m_start, m_end);

  // match each point to the closest code
  u8 closest[16];

  Scr3 error = Scr3(DEVIANCE_BASE);
  for (int i = 0; i < count; ++i) {
    int idx = 0;

    // find the closest code
    Vec3 value = Normalize(scale * (offset + values[i]));
    Scr3 dist; MinDeviance3<true>(dist, idx, value, codes);

    // accumulate the error
    AddDeviance(dist, error, freq[i]);

    // save the index
    closest[i] = (u8)idx;
  }

  // save this scheme if it wins
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16]; m_colours->RemapIndices(closest, indices);

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
  Scr3 const* freq = m_colours->GetWeights();
  
  // use a fitting algorithm
  m_start = m_start_candidate;
  m_end   = m_end_candidate;
  if (m_flags & kColourIterativeClusterFits)
    kMeans4();
//if ((m_flags & kColourIterativeClusterFits) >= (kColourIterativeClusterFit))
//  Permute4();
  
  // create a codebook
  Vec3 codes[4]; Codebook4n(codes, m_start, m_end);

  // match each point to the closest code
  u8 closest[16];

  Scr3 error = Scr3(DEVIANCE_BASE);
  for (int i = 0; i < count; ++i) {
    int idx = 0;

    // find the closest code
    Vec3 value = Normalize(scale * (offset + values[i]));
    Scr3 dist; MinDeviance4<true>(dist, idx, value, codes);

    // accumulate the error
    AddDeviance(dist, error, freq[i]);

    // save the index
    closest[i] = (u8)idx;
  }

  // save this scheme if it wins
  if (error < m_besterror) {
    // save the error
    m_besterror = error;

    // remap the indices
    u8 indices[16]; m_colours->RemapIndices(closest, indices);

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
