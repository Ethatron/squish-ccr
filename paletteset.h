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

#ifndef SQUISH_PALETTESET_H
#define SQUISH_PALETTESET_H

#include <squish.h>
#include <memory.h>
#include "maths.h"

namespace squish {

// -----------------------------------------------------------------------------
#if	!defined(USE_PRE)
/*! @brief Represents a set of block palettes
*/
class PaletteSet
{
public:
  static void GetMasks(int flags, int partition, int (&masks)[4]);

  // maximum number of different sets, aligned, the real limit is 3
#define	PS_MAX	4

public:
  // constructor for regular operation
  PaletteSet(u8 const* rgba, int mask, int flags, int partition, int rotation);
  // constructors for managing backups of palette-sets
  PaletteSet() {};
  PaletteSet(PaletteSet const &palette) { memcpy(this, &palette, sizeof(this)); };

  // active attributes based the parameters passed on initializaton
  int GetSets() const { return m_numsets; }
  int GetRotation() const { return m_rotid; }
  int GetPartition() const { return m_partid; }

  // information determined when the palette-set has been formed
  bool IsTransparent() const { return m_transparent; }
  bool IsSeperateAlpha() const { return /*m_transparent &&*/ m_seperatealpha; }
  Vec4 const* GetPoints(int idx) const { return m_points[idx]; }
  float const* GetWeights(int idx) const { return m_weights[idx]; }
  int GetCount(int idx) const { return m_count[idx]; }
  int GetCount() const {
    return             m_count[0]      +
    /*(m_numsets > 0 ? m_count[0] : 0)*/ (m_seperatealpha ? m_count[m_numsets + 0] : 0) +
      (m_numsets > 1 ? m_count[1] : 0) /*(m_seperatealpha ? m_count[m_numsets + 1] : 0)*/ +
      (m_numsets > 2 ? m_count[2] : 0) /*(m_seperatealpha ? m_count[m_numsets + 2] : 0)*/; }
  u8 const* GetFrequencies(int idx) const { return m_frequencies[idx]; }
  u8 GetMaxFrequency(int idx) const {
    u8 mx = 1; for (int i = 0; i < m_count[idx]; ++i) mx = std::max(mx, m_frequencies[idx][i]); return mx; }

  // map from the set to indices and back to colours
  void RemapIndices(u8 const* source, u8* target, int set) const;
  void UnmapIndices(u8 const* source, u8* rgba, int set, int *codes, int cmask) const;

private:
  int   m_numsets;
  int   m_rotid;
  int   m_partid;
  int   m_partmask;
  bool  m_seperatealpha;

  int   m_mask[4];
  int   m_count[4];
  Vec4  m_points[4][16];
  float m_weights[4][16];
  u8    m_frequencies[4][16];
  int   m_remap[4][16];
  bool  m_transparent;
};
#endif

// -----------------------------------------------------------------------------
#if	defined(USE_AMP) || defined(USE_COMPUTE)
struct PaletteSet_CCR
{
public_hlsl
  int GetCount() amp_restricted;
  point16 GetPoints() amp_restricted;
  weight16 GetWeights() amp_restricted;
  bool IsTransparent() amp_restricted;

  void CountSet    (tile_barrier barrier, const int thread,
		    pixel16 rgba, int mask, const bool tresh, const bool trans) amp_restricted;
  void WriteSet    (tile_barrier barrier, const int thread,
		    lineC2 cline, inout index16x2 source, out code64 block, const int is4,
		    IndexBlockLUT yArr) amp_restricted;

protected_hlsl
  void RemapIndices(tile_barrier barrier, const int thread,
		    inout index16x2 source) amp_restricted;

#if	!defined(USE_COMPUTE)
private_hlsl
  int m_transparent;
  int     m_count;
  float4 m_points[16];
  float m_weights[16];
  int     m_remap[16];
  ccr8  m_indices[16];
#endif
};

#if	defined(USE_COMPUTE)
  tile_static int m_transparent;
  tile_static int     m_count;
  tile_static float4 m_points[16];
  tile_static float m_weights[16];
  tile_static int     m_remap[16];
  tile_static ccr8  m_indices[16];
#endif

#if	!defined(USE_COMPUTE)
  typedef PaletteSet_CCR &PaletteSet_CCRr;
#else
  typedef PaletteSet_CCR  PaletteSet_CCRr;
#endif

#endif

} // namespace sqish

#endif // ndef SQUISH_PALETTESET_H