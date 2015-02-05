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

#include "palettefit.h"
#include "paletteset.h"

#if 1 //ndef NDEBUG
#include "inlineables.cpp"
#endif

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
extern Vec4 g_metric[8];

static const int skip[2][4] = {
  {0, 1 /* 1 to 1 */,  3 /* 3 to  3 */,  7 /*  7 to  7 */},
  {0, 3 /* 1 to 3 */, 15 /* 7 to 15 */, 63 /* 37 to 63 */}
};

static const int maps[2][64] = {
  { 0 << SBSTART | 0 << SBEND, 1 << SBSTART | 1 << SBEND, // 000|000 + 001|001
    2 << SBSTART | 2 << SBEND, 3 << SBSTART | 3 << SBEND, // 010|010 + 011|011
    4 << SBSTART | 4 << SBEND, 5 << SBSTART | 5 << SBEND, // 100|100 + 101|101
    6 << SBSTART | 6 << SBEND, 7 << SBSTART | 7 << SBEND} // 110|110 + 111|111
  ,
  { 0 << SBSTART | 0 << SBEND,
			       1 << SBSTART | 0 << SBEND, 0 << SBSTART | 1 << SBEND, 1 << SBSTART | 1 << SBEND,

    0 << SBSTART | 2 << SBEND,
    2 << SBSTART | 0 << SBEND,
    2 << SBSTART | 2 << SBEND,
			       1 << SBSTART | 2 << SBEND, 0 << SBSTART | 3 << SBEND, 1 << SBSTART | 3 << SBEND,
                               3 << SBSTART | 0 << SBEND, 2 << SBSTART | 1 << SBEND, 3 << SBSTART | 1 << SBEND,
                               3 << SBSTART | 2 << SBEND, 2 << SBSTART | 3 << SBEND, 3 << SBSTART | 3 << SBEND,

    0 << SBSTART | 4 << SBEND,
    0 << SBSTART | 6 << SBEND,
    2 << SBSTART | 4 << SBEND,
    2 << SBSTART | 6 << SBEND,
    4 << SBSTART | 0 << SBEND,
    4 << SBSTART | 2 << SBEND,
    4 << SBSTART | 4 << SBEND,
    4 << SBSTART | 6 << SBEND,
    6 << SBSTART | 0 << SBEND,
    6 << SBSTART | 2 << SBEND,
    6 << SBSTART | 4 << SBEND,
    6 << SBSTART | 6 << SBEND,
                               5 << SBSTART | 4 << SBEND, 4 << SBSTART | 5 << SBEND, 5 << SBSTART | 5 << SBEND,
                               5 << SBSTART | 0 << SBEND, 4 << SBSTART | 1 << SBEND, 5 << SBSTART | 1 << SBEND,
			       1 << SBSTART | 4 << SBEND, 0 << SBSTART | 5 << SBEND, 1 << SBSTART | 5 << SBEND,

                               1 << SBSTART | 6 << SBEND, 0 << SBSTART | 7 << SBEND, 1 << SBSTART | 7 << SBEND,

                               3 << SBSTART | 4 << SBEND, 2 << SBSTART | 5 << SBEND, 3 << SBSTART | 5 << SBEND,
                               3 << SBSTART | 6 << SBEND, 2 << SBSTART | 7 << SBEND, 3 << SBSTART | 7 << SBEND,
                               5 << SBSTART | 2 << SBEND, 4 << SBSTART | 3 << SBEND, 5 << SBSTART | 3 << SBEND,
                               5 << SBSTART | 6 << SBEND, 4 << SBSTART | 7 << SBEND, 5 << SBSTART | 7 << SBEND,
                               7 << SBSTART | 0 << SBEND, 6 << SBSTART | 1 << SBEND, 7 << SBSTART | 1 << SBEND,
                               7 << SBSTART | 2 << SBEND, 6 << SBSTART | 3 << SBEND, 7 << SBSTART | 3 << SBEND,
                               7 << SBSTART | 4 << SBEND, 6 << SBSTART | 5 << SBEND, 7 << SBSTART | 5 << SBEND,
                               7 << SBSTART | 6 << SBEND, 6 << SBSTART | 7 << SBEND, 7 << SBSTART | 7 << SBEND}
};

/*  Mode NS PB RB ISB CB AB EPB SPB IB IB2	NS: Number of subsets in each partition
 *  ---- -- -- -- --- -- -- --- --- -- ---	PB: Partition bits
 *  0    3  4  0  0   4  0  1   0   3  0	RB: Rotation bits
 *  1    2  6  0  0   6  0  0   1   3  0	ISB: Index selection bits
 *  2    3  6  0  0   5  0  0   0   2  0	CB: Color bits
 *  3    2  6  0  0   7  0  1   0   2  0	AB: Alpha bits
 *  4    1  0  2  1   5  6  0   0   2  3	EPB: Endpoint P-bits
 *  5    1  0  2  0   7  8  0   0   2  2	SPB: Shared P-bits
 *  6    1  0  0  0   7  7  1   0   4  0	IB: Index bits per element
 *  7    2  6  0  0   5  5  1   0   2  0	IB2: Secondary index bits per element
 */
static const struct {
  char
   NS,PB,RB,ISB,CB,AB,EPB,SPB,IB,IB2;
} PBcfg[8] = {
  { 3, 4, 0, 0,  4, 0, 1,  0,  3, 0 },
  { 2, 6, 0, 0,  6, 0, 0,  1,  3, 0 },
  { 3, 6, 0, 0,  5, 0, 0,  0,  2, 0 },
  { 2, 6, 0, 0,  7, 0, 1,  0,  2, 0 },
  { 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
  { 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },
  { 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
  { 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 }
};

int PaletteFit::GetNumSets(int mode) {
  return
    PBcfg[mode].NS;
}

int PaletteFit::GetPartitionBits(int mode) {
  return
    PBcfg[mode].PB;
}

int PaletteFit::GetIndexBits(int mode) {
  return
    PBcfg[mode].IB + (PBcfg[mode].IB2 << 16);
}

int PaletteFit::GetRotationBits(int mode) {
  return
    PBcfg[mode].RB;
}

int PaletteFit::GetSelectionBits(int mode) {
  return
    PBcfg[mode].ISB;
}

int PaletteFit::GetSharedBits(int mode) {
  if (PBcfg[mode].EPB) return (1 << (PBcfg[mode].NS * PBcfg[mode].EPB * 2)) - 1;
  if (PBcfg[mode].SPB) return (1 << (PBcfg[mode].NS * PBcfg[mode].SPB * 1)) - 1;
  return 0;
}

const int *PaletteFit::GetSharedMap(int mode) {
  if (PBcfg[mode].EPB) return maps[1];
  if (PBcfg[mode].SPB) return maps[0];
  return NULL;
}

int PaletteFit::GetSharedSkip(int mode) {
  if (PBcfg[mode].EPB) return skip[1][PBcfg[mode].NS];
  if (PBcfg[mode].SPB) return skip[0][PBcfg[mode].NS];
  return NULL;
}

int PaletteFit::GetPrecisionBits(int mode) {
  return
    ((PBcfg[mode].CB + (PBcfg[mode].CB ? PBcfg[mode].EPB + PBcfg[mode].SPB : 0)) <<  0) |
    ((PBcfg[mode].AB + (PBcfg[mode].AB ? PBcfg[mode].EPB + PBcfg[mode].SPB : 0)) << 16);
}

PaletteFit::PaletteFit(PaletteSet const* palette, int flags, int swap, int shared)
  : m_palette(palette), m_swapindex(-1), m_flags(flags), m_sharedbits(-1)
{
  int const ix = m_palette->GetRotation();

  // the weighting of colour against alpha is 1:1, which
  // is a middle-ground of quality. increasing the alpha-weight
  // will start showing false colour assignments, decreasing
  // alpha-weight will start showing flatness of moderately
  // transparent areas
#ifdef FEATURE_METRIC_ROOTED
  const Vec4 onehalf(sqrtf(0.5f));
  const Vec4 onethird(sqrtf(0.3333f));
#else
  const Vec4 onehalf(0.5f);
  const Vec4 onethird(0.3333f);
#endif

  // initialize the metric
  m_metric[0] = g_metric[(flags & kColourMetrics) >> 4];

  // sum is 1.0f
  if (!m_palette->IsTransparent())
    m_metric[0] = KillW(m_metric[0]);
  else
    m_metric[0] = m_metric[0] * onehalf;

  // swap channel-weights
  switch (ix) {
    case 1: m_metric[0] = Exchange<0,3>(m_metric[0]); break;
    case 2: m_metric[0] = Exchange<1,3>(m_metric[0]); break;
    case 3: m_metric[0] = Exchange<2,3>(m_metric[0]); break;
  }

  // split metric into two, double sum is 1.0f
  if (m_palette->IsSeperateAlpha()) {
    // alpha-metric in xyz
    m_metric[2] = KillW(m_metric[0].SplatW()) * onethird;

    // alpha-metric in w
    m_metric[1] = OnlyW(m_metric[0]);
    // color-metric in xyz
    m_metric[0] = KillW(m_metric[0]);
  }
  else
    m_metric[2] =
    m_metric[1] =
    m_metric[0];

  // initialize the best error
  m_besterror = Scr4(FLT_MAX);
  m_best = false;

  m_mode = ((m_flags & kVariableCodingModes) >> 24) - 1;
  m_sharedmap = GetSharedMap(m_mode);
  m_swapindex = swap;
  m_sharedbits = SR(shared);
}

void PaletteFit::Compress(void* block, vQuantizer &q)
{
  if (m_mode < 0) {
    /*
    for (int m = 0; m < 7; m++) {
      if (m_palette->GetSets() == GetNumSets(m))
      if (m_palette->GetRotation() < GetRotationBits(m))
      if (m_palette->GetPartition() < GetPartitionBits(m))
        ...
    }
    */

    switch (m_palette->GetSets()) {
      case 1:
	if (m_palette->GetPartition() < 1) {
	  // 1 partition
	  if (!m_palette->GetRotation())
	    m_mode = 4, Compress(block, q);

	  // 1 swap-bit
	  if (m_swapindex == -1) {
	    m_mode = 5, m_swapindex = 0, Compress(block, q);
	    m_mode = 5, m_swapindex = 1, Compress(block, q);
	    m_swapindex = -1;
	  }

	  // 0 swap-bit
	  m_mode = 6, Compress(block, q);
	}
	break;
      case 2:
	if (m_palette->GetRotation())
	  break;

	if (!m_palette->IsTransparent()) {
	  // 64 partitions
	  m_mode = 1, Compress(block, q);
	  m_mode = 3, Compress(block, q);
	}

	// 64 partitions
	m_mode = 7, Compress(block, q);
	break;
      case 3:
	if (m_palette->GetRotation())
	  break;
	if (m_palette->IsTransparent())
	  break;

	if (m_palette->GetPartition() < (1 << GetPartitionBits(0)))
	  // mode 0 has only 16 partitions
	  m_mode = 0, Compress(block, q);

	// 64 partitions
	m_mode = 2, Compress(block, q);
	break;
    }
  }
  else
    Compress(block, q, m_mode);
}

#if 1 //ndef NDEBUG
void PaletteFit::SumError(u8 (&closest)[4][16], vQuantizer &q, int mode, Scr4 &error) {
  int ib = GetIndexBits(mode);
  int jb = ib >> 16; ib = ib & 0xFF;
  int cb = GetPrecisionBits(mode);
  int ab = cb >> 16; cb = cb & 0xFF;
  int zb = GetSharedField();

  q.ChangeShared(cb, cb, cb, ab, zb);

  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;

  assume((isets >  0) && (isets <= 3));
  assume((asets >= 0) && (asets <= 3));
  assume(((isets    +    asets) <= 3));

  // create a codebook
  Vec4 codes[1 << 4];

  // loop over all sets
  for (int s = 0, sb = zb; s < (isets + asets); s++, sb >>= 1) {
    // how big is the codebook for the current set
    // swap the code-book when the swap-index bit is set
    int kb = ((s < isets) ^ (!!m_swapindex)) ? ib : jb;

    // in case of separate alpha the colors of the alpha-set have all been set to alpha
    Vec4 metric = m_metric[s < isets ? 0 : 1];

    // cache some values
    int const count = m_palette->GetCount(s);
    Vec4 const* values = m_palette->GetPoints(s);
    Scr4 const* freq = m_palette->GetWeights(s);

    // snap floating-point-values to the integer-lattice
    Vec4 start = q.SnapToLattice(m_start[s], sb, 1 << SBSTART);
    Vec4 end   = q.SnapToLattice(m_end  [s], sb, 1 << SBEND);

    // resolve "metric * (value - code)" to "metric * value - metric * code"
    CodebookP(codes, kb, metric * start, metric * end);

    // error-sum
    for (int i = 0; i < count; ++i)
      error += LengthSquared(metric * values[i] - codes[closest[s][i]]) * freq[i];
  }
}

void PaletteFit::Decompress(u8 *rgba, vQuantizer &q, int mode)
{
  int ib = GetIndexBits(mode);
  int jb = ib >> 16; ib = ib & 0xFF;
  int cb = GetPrecisionBits(mode);
  int ab = cb >> 16; cb = cb & 0xFF;
  int zb = GetSharedField();

  q.ChangeShared(cb, cb, cb, ab, zb);

  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;

  assume((isets >  0) && (isets <= 3));
  assume((asets >= 0) && (asets <= 3));
  assume(((isets    +    asets) <= 3));

  // create a codebook
  unsigned int codes[1 << 4];

  // loop over all sets
  for (int s = 0, sb = zb; s < (isets + asets); s++, sb >>= 1) {
    // how big is the codebook for the current set
    int kb = ((s < isets) ^ (!!m_swapindex)) ? ib : jb;
    int mk = !asets ? 0xFFFFFFFF : (s < isets ? 0x00FFFFFF : 0xFF000000);

    // snap floating-point-values to the integer-lattice
    Vec4 fstart = m_start[s];
    Vec4 fend   = m_end  [s];

    // because the original alpha-channel's weight was killed it is completely random and need to be set to 1.0f
    if (!m_palette->IsTransparent()) {
      switch (m_palette->GetRotation()) {
	default: if (s < isets + asets) fstart.Set<3>(1.0f), fend.Set<3>(1.0f); break;
	case  1: if (s <         isets) fstart.Set<0>(1.0f), fend.Set<0>(1.0f); break;
	case  2: if (s <         isets) fstart.Set<1>(1.0f), fend.Set<1>(1.0f); break;
	case  3: if (s <         isets) fstart.Set<2>(1.0f), fend.Set<2>(1.0f); break;
      }
    }

    // snap floating-point-values to the integer-lattice
    Col4 istart = q.QuantizeToInt(fstart, sb, 1 << SBSTART);
    Col4 iend   = q.QuantizeToInt(fend  , sb, 1 << SBEND);

    istart = Col4(
      (istart.R() << (8 - cb)) + (istart.R() >> (cb - (8 - cb))),
      (istart.G() << (8 - cb)) + (istart.G() >> (cb - (8 - cb))),
      (istart.B() << (8 - cb)) + (istart.B() >> (cb - (8 - cb))),
      (istart.A() << (8 - ab)) + (istart.A() >> (ab - (8 - ab)))
    );

    iend = Col4(
      (iend.R() << (8 - cb)) + (iend.R() >> (cb - (8 - cb))),
      (iend.G() << (8 - cb)) + (iend.G() >> (cb - (8 - cb))),
      (iend.B() << (8 - cb)) + (iend.B() >> (cb - (8 - cb))),
      (iend.A() << (8 - ab)) + (iend.A() >> (ab - (8 - ab)))
    );

    int ccs;
    switch (kb) {
      case 2: ccs = CodebookP<2>(codes, istart, iend); break;
      case 3: ccs = CodebookP<3>(codes, istart, iend); break;
      case 4: ccs = CodebookP<4>(codes, istart, iend); break;
    }

    m_palette->UnmapIndices(m_indices[s < isets ? 0 : 1], rgba, s, codes, mk);
  }

  switch (m_palette->GetRotation()) {
    case 1: for (int i = 0; i < 16; ++i) std::swap(rgba[4 * i + 0], rgba[4 * i + 3]); break;
    case 2: for (int i = 0; i < 16; ++i) std::swap(rgba[4 * i + 1], rgba[4 * i + 3]); break;
    case 3: for (int i = 0; i < 16; ++i) std::swap(rgba[4 * i + 2], rgba[4 * i + 3]); break;
  }

  if (!m_palette->IsTransparent()) {
    for (int i = 0; i < 16; ++i) {
      assert(rgba[4 * i + 3] == 0xFF);
    }
  }
}
#endif
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
void PaletteFit_CCR::AssignSet(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, int metric, int fit) amp_restricted
{
  // -1 remap goes to source[0][15]
  wavefrnt_for(miscan, 16) {
    m_matches[0][miscan] = 3;
//  m_matches[1][miscan] = 0;
  }
}

#if 0
void PaletteFit_CCR::Compress(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, out code64 block, const bool trans,
			     IndexBlockLUT yArr) amp_restricted
{
  /* all or nothing branches, OK, same for all threads */
  bool isBtc1 = (trans);
  if (isBtc1 && m_palette.IsTransparent())
    Compress3 (barrier, thread, m_palette, block, yArr);
  else if (!isBtc1)
    Compress4 (barrier, thread, m_palette, block, yArr);
  else
    Compress34(barrier, thread, m_palette, block, yArr);
}
#endif

void PaletteFit_CCR::Compress3(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, out code64 block,
			      IndexBlockLUT yArr) amp_restricted
{
  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    m_palette.WriteSet(barrier, thread, m_line, m_matches, block, 0, yArr);
  }
}

void PaletteFit_CCR::Compress4(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, out code64 block,
			      IndexBlockLUT yArr) amp_restricted
{
  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    m_palette.WriteSet(barrier, thread, m_line, m_matches, block, 1, yArr);
  }
}

void PaletteFit_CCR::Compress34(tile_barrier barrier, const int thread, PaletteSet_CCRr m_palette, out code64 block, const int is4,
			       IndexBlockLUT yArr) amp_restricted
{
  // save this scheme if it wins
  {
    // remap the indices
    // save the block
    m_palette.WriteSet(barrier, thread, m_line, m_matches, block, is4, yArr);
  }
}
#endif

} // namespace squish
