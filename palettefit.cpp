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
#if	!defined(USE_PRE)
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
    PBcfg[mode].IB | (PBcfg[mode].IB2 << 16);
}

int PaletteFit::GetRotationBits(int mode) {
  return
    PBcfg[mode].RB;
}

int PaletteFit::GetSelectionBits(int mode) {
  return
    PBcfg[mode].ISB;
}

int PaletteFit::GetPrecisionBits(int mode) {
  return
    ((PBcfg[mode].CB + (PBcfg[mode].CB ? PBcfg[mode].EPB + PBcfg[mode].SPB : 0)) <<  0) |
    ((PBcfg[mode].AB + (PBcfg[mode].AB ? PBcfg[mode].EPB + PBcfg[mode].SPB : 0)) << 16);
}

PaletteFit::PaletteFit(PaletteSet const* palette, int flags, int swap)
  : m_palette(palette), m_swapindex(-1), m_flags(flags)
{
  int const ix = m_palette->GetRotation();

  // initialize the metric
  bool perceptual = ((m_flags & kColourMetricPerceptual) != 0);
  bool unit       = ((m_flags & kColourMetricUnit      ) != 0);

  // sum is 2.0f
  if (unit)
    m_metric[0] = Vec4(0.5000f, 0.5000f, 0.0000f, 1.0000f);
  else if (perceptual)	// linear RGB luminance
    m_metric[0] = Vec4(0.2126f, 0.7152f, 0.0722f, 1.0000f);
  else
    m_metric[0] = Vec4(0.3333f, 0.3334f, 0.3333f, 1.0000f);

  // sum is 1.0f
  if (!m_palette->IsTransparent())
    m_metric[0] = KillW(m_metric[0]);
  else
    m_metric[0] = m_metric[0] * 0.5f;

  // swap channel-weights
  switch (ix) {
    case 1: m_metric[0] = Exchange<0, 3>(m_metric[0]); break;
    case 2: m_metric[0] = Exchange<1, 3>(m_metric[0]); break;
    case 3: m_metric[0] = Exchange<2, 3>(m_metric[0]); break;
  }

  // split metric into two, double sum is 1.0f
  if (m_palette->IsSeperateAlpha()) {
    // alpha-metric in xyz
    m_metric[2] = KillW(m_metric[0].SplatW()) * Vec4(0.3333f, 0.3334f, 0.3333f, 0.0000f);

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
  m_besterror = Vec4(FLT_MAX);
  m_best = false;

  m_swapindex = swap;
}

void PaletteFit::Compress(void* block)
{
  int varMode = ((m_flags & kVariableCodingModes) >> 24) - 1;
  if (varMode < 0) {
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
	    Compress(block, 4);

	  // 1 swap-bit
	  if (m_swapindex == -1) {
	    m_swapindex = 0;
	    Compress(block, 5);
	    m_swapindex = 1;
	    Compress(block, 5);
	    m_swapindex = -1;
	  }

	  // 0 swap-bit
	  Compress(block, 6);
	}
	break;
      case 2:
	if (m_palette->GetRotation())
	  break;

	if (!m_palette->IsTransparent()) {
	  // 64 partitions
	  Compress(block, 1);
	  Compress(block, 3);
	}

	// 64 partitions
	Compress(block, 7);
	break;
      case 3:
	if (m_palette->GetRotation())
	  break;
	if (m_palette->IsTransparent())
	  break;

	if (m_palette->GetPartition() < (1 << GetPartitionBits(0)))
	  // mode 0 has only 16 partitions
	  Compress(block, 0);

	// 64 partitions
	Compress(block, 2);
	break;
    }
  }
  else
    Compress(block, varMode);
}

#if 1 //ndef NDEBUG
void PaletteFit::SumError(u8 (&closest)[4][16], int mode, Vec4 &error) {
  int ib = GetIndexBits(mode);
  int jb = ib >> 16; ib = ib & 0xFF;
  int cb = GetPrecisionBits(mode);
  int ab = cb >> 16; cb = cb & 0xFF;
  
  vQuantizer q = vQuantizer(cb, cb, cb, ab);

  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;

  // create a codebook
  Vec4 codes[1 << 4];

  // loop over all sets
  for (int s = 0; s < (isets + asets); s++) {
    // how big is the codebook for the current set
    int kb = ((s < isets) ^ (!!m_swapindex)) ? ib : jb;

    // in case of separate alpha the colors of the alpha-set have all been set to alpha
    Vec4 metric = m_metric[s < isets ? 0 : 1];

    // cache some values
    int const count = m_palette->GetCount(s);
    Vec4 const* values = m_palette->GetPoints(s);
    u8 const* freq = m_palette->GetFrequencies(s);
    
    // snap floating-point-values to the integer-lattice
    Vec4 start = q.SnapToLattice(m_start[s]);
    Vec4 end   = q.SnapToLattice(m_end  [s]);

    // swap the code-book when the swap-index bit is set
    CodebookP(codes, kb, start, end);

    // error-sum
    for (int i = 0; i < count; ++i)
      error += LengthSquared(metric * (values[i] - codes[closest[s][i]])) * freq[i];
  }
}

void PaletteFit::Decompress(u8 *rgba, int mode)
{
  int ib = GetIndexBits(mode);
  int jb = ib >> 16; ib = ib & 0xFF;
  int cb = GetPrecisionBits(mode);
  int ab = cb >> 16; cb = cb & 0xFF;
  
  vQuantizer q = vQuantizer(cb, cb, cb, ab);

  // snap floating-point-values to the integer-lattice and save
#define gridp(b)  ((1 << b) - 1)
  Vec4 const grid   (1.0f * gridp(cb), 1.0f * gridp(cb), 1.0f * gridp(cb), ab ? 1.0f * gridp(ab) : 1.0f);
  Vec4 const gridrcp(1.0f / gridp(cb), 1.0f / gridp(cb), 1.0f / gridp(cb), ab ? 1.0f / gridp(ab) : 1.0f);
  Vec4 const half   (0.5f);

  // the alpha-set (in theory we can do separate alpha + separate partitioning, but's not codeable)
  int const isets = m_palette->GetSets();
  int const asets = m_palette->IsSeperateAlpha() ? isets : 0;

  // create a codebook
  int codes[1 << 4];

  // loop over all sets
  for (int s = 0; s < (isets + asets); s++) {
    // how big is the codebook for the current set
    int kb = ((s < isets) ^ (!!m_swapindex)) ? ib : jb;
    int mk = !asets ? 0xFFFFFFFF : (s < isets ? 0x00FFFFFF : 0xFF000000);

    // old original quantizer
    Vec4 zstart = Truncate(MultiplyAdd(grid, m_start[s].Clamp(), half)) * gridrcp;
    Vec4 zend   = Truncate(MultiplyAdd(grid, m_end  [s].Clamp(), half)) * gridrcp;

    Col4 xstart = FloatToInt<true>(zstart * 255.0f);
    Col4 xend   = FloatToInt<true>(zend * 255.0f);
    
    // new exact lattice quantizer
    Vec4 fstart = q.SnapToLattice(m_start[s]);
    Vec4 fend   = q.SnapToLattice(m_end  [s]);
    
    Col4 jstart = q.QuantizeToInt(m_start[s]);
    Col4 jend   = q.QuantizeToInt(m_end  [s]);
    
    Col4 estart = jstart * Col4(1 << (8 - cb), 1 << (8 - cb), 1 << (8 - cb), 1 << (8 - ab));
    Col4 eend   = jend   * Col4(1 << (8 - cb), 1 << (8 - cb), 1 << (8 - cb), 1 << (8 - ab));
    
         estart = estart * Col4((1 << 8) + (1 << (8 - cb)), (1 << 8) + (1 << (8 - cb)), (1 << 8) + (1 << (8 - cb)), (1 << 8) + (1 << (8 - ab)));
         eend   = eend   * Col4((1 << 8) + (1 << (8 - cb)), (1 << 8) + (1 << (8 - cb)), (1 << 8) + (1 << (8 - cb)), (1 << 8) + (1 << (8 - ab)));
    
         estart = estart >> 8;
         eend   = eend   >> 8;
    
    Col4 istart = FloatToInt<true>(fstart * 255.0f);
    Col4 iend   = FloatToInt<true>(fend * 255.0f);

    int ccs;
    switch (kb) {
      case 2: ccs = CodebookP<2>(codes, istart, iend); break;
      case 3: ccs = CodebookP<3>(codes, istart, iend); break;
      case 4: ccs = CodebookP<4>(codes, istart, iend); break;
    }

    m_palette->UnmapIndices(m_indices[s < isets ? 0 : 1], rgba, s, codes, mk);
  }

  switch (m_palette->GetRotation()) {
    case 1:
      for (int i = 0; i < 16; ++i)
	std::swap(rgba[4 * i + 0], rgba[4 * i + 3]);
      break;
    case 2:
      for (int i = 0; i < 16; ++i)
	std::swap(rgba[4 * i + 1], rgba[4 * i + 3]);
      break;
    case 3:
      for (int i = 0; i < 16; ++i)
	std::swap(rgba[4 * i + 2], rgba[4 * i + 3]);
      break;
  }

  for (int i = 0; i < 16; ++i) {
    assert(rgba[4 * i + 3] == 0xFF);
  }
}
#endif
#endif

/* *****************************************************************************
 */
#if	defined(USE_AMP) || defined(USE_COMPUTE)
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
