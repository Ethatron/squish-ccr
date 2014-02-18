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

#include "palettesinglefit.h"
#include "paletteset.h"
#include "paletteblock.h"

#include "inlineables.cpp"

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
struct SP_SourceBlock
{
  u8 start;
  u8 end;
  u8 error;
};

struct PaletteSingleLookup2
{
  SP_SourceBlock sources[2];
};

struct PaletteSingleLookup4
{
  SP_SourceBlock sources[4];
};

struct PaletteSingleLookup8
{
  SP_SourceBlock sources[8];
};

#define	SCL_ITERATIVE
#include "palettesinglelookup.inl"

PaletteSingleFit::PaletteSingleFit(PaletteSet const* palette, int flags, int swap, int shared)
  : PaletteFit(palette, flags, swap, shared)
{
}

Scr4 PaletteSingleFit::ComputeEndPoints(int set, Vec4 const &metric, int cb, int ab, int sb, int ib, u8 cmask)
{
#if	!defined(FEATURE_SHAREDBITS_TRIALS)
  // silence the compiler
  bool hb = !!sb; hb = false;

#define sp_lookup_5u1_4_sb_    sp_lookup_6_4
#define sp_lookup_7u1_4_sb_    sp_lookup_8_4
#define sp_lookup_4u1_8_sb_    sp_lookup_5_8
#define sp_lookup_6s1_8_sb_    sp_lookup_7_8
#define sp_lookup_7u1_16_sb_   sp_lookup_8_16

#define sp_lookup_5u1_4_ck_    sp_lookup_6_4
#define sp_lookup_4u1_8_ck_    sp_lookup_5_8
#else
  // merge start and end shared bits
  int lb = (sb & 1) + ((sb >> (SBEND - 1)) & 2);

  // allow bailout if whole == -1
#define sp_lookup_5u1_4_sb_    (SK(sb) ? sp_lookup_6_4  : sp_lookup_5u1_4[lb])
#define sp_lookup_7u1_4_sb_    (SK(sb) ? sp_lookup_8_4  : sp_lookup_7u1_4[lb])
#define sp_lookup_4u1_8_sb_    (SK(sb) ? sp_lookup_5_8  : sp_lookup_4u1_8[lb])
#define sp_lookup_6s1_8_sb_    (SK(sb) ? sp_lookup_7_8  : sp_lookup_6s1_8[lb&1])
#define sp_lookup_7u1_16_sb_   (SK(sb) ? sp_lookup_8_16 : sp_lookup_7u1_16[lb])

#define sp_lookup_5u1_4_ck_    (SK(sb) ? sp_lookup_6_4  : sp_lookup_5u1_4_sb_)
#define sp_lookup_4u1_8_ck_    (SK(sb) ? sp_lookup_5_8  : sp_lookup_4u1_8_sb_)
#endif

  assume(ib >= 2 && ib <= 4);
  switch (ib) {
    case 2: {
      PaletteSingleLookup2 const* cl;
      PaletteSingleLookup2 const* al;

      assume(cb >= 5 && cb <= 8);
      switch (cb) {
	case  5: cl = sp_lookup_5_4; break;		//{ 3, 6, 0, 0,  5, 0, 0,  0,  2, 0 },
	case  6: cl = sp_lookup_5u1_4_sb_; break;	//{ 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 },
	case  7: cl = sp_lookup_7_4; break;		//{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },
	case  8: cl = sp_lookup_7u1_4_sb_; break;	//{ 2, 6, 0, 0,  7, 0, 1,  0,  2, 0 },
	default: assume(0); break;
      }

      assume(ab == 0 || ab == 6 || ab == 8);
      switch (ab) {
	case  6: al = sp_lookup_5u1_4_ck_; break;	//{ 2, 6, 0, 0,  5, 5, 1,  0,  2, 0 }, / { 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
	case  8: al = sp_lookup_8_4; break;		//{ 1, 0, 2, 0,  7, 8, 0,  0,  2, 2 },
	default: al = NULL; break;			//{ 3, 6, 0, 0,  5, 0, 0,  0,  2, 0 }, / { 2, 6, 0, 0,  7, 0, 1,  0,  2, 0 },
      }

      PaletteSingleLookup2 const* const lookups[] =
      { cl, cl, cl, al };

      return ComputeEndPoints(set, metric, lookups, cmask);
    } break;
    case 3: {
      PaletteSingleLookup4 const* cl;
      PaletteSingleLookup4 const* al;

      assume(cb == 5 || cb == 7);
      switch (cb) {
	case  5: cl = sp_lookup_4u1_8_ck_; break;	//{ 3, 4, 0, 0,  4, 0, 1,  0,  3, 0 }, / { 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
	case  7: cl = sp_lookup_6s1_8_sb_; break;	//{ 2, 6, 0, 0,  6, 0, 0,  1,  3, 0 },
	default: assume(0); break;
      }

      assume(ab == 0 || ab == 6);
      switch (ab) {
	case  6: al = sp_lookup_6_8; break;		//{ 1, 0, 2, 1,  5, 6, 0,  0,  2, 3 },
	default: assume(0); break;
      }

      PaletteSingleLookup4 const* const lookups[] =
      { cl, cl, cl, al };

      return ComputeEndPoints(set, metric, lookups, cmask);
    } break;
    case 4: {
      PaletteSingleLookup8 const* cl;
      PaletteSingleLookup8 const* al;

      assume(cb == 8);
      switch (cb) {
	case  8: cl = sp_lookup_7u1_16_sb_; break;	//{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
	default: assume(0); break;
      }

      assume(ab == 8);
      switch (ab) {
	case  8: al = sp_lookup_7u1_16_sb_; break;	//{ 1, 0, 0, 0,  7, 7, 1,  0,  4, 0 },
	default: assume(0); break;
      }

      PaletteSingleLookup8 const* const lookups[] =
      { cl, cl, cl, al };

      return ComputeEndPoints(set, metric, lookups, cmask);
    } break;

    default:
      return Scr4(INT_MAX);
      break;
  }
}

Scr4 PaletteSingleFit::ComputeEndPoints(int set, Vec4 const &metric, PaletteSingleLookup2 const* const* lookups, u8 mask)
{
  // check each index combination (endpoint or intermediate)
  Scr4 besterror = Scr4(FLT_MAX);

  // grab the single entry
  Vec4 const* values = m_palette->GetPoints(set);
  const float *eLUT = ComputeGammaLUT(false);

  // values are directly out of the codebook and
  // natural numbers / 255, no need to round
  PackBytes(FloatToInt<true>((*values) * Vec4(255.0f)), (unsigned int &)(m_entry[set]));

  /*
  assert(m_entry[set][0] == (u8)FloatToInt<true,false>(255.0f * values->X(), 255));
  assert(m_entry[set][1] == (u8)FloatToInt<true,false>(255.0f * values->Y(), 255));
  assert(m_entry[set][2] == (u8)FloatToInt<true,false>(255.0f * values->Z(), 255));
  assert(m_entry[set][3] == (u8)FloatToInt<true,false>(255.0f * values->W(), 255));
   */

  for (int index = 0; index < 2; ++index) {
    // check the error for this codebook index
    SP_SourceBlock const* sources[4] = {NULL};
    Vec4 cerror(0.0f);
    
    for (int channel = 0, cmask = mask; channel < 4; ++channel, cmask >>= 1) {
      // skip if it's completely irrelevant what's in a specific channel
      if (cmask & 1) {
	// grab the lookup table and index for this channel
	PaletteSingleLookup2 const* lookup = lookups[channel];
	int target = m_entry[set][channel];

	// store a pointer to the source for this channel
	sources[channel] = lookup[target].sources + index;

	// accumulate the error
	cerror.GetO(channel) = eLUT[sources[channel]->error];
      }
    }

    // calculate the error
    Scr4 error = LengthSquared(metric * cerror);

    // keep it if the error is lower
    if (error < besterror) {
      // save the error
      besterror = error;

      Col4 s = Col4(
	sources[0] ? sources[0]->start : 0x00,
	sources[1] ? sources[1]->start : 0x00,
	sources[2] ? sources[2]->start : 0x00,
	sources[3] ? sources[3]->start : 0xFF);
      Col4 e = Col4(
	sources[0] ? sources[0]->end   : 0x00,
	sources[1] ? sources[1]->end   : 0x00,
	sources[2] ? sources[2]->end   : 0x00,
	sources[3] ? sources[3]->end   : 0xFF);

      m_start[set] = Vec4(s) * (1.0f / 255.0f);
      m_end  [set] = Vec4(e) * (1.0f / 255.0f);

      m_index = (u8)(1 * index);

      // early out
      if (!(besterror > Scr4(0.0f)))
	return besterror;
    }
  }

  return besterror;
}

Scr4 PaletteSingleFit::ComputeEndPoints(int set, Vec4 const &metric, PaletteSingleLookup4 const* const* lookups, u8 mask)
{
  // check each index combination (endpoint or intermediate)
  Scr4 besterror = Scr4(FLT_MAX);

  // grab the single entry
  Vec4 const* values = m_palette->GetPoints(set);
  const float *eLUT = ComputeGammaLUT(false);

  // values are directly out of the codebook and
  // natural numbers / 255, no need to round
  PackBytes(FloatToInt<true>((*values) * Vec4(255.0f)), (unsigned int &)(m_entry[set]));

  /*
  assert(m_entry[set][0] == (u8)FloatToInt<true,false>(255.0f * values->X(), 255));
  assert(m_entry[set][1] == (u8)FloatToInt<true,false>(255.0f * values->Y(), 255));
  assert(m_entry[set][2] == (u8)FloatToInt<true,false>(255.0f * values->Z(), 255));
  assert(m_entry[set][3] == (u8)FloatToInt<true,false>(255.0f * values->W(), 255));
   */

  for (int index = 0; index < 4; ++index) {
    // check the error for this codebook index
    SP_SourceBlock const* sources[4] = {NULL};
    Vec4 cerror(0.0f);
    
    for (int channel = 0, cmask = mask; channel < 4; ++channel, cmask >>= 1) {
      // skip if it's completely irrelevant what's in a specific channel
      if (cmask & 1) {
	// grab the lookup table and index for this channel
	PaletteSingleLookup4 const* lookup = lookups[channel];
	int target = m_entry[set][channel];

	// store a pointer to the source for this channel
	sources[channel] = lookup[target].sources + index;

	// accumulate the error
	cerror.GetO(channel) = eLUT[sources[channel]->error];
      }
    }

    // calculate the error
    Scr4 error = LengthSquared(metric * cerror);

    // keep it if the error is lower
    if (error < besterror) {
      besterror = error;

      Col4 s = Col4(
	sources[0] ? sources[0]->start : 0x00,
	sources[1] ? sources[1]->start : 0x00,
	sources[2] ? sources[2]->start : 0x00,
	sources[3] ? sources[3]->start : 0xFF);
      Col4 e = Col4(
	sources[0] ? sources[0]->end   : 0x00,
	sources[1] ? sources[1]->end   : 0x00,
	sources[2] ? sources[2]->end   : 0x00,
	sources[3] ? sources[3]->end   : 0xFF);

      m_start[set] = Vec4(s) * (1.0f / 255.0f);
      m_end  [set] = Vec4(e) * (1.0f / 255.0f);

      m_index = (u8)(1 * index);

      // early out
      if (!(besterror > Scr4(0.0f)))
	return besterror;
    }
  }

  return besterror;
}

Scr4 PaletteSingleFit::ComputeEndPoints(int set, Vec4 const &metric, PaletteSingleLookup8 const* const* lookups, u8 mask)
{
  // check each index combination (endpoint or intermediate)
  Scr4 besterror = Scr4(FLT_MAX);

  // grab the single entry
  Vec4 const* values = m_palette->GetPoints(set);
  const float *eLUT = ComputeGammaLUT(false);

  // values are directly out of the codebook and
  // natural numbers / 255, no need to round
  PackBytes(FloatToInt<true>((*values) * Vec4(255.0f)), (unsigned int &)(m_entry[set]));

  /*
  assert(m_entry[set][0] == (u8)FloatToInt<true,false>(255.0f * values->X(), 255));
  assert(m_entry[set][1] == (u8)FloatToInt<true,false>(255.0f * values->Y(), 255));
  assert(m_entry[set][2] == (u8)FloatToInt<true,false>(255.0f * values->Z(), 255));
  assert(m_entry[set][3] == (u8)FloatToInt<true,false>(255.0f * values->W(), 255));
  */ 

  for (int index = 0; index < 8; ++index) {
    // check the error for this codebook index
    SP_SourceBlock const* sources[4] = {NULL};
    Vec4 cerror(0.0f);

    for (int channel = 0, cmask = mask; channel < 4; ++channel, cmask >>= 1) {
      // skip if it's completely irrelevant what's in a specific channel
      if (cmask & 1) {
	// grab the lookup table and index for this channel
	PaletteSingleLookup8 const* lookup = lookups[channel];
	int target = m_entry[set][channel];

	// store a pointer to the source for this channel
	sources[channel] = lookup[target].sources + index;

	// accumulate the error
	cerror.GetO(channel) = eLUT[sources[channel]->error];
      }
    }

    // calculate the error
    Scr4 error = LengthSquared(metric * cerror);

    // keep it if the error is lower
    if (error < besterror) {
      besterror = error;

      Col4 s = Col4(
	sources[0] ? sources[0]->start : 0x00,
	sources[1] ? sources[1]->start : 0x00,
	sources[2] ? sources[2]->start : 0x00,
	sources[3] ? sources[3]->start : 0xFF);
      Col4 e = Col4(
	sources[0] ? sources[0]->end   : 0x00,
	sources[1] ? sources[1]->end   : 0x00,
	sources[2] ? sources[2]->end   : 0x00,
	sources[3] ? sources[3]->end   : 0xFF);

      m_start[set] = Vec4(s) * (1.0f / 255.0f);
      m_end  [set] = Vec4(e) * (1.0f / 255.0f);

      m_index = (u8)(1 * index);

      // early out
      if (!(besterror > Scr4(0.0f)))
	return besterror;
    }
  }

  return besterror;
}
#endif

} // namespace squish

#if	defined(SBL_FLAT)
#include "palettesinglefit_ccr_flat.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 1)
#include "palettesinglefit_ccr_packed.cpp"
#elif	defined(SBL_PACKED) && (SBL_PACKED == 2)
#include "palettesinglefit_ccr_packed_copy.cpp"
#elif	defined(SBL_VECTOR)
#include "palettesinglefit_ccr_vector.cpp"
#endif
