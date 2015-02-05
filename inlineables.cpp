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

#ifndef SQUISH_INLINEABLES_CPP
#define SQUISH_INLINEABLES_CPP

#include <assert.h>
#include "config.h"
#include "math.h"
#include "simd.h"

#pragma warning(disable: 4505)
#pragma warning(disable: 4127)

namespace squish {

extern const u16 weights_u16[5][16];
extern const Vec4 weights_V4[5][16];
extern const Col4 weights_C4[5][16];

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
template<const bool round, const bool clamp>
static doinline int passreg FloatToInt(float a, int limit) ccr_restricted
{
  // use ANSI round-to-zero behaviour to get round-to-nearest
  assert((a >= 0.0f) || !round);
  int i = (int)(a + (round ? 0.5f : 0.0f));

  // clamp to the limit
  if (clamp) {
    if (i < 0)
      i = 0;
    else if (i > limit)
      i = limit;
  }

  // done
  return i;
}

static doinline int passreg FloatTo565(Vec3::Arg colour) ccr_restricted
{
  // get the components in the correct range
  cQuantizer3<5,6,5> q = cQuantizer3<5,6,5>();
  Col3 rgb = q.LatticeToIntClamped(colour);

  int r = rgb.R();
  int g = rgb.G();
  int b = rgb.B();

  /* not necessarily true
  assert(r == FloatToInt<true,false>(31.0f * colour.X(), 31));
  assert(g == FloatToInt<true,false>(63.0f * colour.Y(), 63));
  assert(b == FloatToInt<true,false>(31.0f * colour.Z(), 31));
   */

  // pack into a single value
  return (r << 11) + (g << 5) + b;
}

static doinline int passreg Unpack565(u8 const* packed, u8* colour) ccr_restricted
{
  // build the packed value
  int value = ((int)packed[0] << 0) + ((int)packed[1] << 8);

  // get the components in the stored range
  u8 red   = (u8)((value >> 11) & 0x1F);
  u8 green = (u8)((value >>  5) & 0x3F);
  u8 blue  = (u8)( value        & 0x1F);

  // scale up to 8 bits
  colour[0] = (red   << 3) + (red   >> 2);
  colour[1] = (green << 2) + (green >> 4);
  colour[2] = (blue  << 3) + (blue  >> 2);
  colour[3] = 255;

  // return the value
  return value;
}

static doinline int passreg FloatTo88(Vec3::Arg colour) ccr_restricted
{
  // get the components in the correct range
  Col3 rgb = FloatToInt<true>(colour * 255.0f);

  int r = rgb.R();
  int g = rgb.G();

  /* not necessarily true
  assert(r == FloatToInt(255.0f * colour.X(), 255));
  assert(g == FloatToInt(255.0f * colour.Y(), 255));
   */

  // pack into a single value
  return (r << 8) + g;
}

static doinline int passreg Unpack88(u8 const* packed, u8* colour) ccr_restricted
{
  // build the packed value
  int value = ((int)packed[0] << 0) + ((int)packed[1] << 8);

  // get the components in the stored range
  u8 red   = (u8)((value >> 8) & 0xFF);
  u8 green = (u8)( value       & 0xFF);

  // scale up to 8 bits
  colour[0] = (red  );
  colour[1] = (green);
  colour[2] = 0;
  colour[3] = 255;

  // return the value
  return value;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
#define	COLORA	0
#define	UNIQUE	1
#define	SHARED	1	// 2
#define	FIELDN	2	// 3

static const vQuantizer q7778s0(7, 7, 7, 8,  0);
static const vQuantizer q5556s0(5, 5, 5, 6,  0);
static const vQuantizer q8888s0(8, 8, 8, 8,  0);
static const vQuantizer q6666s0(6, 6, 6, 6,  0);
static const vQuantizer q8880s0(8, 8, 8, 0,  0);
static const vQuantizer q7770s0(7, 7, 7, 0,  0);
static const vQuantizer q5550s0(5, 5, 5, 0,  0);

static const vQuantizer q7778s1(7, 7, 7, 8, ~0);
static const vQuantizer q5556s1(5, 5, 5, 6, ~0);
static const vQuantizer q8888s1(8, 8, 8, 8, ~0);
static const vQuantizer q6666s1(6, 6, 6, 6, ~0);
static const vQuantizer q8880s1(8, 8, 8, 0, ~0);
static const vQuantizer q7770s1(7, 7, 7, 0, ~0);
static const vQuantizer q5550s1(5, 5, 5, 0, ~0);

#define vGetQuantizer(r, g, b, a)					\
	(((r) == 7) && ((a) == 8)                ? q7778s1 :		\
	(((r) == 5) && ((a) == 6)                ? q5556s1 :		\
	(((r) == 5) && ((a) == 0)                ? q5550s1 :		\
	(((r) == 8) && ((a) == 8)                ? q8888s1 :		\
	(((r) == 6) && ((a) == 6)                ? q6666s1 :		\
	(((r) == 8) && ((a) == 1)                ? q8880s1 :		\
	(((r) == 7) && ((a) == 1)                ? q7770s1 :		\
	(((r) == 5) && ((a) == 1)                ? q5550s1 :		\
	(vQuantizer&)*(vQuantizer*)nullptr))))))))

#define eGetQuantizer(r, g, b, a, e)					\
	(((r) == 7) && ((a) == 8) && ((e) == ~0) ? q7778s1 :		\
	(((r) == 5) && ((a) == 6) && ((e) == ~0) ? q5556s1 :		\
	(((r) == 5) && ((a) == 0) && ((e) == ~0) ? q5550s1 :		\
	(((r) == 8) && ((a) == 8) && ((e) ==  0) ? q8888s0 :		\
	(((r) == 6) && ((a) == 6) && ((e) ==  0) ? q6666s0 :		\
	(((r) == 8) && ((a) == 1) && ((e) ==  0) ? q8880s0 :		\
	(((r) == 7) && ((a) == 1) && ((e) ==  0) ? q7770s0 :		\
	(((r) == 5) && ((a) == 1) && ((e) ==  0) ? q5550s0 :		\
	(vQuantizer&)*(vQuantizer*)nullptr))))))))

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Vec4 (&colour)[1], Col4 (&field)[1][FIELDN], int bitset) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);
  assert((!eb && !sb && !(~bitset)) || ((eb || sb) && (~bitset)));

  // we can't just drop the eb/sb bits in fp-representation, we have to use the exact quantizer
  const vQuantizer &q = eGetQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb,
    eb + sb ? 0 : ~0
  );

  // pack into a single value
          field[0][COLORA] = q.QuantizeToInt(colour[0], bitset, 1 << 0);
          field[0][COLORA] = ShiftRight<eb + sb>(field[0][COLORA]);

  if (eb) field[0][UNIQUE] = (Col4(bitset) >> 0) & Col4(1);
  if (sb) field[0][SHARED] = (Col4(bitset) >> 0) & Col4(1);
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Vec4 (&colour)[1], Col4 (&field)[1][FIELDN]) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);

//const int rm = (1 << (rb + eb + sb)) - 1;
//const int gm = (1 << (gb + eb + sb)) - 1;
//const int bm = (1 << (bb + eb + sb)) - 1;
//const int am = (1 << (ab + eb + sb)) - 1;
  const int em = (1 << (     eb + sb)) - 1;
  const int sm = (1 << (          sb)) - 1;

  const vQuantizer &q = vGetQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb
  );

  Col4 const umask(em);
  Col4 const smask(sm);

  Col4 rexplicit[1];
  Col4 rcomplete[1];
  Col4 runique[1];
  Col4 rshared[1];

  Col4 _e[1];
  Col4 _s[1];

  // get the components in the complete range
  rcomplete[0] = q.QuantizeToInt(colour[0]);

  // get the components in the explicit range
  rexplicit[0] = ShiftRight<eb + sb>(rcomplete[0]);

  if (eb + sb) {
    /* stick g in a for skewing the rounding */
    if (!ab) {
      rcomplete[0] = Shuffle<1, 3>(rcomplete[0]);
    }
#if 0
    else {
      // TODO: this doesn't consider rotations!
      Col4 z0 = IsNotZero(rcomplete[0]).SplatA();
      Col4 o0 = IsOne    (rcomplete[0]).SplatA();

      /* if alpha is black the shared bit must be 0 */
      rcomplete[0] &= z0 | Col4(~em);

      /* if alpha is white the shared bit must be 1 */
      rcomplete[0] |= ShiftRight<32 - (eb + sb)>(o0);
    }
#endif

    if (eb) {
      // get the components in the unique range
      runique[0] = ShiftRight<sb>(rcomplete[0] & umask);

      _e[0] = HorizontalAddTiny(runique[0], Col4(2, 0, 0, 0)) >> 2;
    }

    if (sb) {
      // get the components in the shared range
      rshared[0] = (rcomplete[0] & smask);

      _s[0] = rshared[0];
    }
  }

  // pack into a single value
  field[0][COLORA] = rexplicit[0]; if (eb) field[0][UNIQUE] = _e[0]; if (sb) field[0][SHARED] = _s[0];
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Col4 (&fielda)[1][FIELDN], Col4 (&fieldb)[1][FIELDN]) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);

  if (eb + sb) {
    if (sb) {
      Col4 _s[1];

      _s[0] = (HorizontalAddTiny(fielda[0][SHARED], fieldb[0][SHARED]) + Col4(2 + 2)) >> 3;

      // pack into a single value
      fielda[0][SHARED] = fieldb[0][SHARED] = _s[0];
    }
  }
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Vec4 (&colour)[2], Col4 (&field)[2][FIELDN], int bitset) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);
  assert((!eb && !sb && !(~bitset)) || ((eb || sb) && (~bitset)));

  // we can't just drop the eb/sb bits in fp-representation, we have to use the exact quantizer
  const vQuantizer &q = eGetQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb,
    eb + sb ? 0 : ~0
  );

  // pack into a single value
          field[0][COLORA] = q.QuantizeToInt(colour[0], bitset, 1 << 0);
          field[1][COLORA] = q.QuantizeToInt(colour[1], bitset, 1 << 1);
          field[0][COLORA] = ShiftRight<eb + sb>(field[0][COLORA]);
          field[1][COLORA] = ShiftRight<eb + sb>(field[1][COLORA]);

  if (eb) field[0][UNIQUE] = (Col4(bitset) >> 0) & Col4(1);
  if (eb) field[1][UNIQUE] = (Col4(bitset) >> 1) & Col4(1);
  if (sb) field[0][SHARED] = (Col4(bitset) >> 0) & Col4(1);
  if (sb) field[1][SHARED] = (Col4(bitset) >> 1) & Col4(1);
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Vec4 (&colour)[2], Col4 (&field)[2][FIELDN]) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);

//const int rm = (1 << (rb + eb + sb)) - 1;
//const int gm = (1 << (gb + eb + sb)) - 1;
//const int bm = (1 << (bb + eb + sb)) - 1;
//const int am = (1 << (ab + eb + sb)) - 1;
  const int em = (1 << (     eb + sb)) - 1;
  const int sm = (1 << (          sb)) - 1;

  const vQuantizer &q = vGetQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb
  );

  Col4 const umask(em);
  Col4 const smask(sm);

  Col4 rcomplete[2];
  Col4 rexplicit[2];
  Col4 runique[2];
  Col4 rshared[2];

  Col4 _e[2];
  Col4 _s[2];

  // get the components in the complete range
  rcomplete[0] = q.QuantizeToInt(colour[0]);
  rcomplete[1] = q.QuantizeToInt(colour[1]);

  // get the components in the explicit range
	rexplicit[0] = ShiftRight<eb + sb>(rcomplete[0]);
	rexplicit[1] = ShiftRight<eb + sb>(rcomplete[1]);

  if (eb + sb) {
    /* stick g in a for skewing the rounding */
    if (!ab) {
      rcomplete[0] = Shuffle<1, 3>(rcomplete[0]);
      rcomplete[1] = Shuffle<1, 3>(rcomplete[1]);
    }
#if 0
    else {
      // TODO: this doesn't consider rotations!
      Col4 z0 = IsNotZero(rcomplete[0]).SplatA();
      Col4 z1 = IsNotZero(rcomplete[1]).SplatA();
      Col4 o0 = IsOne    (rcomplete[0]).SplatA();
      Col4 o1 = IsOne    (rcomplete[1]).SplatA();

      /* if alpha is black the shared bit must be 0 */
      rcomplete[0] &= z0 | Col4(~em);
      rcomplete[1] &= z1 | Col4(~em);

      /* if alpha is white the shared bit must be 1 */
			rcomplete[0] |= ShiftRight<32 - (eb + sb)>(o0);
			rcomplete[1] |= ShiftRight<32 - (eb + sb)>(o1);
    }
#endif

    if (eb) {
      // get the components in the unique range
			runique[0] = ShiftRight<sb>(rcomplete[0] & umask);
			runique[1] = ShiftRight<sb>(rcomplete[1] & umask);

      _e[0] = HorizontalAddTiny(runique[0], Col4(2, 0, 0, 0)) >> 2;
      _e[1] = HorizontalAddTiny(runique[1], Col4(2, 0, 0, 0)) >> 2;
    }

    if (sb) {
      // get the components in the shared range
      rshared[0] = (rcomplete[0] & smask);
      rshared[1] = (rcomplete[1] & smask);

      _s[0] = rshared[0];
      _s[1] = rshared[1];
    }
  }

  // pack into a single value
  field[0][COLORA] = rexplicit[0]; if (eb) field[0][UNIQUE] = _e[0]; if (sb) field[0][SHARED] = _s[0];
  field[1][COLORA] = rexplicit[1]; if (eb) field[1][UNIQUE] = _e[1]; if (sb) field[1][SHARED] = _s[1];
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Col4 (&fielda)[2][FIELDN], Col4 (&fieldb)[2][FIELDN]) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);

  if (eb + sb) {
    if (sb) {
      Col4 _s[2];

      _s[0] = (HorizontalAddTiny(fielda[0][SHARED], fieldb[0][SHARED]) + Col4(2 + 2)) >> 3;
      _s[1] = (HorizontalAddTiny(fielda[1][SHARED], fieldb[1][SHARED]) + Col4(2 + 2)) >> 3;

      // pack into a single value
      fielda[0][SHARED] = fieldb[0][SHARED] = _s[0];
      fielda[1][SHARED] = fieldb[1][SHARED] = _s[1];
    }
  }
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Vec4 (&colour)[3], Col4 (&field)[3][FIELDN], int bitset) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);
  assert((!eb && !sb && !(~bitset)) || ((eb || sb) && (~bitset)));

  // we can't just drop the eb/sb bits in fp-representation, we have to use the exact quantizer
  const vQuantizer &q = eGetQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb,
    eb + sb ? 0 : ~0
  );

  // pack into a single value
          field[0][COLORA] = q.QuantizeToInt(colour[0], bitset, 1 << 0);
          field[1][COLORA] = q.QuantizeToInt(colour[1], bitset, 1 << 1);
          field[2][COLORA] = q.QuantizeToInt(colour[2], bitset, 1 << 2);
          field[0][COLORA] = ShiftRight<eb + sb>(field[0][COLORA]);
          field[1][COLORA] = ShiftRight<eb + sb>(field[1][COLORA]);
          field[2][COLORA] = ShiftRight<eb + sb>(field[2][COLORA]);

  if (eb) field[0][UNIQUE] = (Col4(bitset) >> 0) & Col4(1);
  if (eb) field[1][UNIQUE] = (Col4(bitset) >> 1) & Col4(1);
  if (eb) field[2][UNIQUE] = (Col4(bitset) >> 2) & Col4(1);
  if (sb) field[0][SHARED] = (Col4(bitset) >> 0) & Col4(1);
  if (sb) field[1][SHARED] = (Col4(bitset) >> 1) & Col4(1);
  if (sb) field[2][SHARED] = (Col4(bitset) >> 2) & Col4(1);
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Vec4 (&colour)[3], Col4 (&field)[3][FIELDN]) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);

//const int rm = (1 << (rb + eb + sb)) - 1;
//const int gm = (1 << (gb + eb + sb)) - 1;
//const int bm = (1 << (bb + eb + sb)) - 1;
//const int am = (1 << (ab + eb + sb)) - 1;
  const int em = (1 << (     eb + sb)) - 1;
  const int sm = (1 << (          sb)) - 1;

  const vQuantizer &q = vGetQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb
  );

  Col4 const umask(em);
  Col4 const smask(sm);

  Col4 rcomplete[3];
  Col4 rexplicit[3];
  Col4 runique[3];
  Col4 rshared[3];

  Col4 _e[3];
  Col4 _s[3];

  // get the components in the complete range
  rcomplete[0] = q.QuantizeToInt(colour[0]);
  rcomplete[1] = q.QuantizeToInt(colour[1]);
  rcomplete[2] = q.QuantizeToInt(colour[2]);

  // get the components in the explicit range
  rexplicit[0] = ShiftRight<eb + sb>(rcomplete[0]);
  rexplicit[1] = ShiftRight<eb + sb>(rcomplete[1]);
  rexplicit[2] = ShiftRight<eb + sb>(rcomplete[2]);

  if (eb + sb) {
    /* stick g in a for skewing the rounding */
    if (!ab) {
      rcomplete[0] = Shuffle<1, 3>(rcomplete[0]);
      rcomplete[1] = Shuffle<1, 3>(rcomplete[1]);
      rcomplete[2] = Shuffle<1, 3>(rcomplete[2]);
    }
#if 0
    else {
      // TODO: this doesn't consider rotations!
      Col4 z0 = IsNotZero(rcomplete[0]).SplatA();
      Col4 z1 = IsNotZero(rcomplete[1]).SplatA();
      Col4 z2 = IsNotZero(rcomplete[2]).SplatA();
      Col4 o0 = IsOne    (rcomplete[0]).SplatA();
      Col4 o1 = IsOne    (rcomplete[1]).SplatA();
      Col4 o2 = IsOne    (rcomplete[2]).SplatA();

      /* if alpha is black the shared bit must be 0 */
      rcomplete[0] &= z0 | Col4(~em);
      rcomplete[1] &= z1 | Col4(~em);
      rcomplete[2] &= z2 | Col4(~em);

      /* if alpha is white the shared bit must be 1 */
			rcomplete[0] |= ShiftRight<32 - (eb + sb)>(o0);
			rcomplete[1] |= ShiftRight<32 - (eb + sb)>(o1);
			rcomplete[2] |= ShiftRight<32 - (eb + sb)>(o2);
    }
#endif

    if (eb) {
      // get the components in the unique range
			runique[0] = ShiftRight<sb>(rcomplete[0] & umask);
			runique[1] = ShiftRight<sb>(rcomplete[1] & umask);
			runique[2] = ShiftRight<sb>(rcomplete[2] & umask);

      _e[0] = HorizontalAddTiny(runique[0], Col4(2, 0, 0, 0)) >> 2;
      _e[1] = HorizontalAddTiny(runique[1], Col4(2, 0, 0, 0)) >> 2;
      _e[2] = HorizontalAddTiny(runique[2], Col4(2, 0, 0, 0)) >> 2;
    }

    if (sb) {
      // get the components in the shared range
      rshared[0] = (rcomplete[0] & smask);
      rshared[1] = (rcomplete[1] & smask);
      rshared[2] = (rcomplete[2] & smask);

      _s[0] = rshared[0];
      _s[1] = rshared[1];
      _s[2] = rshared[2];
    }
  }

  // pack into a single value
  field[0][COLORA] = rexplicit[0]; if (eb) field[0][UNIQUE] = _e[0]; if (sb) field[0][SHARED] = _s[0];
  field[1][COLORA] = rexplicit[1]; if (eb) field[1][UNIQUE] = _e[1]; if (sb) field[1][SHARED] = _s[1];
  field[2][COLORA] = rexplicit[2]; if (eb) field[2][UNIQUE] = _e[2]; if (sb) field[2][SHARED] = _s[2];
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Col4 (&fielda)[3][FIELDN], Col4 (&fieldb)[3][FIELDN]) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);

  if (eb + sb) {
    if (sb) {
      Col4 _s[3];

      _s[0] = (HorizontalAddTiny(fielda[0][SHARED], fieldb[0][SHARED]) + Col4(2 + 2)) >> 3;
      _s[1] = (HorizontalAddTiny(fielda[1][SHARED], fieldb[1][SHARED]) + Col4(2 + 2)) >> 3;
      _s[2] = (HorizontalAddTiny(fielda[2][SHARED], fieldb[2][SHARED]) + Col4(2 + 2)) >> 3;

      // pack into a single value
      fielda[0][SHARED] = fieldb[0][SHARED] = _s[0];
      fielda[1][SHARED] = fieldb[1][SHARED] = _s[1];
      fielda[2][SHARED] = fieldb[2][SHARED] = _s[2];
    }
  }
}

/* -----------------------------------------------------------------------------
 */
template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg UnpackFrom(Col4 (&field)[1][FIELDN]) ccr_restricted {
  /* not both yet */
  assert(!eb || !sb);

  static const Col4 scale = Col4(
    // rb = gb = bb = 3, ab = 4, eb = 1, sb = 0 -> 4/5
    // 1 * (rb = gb = bb) + eb + sb ->  4 => (1 << 4) + 1 => 17
    // 1 * (     ab     ) + eb + sb ->  5 => (1 << 5) + 1 => 33
    // 2 * (rb = gb = bb) + eb + sb ->  8 => 16 -  8 = 8 => 17 << 8 => 4352
    // 2 * (     ab     ) + eb + sb -> 10 => 16 - 10 = 6 => 33 << 6 => 2112
    // 00000000 00001111, 00000000 00001111, 00000000 00001111, 00000000 00011111
    // 00000000 11111111, 00000000 11111111, 00000000 11111111, 00000011 11111111
    // 11111111 00000000, 11111111 00000000, 11111111 00000000, 11111111 11000000

    ((1 << (rb + eb + sb)) + 1) << (16 - 2 * (rb + eb + sb)),
    ((1 << (gb + eb + sb)) + 1) << (16 - 2 * (gb + eb + sb)),
    ((1 << (bb + eb + sb)) + 1) << (16 - 2 * (bb + eb + sb)),
    ((1 << (ab + eb + sb)) + 1) << (16 - 2 * (ab + eb + sb))
  );

  // insert the 1 unique bit
  if (eb) {
    field[0][COLORA] <<= eb;

    field[0][COLORA] |= field[0][UNIQUE];
  }

  // insert the 1 shared bit
  if (sb) {
    field[0][COLORA] <<= sb;

    field[0][COLORA] |= field[0][SHARED];
  }

  if ((((rb + eb + sb) != 8)) ||
      (((gb + eb + sb) != 8)) ||
      (((bb + eb + sb) != 8)) ||
      (((ab + eb + sb) != 8) && ab)) {
    // extend X+Y+Z bits to 8 bits
    field[0][COLORA] = (field[0][COLORA] * scale) >> 8;
  }

  // set A to opaque
  if (!ab) {
    field[0][COLORA] = KillA(field[0][COLORA]);
  }
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg UnpackFrom(Col4 (&field)[2][FIELDN]) ccr_restricted {
  /* not both yet */
  assert(!eb || !sb);

  static const Col4 scale = Col4(
    // rb = gb = bb = 3, ab = 4, eb = 1, sb = 0 -> 4/5
    // 1 * (rb = gb = bb) + eb + sb ->  4 => (1 << 4) + 1 => 17
    // 1 * (     ab     ) + eb + sb ->  5 => (1 << 5) + 1 => 33
    // 2 * (rb = gb = bb) + eb + sb ->  8 => 16 -  8 = 8 => 17 << 8 => 4352
    // 2 * (     ab     ) + eb + sb -> 10 => 16 - 10 = 6 => 33 << 6 => 2112
    // 00000000 00001111, 00000000 00001111, 00000000 00001111, 00000000 00011111
    // 00000000 11111111, 00000000 11111111, 00000000 11111111, 00000011 11111111
    // 11111111 00000000, 11111111 00000000, 11111111 00000000, 11111111 11000000

    ((1 << (rb + eb + sb)) + 1) << (16 - 2 * (rb + eb + sb)),
    ((1 << (gb + eb + sb)) + 1) << (16 - 2 * (gb + eb + sb)),
    ((1 << (bb + eb + sb)) + 1) << (16 - 2 * (bb + eb + sb)),
    ((1 << (ab + eb + sb)) + 1) << (16 - 2 * (ab + eb + sb))
  );

  // insert the 1 unique bit
  if (eb) {
    field[1][COLORA] <<= eb;
    field[0][COLORA] <<= eb;

    field[1][COLORA] |= field[1][UNIQUE];
    field[0][COLORA] |= field[0][UNIQUE];
  }

  // insert the 1 shared bit
  if (sb) {
    field[1][COLORA] <<= sb;
    field[0][COLORA] <<= sb;

    field[1][COLORA] |= field[1][SHARED];
    field[0][COLORA] |= field[0][SHARED];
  }

  if ((((rb + eb + sb) != 8)) ||
      (((gb + eb + sb) != 8)) ||
      (((bb + eb + sb) != 8)) ||
      (((ab + eb + sb) != 8) && ab)) {
    // extend X+Y+Z bits to 8 bits
    field[1][COLORA] = (field[1][COLORA] * scale) >> 8;
    field[0][COLORA] = (field[0][COLORA] * scale) >> 8;
  }

  // set A to opaque
  if (!ab) {
    field[1][COLORA] = KillA(field[1][COLORA]);
    field[0][COLORA] = KillA(field[0][COLORA]);
  }
}

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg UnpackFrom(Col4 (&field)[3][FIELDN]) ccr_restricted {
  /* not both yet */
  assert(!eb || !sb);

  static const Col4 scale = Col4(
    // rb = gb = bb = 3, ab = 4, eb = 1, sb = 0 -> 4/5
    // 1 * (rb = gb = bb) + eb + sb ->  4 => (1 << 4) + 1 => 17
    // 1 * (     ab     ) + eb + sb ->  5 => (1 << 5) + 1 => 33
    // 2 * (rb = gb = bb) + eb + sb ->  8 => 16 -  8 = 8 => 17 << 8 => 4352
    // 2 * (     ab     ) + eb + sb -> 10 => 16 - 10 = 6 => 33 << 6 => 2112
    // 00000000 00001111, 00000000 00001111, 00000000 00001111, 00000000 00011111
    // 00000000 11111111, 00000000 11111111, 00000000 11111111, 00000011 11111111
    // 11111111 00000000, 11111111 00000000, 11111111 00000000, 11111111 11000000

    ((1 << (rb + eb + sb)) + 1) << (16 - 2 * (rb + eb + sb)),
    ((1 << (gb + eb + sb)) + 1) << (16 - 2 * (gb + eb + sb)),
    ((1 << (bb + eb + sb)) + 1) << (16 - 2 * (bb + eb + sb)),
    ((1 << (ab + eb + sb)) + 1) << (16 - 2 * (ab + eb + sb))
  );

  // insert the 1 unique bit
  if (eb) {
    field[2][COLORA] <<= eb;
    field[1][COLORA] <<= eb;
    field[0][COLORA] <<= eb;

    field[2][COLORA] |= field[2][UNIQUE];
    field[1][COLORA] |= field[1][UNIQUE];
    field[0][COLORA] |= field[0][UNIQUE];
  }

  // insert the 1 shared bit
  if (sb) {
    field[2][COLORA] <<= sb;
    field[1][COLORA] <<= sb;
    field[0][COLORA] <<= sb;

    field[2][COLORA] |= field[2][SHARED];
    field[1][COLORA] |= field[1][SHARED];
    field[0][COLORA] |= field[0][SHARED];
  }

  if ((((rb + eb + sb) != 8)) ||
      (((gb + eb + sb) != 8)) ||
      (((bb + eb + sb) != 8)) ||
      (((ab + eb + sb) != 8) && ab)) {
    // extend X+Y+Z bits to 8 bits
    field[2][COLORA] = (field[2][COLORA] * scale) >> 8;
    field[1][COLORA] = (field[1][COLORA] * scale) >> 8;
    field[0][COLORA] = (field[0][COLORA] * scale) >> 8;
  }

  // set A to opaque
  if (!ab) {
    field[2][COLORA] = KillA(field[2][COLORA]);
    field[1][COLORA] = KillA(field[1][COLORA]);
    field[0][COLORA] = KillA(field[0][COLORA]);
  }
}

/* *****************************************************************************
 * http://embeddedgurus.com/stack-overflow/2009/06/division-of-integers-by-constants/
 * Divide by 3:  (((uint32_t)A * (uint32_t)0xAAAB) >> 16) >> 1
 * Divide by 5:  (((uint32_t)A * (uint32_t)0xCCCD) >> 16) >> 2
 * Divide by 7: ((((uint32_t)A * (uint32_t)0x2493) >> 16) + A) >> 1) >> 2
 */
static doinline void passreg Codebook3or4(u8 (&codes)[4*4], bool bw) ccr_restricted
{
  // generate the midpoints
  for (int i = 0; i < 3; ++i) {
    const int c = codes[0 + i];
    const int d = codes[4 + i];

    if (bw) {
      codes[ 8 + i] = (u8)(((1 * c + 1 * d)         ) >>  1);
      codes[12 + i] = 0;
    }
    else {
      codes[ 8 + i] = (u8)(((2 * c + 1 * d) * 0xAAAB) >> 17);
      codes[12 + i] = (u8)(((1 * c + 2 * d) * 0xAAAB) >> 17);
    }
  }

  // fill in alpha for the intermediate values
  codes[ 8 + 3] =          255;
  codes[12 + 3] = bw ? 0 : 255;
}

template<const int prc>
static doinline void passreg Codebook6or8(u8 (&codes)[8*1], bool bw) ccr_restricted
{
  // generate the midpoints
  for (int i = 0; i < 1; ++i) {
    const int c = codes[0 + i];
    const int d = codes[1 + i];
    int cd;
    
    if (bw) {
      cd = (4 * c + 1 * d); codes[2 + i] = (u8)((cd * 0xCCCD) >> 18);
      cd = (3 * c + 2 * d); codes[3 + i] = (u8)((cd * 0xCCCD) >> 18);
      cd = (2 * c + 3 * d); codes[4 + i] = (u8)((cd * 0xCCCD) >> 18);
      cd = (1 * c + 4 * d); codes[5 + i] = (u8)((cd * 0xCCCD) >> 18);

      codes[6 + i] = (u8)0;
      codes[7 + i] = (u8)255;
    }
    else {
      cd = (6 * c + 1 * d); codes[2 + i] = (u8)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (5 * c + 2 * d); codes[3 + i] = (u8)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (4 * c + 3 * d); codes[4 + i] = (u8)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (3 * c + 4 * d); codes[5 + i] = (u8)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (2 * c + 5 * d); codes[6 + i] = (u8)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (1 * c + 6 * d); codes[7 + i] = (u8)((((cd * 0x2493) >> 16) + cd) >> 3);
    }
  }
}

template<const int prc>
static doinline void passreg Codebook6or8(s8 (&codes)[8*1], bool bw) ccr_restricted
{
  // generate the midpoints
  for (int i = 0; i < 1; ++i) {
    const int c = codes[0 + i];
    const int d = codes[1 + i];
    int cd;
    
    if (bw) {
      cd = (4 * c + 1 * d); codes[2 + i] = (s8)((cd * 0x3334) >> 16) + (cd < 0);
      cd = (3 * c + 2 * d); codes[3 + i] = (s8)((cd * 0x3334) >> 16) + (cd < 0);
      cd = (2 * c + 3 * d); codes[4 + i] = (s8)((cd * 0x3334) >> 16) + (cd < 0);
      cd = (1 * c + 4 * d); codes[5 + i] = (s8)((cd * 0x3334) >> 16) + (cd < 0);

      codes[6 + i] = (s8)-127;
      codes[7 + i] = (s8) 127;

      assert(s8(codes[2]) == (((s8(4) * s8(codes[0])) + (s8(1) * s8(codes[1]))) / 5));
      assert(s8(codes[3]) == (((s8(3) * s8(codes[0])) + (s8(2) * s8(codes[1]))) / 5));
      assert(s8(codes[4]) == (((s8(2) * s8(codes[0])) + (s8(3) * s8(codes[1]))) / 5));
      assert(s8(codes[5]) == (((s8(1) * s8(codes[0])) + (s8(4) * s8(codes[1]))) / 5));
      assert(s8(codes[6]) == (-127));
      assert(s8(codes[7]) == ( 127));
    }
    else {
      cd = (6 * c + 1 * d); codes[2 + i] = (s8)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (5 * c + 2 * d); codes[3 + i] = (s8)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (4 * c + 3 * d); codes[4 + i] = (s8)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (3 * c + 4 * d); codes[5 + i] = (s8)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (2 * c + 5 * d); codes[6 + i] = (s8)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (1 * c + 6 * d); codes[7 + i] = (s8)((cd * 0x4925) >> 17) + (cd < 0);
  
      assert(s8(codes[2]) == (((s8(6) * s8(codes[0])) + (s8(1) * s8(codes[1]))) / 7));
      assert(s8(codes[3]) == (((s8(5) * s8(codes[0])) + (s8(2) * s8(codes[1]))) / 7));
      assert(s8(codes[4]) == (((s8(4) * s8(codes[0])) + (s8(3) * s8(codes[1]))) / 7));
      assert(s8(codes[5]) == (((s8(3) * s8(codes[0])) + (s8(4) * s8(codes[1]))) / 7));
      assert(s8(codes[6]) == (((s8(2) * s8(codes[0])) + (s8(5) * s8(codes[1]))) / 7));
      assert(s8(codes[7]) == (((s8(1) * s8(codes[0])) + (s8(6) * s8(codes[1]))) / 7));
    }
  }
}

template<const int prc>
static doinline void passreg Codebook6or8(u16 (&codes)[8*1], bool bw) ccr_restricted
{
  // generate the midpoints
  for (int i = 0; i < 1; ++i) {
    const int c = codes[0 + i];
    const int d = codes[1 + i];
    int cd;
    
    if (bw) {
      cd = (4 * c + 1 * d); codes[2 + i] = (u16)((cd * 0xCCCD) >> 18);
      cd = (3 * c + 2 * d); codes[3 + i] = (u16)((cd * 0xCCCD) >> 18);
      cd = (2 * c + 3 * d); codes[4 + i] = (u16)((cd * 0xCCCD) >> 18);
      cd = (1 * c + 4 * d); codes[5 + i] = (u16)((cd * 0xCCCD) >> 18);

      codes[6 + i] = (u16)  0 << prc;
      codes[7 + i] = (u16)255 << prc;
    }
    else {
      cd = (6 * c + 1 * d); codes[2 + i] = (u16)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (5 * c + 2 * d); codes[3 + i] = (u16)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (4 * c + 3 * d); codes[4 + i] = (u16)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (3 * c + 4 * d); codes[5 + i] = (u16)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (2 * c + 5 * d); codes[6 + i] = (u16)((((cd * 0x2493) >> 16) + cd) >> 3);
      cd = (1 * c + 6 * d); codes[7 + i] = (u16)((((cd * 0x2493) >> 16) + cd) >> 3);
    }
  }
}

template<const int prc>
static doinline void passreg Codebook6or8(s16 (&codes)[8*1], bool bw) ccr_restricted
{
  // generate the midpoints
  for (int i = 0; i < 1; ++i) {
    const int c = codes[0 + i];
    const int d = codes[1 + i];
    int cd;
    
    if (bw) {
      cd = (4 * c + 1 * d); codes[2 + i] = (s16)((cd * 0x3334) >> 16) + (cd < 0);
      cd = (3 * c + 2 * d); codes[3 + i] = (s16)((cd * 0x3334) >> 16) + (cd < 0);
      cd = (2 * c + 3 * d); codes[4 + i] = (s16)((cd * 0x3334) >> 16) + (cd < 0);
      cd = (1 * c + 4 * d); codes[5 + i] = (s16)((cd * 0x3334) >> 16) + (cd < 0);

      codes[6 + i] = (s16)-127 << prc;
      codes[7 + i] = (s16) 127 << prc;

      assert(s16(codes[2]) == (((s16(4) * s16(codes[0])) + (s16(1) * s16(codes[1]))) / 5));
      assert(s16(codes[3]) == (((s16(3) * s16(codes[0])) + (s16(2) * s16(codes[1]))) / 5));
      assert(s16(codes[4]) == (((s16(2) * s16(codes[0])) + (s16(3) * s16(codes[1]))) / 5));
      assert(s16(codes[5]) == (((s16(1) * s16(codes[0])) + (s16(4) * s16(codes[1]))) / 5));
      assert(s16(codes[6]) == (-127 << prc));
      assert(s16(codes[7]) == ( 127 << prc));
    }
    else {
      cd = (6 * c + 1 * d); codes[2 + i] = (s16)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (5 * c + 2 * d); codes[3 + i] = (s16)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (4 * c + 3 * d); codes[4 + i] = (s16)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (3 * c + 4 * d); codes[5 + i] = (s16)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (2 * c + 5 * d); codes[6 + i] = (s16)((cd * 0x4925) >> 17) + (cd < 0);
      cd = (1 * c + 6 * d); codes[7 + i] = (s16)((cd * 0x4925) >> 17) + (cd < 0);
  
      assert(s16(codes[2]) == (((s16(6) * s16(codes[0])) + (s16(1) * s16(codes[1]))) / 7));
      assert(s16(codes[3]) == (((s16(5) * s16(codes[0])) + (s16(2) * s16(codes[1]))) / 7));
      assert(s16(codes[4]) == (((s16(4) * s16(codes[0])) + (s16(3) * s16(codes[1]))) / 7));
      assert(s16(codes[5]) == (((s16(3) * s16(codes[0])) + (s16(4) * s16(codes[1]))) / 7));
      assert(s16(codes[6]) == (((s16(2) * s16(codes[0])) + (s16(5) * s16(codes[1]))) / 7));
      assert(s16(codes[7]) == (((s16(1) * s16(codes[0])) + (s16(6) * s16(codes[1]))) / 7));
    }
  }
}

static int passreg CodebookP(u8 *codes, int bits) ccr_restricted
{
  // generate the midpoints
  for (int m = 0; m < 4; ++m) {
    const int j = (1 << bits) - 1;

    const int c = codes[0 * 4 + m];
    const int d = codes[j * 4 + m];

    // the quantizer is not equi-distant, but it is symmetric
    for (int i = 1; i < j; i++) {
      int s = (weights_u16[bits][j - i]) * c;
      int e = (weights_u16[bits][i + 0]) * d;

      codes[i * 4 + m] = (u8)((s + e + 32) >> 6);
    }
  }

  return (1 << bits);
}

/* -----------------------------------------------------------------------------
 */
static doinline void passreg Codebook3(Vec3 (&codes)[3], Vec3::Arg start, Vec3::Arg end) ccr_restricted
{
  codes[0] = start;
  codes[1] = end;
  codes[2] = (0.5f * start) + (0.5f * end);
//codes[3] = 0;
}

static doinline void passreg Codebook4(Vec3 (&codes)[4], Vec3::Arg start, Vec3::Arg end) ccr_restricted
{
  codes[0] = start;
  codes[1] = end;
  codes[2] = (2.0f / 3.0f) * start + (1.0f / 3.0f) * end;
  codes[3] = (1.0f / 3.0f) * start + (2.0f / 3.0f) * end;
}

static doinline void passreg Codebook6(Vec4 (&codes)[8], Vec4::Arg start, Vec4::Arg end) ccr_restricted
{
  codes[0] = start;
  codes[1] = end;
  codes[2] = (4.0f / 5.0f) * start + (1.0f / 5.0f) * end;
  codes[3] = (3.0f / 5.0f) * start + (2.0f / 5.0f) * end;
  codes[4] = (2.0f / 5.0f) * start + (3.0f / 5.0f) * end;
  codes[5] = (1.0f / 5.0f) * start + (4.0f / 5.0f) * end;
  codes[6] = Vec4(  0.0f);
  codes[7] = Vec4(255.0f);
}

static doinline void passreg Codebook8(Vec4 (&codes)[8], Vec4::Arg start, Vec4::Arg end) ccr_restricted
{
  codes[0] = start;
  codes[1] = end;
  codes[2] = (6.0f / 7.0f) * start + (1.0f / 7.0f) * end;
  codes[3] = (5.0f / 7.0f) * start + (2.0f / 7.0f) * end;
  codes[4] = (4.0f / 7.0f) * start + (3.0f / 7.0f) * end;
  codes[5] = (3.0f / 7.0f) * start + (4.0f / 7.0f) * end;
  codes[4] = (2.0f / 7.0f) * start + (5.0f / 7.0f) * end;
  codes[5] = (1.0f / 7.0f) * start + (6.0f / 7.0f) * end;
}

/* -----------------------------------------------------------------------------
 */
#ifdef	FEATURE_NORMALFIT_UNITGUARANTEE
#define DISARM	true
#else
#define DISARM	false
#endif

static doinline void passreg Codebook3n(Vec3 (&codes)[3], Vec3::Arg start, Vec3::Arg end) ccr_restricted
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  Codebook3(codes, start, end);

  codes[0] = Normalize(scale * (offset + codes[0]));
  codes[1] = Normalize(scale * (offset + codes[1]));
  codes[2] = Normalize(scale * (offset + codes[2]));
}

static doinline void passreg Codebook4n(Vec3 (&codes)[4], Vec3::Arg start, Vec3::Arg end) ccr_restricted
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  Codebook4(codes, start, end);

  codes[0] = Normalize(scale * (offset + codes[0]));
  codes[1] = Normalize(scale * (offset + codes[1]));
  codes[2] = Normalize(scale * (offset + codes[2]));
  codes[3] = Normalize(scale * (offset + codes[3]));
}

static doinline void passreg Codebook3nc(Vec3 (&codes)[3], Vec3::Arg start, Vec3::Arg end) ccr_restricted
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  Codebook3(codes, start, end);
  
  codes[0] = Complement<DISARM>(scale * (offset + codes[0]));
  codes[1] = Complement<DISARM>(scale * (offset + codes[1]));
  codes[2] = Complement<DISARM>(scale * (offset + codes[2]));
}

static doinline void passreg Codebook4nc(Vec3 (&codes)[4], Vec3::Arg start, Vec3::Arg end) ccr_restricted
{
  const Vec3 scale  = Vec3( 1.0f / 0.5f);
  const Vec3 offset = Vec3(-1.0f * 0.5f);

  Codebook4(codes, start, end);

  codes[0] = Complement<DISARM>(scale * (offset + codes[0]));
  codes[1] = Complement<DISARM>(scale * (offset + codes[1]));
  codes[2] = Complement<DISARM>(scale * (offset + codes[2]));
  codes[3] = Complement<DISARM>(scale * (offset + codes[3]));
}

/* -----------------------------------------------------------------------------
 */
#define CODEBOOKLQ_PRECISIONBITS	(0)
#define CODEBOOKHQ_PRECISIONBITS	(5)
#define CODEBOOKLQ_MULTIPLIER		(1 << CODEBOOKLQ_PRECISIONBITS)
#define CODEBOOKHQ_MULTIPLIER		(1 << CODEBOOKHQ_PRECISIONBITS)

template<const int min, const int max, const int pb>
static doinline void passreg Codebook6(Col8 &codes, Col8::Arg start, Col8::Arg end) ccr_restricted
{
  // max unsigned: (5 * 255) << 5 = 40800 / 0x9F60 fits unsigned short
  // max   signed: (5 * 127) << 5 = 20320 / 0x4F60 fits   signed short
  const Col8 smul = Col8(0x05 << pb, 0x00 << pb, 0x04 << pb, 0x03 << pb, 0x02 << pb, 0x01 << pb, 0x00 << pb, 0x00 << pb);
  const Col8 emul = Col8(0x00 << pb, 0x05 << pb, 0x01 << pb, 0x02 << pb, 0x03 << pb, 0x04 << pb, 0x00 << pb, 0x00 << pb);
  const Col8 mask = Col8(0x00 << pb, 0x00 << pb, 0x00 << pb, 0x00 << pb, 0x00 << pb, 0x00 << pb, min  << pb, max  << pb);

  // range [0,2*5*255]
  Col8 ipol = (smul * start) + (emul * end);
  
  if (min >= 0)
    // max unsigned:  0x9F60 * 0xCCCD = 0x7F801FE0 = 0x7F80		255 << 7
    codes = ((ipol * 0xCCCDU) >> 2U) + mask;
  else
    // max   signed:  0x4F60 * 0x3334 = 0x0FE03F80 = 0x0FE0		127 << 5
    codes = ((ipol * 0x6667 ) >> 1 ) + mask - CompareAllLessThan(ipol, Col8(0,0,0,0,0,0, 0x8000, 0x8000));

  assert(s16(codes[0]) == (((s16(smul[0]) * s16(start[0])) + (s16(emul[0]) * s16(end[0]))) / 5 + s16(mask[0])));
  assert(s16(codes[1]) == (((s16(smul[1]) * s16(start[1])) + (s16(emul[1]) * s16(end[1]))) / 5 + s16(mask[1])));
  assert(s16(codes[2]) == (((s16(smul[2]) * s16(start[2])) + (s16(emul[2]) * s16(end[2]))) / 5 + s16(mask[2])));
  assert(s16(codes[3]) == (((s16(smul[3]) * s16(start[3])) + (s16(emul[3]) * s16(end[3]))) / 5 + s16(mask[3])));
  assert(s16(codes[4]) == (((s16(smul[4]) * s16(start[4])) + (s16(emul[4]) * s16(end[4]))) / 5 + s16(mask[4])));
  assert(s16(codes[5]) == (((s16(smul[5]) * s16(start[5])) + (s16(emul[5]) * s16(end[5]))) / 5 + s16(mask[5])));
  assert(s16(codes[6]) ==                                                                        s16(mask[6]) );
  assert(s16(codes[7]) ==                                                                        s16(mask[7]) );
}

template<const int min, const int max, const int pb>
static doinline void passreg Codebook8(Col8 &codes, Col8::Arg start, Col8::Arg end) ccr_restricted
{
  // max unsigned: (7 * 255) << 5 = 57120 / 0xDF20 fits unsigned short
  // max   signed: (7 * 127) << 5 = 28448 / 0x6F20 fits   signed short
  const Col8 smul = Col8(0x07 << pb, 0x00 << pb, 0x06 << pb, 0x05 << pb, 0x04 << pb, 0x03 << pb, 0x02 << pb, 0x01 << pb);
  const Col8 emul = Col8(0x00 << pb, 0x07 << pb, 0x01 << pb, 0x02 << pb, 0x03 << pb, 0x04 << pb, 0x05 << pb, 0x06 << pb);

  // range [0,2*7*255]
  Col8 ipol = (smul * start) + (emul * end);
  
  if (min >= 0)
    // max unsigned:  0xDF20 * 0x2493 = 0x1FE09F60 + 0xDF20 = FF00	255 << 8
    codes = (((ipol * 0x2493U) + ipol) >> 3U);
  else
    // max unsigned:  0x6F20 * 0x4925 = 0x1FC02FA0          = 1FC0	127 << 6
    codes = (((ipol * 0x4925 )       ) >> 1 ) - CompareAllLessThan(ipol, Col8(0,0,0,0,0,0,0,0));
  
  assert(s16(codes[0]) == (((s16(smul[0]) * s16(start[0])) + (s16(emul[0]) * s16(end[0]))) / 7));
  assert(s16(codes[1]) == (((s16(smul[1]) * s16(start[1])) + (s16(emul[1]) * s16(end[1]))) / 7));
  assert(s16(codes[2]) == (((s16(smul[2]) * s16(start[2])) + (s16(emul[2]) * s16(end[2]))) / 7));
  assert(s16(codes[3]) == (((s16(smul[3]) * s16(start[3])) + (s16(emul[3]) * s16(end[3]))) / 7));
  assert(s16(codes[4]) == (((s16(smul[4]) * s16(start[4])) + (s16(emul[4]) * s16(end[4]))) / 7));
  assert(s16(codes[5]) == (((s16(smul[5]) * s16(start[5])) + (s16(emul[5]) * s16(end[5]))) / 7));
  assert(s16(codes[6]) == (((s16(smul[6]) * s16(start[6])) + (s16(emul[6]) * s16(end[6]))) / 7));
  assert(s16(codes[7]) == (((s16(smul[7]) * s16(start[7])) + (s16(emul[7]) * s16(end[7]))) / 7));
}

/* -----------------------------------------------------------------------------
 */
static int passreg CodebookP(Vec4 *codes, int bits, Vec4::Arg start, Vec4::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  codes[0] = start;
  codes[j] = end;

  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    const Vec4 s = weights_V4[bits][j - i] * start;
    const Vec4 e = weights_V4[bits][i + 0] * end;

    codes[i] = s + e;
  }

  return (1 << bits);
}

static int passreg CodebookP(Col4 *codes, int bits, Col4::Arg start, Col4::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  codes[0] = start;
  codes[j] = end;

  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    const Col4 s = (weights_C4[bits][j - i]) * start;
    const Col4 e = (weights_C4[bits][i + 0]) * end;

    codes[i] = (s + e + Col4(32)) >> 6;
  }

  return (1 << bits);
}

template<const int bits>
static int passreg CodebookP(unsigned int *codes, Col4::Arg start, Col4::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  PackBytes(start, codes[0]);
  PackBytes(end  , codes[j]);

  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    const Col4 s = (weights_C4[bits][j - i]) * start;
    const Col4 e = (weights_C4[bits][i + 0]) * end;

    PackBytes((s + e + Col4(32)) >> 6, codes[i]);
  }

  return (1 << bits);
}

static int passreg CodebookP(Col3 *codes, int bits, Col3::Arg start, Col3::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  codes[0] = start;
  codes[j] = end;

  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    Col3 s = Mul16x16u(weights_C4[bits][j - i].GetCol3(), start);
    Col3 e = Mul16x16u(weights_C4[bits][i + 0].GetCol3(), end);

    codes[i] = (s + e + Col3(32)) >> 6;
  }

  return (1 << bits);
}

template<const int bits>
static int passreg CodebookP(unsigned__int64 *codes, Col3::Arg start, Col3::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  PackWords(Col4(start), codes[0]);
  PackWords(Col4(end  ), codes[j]);

  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    Col3 s = Mul16x16u(weights_C4[bits][j - i].GetCol3(), start);
    Col3 e = Mul16x16u(weights_C4[bits][i + 0].GetCol3(), end  );

    PackWords(Col4((s + e + Col3(32)) >> 6), codes[i]);
  }

  return (1 << bits);
}
/* -----------------------------------------------------------------------------
 */
static int passreg CodebookPn(Vec4 *codes, int bits, Vec4::Arg start, Vec4::Arg end) ccr_restricted
{
  const Vec4 scale  = Vec4( 1.0f / 0.5f);
  const Vec4 offset = Vec4(-1.0f * 0.5f);

  CodebookP(codes, bits, start, end);
  
  const int j = (1 << bits) - 1;
  for (int i = 0; i <= j; i++)
    codes[i] = TransferW(Normalize(KillW(scale * (offset + codes[i]))), codes[i]);

  return (1 << bits);
}
#endif

/* *****************************************************************************
 */
#define	DISTANCE_BASE	0.0f

template<typename dtyp>
static doinline void AddDistance(dtyp const &dist, dtyp &error) {
  error += dist;
}

template<typename dtyp>
static doinline void AddDistance(dtyp const &dist, dtyp &error, dtyp const &freq) {
  error += dist * freq;
}

template<const bool which, const int elements>
static doinline void MinDistance3(Scr3 &dist, int &index, Vec3 &value, Vec3 (&codes)[elements]) {
  Scr3 d0 = LengthSquared(value - codes[0]);
  Scr3 d1 = LengthSquared(value - codes[1]);
  Scr3 d2 = LengthSquared(value - codes[2]);

  // encourage OoO
  Scr3 da = Min(d0, d1);
  Scr3 db =    (d2    );
  dist    = Min(da, db);

  if (which) {
    // will cause VS to make them all cmovs
    if (d2 == dist) { index = 2; }
    if (d1 == dist) { index = 1; }
    if (d0 == dist) { index = 0; }
  }
}

template<const bool which, const int elements>
static doinline void MinDistance4(Scr3 &dist, int &index, Vec3 &value, Vec3 (&codes)[elements]) {
  Scr3 d0 = LengthSquared(value - codes[0]);
  Scr3 d1 = LengthSquared(value - codes[1]);
  Scr3 d2 = LengthSquared(value - codes[2]);
  Scr3 d3 = LengthSquared(value - codes[3]);

  // encourage OoO
  Scr3 da = Min(d0, d1);
  Scr3 db = Min(d2, d3);
  dist    = Min(da, db);

  if (which) {
    // will cause VS to make them all cmovs
    if (d3 == dist) { index = 3; }
    if (d2 == dist) { index = 2; }
    if (d1 == dist) { index = 1; }
    if (d0 == dist) { index = 0; }
  }
}

template<const bool which, const int elements>
static doinline void MinDistance4(Scr4 &dist, int &index, Vec4 &value, Vec4 (&codes)[elements], int &offset) {
  Scr4 d0 = LengthSquared(value - codes[offset + 0]);
  Scr4 d1 = LengthSquared(value - codes[offset + 1]);
  Scr4 d2 = LengthSquared(value - codes[offset + 2]);
  Scr4 d3 = LengthSquared(value - codes[offset + 3]);

  // encourage OoO
  Scr4 da = Min(d0, d1);
  Scr4 db = Min(d2, d3);
  dist    = Min(da, dist);
  dist    = Min(db, dist);

  if (which) {
    // will cause VS to make them all cmovs
    if (d0 == dist) { index = offset; } offset++;
    if (d1 == dist) { index = offset; } offset++;
    if (d2 == dist) { index = offset; } offset++;
    if (d3 == dist) { index = offset; } offset++;
  }
}

/* -----------------------------------------------------------------------------
 */
#define	DEVIANCE_SQUARE
#ifndef	DEVIANCE_SQUARE
#define	DEVIANCE_BASE	16.0f
#define	DEVIANCE_MAX	-1.0f	// largest angle
#define	DEVIANCE_MAXSUM	32.0f	// 16 * [-1,+1]
#else
#define	DEVIANCE_BASE	0.0f
#define	DEVIANCE_MAX	-1.0f	// largest angle
#define	DEVIANCE_MAXSUM	64.0f	// 16 * [-2, 0]²
#endif

template<typename dtyp>
static doinline void AddDeviance(dtyp const &dist, dtyp &error) {
#ifndef	DEVIANCE_SQUARE
  error -= dist;
#else
  // convert range [-1,+1] -> [-2,0]
  dtyp ang = dist + dtyp(DEVIANCE_MAX);
  dtyp sqr = ang * ang;

  error += sqr;
#endif
}

template<typename dtyp>
static doinline void AddDeviance(dtyp const &dist, dtyp &error, dtyp const &freq) {
#ifndef	DEVIANCE_SQUARE
  error -= dist * freq;
#else
  // convert range [-1,+1] -> [-2,0]
  dtyp ang = dist + dtyp(DEVIANCE_MAX);
  dtyp sqr = ang * ang;

  error += sqr * freq;
#endif
}

template<const bool which, const int elements>
static doinline void MinDeviance3(Scr3 &dist, int &index, Vec3 const &value, Vec3 const (&codes)[elements]) {
  Scr3 d0 = Dot(value, codes[0]);
  Scr3 d1 = Dot(value, codes[1]);
  Scr3 d2 = Dot(value, codes[2]);

  // encourage OoO
  Scr3 da = Max(d0, d1);
  Scr3 db =    (d2    );
  dist    = Max(da, db);

  if (which) {
    // will cause VS to make them all cmovs
    if (d2 == dist) { index = 2; }
    if (d1 == dist) { index = 1; }
    if (d0 == dist) { index = 0; }
  }
}

template<const bool which, const int elements>
static doinline void MinDeviance4(Scr3 &dist, int &index, Vec3 const &value, Vec3 const (&codes)[elements]) {
  Scr3 d0 = Dot(value, codes[0]);
  Scr3 d1 = Dot(value, codes[1]);
  Scr3 d2 = Dot(value, codes[2]);
  Scr3 d3 = Dot(value, codes[3]);

  // encourage OoO
  Scr3 da = Max(d0, d1);
  Scr3 db = Max(d2, d3);
  dist    = Max(da, db);

  if (which) {
    // will cause VS to make them all cmovs
    if (d3 == dist) { index = 3; }
    if (d2 == dist) { index = 2; }
    if (d1 == dist) { index = 1; }
    if (d0 == dist) { index = 0; }
  }
}

template<const bool which, const int elements>
static doinline void MinDeviance4(Scr4 &dist, int &index, Vec4 const &value, Vec4 const (&codes)[elements], int &offset) {
  Scr4 d0 = Dot(value, codes[offset + 0]);
  Scr4 d1 = Dot(value, codes[offset + 1]);
  Scr4 d2 = Dot(value, codes[offset + 2]);
  Scr4 d3 = Dot(value, codes[offset + 3]);

  // encourage OoO
  Scr4 da = Max(d0, d1);
  Scr4 db = Max(d2, d3);
  dist = Max(da, dist);
  dist = Max(db, dist);

  if (which) {
    // will cause VS to make them all cmovs
    if (d0 == dist) { index = offset; } offset++;
    if (d1 == dist) { index = offset; } offset++;
    if (d2 == dist) { index = offset; } offset++;
    if (d3 == dist) { index = offset; } offset++;
  }
}

/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
 */
template<const bool which, const int elements>
static doinline void MinDeviance3c(Scr3 &dist, int &index, Vec3 const &value, Vec3 const (&codes)[elements]) {
  Scr3 d0 = Dot(value, codes[0]);
  Scr3 d1 = Dot(value, codes[1]);
  Scr3 d2 = Dot(value, codes[2]);
  
  // select the smallest deviation (NaN as first arg is ignored!)
  dist = Max(d0, Scr3(DEVIANCE_MAX));
  dist = Max(d1, dist);
  dist = Max(d2, dist);

  if (which) {
    // will cause VS to make them all cmovs
    if (d2 == dist) { index = 2; }
    if (d1 == dist) { index = 1; }
    if (d0 == dist) { index = 0; }
  }
}

template<const bool which, const int elements>
static doinline void MinDeviance4c(Scr3 &dist, int &index, Vec3 const &value, Vec3 const (&codes)[elements]) {
  Scr3 d0 = Dot(value, codes[0]);
  Scr3 d1 = Dot(value, codes[1]);
  Scr3 d2 = Dot(value, codes[2]);
  Scr3 d3 = Dot(value, codes[3]);

  // select the smallest deviation (NaN as first arg is ignored!)
  dist = Max(d0, Scr3(DEVIANCE_MAX));
  dist = Max(d1, dist);
  dist = Max(d2, dist);
  dist = Max(d3, dist);

  if (which) {
    // will cause VS to make them all cmovs
    if (d3 == dist) { index = 3; }
    if (d2 == dist) { index = 2; }
    if (d1 == dist) { index = 1; }
    if (d0 == dist) { index = 0; }
  }
}

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
static int3 FloatToInt(float3 a, int3 limit) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // use ANSI round-to-zero behaviour to get round-to-nearest
  int3 i = (int3)round(a);

  // clamp to the limit
  return minimax(i, 0, limit);
}

static int4 FloatToInt(float4 a, int4 limit) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // use ANSI round-to-zero behaviour to get round-to-nearest
  int4 i = (int4)round(a);

  // clamp to the limit
  return minimax(i, 0, limit);
}

static int2 QuantizeFloatToInt(float2 a, int2 limit) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // clamp to [0,1], multiply, round
  // AMP: prevents the use of imin/imax because of the missing
  //      saturated integer mov
  return (int2)round(saturate(a) * limit);
}

static int3 QuantizeFloatToInt(float3 a, int3 limit) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // clamp to [0,1], multiply, round
  // AMP: prevents the use of imin/imax because of the missing
  //      saturated integer mov
  return (int3)round(saturate(a) * limit);
}

static int4 QuantizeFloatToInt(float4 a, int4 limit) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // clamp to [0,1], multiply, round
  return (int4)round(saturate(a) * limit);
}

static int FloatTo565(float3 colour) amp_restricted
{
  // get the components in the correct range
  int3 quant = QuantizeFloatToInt(colour, int3(31, 63, 31));

  // pack into a single value
  return (quant.r << 11) + (quant.g << 5) + (quant.b << 0);
}

static void FloatTo565(lineC2 colour, out lineI2 values) amp_restricted
{
  // get the components in the correct range
#if 0
  int3 quant[CVALS];

  quant[CSTRT] = QuantizeFloatToInt(colour[CSTRT], int3(31, 63, 31));
  quant[CSTOP] = QuantizeFloatToInt(colour[CSTOP], int3(31, 63, 31));

  // pack into a single value
  values[CSTRT] = (quant[CSTRT].r << 11) + (quant[CSTRT].g << 5) + quant[CSTRT].b;
  values[CSTOP] = (quant[CSTOP].r << 11) + (quant[CSTOP].g << 5) + quant[CSTOP].b;
#else
  // AMP: prefer vectorized (saves 1 instruction)
  int4 quantr;
  int2 quantg;
  int4 quantb;

  quantb = QuantizeFloatToInt(float4(
    colour[CSTRT].r, colour[CSTOP].r,
    colour[CSTRT].b, colour[CSTOP].b
  ), 31);

  quantg = QuantizeFloatToInt(float2(
    colour[CSTRT].g, colour[CSTOP].g
  ), 63);

  quantr = quantb << 11;
  quantg = quantg <<  5;
  quantb = quantb <<  0;

  // pack into a single value
  values[CSTRT] = (quantr.x) + (quantg.x) + (quantb.z);
  values[CSTOP] = (quantr.y) + (quantg.y) + (quantb.w);
#endif
}

static int3 Unpack565(int value) amp_restricted
{
  // get the components in the stored range
  ccr8 red   = (ccr8)((value >> 11) & 0x1F);
  ccr8 green = (ccr8)((value >>  5) & 0x3F);
  ccr8 blue  = (ccr8)( value        & 0x1F);

  // scale up to 8 bits & return the value
  return int3(
    (red   << 3) + (red   >> 2),
    (green << 2) + (green >> 4),
    (blue  << 3) + (blue  >> 2)
  );
}
#endif

}

#endif
