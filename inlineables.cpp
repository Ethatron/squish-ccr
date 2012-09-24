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
#if	!defined(USE_PRE)
static doinline int passreg FloatToInt(float a, int limit) ccr_restricted
{
  // use ANSI round-to-zero behaviour to get round-to-nearest
  int i = (int)(a + 0.5f);

  // clamp to the limit
  if (i < 0)
    i = 0;
  else if (i > limit)
    i = limit;

  // done
  return i;
}

static doinline int passreg FloatTo565(Vec3::Arg colour) ccr_restricted
{
  // get the components in the correct range
  int r = FloatToInt(31.0f * colour.X(), 31);
  int g = FloatToInt(63.0f * colour.Y(), 63);
  int b = FloatToInt(31.0f * colour.Z(), 31);

  // pack into a single value
  return (r << 11) | (g << 5) | b;
}

static doinline int passreg Unpack565(u8 const* packed, u8* colour) ccr_restricted
{
  // build the packed value
  int value = (int)packed[0] | ((int)packed[1] << 8);

  // get the components in the stored range
  u8 red   = (u8)((value >> 11) & 0x1F);
  u8 green = (u8)((value >>  5) & 0x3F);
  u8 blue  = (u8)( value        & 0x1F);

  // scale up to 8 bits
  colour[0] = (red   << 3) | (red   >> 2);
  colour[1] = (green << 2) | (green >> 4);
  colour[2] = (blue  << 3) | (blue  >> 2);
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

template<const int rb, const int gb, const int bb, const int ab, const int eb, const int sb>
static doinline void passreg FloatTo(Vec4 (&colour)[1], Col4 (&field)[1][FIELDN], int bitset) ccr_restricted
{
  /* not both yet */
  assert(!eb || !sb);

  // we can't just drop the eb/sb bits in fp-representation, we have to use the exact quantizer
  vQuantizer q = vQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb
  );
  
  // pack into a single value
          field[0][COLORA] = q.QuantizeToInt(colour[0]);
          field[0][COLORA] = field[0][COLORA] >> (eb + sb);

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

  vQuantizer q = vQuantizer(
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
  rexplicit[0] = rcomplete[0] >> (eb + sb);

  if (eb + sb) {
    /* stick g in a for skewing the rounding */
    if (!ab) {
      rcomplete[0] = Shuffle<1, 3>(rcomplete[0]);
    }
    else {
      Col4 z0 = IsNotZero(rcomplete[0]).SplatA();
      Col4 o0 = IsOne    (rcomplete[0]).SplatA();

      /* if alpha is black the shared bit must be 0 */
      rcomplete[0] &= z0 | Col4(~em);

      /* if alpha is white the shared bit must be 1 */
      rcomplete[0] |= o0 >> (32 - (eb + sb));
    }

    if (eb) {
      // get the components in the unique range
      runique[0] = (rcomplete[0] & umask) >> sb;

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

  // we can't just drop the eb/sb bits in fp-representation, we have to use the exact quantizer
  vQuantizer q = vQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb
  );
  
  // pack into a single value
          field[0][COLORA] = q.QuantizeToInt(colour[0]);
          field[1][COLORA] = q.QuantizeToInt(colour[1]);
          field[0][COLORA] = field[0][COLORA] >> (eb + sb);
          field[1][COLORA] = field[1][COLORA] >> (eb + sb);

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

  vQuantizer q = vQuantizer(
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
  rexplicit[0] = rcomplete[0] >> (eb + sb);
  rexplicit[1] = rcomplete[1] >> (eb + sb);

  if (eb + sb) {
    /* stick g in a for skewing the rounding */
    if (!ab) {
      rcomplete[0] = Shuffle<1, 3>(rcomplete[0]);
      rcomplete[1] = Shuffle<1, 3>(rcomplete[1]);
    }
    else {
      Col4 z0 = IsNotZero(rcomplete[0]).SplatA();
      Col4 z1 = IsNotZero(rcomplete[1]).SplatA();
      Col4 o0 = IsOne    (rcomplete[0]).SplatA();
      Col4 o1 = IsOne    (rcomplete[1]).SplatA();

      /* if alpha is black the shared bit must be 0 */
      rcomplete[0] &= z0 | Col4(~em);
      rcomplete[1] &= z1 | Col4(~em);

      /* if alpha is white the shared bit must be 1 */
      rcomplete[0] |= o0 >> (32 - (eb + sb));
      rcomplete[1] |= o1 >> (32 - (eb + sb));
    }

    if (eb) {
      // get the components in the unique range
      runique[0] = (rcomplete[0] & umask) >> sb;
      runique[1] = (rcomplete[1] & umask) >> sb;

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

  // we can't just drop the eb/sb bits in fp-representation, we have to use the exact quantizer
  vQuantizer q = vQuantizer(
    rb + eb + sb,
    gb + eb + sb,
    bb + eb + sb,
    ab + eb + sb
  );
  
  // pack into a single value
          field[0][COLORA] = q.QuantizeToInt(colour[0]);
          field[1][COLORA] = q.QuantizeToInt(colour[1]);
          field[2][COLORA] = q.QuantizeToInt(colour[2]);
          field[0][COLORA] = field[0][COLORA] >> (eb + sb);
          field[1][COLORA] = field[1][COLORA] >> (eb + sb);
          field[2][COLORA] = field[2][COLORA] >> (eb + sb);

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

  vQuantizer q = vQuantizer(
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
  rexplicit[0] = rcomplete[0] >> (eb + sb);
  rexplicit[1] = rcomplete[1] >> (eb + sb);
  rexplicit[2] = rcomplete[2] >> (eb + sb);

  if (eb + sb) {
    /* stick g in a for skewing the rounding */
    if (!ab) {
      rcomplete[0] = Shuffle<1, 3>(rcomplete[0]);
      rcomplete[1] = Shuffle<1, 3>(rcomplete[1]);
      rcomplete[2] = Shuffle<1, 3>(rcomplete[2]);
    }
    else {
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
      rcomplete[0] |= o0 >> (32 - (eb + sb));
      rcomplete[1] |= o1 >> (32 - (eb + sb));
      rcomplete[2] |= o2 >> (32 - (eb + sb));
    }

    if (eb) {
      // get the components in the unique range
      runique[0] = (rcomplete[0] & umask) >> sb;
      runique[1] = (rcomplete[1] & umask) >> sb;
      runique[2] = (rcomplete[2] & umask) >> sb;

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

/* -----------------------------------------------------------------------------
 */
static doinline void passreg Codebook3(u8 (&codes)[4*4], bool bw) ccr_restricted
{
  // generate the midpoints
  for (int i = 0; i < 3; ++i) {
    int c = codes[0 + i];
    int d = codes[4 + i];

    if (bw) {
      codes[ 8 + i] = (u8)((1*c + 1*d) / 2);
      codes[12 + i] = 0;
    }
    else {
      codes[ 8 + i] = (u8)((2*c + 1*d) / 3);
      codes[12 + i] = (u8)((1*c + 2*d) / 3);
    }
  }

  // fill in alpha for the intermediate values
  codes[ 8 + 3] = 255;
  codes[12 + 3] = bw ? 0 : 255;
}

static doinline void passreg Codebook4(u8 (&codes)[4*4]) ccr_restricted
{
  // generate the midpoints
  for (int i = 0; i < 3; ++i) {
    int c = codes[0 + i];
    int d = codes[4 + i];

    {
      codes[ 8 + i] = (u8)((2*c + 1*d) / 3);
      codes[12 + i] = (u8)((1*c + 2*d) / 3);
    }
  }

  // fill in alpha for the intermediate values
  codes[ 8 + 3] = 255;
  codes[12 + 3] = 255;
}

static int passreg CodebookP(u8 *codes, int bits) ccr_restricted
{
  // generate the midpoints
  for (int m = 0; m < 4; ++m) {
    const int j = (1 << bits) - 1;

    int c = codes[0 * 4 + m];
    int d = codes[j * 4 + m];

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

static int passreg CodebookP(Vec4 *codes, int bits, Vec4::Arg start, Vec4::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  codes[0] = start;
  codes[j] = end;
  
  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    Vec4 s = weights_V4[bits][j - i] * start;
    Vec4 e = weights_V4[bits][i + 0] * end;

    codes[i] = s + e;
  }

  return (1 << bits);
}

/* -----------------------------------------------------------------------------
 */
static int passreg CodebookP(Col4 *codes, int bits, Col4::Arg start, Col4::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  codes[0] = start;
  codes[j] = end;

  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    Col4 s = (weights_C4[bits][j - i]) * start;
    Col4 e = (weights_C4[bits][i + 0]) * end;

    codes[i] = (s + e + Col4(32)) >> 6;
  }

  return (1 << bits);
}

template<const int bits>
static int passreg CodebookP(int *codes, Col4::Arg start, Col4::Arg end) ccr_restricted
{
  const int j = (1 << bits) - 1;

  PackBytes(start, codes[0]);
  PackBytes(end  , codes[j]);

  // the quantizer is not equi-distant, but it is symmetric
  for (int i = 1; i < j; i++) {
    Col4 s = (weights_C4[bits][j - i]) * start;
    Col4 e = (weights_C4[bits][i + 0]) * end;

    PackBytes((s + e + Col4(32)) >> 6, codes[i]);
  }

  return (1 << bits);
}

#endif

/* *****************************************************************************
 */
#if	defined(USE_AMP) || defined(USE_COMPUTE)
static int3 FloatToInt(float3 a, int3 limit) amp_restricted
{
#if	!defined(USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // use ANSI round-to-zero behaviour to get round-to-nearest
  int3 i = (int3)round(a);

  // clamp to the limit
  return minimax(i, 0, limit);
}

static int4 FloatToInt(float4 a, int4 limit) amp_restricted
{
#if	!defined(USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // use ANSI round-to-zero behaviour to get round-to-nearest
  int4 i = (int4)round(a);

  // clamp to the limit
  return minimax(i, 0, limit);
}

static int2 QuantizeFloatToInt(float2 a, int2 limit) amp_restricted
{
#if	!defined(USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // clamp to [0,1], multiply, round
  // AMP: prevents the use of imin/imax because of the missing
  //      saturated integer mov
  return (int2)round(saturate(a) * limit);
}

static int3 QuantizeFloatToInt(float3 a, int3 limit) amp_restricted
{
#if	!defined(USE_COMPUTE)
  using namespace Concurrency::vector_math;
#endif

  // clamp to [0,1], multiply, round
  // AMP: prevents the use of imin/imax because of the missing
  //      saturated integer mov
  return (int3)round(saturate(a) * limit);
}

static int4 QuantizeFloatToInt(float4 a, int4 limit) amp_restricted
{
#if	!defined(USE_COMPUTE)
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
  return (quant.r << 11) | (quant.g << 5) | quant.b;
}

static void FloatTo565(lineC2 colour, out lineI2 values) amp_restricted
{
  // get the components in the correct range
#if 0
  int3 quant[CVALS];

  quant[CSTRT] = QuantizeFloatToInt(colour[CSTRT], int3(31, 63, 31));
  quant[CSTOP] = QuantizeFloatToInt(colour[CSTOP], int3(31, 63, 31));

  // pack into a single value
  values[CSTRT] = (quant[CSTRT].r << 11) | (quant[CSTRT].g << 5) | quant[CSTRT].b;
  values[CSTOP] = (quant[CSTOP].r << 11) | (quant[CSTOP].g << 5) | quant[CSTOP].b;
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
  values[CSTRT] = (quantr.x) | (quantg.x) | (quantb.z);
  values[CSTOP] = (quantr.y) | (quantg.y) | (quantb.w);
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
    (red   << 3) | (red   >> 2),
    (green << 2) | (green >> 4),
    (blue  << 3) | (blue  >> 2)
  );
}
#endif

}

#endif
