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

#if defined(_MSC_VER) && (_MSC_VER > 1300)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>
#endif

#include "paletteblock.h"

#define SBSTART	0
#define SBEND	3

#include "inlineables.cpp"

#pragma warning(disable: 4100)

namespace squish {

/* *****************************************************************************
 */
#define S1(x)  u16(x)
  // http://msdn.microsoft.com/en-us/library/hh308953%28v=vs.85%29
  const u16 weights_u16[5][16] = {
    {S1(0)                                                                                                                        },  // 0
    {S1(0),                                                                                                                 S1(64)},  // 1
    {S1(0),                                 S1(21),                                 S1(43),                                 S1(64)},  // 2
    {S1(0),          S1(9),         S1(18),         S1(27),                 S1(37),         S1(46),         S1(55),         S1(64)},  // 3
    {S1(0),  S1(4),  S1(9), S1(13), S1(17), S1(21), S1(26), S1(30), S1(34), S1(38), S1(43), S1(47), S1(51), S1(55), S1(60), S1(64)}   // 4
  };
#undef S1

#define V4(x)  Vec4((1.0f / 64.0f) * x)
  // http://msdn.microsoft.com/en-us/library/hh308953%28v=vs.85%29
  const Vec4 weights_V4[5][16] = {
    {V4(0)                                                                                                                        },  // 0
    {V4(0),                                                                                                                 V4(64)},  // 1
    {V4(0),                                 V4(21),                                 V4(43),                                 V4(64)},  // 2
    {V4(0),          V4(9),         V4(18),         V4(27),                 V4(37),         V4(46),         V4(55),         V4(64)},  // 3
    {V4(0),  V4(4),  V4(9), V4(13), V4(17), V4(21), V4(26), V4(30), V4(34), V4(38), V4(43), V4(47), V4(51), V4(55), V4(60), V4(64)}   // 4
  };
#undef V4

#define C4(x)  Col4(x)
  // http://msdn.microsoft.com/en-us/library/hh308953%28v=vs.85%29
  const Col4 weights_C4[5][16] = {
    {C4(0)                                                                                                                        },  // 0
    {C4(0),                                                                                                                 C4(64)},  // 1
    {C4(0),                                 C4(21),                                 C4(43),                                 C4(64)},  // 2
    {C4(0),          C4(9),         C4(18),         C4(27),                 C4(37),         C4(46),         C4(55),         C4(64)},  // 3
    {C4(0),  C4(4),  C4(9), C4(13), C4(17), C4(21), C4(26), C4(30), C4(34), C4(38), C4(43), C4(47), C4(51), C4(55), C4(60), C4(64)}   // 4
  };
#undef C4

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
#ifdef __GNUC__
  __attribute__ ((__aligned__ (16)))
#else
  __declspec(align(16))
#endif

extern const unsigned int blockxor[64][/*6*/5][4] = {
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x0000ffff,0x0000ffff,0x0000ffff}, {0xffff0000,0xffff0000,0xffff0000,0xffff0000}, {0x0000ffff,0x0000ffff,0x000000ff,0x00000000}, {0xffff0000,0xffff0000,0xff000000,0x00000000}, {0x00000000,0x00000000,0x00ffff00,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ffffff,0x00ffffff,0x00ffffff,0x00ffffff}, {0xff000000,0xff000000,0xff000000,0xff000000}, {0x00ffffff,0x0000ffff,0x00000000,0x00000000}, {0xff000000,0xffff0000,0xffff0000,0xff000000}, {0x00000000,0x00000000,0x0000ffff,0x00ffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x000000ff,0x000000ff,0x000000ff,0x000000ff}, {0xffffff00,0xffffff00,0xffffff00,0xffffff00}, {0xffffffff,0x00ffff00,0x00000000,0x00000000}, {0x00000000,0xff000000,0xffff0000,0xffff0000}, {0x00000000,0x000000ff,0x0000ffff,0x0000ffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ffffff,0x0000ffff,0x0000ffff,0x000000ff}, {0xff000000,0xffff0000,0xffff0000,0xffffff00}, {0x000000ff,0x0000ffff,0x0000ffff,0x000000ff}, {0x00000000,0x00000000,0xffff0000,0xffffff00}, {0xffffff00,0xffff0000,0x00000000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0x00ffffff,0x00ffffff,0x0000ffff}, {0x00000000,0xff000000,0xff000000,0xffff0000}, {0xffffffff,0xffffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0x0000ffff,0x0000ffff}, {0x00000000,0x00000000,0xffff0000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x000000ff,0x000000ff,0x00000000}, {0xffff0000,0xffffff00,0xffffff00,0xffffffff}, {0x0000ffff,0x0000ffff,0x0000ffff,0x0000ffff}, {0xffff0000,0xffff0000,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffff0000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ffffff,0x0000ffff,0x000000ff,0x00000000}, {0xff000000,0xffff0000,0xffffff00,0xffffffff}, {0x0000ffff,0x0000ffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffffff,0xffffffff}, {0xffff0000,0xffff0000,0x00000000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0x00ffffff,0x0000ffff,0x000000ff}, {0x00000000,0xff000000,0xffff0000,0xffffff00}, {0x0000ffff,0x0000ffff,0x00000000,0x00000000}, {0xffff0000,0xffff0000,0xffff0000,0xffff0000}, {0x00000000,0x00000000,0x0000ffff,0x0000ffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffffff,0x00ffffff,0x0000ffff}, {0x00000000,0x00000000,0xff000000,0xffff0000}, {0xffffffff,0xffffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffffff,0x00000000}, {0x00000000,0x00000000,0x00000000,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x000000ff,0x00000000,0x00000000}, {0xffff0000,0xffffff00,0xffffffff,0xffffffff}, {0xffffffff,0x00000000,0x00000000,0x00000000}, {0x00000000,0xffffffff,0xffffffff,0x00000000}, {0x00000000,0x00000000,0x00000000,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0x00ffffff,0x000000ff,0x00000000}, {0x00000000,0xff000000,0xffffff00,0xffffffff}, {0xffffffff,0x00000000,0x00000000,0x00000000}, {0x00000000,0xffffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffffff,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffffff,0x00ffffff,0x000000ff}, {0x00000000,0x00000000,0xff000000,0xffffff00}, {0x0000ffff,0x0000ffff,0x0000ffff,0x0000ffff}, {0x00ff0000,0x00ff0000,0x00ff0000,0x00ff0000}, {0xff000000,0xff000000,0xff000000,0xff000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ffffff,0x000000ff,0x00000000,0x00000000}, {0xff000000,0xffffff00,0xffffffff,0xffffffff}, {0x000000ff,0x000000ff,0x000000ff,0x000000ff}, {0x00ffff00,0x00ffff00,0x00ffff00,0x00ffff00}, {0xff000000,0xff000000,0xff000000,0xff000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffffff,0xffffffff}, {0x000000ff,0x000000ff,0x000000ff,0x000000ff}, {0x0000ff00,0x0000ff00,0x0000ff00,0x0000ff00}, {0xffff0000,0xffff0000,0xffff0000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0x00000000,0x00000000,0x00000000}, {0x00000000,0xffffffff,0xffffffff,0xffffffff}, {0x0000ffff,0x000000ff,0x00000000,0x00000000}, {0xffff0000,0x00ffff00,0x0000ffff,0x000000ff}, {0x00000000,0xff000000,0xffff0000,0xffffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffffff,0xffffffff,0x00000000}, {0x00000000,0x00000000,0x00000000,0xffffffff}, {0x0000ffff,0x00ffff00,0xffff0000,0xff000000}, {0xffff0000,0xff000000,0x00000000,0x00000000}, {0x00000000,0x000000ff,0x0000ffff,0x00ffffff} },

  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffff00,0xff000000,0x00000000}, {0x00000000,0x000000ff,0x00ffffff,0xffffffff}, {0x00ffffff,0x0000ffff,0x000000ff,0x00000000}, {0xff000000,0xffff0000,0x00ffff00,0x0000ffff}, {0x00000000,0x00000000,0xff000000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x000000ff,0x00ffffff,0xffffffff,0xffffffff}, {0xffffff00,0xff000000,0x00000000,0x00000000}, {0x000000ff,0x0000ffff,0x00ffff00,0xffff0000}, {0xffffff00,0xffff0000,0xff000000,0x00000000}, {0x00000000,0x00000000,0x000000ff,0x0000ffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffffff,0xffffff00,0xff000000}, {0x00000000,0x00000000,0x000000ff,0x00ffffff}, {0xffffffff,0x00000000,0x00000000,0x00000000}, {0x00000000,0x0000ffff,0x0000ffff,0x0000ffff}, {0x00000000,0xffff0000,0xffff0000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x000000ff,0x0000ffff,0x00ffffff,0xffffffff}, {0xffffff00,0xffff0000,0xff000000,0x00000000}, {0x0000ffff,0x0000ffff,0x0000ffff,0x00000000}, {0x00000000,0x00000000,0x00000000,0xffffffff}, {0xffff0000,0xffff0000,0xffff0000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x00ffffff,0xffffffff,0xffffffff}, {0xffff0000,0xff000000,0x00000000,0x00000000}, {0x000000ff,0x000000ff,0x000000ff,0x000000ff}, {0xffffff00,0xffffff00,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffff00,0xffffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffff00,0xffff0000,0xff000000}, {0x00000000,0x000000ff,0x0000ffff,0x00ffffff}, {0x00ffffff,0x00ffffff,0x00000000,0x00000000}, {0xff000000,0xff000000,0xff000000,0xff000000}, {0x00000000,0x00000000,0x00ffffff,0x00ffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffffff,0xffffff00,0xffff0000}, {0x00000000,0x00000000,0x000000ff,0x0000ffff}, {0xffffffff,0x0000ffff,0x000000ff,0x000000ff}, {0x00000000,0xffff0000,0x0000ff00,0x0000ff00}, {0x00000000,0x00000000,0xffff0000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x000000ff,0x0000ffff,0x0000ffff,0x00ffffff}, {0xffffff00,0xffff0000,0xffff0000,0xff000000}, {0xffffffff,0xffff0000,0xff000000,0xff000000}, {0x00000000,0x0000ffff,0x00ff0000,0x00ff0000}, {0x00000000,0x00000000,0x0000ffff,0x0000ffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x00ffffff,0x00ffffff,0xffffffff}, {0xffff0000,0xff000000,0xff000000,0x00000000}, {0x000000ff,0x000000ff,0x0000ffff,0xffffffff}, {0x0000ff00,0x0000ff00,0xffff0000,0x00000000}, {0xffff0000,0xffff0000,0x00000000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffffff00,0xffffff00,0xffff0000}, {0x00000000,0x000000ff,0x000000ff,0x0000ffff}, {0x0000ffff,0x0000ffff,0x00000000,0x00000000}, {0x00ff0000,0x00ff0000,0x0000ffff,0x00000000}, {0xff000000,0xff000000,0xffff0000,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0xff0000ff,0xff0000ff,0xff0000ff}, {0x00ffff00,0x00ffff00,0x00ffff00,0x00ffff00}, {0xff0000ff,0x00000000,0x00000000,0xff0000ff}, {0x00ffff00,0xff0000ff,0xff0000ff,0x00ffff00}, {0x00000000,0x00ffff00,0x00ffff00,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0xff0000ff,0xff0000ff,0xffff0000}, {0xffff0000,0x00ffff00,0x00ffff00,0x0000ffff}, {0xffffffff,0xff0000ff,0x00000000,0x00000000}, {0x00000000,0x00ffff00,0xff0000ff,0xff0000ff}, {0x00000000,0x00000000,0x00ffff00,0x00ffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ffffff,0x000000ff,0xff000000,0xffffff00}, {0xff000000,0xffffff00,0x00ffffff,0x000000ff}, {0x0000ffff,0x00ff0000,0x00ff0000,0x0000ffff}, {0x00000000,0x0000ffff,0x0000ffff,0x00000000}, {0xffff0000,0xff000000,0xff000000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0x00000000,0x00000000,0xffffffff}, {0x00000000,0xffffffff,0xffffffff,0x00000000}, {0xff0000ff,0xff0000ff,0x00ffff00,0x00000000}, {0x00ffff00,0x00ffff00,0x00000000,0x00000000}, {0x00000000,0x00000000,0xff0000ff,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x000000ff,0x00ffffff,0xffffff00,0xff000000}, {0xffffff00,0xff000000,0x000000ff,0x00ffffff}, {0x0000ffff,0x000000ff,0x000000ff,0x0000ffff}, {0xffff0000,0x0000ff00,0x0000ff00,0xffff0000}, {0x00000000,0xffff0000,0xffff0000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x00ffff00,0x00ffff00,0xffff0000}, {0xffff0000,0xff0000ff,0xff0000ff,0x0000ffff}, {0xffffffff,0xffffff00,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffff0000,0xff000000}, {0x00000000,0x000000ff,0x0000ffff,0x00ffffff} },

  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ff00ff,0x00ff00ff,0x00ff00ff,0x00ff00ff}, {0xff00ff00,0xff00ff00,0xff00ff00,0xff00ff00}, {0xffffffff,0x00ffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0x0000ffff,0x000000ff}, {0x00000000,0xff000000,0xffff0000,0xffffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0x00000000,0xffffffff,0x00000000}, {0x00000000,0xffffffff,0x00000000,0xffffffff}, {0x000000ff,0x0000ffff,0x0000ffff,0x0000ffff}, {0x00000000,0x00000000,0x00ff0000,0xffff0000}, {0xffffff00,0xffff0000,0xff000000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ff00ff,0xff00ff00,0x00ff00ff,0xff00ff00}, {0xff00ff00,0x00ff00ff,0xff00ff00,0x00ff00ff}, {0x0000ffff,0x0000ffff,0x0000ffff,0x000000ff}, {0xffff0000,0x00ff0000,0x00000000,0x00000000}, {0x00000000,0xff000000,0xffff0000,0xffffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x0000ffff,0xffff0000,0xffff0000}, {0xffff0000,0xffff0000,0x0000ffff,0x0000ffff}, {0xff0000ff,0xff0000ff,0xff0000ff,0xff0000ff}, {0x0000ff00,0x0000ff00,0x0000ff00,0x0000ff00}, {0x00ff0000,0x00ff0000,0x00ff0000,0x00ff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0xffff0000,0x0000ffff,0xffff0000}, {0xffff0000,0x0000ffff,0xffff0000,0x0000ffff}, {0xffffffff,0x00000000,0x00000000,0xffffffff}, {0x00000000,0xffffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffffff,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ff00ff,0x00ff00ff,0xff00ff00,0xff00ff00}, {0xff00ff00,0xff00ff00,0x00ff00ff,0x00ff00ff}, {0xff0000ff,0x00ff0000,0x0000ff00,0xff0000ff}, {0x0000ff00,0xff0000ff,0x00ff0000,0x0000ff00}, {0x00ff0000,0x0000ff00,0xff0000ff,0x00ff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0x00ffff00,0xff0000ff,0x00ffff00}, {0x00ffff00,0xff0000ff,0x00ffff00,0xff0000ff}, {0xff0000ff,0x0000ff00,0x00ff0000,0xff0000ff}, {0x0000ff00,0x00ff0000,0xff0000ff,0x0000ff00}, {0x00ff0000,0xff0000ff,0x0000ff00,0x00ff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ff00ff,0xff00ff00,0xff00ff00,0x00ff00ff}, {0xff00ff00,0x00ff00ff,0x00ff00ff,0xff00ff00}, {0x0000ffff,0xffff0000,0x00000000,0x0000ffff}, {0xffff0000,0x00000000,0x0000ffff,0xffff0000}, {0x00000000,0x0000ffff,0xffff0000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x000000ff,0x0000ffff,0xffff0000,0xff000000}, {0xffffff00,0xffff0000,0x0000ffff,0x00ffffff}, {0x0000ffff,0x00000000,0xffff0000,0x0000ffff}, {0xffff0000,0x0000ffff,0x00000000,0xffff0000}, {0x00000000,0xffff0000,0x0000ffff,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ffffff,0x0000ffff,0xffff0000,0xffffff00}, {0xff000000,0xffff0000,0x0000ffff,0x000000ff}, {0x00ff00ff,0x00ff00ff,0x00000000,0x00000000}, {0xff00ff00,0xff00ff00,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffffff,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0xff00ffff,0xffff00ff,0xffff0000}, {0xffff0000,0x00ff0000,0x0000ff00,0x0000ffff}, {0xffffffff,0xffffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0xff00ff00,0xff00ff00}, {0x00000000,0x00000000,0x00ff00ff,0x00ff00ff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x0000ff00,0x00ff0000,0xffff0000}, {0xffff0000,0xffff00ff,0xff00ffff,0x0000ffff}, {0x0000ffff,0x00000000,0x0000ffff,0x00000000}, {0x00000000,0x0000ffff,0x00000000,0x0000ffff}, {0xffff0000,0xffff0000,0xffff0000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0x00ffff00,0x00ffff00,0xff0000ff}, {0x00ffff00,0xff0000ff,0xff0000ff,0x00ffff00}, {0x0000ffff,0x0000ffff,0x0000ffff,0x0000ffff}, {0x00000000,0xffff0000,0x00000000,0xffff0000}, {0xffff0000,0x00000000,0xffff0000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0xffff0000,0xffff0000,0x0000ffff}, {0xffff0000,0x0000ffff,0x0000ffff,0xffff0000}, {0xff0000ff,0x00000000,0xff0000ff,0x00000000}, {0x00000000,0xff0000ff,0x00000000,0xff0000ff}, {0x00ffff00,0x00ffff00,0x00ffff00,0x00ffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0xff0000ff,0x00ffff00,0x00ffff00}, {0x00ffff00,0x00ffff00,0xff0000ff,0xff0000ff}, {0x00ff00ff,0x00000000,0x00000000,0x00ff00ff}, {0xff00ff00,0x00000000,0x00000000,0xff00ff00}, {0x00000000,0xffffffff,0xffffffff,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xff0000ff,0xff0000ff,0xffffffff}, {0x00000000,0x00ffff00,0x00ffff00,0x00000000}, {0xffffffff,0x00000000,0x00000000,0x00000000}, {0x00000000,0xff00ff00,0xff00ff00,0xff00ff00}, {0x00000000,0x00ff00ff,0x00ff00ff,0x00ff00ff} },

  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffff00ff,0xff000000,0xffff00ff,0xffffffff}, {0x0000ff00,0x00ffffff,0x0000ff00,0x00000000}, {0x00ff00ff,0x00ff00ff,0x00ff00ff,0x00000000}, {0xff00ff00,0xff00ff00,0xff00ff00,0x00000000}, {0x00000000,0x00000000,0x00000000,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff00ffff,0x000000ff,0xff00ffff,0xffffffff}, {0x00ff0000,0xffffff00,0x00ff0000,0x00000000}, {0x000000ff,0x000000ff,0x000000ff,0x000000ff}, {0x00000000,0xffffff00,0x00000000,0xffffff00}, {0xffffff00,0x00000000,0xffffff00,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xff00ffff,0x000000ff,0xff00ffff}, {0x00000000,0x00ff0000,0xffffff00,0x00ff0000}, {0x00ffffff,0x00000000,0x00ffffff,0x00000000}, {0x00000000,0x00ffffff,0x00000000,0x00ffffff}, {0xff000000,0xff000000,0xff000000,0xff000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0xffff00ff,0xff000000,0xffff00ff}, {0x00000000,0x0000ff00,0x00ffffff,0x0000ff00}, {0xffffffff,0x00000000,0x00000000,0x00000000}, {0x00000000,0x00ffff00,0x00ffff00,0x00ffff00}, {0x00000000,0xff0000ff,0xff0000ff,0xff0000ff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0xffff0000,0x00ffff00,0x0000ffff}, {0x00ffff00,0x0000ffff,0xff0000ff,0xffff0000}, {0x000000ff,0x000000ff,0x000000ff,0x000000ff}, {0x00000000,0xffffff00,0xffffff00,0x00000000}, {0xffffff00,0x00000000,0x00000000,0xffffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0xff0000ff,0xffff0000,0x00ffff00}, {0xffff0000,0x00ffff00,0x0000ffff,0xff0000ff}, {0x00ffffff,0x00000000,0x00000000,0x00ffffff}, {0x00000000,0x00ffffff,0x00ffffff,0x00000000}, {0xff000000,0xff000000,0xff000000,0xff000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0x0000ffff,0x00ffff00,0xffff0000}, {0x00ffff00,0xffff0000,0xff0000ff,0x0000ffff}, {0xff0000ff,0xff0000ff,0xff0000ff,0x00000000}, {0x00ffff00,0x00ffff00,0x00ffff00,0x00000000}, {0x00000000,0x00000000,0x00000000,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x00ffff00,0xffff0000,0xff0000ff}, {0xffff0000,0xff0000ff,0x0000ffff,0x00ffff00}, {0xffffffff,0xffffffff,0x00000000,0x00000000}, {0x00000000,0x00000000,0x00ffff00,0x00ffff00}, {0x00000000,0x00000000,0xff0000ff,0xff0000ff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0xffff0000,0xffff0000,0x00ffff00}, {0x00ffff00,0x0000ffff,0x0000ffff,0xff0000ff}, {0xff0000ff,0xff0000ff,0x00000000,0x00000000}, {0x00ffff00,0x00ffff00,0x00000000,0x00000000}, {0x00000000,0x00000000,0xffffffff,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff0000ff,0x0000ffff,0x0000ffff,0x00ffff00}, {0x00ffff00,0xffff0000,0xffff0000,0xff0000ff}, {0x0000ffff,0x0000ffff,0x0000ffff,0x0000ffff}, {0x00000000,0xffff0000,0xffff0000,0x00000000}, {0xffff0000,0x00000000,0x00000000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x000000ff,0xff000000,0xffffff00,0x00ffffff}, {0xffffff00,0x00ffffff,0x000000ff,0xff000000}, {0x0000ffff,0x00000000,0x00000000,0x0000ffff}, {0x00000000,0x0000ffff,0x0000ffff,0x00000000}, {0xffff0000,0xffff0000,0xffff0000,0xffff0000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x00ffffff,0xffffff00,0xff000000,0x000000ff}, {0xff000000,0x000000ff,0x00ffffff,0xffffff00}, {0xffffffff,0xffffffff,0xffffffff,0x00000000}, {0x00000000,0x00000000,0x00000000,0x00ffff00}, {0x00000000,0x00000000,0x00000000,0xff0000ff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffffffff,0x00000000,0x0000ffff,0x0000ffff}, {0x00000000,0xffffffff,0xffff0000,0xffff0000}, {0x00ffffff,0x00ffffff,0x00ffffff,0x00ffffff}, {0x00000000,0xff000000,0x00000000,0xff000000}, {0xff000000,0x00000000,0xff000000,0x00000000} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0x0000ffff,0x0000ffff,0x00000000,0xffffffff}, {0xffff0000,0xffff0000,0xffffffff,0x00000000}, {0x000000ff,0x00000000,0x000000ff,0x00000000}, {0x00000000,0x000000ff,0x00000000,0x000000ff}, {0xffffff00,0xffffff00,0xffffff00,0xffffff00} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xff00ffff,0xff00ffff,0xff000000,0xff000000}, {0x00ff0000,0x00ff0000,0x00ffffff,0x00ffffff}, {0x00ff00ff,0x00000000,0x00000000,0x00000000}, {0xff00ff00,0x00000000,0x00000000,0x00000000}, {0x00000000,0xffffffff,0xffffffff,0xffffffff} },
  { /*{0xffffffff,0xffffffff,0xffffffff,0xffffffff},*/ {0xffff00ff,0xffff00ff,0x000000ff,0x000000ff}, {0x0000ff00,0x0000ff00,0xffffff00,0xffffff00}, {0x000000ff,0x0000ff00,0x00ff0000,0xff000000}, {0xffffff00,0xffff0000,0xff000000,0x00000000}, {0x00000000,0x000000ff,0x0000ffff,0x00ffffff} }
};

extern const int shorterindex[64][/*6*/5] = {
 //   1       1 | 2      1 | 2 | 3     1 | 2 | 3 (ordered)
  { /*0,*/  /*0,*/15,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3, 8,  /*0,*/ 3, 8},
  { /*0,*/  /*0,*/15,  /*0,*/15, 8,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 8,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 6,15,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 6,15,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 6,15,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3, 8,  /*0,*/ 3, 8},

  { /*0,*/  /*0,*/15,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 3, 8,  /*0,*/ 3, 8},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 3,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 3, 8,  /*0,*/ 3, 8},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 6,15,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/10, 8,  /*0,*/ 8,10},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 5, 3,  /*0,*/ 3, 5},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 8, 6,  /*0,*/ 6, 8},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 6,10,  /*0,*/ 6,10},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15,10,  /*0,*/10,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 8,  /*0,*/ 8,15},

  { /*0,*/  /*0,*/15,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 6,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 5,10,  /*0,*/ 5,10},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 6,10,  /*0,*/ 6,10},
  { /*0,*/  /*0,*/ 8,  /*0,*/10, 8,  /*0,*/ 8,10},
  { /*0,*/  /*0,*/15,  /*0,*/ 8, 9,  /*0,*/ 8, 9},
  { /*0,*/  /*0,*/15,  /*0,*/15,10,  /*0,*/10,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 6,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 8,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 3,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 6,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 6,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/ 6,  /*0,*/15, 8,  /*0,*/ 8,15},

  { /*0,*/  /*0,*/ 6,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 3,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 6,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/10,15,  /*0,*/10,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 5,15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/10,15,  /*0,*/10,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/13,15,  /*0,*/13,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/12,15,  /*0,*/12,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 3,15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3, 8,  /*0,*/ 3, 8}
};

/* -----------------------------------------------------------------------------
 */
template<const int sets, const int ibits, const int begin>
static void passreg WritePaletteBlock(int partition, Col4 (&idx)[1], Col4 &blkl, Col4 &blkh)
{
  Col4 iidx = idx[0];
  Col4 iblk;

  /* none of the cases straddles the lo to hi border, all go into hi
   *
   * WritePaletteBlock<3, 3, 83>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<2, 3, 82>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<3, 2, 99>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<2, 2, 98>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<1, 4, 66>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<2, 2, 98>(partition, remapped, blkl, blkh);
   */
  assume(partition >= 0 && partition <= 64);

  // max 16*4 -> 64bits
  iblk = CopyBits<ibits, ibits *  0>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  1>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  2>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  3>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  4>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  5>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  6>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  7>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  8>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  9>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 10>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 11>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 12>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 13>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 14>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 15>(iblk, iidx); //dx = ShiftRight<8>(iidx);

  // one bit is omited per set, so it's one instruction per set + 1
  switch (sets) {
    case 1: {
      // always index 0
      blkh = CopyBits<ibits *  1 - 1, begin - 64 +         0>(blkh, iblk);
      iblk = ShiftRightHalf<ibits *  1>(iblk);
      blkh = CopyBits<ibits * 15 + 0, begin - 64 + ibits - 1>(blkh, iblk);
//    iblk = ShiftRightHalf<ibits * 15>(iblk);
    } break;
    case 2: {
      // always index 0
      blkh = CopyBits<ibits *  1 - 1, begin - 64 +         0>(blkh, iblk);
      iblk = ShiftRightHalf<ibits *  1>(iblk);
      blkh = CopyBits<ibits * 15 + 0, begin - 64 + ibits - 1>(blkh, iblk);
//    iblk = ShiftRightHalf<ibits * 15>(iblk);

      // if it's not the last bit which is cut (low probability)
      int len, bgn = begin - 64;
      if ((len = shorterindex[partition][0]) < 15) {
        len = (len * ibits) + ibits;
        bgn = (bgn + len - 2);

      	// subtract the already conducted shifts
        iblk = ShiftRightHalf(iblk, len - ibits);

        // remaining length can be anything, length overflow is silent
        blkh = CopyBits(blkh, iblk, ibits * 14, bgn);
//      iblk = ShiftRightHalf<ibits * 15>(iblk);
      }
    } break;
    case 3: {
      // always index 0
      blkh = CopyBits<ibits *  1 - 1, begin - 64 +         0>(blkh, iblk);
      iblk = ShiftRightHalf<ibits *  1>(iblk);
      blkh = CopyBits<ibits * 15 + 0, begin - 64 + ibits - 1>(blkh, iblk);
//    iblk = ShiftRightHalf<ibits * 15>(iblk);

      int bgn = begin - 64;
      int fln, fbn;
      fln = shorterindex[partition][3]; {
        fln = (fln * ibits) + ibits;
        fbn = (bgn + fln - 2);

      	// subtract the already conducted shifts
        iblk = ShiftRightHalf(iblk, fln - ibits);

        // remaining length can be anything, length overflow is silent
        blkh = CopyBits(blkh, iblk, ibits * 14, fbn);
//      iblk = ShiftRightHalf(iblk, ibits * 15);
      }

      // if it's not the last bit which is cut (low probability)
      int sln, sbn;
      if ((sln = shorterindex[partition][4]) < 15) {
        sln = (sln * ibits) + ibits;
        sbn = (bgn + sln - 3);

      	// subtract the already conducted shifts
        iblk = ShiftRightHalf(iblk, sln - fln);

        // remaining length can be anything, length overflow is silent
        blkh = CopyBits(blkh, iblk, ibits * 14, sbn);
//      iblk = ShiftRightHalf(iblk, ibits * 15);
      }
    } break;
  }
}

template<const int sets, const int ibits, const int abits, const int begin>
static void passreg WritePaletteBlock(int partition, Col4 (&idx)[2], Col4 &blkl, Col4 &blkh)
{
#define	fbits	(ibits < abits ? ibits : abits)	// index-block of smaller number of index-bits leads (first)
#define	sbits	(ibits > abits ? ibits : abits)	// index-block of larger number of index-bits follows (second)
  Col4 fblk, sblk;
  Col4 iidx = idx[0];
  Col4 aidx = idx[1];
  Col4 iblk;
  Col4 ablk;

  /* two of the cases straddle the lo to hi border
   *
   * WritePaletteBlock<1, 2,3, 50>(0, remapped, blkl, blkh);
   * WritePaletteBlock<1, 3,2, 50>(0, remapped, blkl, blkh);
   * WritePaletteBlock<1, 2,2, 66>(0, remapped, blkl, blkh);
   *
   * one set only
   */
  assume(partition >= 0 && partition <= 64);

  // max 16*3 -> 48bits
  iblk = CopyBits<ibits, ibits *  0>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  1>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  2>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  3>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  4>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  5>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  6>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  7>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  8>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits *  9>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 10>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 11>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 12>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 13>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 14>(iblk, iidx); iidx = ShiftRight<8>(iidx);
  iblk = CopyBits<ibits, ibits * 15>(iblk, iidx); //dx = ShiftRight<8>(iidx);

  // max 16*3 -> 48bits
  ablk = CopyBits<abits, abits *  0>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  1>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  2>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  3>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  4>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  5>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  6>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  7>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  8>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits *  9>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits * 10>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits * 11>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits * 12>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits * 13>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits * 14>(ablk, aidx); aidx = ShiftRight<8>(aidx);
  ablk = CopyBits<abits, abits * 15>(ablk, aidx); //dx = ShiftRight<8>(aidx);

  // index-block of smaller number of index-bits leads (first)
  // index-block of larger number of index-bits follows (second)
  fblk = (ibits == fbits ? iblk : ablk);
  sblk = (abits == sbits ? ablk : iblk);

  // one bit is omited per set, so it's one instruction per set + 1
  assume(begin == 50 || begin == 66);
  switch(begin) {
    case 50: {
      // always index 0
      blkl = CopyBits<fbits *  1 -  1, /*begin*/50 -  0 +         0>(blkl, fblk);
      fblk = ShiftRightHalf<fbits *  1>(fblk);
      blkl = CopyBits<fbits * 15 +  0, /*begin*/50 -  0 + fbits - 1>(blkl, fblk);

      // block-straddle
      fblk = ShiftRightHalf<(64 - /*begin*/50) + 1 - fbits>(fblk);
    //blkh = MaskBits<fbits * 15 - 15,          64 - 64            >(      fblk);
      blkh = fblk;
    } break;
    case 66: {
      // always index 0
      blkh = CopyBits<fbits *  1 -  1, /*begin*/66 - 64 +         0>(blkh, fblk);
      fblk = ShiftRightHalf<fbits *  1>(fblk);
      blkh = CopyBits<fbits * 15 +  0, /*begin*/66 - 64 + fbits - 1>(blkh, fblk);
//    fblk = ShiftRightHalf<fbits * 15>(fblk);
    } break;
  }

  // one bit is omited per set, so it's one instruction per set + 1
  assume(begin == 50 || begin == 66);
  switch(begin) {
    case 50:
    case 66: {
      // always index 0
      blkh = CopyBits<sbits *  1 - 1, 64 - sbits * 16 +         1>(blkh, sblk);
      sblk = ShiftRightHalf<sbits *  1>(sblk);
      blkh = CopyBits<sbits * 15 + 0, 64 - sbits * 16 + sbits - 0>(blkh, sblk);
//    sblk = ShiftRightHalf<sbits * 15>(sblk);
    } break;
  }
#undef	fbits
#undef	sbits
}

/* -----------------------------------------------------------------------------
 */
#define	ihibit	(1 << (ibits - 1))
#define	ihimsk	((1 << ibits) - 1)
#define	ihixor	Col4((ihimsk << 24) + (ihimsk << 16) + (ihimsk << 8) + (ihimsk << 0))

#define	ahibit	(1 << (abits - 1))
#define	ahimsk	((1 << abits) - 1)
#define	ahixor	Col4((ahimsk << 24) + (ahimsk << 16) + (ahimsk << 8) + (ahimsk << 0))

template<const int set>
static void passreg ExchangeBits(int &sharedbits) {
  const int setbit = (1 << set);
  int switched = ((sharedbits & setbit) << SBEND) + ((sharedbits >> SBEND) & setbit);
  sharedbits = (sharedbits & (~((setbit << SBEND) + (setbit << SBSTART)))) + switched;
}

template<const int ibits>
static void passreg RemapPaletteBlock(int partition, Vec4 (&start)[3], Vec4 (&end)[3], int &sharedbits, Col4 (&idxs)[1], u8 const (&indices)[1][16])
{
#if 0
  int masks[4], xors[4] = {0,0,0,0};

  PaletteSet::GetMasks(flags, masks);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapXYZW(end[0]); xors[0] = ihimsk; }
  /* set 2 */
  if (indices[0][shorterindex[partition][1]] & ihibit) {
    start[1].SwapXYZW(end[1]); xors[1] = ihimsk; }
  /* set 3 */
  if (indices[0][shorterindex[partition][2]] & ihibit) {
    start[2].SwapXYZW(end[2]); xors[2] = ihimsk; }

  /* TODO: 16 bytes fit into one SSE2 register */
  for (int i = 0; i < 16; ++i) {
    int set =
      (((masks[1] >> i) & 1) * 1) +
      (((masks[2] >> i) & 1) * 2);

    indices[0][i] = indices[0][i] ^ xors[set];
  }
#else
  Col4 xors = Col4(0);

  LoadAligned(idxs[0], indices[0]);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapXYZW(end[0]); ExchangeBits<0>(sharedbits); xors  = Col4(blockxor[partition][2]); }
  /* set 2 */
  if (indices[0][shorterindex[partition][1]] & ihibit) {
    start[1].SwapXYZW(end[1]); ExchangeBits<1>(sharedbits); xors |= Col4(blockxor[partition][3]); }
  /* set 3 */
  if (indices[0][shorterindex[partition][2]] & ihibit) {
    start[2].SwapXYZW(end[2]); ExchangeBits<2>(sharedbits); xors |= Col4(blockxor[partition][4]); }

  idxs[0] ^= (xors & ihixor);
#endif
}

template<const int ibits>
static void passreg RemapPaletteBlock(int partition, Vec4 (&start)[2], Vec4 (&end)[2], int &sharedbits, Col4 (&idxs)[1], u8 const (&indices)[1][16])
{
#if 0
  int masks[4], xors[4] = {0,0,0,0};

  PaletteSet::GetMasks(flags, masks);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapXYZW(end[0]); xors[0] = ihimsk; }
  /* set 2 */
  if (indices[0][shorterindex[partition][0]] & ihibit) {
    start[1].SwapXYZW(end[1]); xors[1] = ihimsk; }

  /* TODO: 16 bytes fit into one SSE2 register */
  for (int i = 0; i < 16; ++i) {
    int set =
      (((masks[1] >> i) & 1) * 1);

    indices[0][i] = indices[0][i] ^ xors[set];
  }
#else
  Col4 xors = Col4(0);

  LoadAligned(idxs[0], indices[0]);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapXYZW(end[0]); ExchangeBits<0>(sharedbits); xors  = Col4(blockxor[partition][0]); }
  /* set 2 */
  if (indices[0][shorterindex[partition][0]] & ihibit) {
    start[1].SwapXYZW(end[1]); ExchangeBits<1>(sharedbits); xors |= Col4(blockxor[partition][1]); }

  idxs[0] ^= (xors & ihixor);
#endif
}

template<const int ibits, const int abits>
static void passreg RemapPaletteBlock(int partition, Vec4 (&start)[1], Vec4 (&end)[1], int &sharedbits, Col4 (&idxs)[2], u8 const (&indices)[2][16])
{
#if 0
  int masks[4], xors[4] = {0,0,0,0};

//PaletteSet::GetMasks(flags, masks);

  /* same for all set 1s */
  if (indices[0][                                 0 ] & ihibit) {
    start[0].SwapXYZ(end[0]); xors[0] = ihibit - 1; }
  if (indices[1][                                 0 ] & ahibit) {
    start[0].SwapW  (end[0]); xors[1] = ahibit - 1; }

  /* TODO: 16 bytes fit into one SSE2 register */
  for (int i = 0; i < 16; ++i) {

    int iset = 0;
    int aset = 1;

    remapped[0][i] = indices[0][i] ^ xors[iset];
    remapped[1][i] = indices[1][i] ^ xors[aset];
  }
#else
  Col4 ixors = Col4(0);
  Col4 axors = Col4(0);

  LoadAligned(idxs[0], indices[0]);
  LoadAligned(idxs[1], indices[1]);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapXYZ (end[0]); ExchangeBits<0>(sharedbits); ixors  = Col4(0xFFFFFFFF); }
  /* same for all set 1s */
  if (indices[1][                        0 ] & ahibit) {
    start[0].SwapW   (end[0]); ExchangeBits<0>(sharedbits); axors  = Col4(0xFFFFFFFF); }

  idxs[0] ^= (ixors & ihixor);
  idxs[1] ^= (axors & ahixor);
#endif
}

template<const int ibits>
static void passreg RemapPaletteBlock(int partition, Vec4 (&start)[1], Vec4 (&end)[1], int &sharedbits, Col4 (&idxs)[1], u8 const (&indices)[1][16])
{
#if 0
  int masks[4], xors[4] = {0,0,0,0};
  u8 const ihibit = 1 << ibits;

//PaletteSet::GetMasks(flags, masks);

  /* same for all set 1s */
  if (indices[0][                         0 ] & ihibit) {
    start[0].SwapXYZW(end[0]); xors[0] = ihimsk; }

  /* TODO: 16 bytes fit into one SSE2 register */
  for (int i = 0; i < 16; ++i) {
    int set = 0;

    indices[0][i] = indices[0][i] ^ xors[set];
  }
#else
  Col4 xors = Col4(0);

  LoadAligned(idxs[0], indices[0]);

  /* same for all set 1s */
  if (indices[0][                        0 ] & ihibit) {
    start[0].SwapXYZW(end[0]); ExchangeBits<0>(sharedbits); xors  = Col4(0xFFFFFFFF); }

  idxs[0] ^= (xors & ihixor);
#endif
}

#undef	ihibit
#undef	ihimsk
#undef	ihixor

#undef	ahibit
#undef	ahimsk
#undef	ahixor

/* -----------------------------------------------------------------------------
 * Remarks
 *
 * Mode 8 (the least significant bit is set to 0x00) is reserved. Do not use it
 * in your encoder. If you pass this mode to the hardware, a block initialized
 * to all zeroes is returned.
 *
 * In BC7, you can encode the alpha component in one of the following ways:
 *
 *   - Block types without explicit alpha component encoding. In these blocks,
 *           the color endpoints have an RGB-only encoding, with the alpha component
 *           decoded to 1.0 for all texels.
 *   - Block types with combined color and alpha components. In these blocks, the
 *           endpoint color values are specified in the RGBA format, and the alpha
 *           component values are interpolated along with the color values.
 *   - Block types with separated color and alpha components. In these blocks the
 *           color and alpha values are specified separately, each with their own set
 *           of indices. As a result, they have an effective vector and a scalar
 *           channel separately encoded, where the vector commonly specifies the color
 *           channels [R, G, B] and the scalar specifies the alpha channel [A]. To
 *           support this approach, a separate 2-bit field is provided in the encoding,
 *           which permits the specification of the separate channel encoding as a
 *           scalar value. As a result, the block can have one of the following four
 *           different representations of this alpha encoding (as indicated by the 2-bit
 *           field):
 *       + RGB|A: alpha channel separate
 *       + AGB|R: "red" color channel separate
 *       + RAB|G: "green" color channel separate
 *       + RGA|B: "blue" color channel separate
 *   - The decoder reorders the channel order back to RGBA after decoding, so the
 *           internal block format is invisible to the developer. Blacks with separate
 *           color and alpha components also have two sets of index data: one for the
 *           vectored set of channels, and one for the scalar channel. (In the case of
 *           Mode 4, these indices are of differing widths [2 or 3 bits]. Mode 4 also
 *           contains a 1-bit selector that specifies whether the vector or the scalar
 *           channel uses the 3-bit indices.)
 *
 */
#define	C	COLORA
#define	U	UNIQUE
#define	S	SHARED

void WritePaletteBlock3_m1(int partition, Vec4 const (&start)[3], Vec4 const (&end)[3], int sharedbits, u8 const (&indices)[1][16], void* block)
{
  Vec4 s[3] = {start[0], start[1], start[2]};
  Vec4 e[3] = {end  [0], end  [1], end  [2]};
  Col4 a[3][FIELDN];
  Col4 b[3][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapPaletteBlock<3>(partition, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_LOWPRC)
  if ((FEATURE_SHAREDBITS_TRIALS == SHAREDBITS_TRIAL_ALL) || (~sharedbits))
  {
    FloatTo<4,4,4,0,1,0>(s, a, sharedbits >> SBSTART);
    FloatTo<4,4,4,0,1,0>(e, b, sharedbits >> SBEND);
  }
  else
#endif
  {
    FloatTo<4,4,4,0,1,0>(s, a);
    FloatTo<4,4,4,0,1,0>(e, b);
    FloatTo<4,4,4,0,1,0>(a, b);		// rounded shared bits
  }

  // 4 bits set 1/2/3 red/green/blue start/stop
  a[0][C] |= (b[0][C] <<=  4);
  a[1][C] |= (b[1][C] <<=  4);
  a[2][C] |= (b[2][C] <<=  4);
  // 8 bits set 1/2/3 red/green/blue start/stop
  a[0][C] |= (a[1][C] <<=  8);
  a[0][C] |= (a[2][C] <<= 16);

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 0 has the following characteristics:
   *
   *     Color components only (no alpha)
   *     3 subsets per block
   *     RGBP 4.4.4.1 endpoints with a unique P-bit per endpoint
   *     3-bit indices
   *     16 partitions
   *
   * [ 1 ]
   * [ 4 bit partition ]
   * [ 4 bits R0 ][ 4 bits R1 }[ 4 bits R2 ][ 4 bits R3 ][ 4 bits R4 ][ 4 bits R5 ]
   * [ 4 bits G0 ][ 4 bits G1 }[ 4 bits G2 ][ 4 bits G3 ][ 4 bits G4 ][ 4 bits G5 ]
   * [ 4 bits B0 ][ 4 bits B1 }[ 4 bits B2 ][ 4 bits B3 ][ 4 bits B4 ][ 4 bits B5 ]
   * [ 1 bit P0 ][ 1 bit P1 ][ 1 bit P2 ][ 1 bit P3 ][ 1 bit P4 ][ 1 bit P5 ]
   * [ 45 bits index ................................................................... ]
   */

  blkl = blkl.SetLong((partition << 1) + (1 << 0));		// 1 mode bit, 4 partition bits

  blkl = CopyBits<24,  5 -  0>(blkl,                a[0][C] );	// 24 bits set 1-3 red   start/stop
  blkl = CopyBits<24, 29 -  0>(blkl, ShiftRight<32>(a[0][C]));	// 24 bits set 1-3 green start/stop
  blkl = CopyBits<11, 53 -  0>(blkl, ShiftRight<64>(a[0][C]));	// 11 bits set 1-3 blue  start/stop
  blkh = MaskBits<13, 64 - 64>(      ShiftRight<75>(a[0][C]));	// 13 bits set 1-3 blue  start/stop

  blkh = CopyBits< 1, 77 - 64>(blkh,                a[0][U] );	// 1 bits set 1 unique start
  blkh = CopyBits< 1, 78 - 64>(blkh,                b[0][U] );	// 1 bits set 1 unique stop
  blkh = CopyBits< 1, 79 - 64>(blkh,                a[1][U] );	// 1 bits set 2 unique start
  blkh = CopyBits< 1, 80 - 64>(blkh,                b[1][U] );	// 1 bits set 2 unique stop
  blkh = CopyBits< 1, 81 - 64>(blkh,                a[2][U] );	// 1 bits set 3 unique start
  blkh = CopyBits< 1, 82 - 64>(blkh,                b[2][U] );	// 1 bits set 3 unique stop

  // 128 - 83 -> 45 index bits + 3 bit from 3 set start/end order -> 16 * 3bit
  WritePaletteBlock<3, 3, 83>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WritePaletteBlock3_m2(int partition, Vec4 const (&start)[2], Vec4 const (&end)[2], int sharedbits, u8 const (&indices)[1][16], void* block)
{
  Vec4 s[2] = {start[0], start[1]};
  Vec4 e[2] = {end  [0], end  [1]};
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapPaletteBlock<3>(partition, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_LOWPRC)
  if ((FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALL) || (~sharedbits))
  {
    FloatTo<6,6,6,0,0,1>(s, a, sharedbits >> SBSTART);
    FloatTo<6,6,6,0,0,1>(e, b, sharedbits >> SBSTART);
  }
  else
#endif
  {
    FloatTo<6,6,6,0,0,1>(s, a);
    FloatTo<6,6,6,0,0,1>(e, b);
    FloatTo<6,6,6,0,0,1>(a, b);		// rounded shared bits
  }

  // 6 bits set 1/2 red/green/blue start/stop
  a[0][C] |= (b[0][C] <<=  6);
  a[1][C] |= (b[1][C] <<=  6);
  // 12 bits set 1/2 red/green/blue start/stop
  a[0][C] |= (a[1][C] <<= 12);

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 1 has the following characteristics:
   *
   *    Color components only (no alpha)
   *    2 subsets per block
   *    RGBP 6.6.6.1 endpoints with a shared P-bit per subset)
   *    3-bit indices
   *    64 partitions
   *
   * [ 01 ]
   * [ 6 bit partition ]
   * [ 6 bits R0 ][ 6 bits R1 ][ 6 bits R2 ][ 6 bits R3 ]
   * [ 6 bits G0 ][ 6 bits G1 ][ 6 bits G2 ][ 6 bits G3 ]
   * [ 6 bits B0 ][ 6 bits B1 ][ 6 bits B2 ][ 6 bits B3 ]
   * [ 1 bit P0 ][ 1 bit P1 ]
   * [ 46 bits index ................................................................... ]
   */

  blkl = blkl.SetLong((partition << 2) + (1 << 1));		// 2 mode bit, 6 partition bits

  blkl = CopyBits<24,  8 -  0>(blkl,                a[0][C] );	// 24 bits set 1-2 red   start/stop
  blkl = CopyBits<24, 32 -  0>(blkl, ShiftRight<32>(a[0][C]));	// 24 bits set 1-2 green start/stop
  blkl = CopyBits< 8, 56 -  0>(blkl, ShiftRight<64>(a[0][C]));	//  8 bits set 1-2 blue  start/stop
  blkh = MaskBits<16, 64 - 64>(      ShiftRight<72>(a[0][C]));	// 16 bits set 1-2 blue  start/stop

  blkh = CopyBits< 1, 80 - 64>(blkh,                a[0][S] );	// 1 bits set 1 shared start/stop
  blkh = CopyBits< 1, 81 - 64>(blkh,                a[1][S] );	// 1 bits set 2 shared start/stop

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  WritePaletteBlock<2, 3, 82>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WritePaletteBlock3_m3(int partition, Vec4 const (&start)[3], Vec4 const (&end)[3], int sharedbits, u8 const (&indices)[1][16], void* block)
{
  Vec4 s[3] = {start[0], start[1], start[2]};
  Vec4 e[3] = {end  [0], end  [1], end  [2]};
  Col4 a[3][FIELDN];
  Col4 b[3][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapPaletteBlock<2>(partition, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALL)
  FloatTo<5,5,5,0,0,0>(s, a, ~0);
  FloatTo<5,5,5,0,0,0>(e, b, ~0);
#else
  FloatTo<5,5,5,0,0,0>(s, a);
  FloatTo<5,5,5,0,0,0>(e, b);
  FloatTo<5,5,5,0,0,0>(a, b);		// rounded shared bits
#endif

  // 5 bits set 1/2/3 red/green/blue start/stop
  a[0][C] |= (b[0][C] <<=  5);
  a[1][C] |= (b[1][C] <<=  5);
  a[2][C] |= (b[2][C] <<=  5);
  // 15 bits set 1/2/3 red/green/blue start/stop
  a[0][C] |= (a[1][C] <<= 10);
  a[0][C] |= (a[2][C] <<= 20);

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 2 has the following characteristics:
   *
   *    Color components only (no alpha)
   *    3 subsets per block
   *    RGB 5.5.5 endpoints
   *    2-bit indices
   *    64 partitions
   *
   * [ 1 ]
   * [ 6 bit partition ]
   * [ 5 bits R0 ][ 5 bits R1 ][ 5 bits R2 ][ 5 bits R3 ][ 5 bits R4 ][ 5 bits R5 ]
   * [ 5 bits G0 ][ 5 bits G1 ][ 5 bits G2 ][ 5 bits G3 ][ 5 bits G4 ][ 5 bits G5 ]
   * [ 5 bits B0 ][ 5 bits B1 ][ 5 bits B2 ][ 5 bits B3 ][ 5 bits B4 ][ 5 bits B5 ]
   * [ 29 bits index ................................................................... ]
   */

  blkl = blkl.SetLong((partition << 3) + (1 << 2));		// 3 mode bit, 6 partition bits

  blkl = CopyBits<30,  9 -  0>(blkl,                a[0][C] );	// 30 bits set 1-3 red   start/stop
  blkl = CopyBits<25, 39 -  0>(blkl, ShiftRight<32>(a[0][C]));	// 25 bits set 1-3 green start/stop
  blkh = MaskBits< 5, 64 - 64>(      ShiftRight<57>(a[0][C]));	//  5 bits set 1-3 green start/stop
  blkh = CopyBits<30, 69 - 64>(blkh, ShiftRight<64>(a[0][C]));	// 30 bits set 1-3 blue start/stop

  // 128 - 99 -> 29 index bits + 3 bit from 3 set start/end order -> 16 * 2bit
  WritePaletteBlock<3, 2, 99>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WritePaletteBlock3_m4(int partition, Vec4 const (&start)[2], Vec4 const (&end)[2], int sharedbits, u8 const (&indices)[1][16], void* block)
{
  Vec4 s[2] = {start[0], start[1]};
  Vec4 e[2] = {end  [0], end  [1]};
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapPaletteBlock<2>(partition, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALL)
  if ((FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALL) || (~sharedbits))
  {
    FloatTo<7,7,7,0,1,0>(s, a, sharedbits >> SBSTART);
    FloatTo<7,7,7,0,1,0>(e, b, sharedbits >> SBEND);
  }
  else
#endif
  {
    FloatTo<7,7,7,0,1,0>(s, a);
    FloatTo<7,7,7,0,1,0>(e, b);
    FloatTo<7,7,7,0,1,0>(a, b);		// rounded shared bits
  }

  // 7 bits set 1/2 red/green/blue start/stop
  a[0][C] |= (b[0][C] <<=  7);
  a[1][C] |= (b[1][C] <<=  7);
  // 14 bits set 1/2 red/green/blue start/stop
  a[0][C] |= (a[1][C] <<= 14);

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 3 has the following characteristics:
   *
   *    Color components only (no alpha)
   *    2 subsets per block
   *    RGBP 7.7.7.1 endpoints with a unique P-bit per subset)
   *    2-bit indices
   *    64 partitions
   *
   * [ 001 ]
   * [ 6 bit partition ]
   * [ 7 bits R0 ][ 7 bits R1 ][ 7 bits R2 ][ 7 bits R3 ]
   * [ 7 bits G0 ][ 7 bits G1 ][ 7 bits G2 ][ 7 bits G3 ]
   * [ 7 bits B0 ][ 7 bits B1 ][ 7 bits B2 ][ 7 bits B3 ]
   * [ 1 bit P0 ][ 1 bit P1 ][ 1 bit P2 ][ 1 bit P3 ]
   * [ 30 bits index ................................................................... ]
   */

  blkl = blkl.SetLong((partition << 4) + (1 << 3));		// 4 mode bit, 6 partition bits

  blkl = CopyBits<28, 10 -  0>(blkl,                a[0][C] );	// 28 bits set 1-2 red   start/stop
  blkl = CopyBits<26, 38 -  0>(blkl, ShiftRight<32>(a[0][C]));	// 26 bits set 1-2 green start/stop
  blkh = MaskBits< 2, 64 - 64>(      ShiftRight<58>(a[0][C]));	//  2 bits set 1-2 green start/stop
  blkh = CopyBits<28, 66 - 64>(blkh, ShiftRight<64>(a[0][C]));	// 28 bits set 1-2 blue  start/stop

  blkh = CopyBits< 1, 94 - 64>(blkh,                a[0][U] );	// 1 bits set 1 unique start
  blkh = CopyBits< 1, 95 - 64>(blkh,                b[0][U] );	// 1 bits set 1 unique stop
  blkh = CopyBits< 1, 96 - 64>(blkh,                a[1][U] );	// 1 bits set 2 unique start
  blkh = CopyBits< 1, 97 - 64>(blkh,                b[1][U] );	// 1 bits set 2 unique stop

  // 128 - 98 -> 30 index bits + 2 bit from 2 set start/end order -> 16 * 2bit
  WritePaletteBlock<2, 2, 98>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WritePaletteBlock4_m5(int r, int ix, Vec4 const (&start)[1], Vec4 const (&end)[1], int sharedbits, u8 const (&indices)[2][16], void* block)
{
  Vec4 s[1] = {start[0]};
  Vec4 e[1] = {end  [0]};
  Col4 a[1][FIELDN];
  Col4 b[1][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[2];

  // remap the indices
  if (!ix)
    RemapPaletteBlock<2,3>(0, s, e, sharedbits, idxs, indices);
  else
    RemapPaletteBlock<3,2>(0, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALPHAONLY)
  FloatTo<5,5,5,6,0,0>(s, a, ~0);
  FloatTo<5,5,5,6,0,0>(e, b, ~0);
#else
  FloatTo<5,5,5,6,0,0>(s, a);
  FloatTo<5,5,5,6,0,0>(e, b);
  FloatTo<5,5,5,6,0,0>(a, b);		// rounded shared bits
#endif

  // 5 bits set 1 red/green/blue start/stop
  // 6 bits set 1 alpha start/stop
  a[0][C] |= (b[0][C] = ShiftLeftLo<5,5,5,6>(b[0][C]));

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 4 has the following characteristics:
   *
   *    Color components with separate alpha component
   *    1 subset per block
   *    RGB 5.5.5 color endpoints
   *    6-bit alpha endpoints
   *    16 x 2-bit indices
   *    16 x 3-bit indices
   *    2-bit component rotation
   *    1-bit index selector (whether the 2- or 3-bit indices are used)
   *
   * [ 00001 ]
   * [ 2 bit rotation ]
   * [ 1 bit idxMode ]
   * [ 5 bits R0 ][ 5 bits R1 ]
   * [ 5 bits G0 ][ 5 bits G1 ]
   * [ 5 bits B0 ][ 5 bits B1 ]
   * [  6 bits A0  ][  6 bits A1  ]
   * [ 31 bits index ....................................................... ]
   * [ 47 bits index ................................................................... ]
   */

  blkl = blkl.SetLong((ix << 7) + (r << 5) + (1 << 4));		// 5 mode bit, 2 rotation bits, 1 index bit

  blkl = CopyBits<10,  8 -  0>(blkl,                a[0][C] );	// 10 bits set 1 red   start/stop
  blkl = CopyBits<10, 18 -  0>(blkl, ShiftRight<32>(a[0][C]));	// 10 bits set 1 green start/stop
  blkl = CopyBits<10, 28 -  0>(blkl, ShiftRight<64>(a[0][C]));	// 10 bits set 1 blue  start/stop
  blkl = CopyBits<12, 38 -  0>(blkl, ShiftRight<96>(a[0][C]));	// 12 bits set 1 alpha start/stop

  // 128 - 50 -> 78 index bits -> 31 + 47 index bits + 2 bit from 2 set start/end order -> 16 * 2bit + 16 * 3bit
  if (!ix)
    WritePaletteBlock<1, 2,3, 50>(0, idxs, blkl, blkh);
  else
    WritePaletteBlock<1, 3,2, 50>(0, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WritePaletteBlock4_m6(int rotation, Vec4 const (&start)[1], Vec4 const (&end)[1], int sharedbits, u8 const (&indices)[2][16], void* block)
{
  Vec4 s[1] = {start[0]};
  Vec4 e[1] = {end  [0]};
  Col4 a[1][FIELDN];
  Col4 b[1][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[2];

  // remap the indices
  RemapPaletteBlock<2,2>(0, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALPHAONLY)
  FloatTo<7,7,7,8,0,0>(s, a, ~0);
  FloatTo<7,7,7,8,0,0>(e, b, ~0);
#else
  FloatTo<7,7,7,8,0,0>(s, a);
  FloatTo<7,7,7,8,0,0>(e, b);
  FloatTo<7,7,7,8,0,0>(a, b);		// rounded shared bits
#endif

  // 7 bits set 1 red/green/blue start/stop
  // 8 bits set 1 alpha start/stop
  a[0][C] |= (b[0][C] = ShiftLeftLo<7,7,7,8>(b[0][C]));

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 5 has the following characteristics:
   *
   *    Color components with separate alpha component
   *    1 subset per block
   *    RGB 7.7.7 color endpoints
   *    6-bit alpha endpoints
   *    16 x 2-bit color indices
   *    16 x 2-bit alpha indices
   *    2-bit component rotation
   *
   * [ 000001 ]
   * [ 2 bit rotation ]
   * [ 7 bits R0 ][ 7 bits R1 ]
   * [ 7 bits G0 ][ 7 bits G1 ]
   * [ 7 bits B0 ][ 7 bits B1 ]
   * [  8 bits A0  ][  8 bits A1  ]
   * [ 31 bits index ................................................................... ]
   * [ 31 bits index ................................................................... ]
   */

  blkl = blkl.SetLong((rotation << 6) + (1 << 5));		// 6 mode bit, 2 rotation bits

  blkl = CopyBits<14,  8 -  0>(blkl,                 a[0][C] );	// 14 bits set 1 red   start/stop
  blkl = CopyBits<14, 22 -  0>(blkl, ShiftRight< 32>(a[0][C]));	// 14 bits set 1 green start/stop
  blkl = CopyBits<14, 36 -  0>(blkl, ShiftRight< 64>(a[0][C]));	// 14 bits set 1 blue  start/stop
  blkl = CopyBits<14, 50 -  0>(blkl, ShiftRight< 96>(a[0][C]));	// 14 bits set 1 alpha start/stop
  blkh = MaskBits< 2, 64 - 64>(      ShiftRight<110>(a[0][C]));	//  2 bits set 1 alpha start/stop

  // 128 - 66 -> 62 index bits -> 31 + 31 index bits + 2 bit from 2 set start/end order -> 16 * 2bit + 16 * 3bit
  WritePaletteBlock<1, 2,2, 66>(0, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WritePaletteBlock4_m7(int partition, Vec4 const (&start)[1], Vec4 const (&end)[1], int sharedbits, u8 const (&indices)[1][16], void* block)
{
  Vec4 s[1] = {start[0]};
  Vec4 e[1] = {end  [0]};
  Col4 a[1][FIELDN];
  Col4 b[1][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapPaletteBlock<4>(partition, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALPHAONLYOPAQUE)
  if ((FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALL) || (~sharedbits)) {
    FloatTo<7,7,7,7,1,0>(s, a, sharedbits >> SBSTART);
    FloatTo<7,7,7,7,1,0>(e, b, sharedbits >> SBEND);
  }
  else
#endif
  {
    FloatTo<7,7,7,7,1,0>(s, a);
    FloatTo<7,7,7,7,1,0>(e, b);
    FloatTo<7,7,7,7,1,0>(a, b);		// rounded shared bits
  }

  // 7 bits set 1 red/green/blue/alpha start/stop
  a[0][C] |= (b[0][C] <<= 7);

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 6 has the following characteristics:
   *
   *    Combined color and alpha components
   *    One subset per block
   *    RGBAP 7.7.7.7.1 color (and alpha) endpoints (unique P-bit per endpoint)
   *    16 x 4-bit indices
   *
   * [ 0000001 ]
   * [ 4 bit partition ]
   * [ 7 bits R0 ][ 7 bits R1 ]
   * [ 7 bits G0 ][ 7 bits G1 ]
   * [ 7 bits B0 ][ 7 bits B1 ]
   * [ 7 bits A0 ][ 7 bits A1 ]
   * [ 1 bit P0 ][ 1 bit P1 ]
   * [ 63 bits index ................................................................... ]
   */

//blkl = blkl.SetLong((partition << 7) + (1 << 6));		// 7 mode bit, 0 partition bits
  blkl = Col4(1 << 6, 0, 0, 0);					// 7 mode bit, 0 partition bits

  blkl = CopyBits<14,  7 -  0>(blkl,                a[0][C] );	// 14 bits set 1 red   start/stop
  blkl = CopyBits<14, 21 -  0>(blkl, ShiftRight<32>(a[0][C]));	// 14 bits set 1 green start/stop
  blkl = CopyBits<14, 35 -  0>(blkl, ShiftRight<64>(a[0][C]));	// 14 bits set 1 blue  start/stop
  blkl = CopyBits<14, 49 -  0>(blkl, ShiftRight<96>(a[0][C]));	// 14 bits set 1 alpha start/stop

  blkl = CopyBits< 1, 63 -  0>(blkl,                a[0][U] );	// 1 bits set 1 unique start
  blkh = MaskBits< 1, 64 - 64>(                     b[0][U] );	// 1 bits set 1 unique stop

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  WritePaletteBlock<1, 4, 65>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

void WritePaletteBlock4_m8(int partition, Vec4 const (&start)[2], Vec4 const (&end)[2], int sharedbits, u8 const (&indices)[1][16], void* block)
{
  Vec4 s[2] = {start[0], start[1]};
  Vec4 e[2] = {end  [0], end  [1]};
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;
  Col4 idxs[1];

  // remap the indices
  RemapPaletteBlock<2>(partition, s, e, sharedbits, idxs, indices);

  // get the packed values
#if	defined(FEATURE_SHAREDBITS_TRIALS) && (FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALPHAONLYOPAQUE)
  if ((FEATURE_SHAREDBITS_TRIALS >= SHAREDBITS_TRIAL_ALL) || (~sharedbits)) {
    FloatTo<5,5,5,5,1,0>(s, a, sharedbits >> SBSTART);
    FloatTo<5,5,5,5,1,0>(e, b, sharedbits >> SBEND);
  }
  else
#endif
  {
    FloatTo<5,5,5,5,1,0>(s, a);
    FloatTo<5,5,5,5,1,0>(e, b);
    FloatTo<5,5,5,5,1,0>(a, b);		// rounded shared bits
  }

  // 5 bits set 1/2 red/green/blue/alpha start/stop
  a[0][C] |= (b[0][C] <<=  5);
  a[1][C] |= (b[1][C] <<=  5);
  // 10 bits set 1/2 red/green/blue/alpha start/stop
  a[0][C] |= (a[1][C] <<= 10);

  /* ..............................................................................................
   * http://msdn.microsoft.com/en-us/library/hh308954%28v=vs.85%29
   *
   * BC7 Mode 7 has the following characteristics:
   *
   *    Combined color and alpha components
   *    2 subsets per block
   *    RGBAP 5.5.5.5.1 color (and alpha) endpoints (unique P-bit per endpoint)
   *    2-bit indices
   *    64 partitions
   *
   * [ 00000001 ]
   * [ 6 bit partition ]
   * [ 5 bits R0 ][ 5 bits R1 ][ 5 bits R2 ][ 5 bits R3 ]
   * [ 5 bits G0 ][ 5 bits G1 ][ 5 bits G2 ][ 5 bits G3 ]
   * [ 5 bits B0 ][ 5 bits B1 ][ 5 bits B2 ][ 5 bits B3 ]
   * [ 5 bits A0 ][ 5 bits A1 ][ 5 bits A2 ][ 5 bits A3 ]
   * [ 1 bit P0 ][ 1 bit P1 ][ 1 bit P2 ][ 1 bit P3 ]
   * [ 30 bits index ................................................................... ]
   */

  blkl = blkl.SetLong((partition << 8) + (1 << 7));		// 8 mode bit, 6 partition bits

  blkl = CopyBits<20, 14 -  0>(blkl,                a[0][C] );	// 20 bits set 1-2 red   start/stop
  blkl = CopyBits<20, 34 -  0>(blkl, ShiftRight<32>(a[0][C]));	// 20 bits set 1-2 green start/stop
  blkl = CopyBits<10, 54 -  0>(blkl, ShiftRight<64>(a[0][C]));	// 10 bits set 1-2 blue  start/stop
  blkh = MaskBits<10, 64 - 64>(      ShiftRight<74>(a[0][C]));	// 10 bits set 1-2 blue  start/stop
  blkh = CopyBits<20, 74 - 64>(blkh, ShiftRight<96>(a[0][C]));	// 20 bits set 1-2 alpha start/stop

  blkh = CopyBits< 1, 94 - 64>(blkh,                a[0][U] );	// 1 bits set 1 unique start
  blkh = CopyBits< 1, 95 - 64>(blkh,                b[0][U] );	// 1 bits set 1 unique stop
  blkh = CopyBits< 1, 96 - 64>(blkh,                a[1][U] );	// 1 bits set 2 unique start
  blkh = CopyBits< 1, 97 - 64>(blkh,                b[1][U] );	// 1 bits set 2 unique stop

  // 128 - 98 -> 30 index bits + 2 bit from 2 set start/end order -> 16 * 2bit
  WritePaletteBlock<2, 2, 98>(partition, idxs, blkl, blkh);

  /* write out */
  StoreUnaligned(blkl, blkh, block);
}

/* *****************************************************************************
 */
extern const u8 whichsetinpartition[64][/*3*/2][16] = {
{ /*{0},*/ {0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1},{0,0,1,1,0,0,1,1,0,2,2,1,2,2,2,2}},
{ /*{0},*/ {0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1},{0,0,0,1,0,0,1,1,2,2,1,1,2,2,2,1}},
{ /*{0},*/ {0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1},{0,0,0,0,2,0,0,1,2,2,1,1,2,2,1,1}},
{ /*{0},*/ {0,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1},{0,2,2,2,0,0,2,2,0,0,1,1,0,1,1,1}},
{ /*{0},*/ {0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1},{0,0,0,0,0,0,0,0,1,1,2,2,1,1,2,2}},
{ /*{0},*/ {0,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1},{0,0,1,1,0,0,1,1,0,0,2,2,0,0,2,2}},
{ /*{0},*/ {0,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1},{0,0,2,2,0,0,2,2,1,1,1,1,1,1,1,1}},
{ /*{0},*/ {0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,1},{0,0,1,1,0,0,1,1,2,2,1,1,2,2,1,1}},
{ /*{0},*/ {0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1},{0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2}},
{ /*{0},*/ {0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1},{0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2}},
{ /*{0},*/ {0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1},{0,0,0,0,1,1,1,1,2,2,2,2,2,2,2,2}},
{ /*{0},*/ {0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1},{0,0,1,2,0,0,1,2,0,0,1,2,0,0,1,2}},
{ /*{0},*/ {0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1},{0,1,1,2,0,1,1,2,0,1,1,2,0,1,1,2}},
{ /*{0},*/ {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1},{0,1,2,2,0,1,2,2,0,1,2,2,0,1,2,2}},
{ /*{0},*/ {0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1},{0,0,1,1,0,1,1,2,1,1,2,2,1,2,2,2}},
{ /*{0},*/ {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1},{0,0,1,1,2,0,0,1,2,2,0,0,2,2,2,0}},

{ /*{0},*/ {0,0,0,0,1,0,0,0,1,1,1,0,1,1,1,1},{0,0,0,1,0,0,1,1,0,1,1,2,1,1,2,2}},
{ /*{0},*/ {0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0},{0,1,1,1,0,0,1,1,2,0,0,1,2,2,0,0}},
{ /*{0},*/ {0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0},{0,0,0,0,1,1,2,2,1,1,2,2,1,1,2,2}},
{ /*{0},*/ {0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0},{0,0,2,2,0,0,2,2,0,0,2,2,1,1,1,1}},
{ /*{0},*/ {0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0},{0,1,1,1,0,1,1,1,0,2,2,2,0,2,2,2}},
{ /*{0},*/ {0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,0},{0,0,0,1,0,0,0,1,2,2,2,1,2,2,2,1}},
{ /*{0},*/ {0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0},{0,0,0,0,0,0,1,1,0,1,2,2,0,1,2,2}},
{ /*{0},*/ {0,1,1,1,0,0,1,1,0,0,1,1,0,0,0,1},{0,0,0,0,1,1,0,0,2,2,1,0,2,2,1,0}},
{ /*{0},*/ {0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0},{0,1,2,2,0,1,2,2,0,0,1,1,0,0,0,0}},
{ /*{0},*/ {0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0},{0,0,1,2,0,0,1,2,1,1,2,2,2,2,2,2}},
{ /*{0},*/ {0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0},{0,1,1,0,1,2,2,1,1,2,2,1,0,1,1,0}},
{ /*{0},*/ {0,0,1,1,0,1,1,0,0,1,1,0,1,1,0,0},{0,0,0,0,0,1,1,0,1,2,2,1,1,2,2,1}},
{ /*{0},*/ {0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0},{0,0,2,2,1,1,0,2,1,1,0,2,0,0,2,2}},
{ /*{0},*/ {0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0},{0,1,1,0,0,1,1,0,2,0,0,2,2,2,2,2}},
{ /*{0},*/ {0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0},{0,0,1,1,0,1,2,2,0,1,2,2,0,0,1,1}},
{ /*{0},*/ {0,0,1,1,1,0,0,1,1,0,0,1,1,1,0,0},{0,0,0,0,2,0,0,0,2,2,1,1,2,2,2,1}},

{ /*{0},*/ {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1},{0,0,0,0,0,0,0,2,1,1,2,2,1,2,2,2}},
{ /*{0},*/ {0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1},{0,2,2,2,0,0,2,2,0,0,1,2,0,0,1,1}},
{ /*{0},*/ {0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0},{0,0,1,1,0,0,1,2,0,0,2,2,0,2,2,2}},
{ /*{0},*/ {0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0},{0,1,2,0,0,1,2,0,0,1,2,0,0,1,2,0}},
{ /*{0},*/ {0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0},{0,0,0,0,1,1,1,1,2,2,2,2,0,0,0,0}},
{ /*{0},*/ {0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0},{0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0}},
{ /*{0},*/ {0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1},{0,1,2,0,2,0,1,2,1,2,0,1,0,1,2,0}},
{ /*{0},*/ {0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1},{0,0,1,1,2,2,0,0,1,1,2,2,0,0,1,1}},
{ /*{0},*/ {0,1,1,1,0,0,1,1,1,1,0,0,1,1,1,0},{0,0,1,1,1,1,2,2,2,2,0,0,0,0,1,1}},
{ /*{0},*/ {0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0},{0,1,0,1,0,1,0,1,2,2,2,2,2,2,2,2}},
{ /*{0},*/ {0,0,1,1,0,0,1,0,0,1,0,0,1,1,0,0},{0,0,0,0,0,0,0,0,2,1,2,1,2,1,2,1}},
{ /*{0},*/ {0,0,1,1,1,0,1,1,1,1,0,1,1,1,0,0},{0,0,2,2,1,1,2,2,0,0,2,2,1,1,2,2}},
{ /*{0},*/ {0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0},{0,0,2,2,0,0,1,1,0,0,2,2,0,0,1,1}},
{ /*{0},*/ {0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1},{0,2,2,0,1,2,2,1,0,2,2,0,1,2,2,1}},
{ /*{0},*/ {0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1},{0,1,0,1,2,2,2,2,2,2,2,2,0,1,0,1}},
{ /*{0},*/ {0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0},{0,0,0,0,2,1,2,1,2,1,2,1,2,1,2,1}},

{ /*{0},*/ {0,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0},{0,1,0,1,0,1,0,1,0,1,0,1,2,2,2,2}},
{ /*{0},*/ {0,0,1,0,0,1,1,1,0,0,1,0,0,0,0,0},{0,2,2,2,0,1,1,1,0,2,2,2,0,1,1,1}},
{ /*{0},*/ {0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,0},{0,0,0,2,1,1,1,2,0,0,0,2,1,1,1,2}},
{ /*{0},*/ {0,0,0,0,0,1,0,0,1,1,1,0,0,1,0,0},{0,0,0,0,2,1,1,2,2,1,1,2,2,1,1,2}},
{ /*{0},*/ {0,1,1,0,1,1,0,0,1,0,0,1,0,0,1,1},{0,2,2,2,0,1,1,1,0,1,1,1,0,2,2,2}},
{ /*{0},*/ {0,0,1,1,0,1,1,0,1,1,0,0,1,0,0,1},{0,0,0,2,1,1,1,2,1,1,1,2,0,0,0,2}},
{ /*{0},*/ {0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0},{0,1,1,0,0,1,1,0,0,1,1,0,2,2,2,2}},
{ /*{0},*/ {0,0,1,1,1,0,0,1,1,1,0,0,0,1,1,0},{0,0,0,0,0,0,0,0,2,1,1,2,2,1,1,2}},
{ /*{0},*/ {0,1,1,0,1,1,0,0,1,1,0,0,1,0,0,1},{0,1,1,0,0,1,1,0,2,2,2,2,2,2,2,2}},
{ /*{0},*/ {0,1,1,0,0,0,1,1,0,0,1,1,1,0,0,1},{0,0,2,2,0,0,1,1,0,0,1,1,0,0,2,2}},
{ /*{0},*/ {0,1,1,1,1,1,1,0,1,0,0,0,0,0,0,1},{0,0,2,2,1,1,2,2,1,1,2,2,0,0,2,2}},
{ /*{0},*/ {0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1},{0,0,0,0,0,0,0,0,0,0,0,0,2,1,1,2}},
{ /*{0},*/ {0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1},{0,0,0,2,0,0,0,1,0,0,0,2,0,0,0,1}},
{ /*{0},*/ {0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0},{0,2,2,2,1,2,2,2,0,2,2,2,1,2,2,2}},
{ /*{0},*/ {0,0,1,0,0,0,1,0,1,1,1,0,1,1,1,0},{0,1,0,1,2,2,2,2,2,2,2,2,2,2,2,2}},
{ /*{0},*/ {0,1,0,0,0,1,0,0,0,1,1,1,0,1,1,1},{0,1,1,1,2,0,1,1,2,2,0,1,2,2,2,0}}
};

/* -----------------------------------------------------------------------------
 */
template<const int sets, const int ibits, const int begin>
void ReadPaletteBlock(int partition, unsigned int *codes, Col4 &blkl, Col4 &blkh, int *out)
{
  Col4 iblk;

  /* none of the cases straddles the lo to hi border, all go into hi
   *
   * WritePaletteBlock<3, 3, 83>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<2, 3, 82>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<3, 2, 99>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<2, 2, 98>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<1, 4, 65>(partition, remapped, blkl, blkh);
   * WritePaletteBlock<2, 2, 98>(partition, remapped, blkl, blkh);
   */

  // one bit is omited per set, so it's one instruction per set + 1
  assume(sets >= 1 && sets <= 3);
  switch(sets) {
    case 1: {
      // always index 0
      iblk = ShiftLeftHalf<ibits>(
        ExtrBits<ibits * 15 - 0, begin - 64 + ibits *  1 - 1>(blkh)) |
        ExtrBits<ibits *  1 - 1, begin - 64 +              0>(blkh)
      ;
    } break;
    case 2: {
      // always index 0
      iblk = ShiftLeftHalf<ibits>(
        ExtrBits<ibits * 15 - 0, begin - 64 + ibits *  1 - 1>(blkh)) |
        ExtrBits<ibits *  1 - 1, begin - 64 +              0>(blkh)
      ;

      // if it's not the last bit which is cut (low probability)
      int len, bgn = begin - 64;
      if ((len = shorterindex[partition][0]) < 15) {
        len = (len * ibits) + ibits;
        bgn = (bgn + len - 2);

        /* no obvious mask-pattern, set of 9 distinct lens
         *
         * bgn=83-64 len={2,6,8}*3={6,18,24}+b={ 9,21,27}-1={ 8,20,26};
         * bgn=82-64 len={2,6,8}*3={6,18,24}+b={ 9,21,27}-1={ 8,20,26};
         * bgn=99-64 len={2,6,8}*2={4,12,16}+b={ 6,14,18}-1={ 5,13,17};
         * bgn=98-64 len={2,6,8}*2={4,12,16}+b={ 6,14,18}-1={ 5,13,17};
         * bgn=65-64 len={2,6,8}*4={8,24,32}+b={12,32,40}-1={11,31,39};
         * bgn=98-64 len={2,6,8}*2={4,12,16}+b={ 6,14,18}-1={ 5,13,17};
         */

      	// remaining length can be anything, length overflow is silent
        blkh = ExtrBits(blkh, ibits * 14, bgn);

        // do the whole shift
	iblk = ShiftLeftHalf(blkh, len) |
	  MaskBits(iblk, len - 1, 0)
	;
      }
    } break;
    case 3: {
      iblk = ShiftLeftHalf<ibits>(
        ExtrBits<ibits * 15 - 0, begin - 64 + ibits *  1 - 1>(blkh)) |
        ExtrBits<ibits *  1 - 1, begin - 64 +              0>(blkh)
      ;

      int bgn = begin - 64;
      int fln, fbn;
      fln = shorterindex[partition][3]; {
        fln = (fln * ibits) + ibits;
        fbn = (bgn + fln - 2);

      	// remaining length can be anything, length overflow is silent
        blkh = ExtrBits(blkh, ibits * 14, fbn);

        // do the whole shift
	iblk = ShiftLeftHalf(blkh, fln) |
	  MaskBits(iblk, fln - 1, 0)
	;
      }

      // if it's not the last bit which is cut (low probability)
      int sln, sbn;
      if ((sln = shorterindex[partition][4]) < 15) {
        sln = (sln * ibits) + ibits;
        sbn = (bgn + sln - 3);

      	// subtract the already conducted shifts
        // remaining length can be anything, length overflow is silent
        blkh = ExtrBits(blkh, ibits * 14, sbn - fbn);

        // do the whole shift
	iblk = ShiftLeftHalf(blkh, sln) |
	  MaskBits(iblk, sln - 1, 0)
	;
      }
    } break;
  }

  // store out the palettes
  u8 const *ios = (sets >= 2 ? whichsetinpartition[partition][sets - 2] : NULL);

  // max 16*4 -> 64bits
  // no scattered writes on SSE
  // movd: 3clk 1/1, extrq: 2clk 1/2
  Col4 c0, c1, c2, c3,
       c4, c5, c6, c7;

  // m1:
  //  c00 = 0 -> [0][0] -> 95 cd e1 ff, c08 = 6 -> [0][6] -> 80 c3 d0 ff
  //  c01 = 4 -> [0][4] -> 87 c6 d5 ff, c09 = 4 -> [1][4] -> 79 c3 cb ff
  //  c02 = 7 -> [0][7] -> 7c c1 cd ff, c10 = 2 -> [1][2] -> 80 c3 cb ff
  //  c03 = 7 -> [0][7] -> 7c c1 cd ff, c11 = 0 -> [1][0] -> 87 c3 cb ff
  //  c04 = 2 -> [0][2] -> 8e ca db ff, c12 = 3 -> [1][3] -> 7c c3 cb ff
  //  c05 = 7 -> [0][7] -> 7c c1 cd ff, c13 = 4 -> [1][4] -> 79 c3 cb ff
  //  c06 = 6 -> [0][6] -> 80 c3 d0 ff, c14 = 3 -> [1][3] -> 7c c3 cb ff
  //  c07 = 0 -> [1][0] -> 87 c3 cb ff, c15 = 2 -> [1][2] -> 80 c3 cb ff

  // m3:
  //  c00 = 0 -> [0][0] -> 92 c8 dc ff, c08 = 3 -> [1][3] -> 8c c8 dc ff
  //  c01 = 2 -> [1][2] -> 8b c6 d9 ff, c09 = 3 -> [1][3] -> 8c c8 dc ff
  //  c02 = 0 -> [1][0] -> 8a c2 d2 ff, c10 = 2 -> [0][2] -> 87 c4 d4 ff
  //  c03 = 1 -> [1][1] -> 8b c4 d5 ff, c11 = 3 -> [0][3] -> 82 c2 d0 ff
  //  c04 = 1 -> [0][1] -> 8d c6 d8 ff, c12 = 3 -> [1][3] -> 8c c8 dc ff
  //  c05 = 1 -> [0][1] -> 8d c6 d8 ff, c13 = 3 -> [1][3] -> 8c c8 dc ff
  //  c06 = 1 -> [1][1] -> 8b c4 d5 ff, c14 = 3 -> [1][3] -> 8c c8 dc ff
  //  c07 = 0 -> [1][0] -> 8a c2 d2 ff, c15 = 3 -> [0][3] -> 82 c2 d0 ff

  // m7:
  //  c00 = 6 -> [-][6] -> cb f1 fa ff, c08 = b -> [-][b] -> d4 f4 fb fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c01 = 5 -> [-][5] -> c9 f0 fa ff, c09 = 8 -> [-][8] -> ce f2 fb fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c02 = 0 -> [-][0] -> bf ed f9 ff, c10 = 3 -> [-][3] -> c5 ef fa ff
  //  c03 = 1 -> [-][1] -> c1 ee f9 ff, c11 = 0 -> [-][0] -> bf ed f9 ff
  //  c04 = 7 -> [-][7] -> cd f1 fa ff, c12 = e -> [-][e] -> da f5 fc fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c05 = 5 -> [-][5] -> c9 f0 fa ff, c13 = b -> [-][b] -> d4 f4 fb fe DirextXTex seems to be buggy, no transparent alpha allowed!
  //  c06 = 4 -> [-][4] -> c7 ef fa ff, c14 = 6 -> [-][6] -> cb f1 fa ff
  //  c07 = 2 -> [-][2] -> c3 ee f9 ff, c15 = 0 -> [-][0] -> bf ed f9 ff

  // m0:
  //  c00 = 3 -> [0][3] -> 8c d7 e1 ff, c08 = 5 -> [0][5] -> 8c d3 d7 ff
  //  c01 = 3 -> [0][3] -> 8c d7 e1 ff, c09 = 6 -> [0][6] -> 8c d0 d3 ff
  //  c02 = 4 -> [1][4] -> 84 c6 d9 ff, c10 = 0 -> [1][0] -> 84 c6 c6 ff
  //  c03 = 3 -> [2][3] -> 95 c9 d0 ff, c11 = 0 -> [2][0] -> 84 c6 c6 ff
  //  c04 = 4 -> [0][4] -> 8c d5 dc ff, c12 = 6 -> [0][6] -> 8c d0 d3 ff
  //  c05 = 4 -> [0][4] -> 8c d5 dc ff, c13 = 7 -> [0][7] -> 8c ce ce ff
  //  c06 = 2 -> [1][2] -> 84 c6 cf ff, c14 = 0 -> [1][0] -> 84 c6 c6 ff
  //  c07 = 0 -> [2][0] -> 84 c6 c6 ff, c15 = 1 -> [2][1] -> 8a c7 c9 ff

  // m3:
  //  c00 = 0 -> [0][0] -> 9c b5 94 ff, c08 = 1 -> [0][1] -> a4 b8 91 ff
  //  c01 = 1 -> [0][1] -> a4 b8 91 ff, c09 = 3 -> [0][3] -> b5 bd 8c ff
  //  c02 = 3 -> [2][3] -> ad b5 8c ff, c10 = 3 -> [2][3] -> ad b5 8c ff
  //  c03 = 1 -> [2][1] -> b2 c0 9d ff, c11 = 2 -> [2][2] -> b0 bb 94 ff
  //  c04 = 0 -> [0][0] -> 9c b5 94 ff, c12 = 2 -> [1][2] -> aa b2 94 ff
  //  c05 = 1 -> [0][1] -> a4 b8 91 ff, c13 = 2 -> [1][2] -> aa b2 94 ff
  //  c06 = 2 -> [2][2] -> b0 bb 94 ff, c14 = 0 -> [1][0] -> 94 ad 94 ff
  //  c07 = 2 -> [2][2] -> b0 bb 94 ff, c15 = 1 -> [1][1] -> 9f b0 94 ff

  // m8:
  //  c00 = 1 -> [0][1] -> c3 3c 10 ee, c08 = 3 -> [0][3] -> 51 00 00 d3
  //  c01 = 0 -> [0][0] -> fb 59 18 fb, c09 = 2 -> [1][2] -> 7b 28 0d aa
  //  c02 = 0 -> [0][0] -> fb 59 18 fb, c10 = 1 -> [1][1] -> b5 4d 17 cf
  //  c03 = 1 -> [0][1] -> c3 3c 10 ee, c11 = 0 -> [1][0] -> eb 71 20 f3
  //  c04 = 2 -> [0][2] -> 89 1d 08 e0, c12 = 3 -> [1][3] -> 45 04 04 86
  //  c05 = 1 -> [0][1] -> c3 3c 10 ee, c13 = 3 -> [1][3] -> 45 04 04 86
  //  c06 = 0 -> [0][0] -> fb 59 18 fb, c14 = 3 -> [1][3] -> 45 04 04 86
  //  c07 = 0 -> [1][0] -> eb 71 20 f3, c15 = 1 -> [1][1] -> b5 4d 17 cf

  c0 = ExtrBits<ibits, ibits *  0>(iblk);
  c1 = ExtrBits<ibits, ibits *  1>(iblk);
  c2 = ExtrBits<ibits, ibits *  2>(iblk);
  c3 = ExtrBits<ibits, ibits *  3>(iblk);
  c4 = ExtrBits<ibits, ibits *  4>(iblk);
  c5 = ExtrBits<ibits, ibits *  5>(iblk);
  c6 = ExtrBits<ibits, ibits *  6>(iblk);
  c7 = ExtrBits<ibits, ibits *  7>(iblk);

  out[ 0] = codes[(sets >= 2 ? ios[ 0] << ibits : 0) + c0.GetLong()];
  out[ 1] = codes[(sets >= 2 ? ios[ 1] << ibits : 0) + c1.GetLong()];
  out[ 2] = codes[(sets >= 2 ? ios[ 2] << ibits : 0) + c2.GetLong()];
  out[ 3] = codes[(sets >= 2 ? ios[ 3] << ibits : 0) + c3.GetLong()];
  out[ 4] = codes[(sets >= 2 ? ios[ 4] << ibits : 0) + c4.GetLong()];
  out[ 5] = codes[(sets >= 2 ? ios[ 5] << ibits : 0) + c5.GetLong()];
  out[ 6] = codes[(sets >= 2 ? ios[ 6] << ibits : 0) + c6.GetLong()];
  out[ 7] = codes[(sets >= 2 ? ios[ 7] << ibits : 0) + c7.GetLong()];

  c0 = ExtrBits<ibits, ibits *  8>(iblk);
  c1 = ExtrBits<ibits, ibits *  9>(iblk);
  c2 = ExtrBits<ibits, ibits * 10>(iblk);
  c3 = ExtrBits<ibits, ibits * 11>(iblk);
  c4 = ExtrBits<ibits, ibits * 12>(iblk);
  c5 = ExtrBits<ibits, ibits * 13>(iblk);
  c6 = ExtrBits<ibits, ibits * 14>(iblk);
  c7 = ExtrBits<ibits, ibits * 15>(iblk);

  out[ 8] = codes[(sets >= 2 ? ios[ 8] << ibits : 0) + c0.GetLong()];
  out[ 9] = codes[(sets >= 2 ? ios[ 9] << ibits : 0) + c1.GetLong()];
  out[10] = codes[(sets >= 2 ? ios[10] << ibits : 0) + c2.GetLong()];
  out[11] = codes[(sets >= 2 ? ios[11] << ibits : 0) + c3.GetLong()];
  out[12] = codes[(sets >= 2 ? ios[12] << ibits : 0) + c4.GetLong()];
  out[13] = codes[(sets >= 2 ? ios[13] << ibits : 0) + c5.GetLong()];
  out[14] = codes[(sets >= 2 ? ios[14] << ibits : 0) + c6.GetLong()];
  out[15] = codes[(sets >= 2 ? ios[15] << ibits : 0) + c7.GetLong()];
}

template<const int sets, const int ibits, const int abits, const int begin>
void ReadPaletteBlock(int partition, unsigned int *icodes, unsigned int *acodes, Col4 &blkl, Col4 &blkh, int *out)
{
#define	fbits	(ibits < abits ? ibits : abits)	// index-block of smaller number of index-bits leads (first)
#define	sbits	(ibits > abits ? ibits : abits)	// index-block of larger number of index-bits follows (second)
  Col4 fblk, sblk;
  Col4 iblk, ablk;

  /* none of the cases straddles the lo to hi border, all go into hi
   *
   * WritePaletteBlock<1, 2,3, 50>(0, remapped, blkl, blkh);
   * WritePaletteBlock<1, 3,2, 50>(0, remapped, blkl, blkh);
   * WritePaletteBlock<1, 2,2, 66>(0, remapped, blkl, blkh);
   *
   */

  // one bit is omited per set, so it's one instruction per set + 1
  switch (begin) {
    case 50: {
      // always index 0
      fblk = ShiftLeftHalf<fbits>(
        ExtrBits<fbits * 15 - 0, /*begin*/50 -  0 + fbits *  1 - 1>(blkl)) |
        ExtrBits<fbits *  1 - 1, /*begin*/50 -  0 +              0>(blkl)
      ;
      sblk = ShiftLeftHalf<sbits>(
        ExtrBits<sbits * 15 - 0, /*begin*/50 - 64 + fbits * 16 - 1 + sbits *  1 - 1>(blkh)) |
        ExtrBits<sbits *  1 - 1, /*begin*/50 - 64 + fbits * 16 - 1 +              0>(blkh)
      ;

      // i-index block straddles 64bit border
      fblk = CopyBits<fbits * 15 - 0, (64 - /*begin*/50) + 1>(fblk, blkh);
    } break;
    case 66: {
      // always index 0
      fblk = ShiftLeftHalf<fbits>(
        ExtrBits<fbits * 15 - 0, /*begin*/66 - 64 + fbits *  1 - 1>(blkh)) |
        ExtrBits<fbits *  1 - 1, /*begin*/66 - 64 +              0>(blkh)
      ;
      sblk = ShiftLeftHalf<sbits>(
        ExtrBits<sbits * 15 - 0, /*begin*/66 - 64 + fbits * 16 - 1 + sbits *  1 - 1>(blkh)) |
        ExtrBits<sbits *  1 - 1, /*begin*/66 - 64 + fbits * 16 - 1 +              0>(blkh)
      ;
    } break;
  }

  // index-block of smaller number of index-bits leads (first)
  // index-block of larger number of index-bits follows (second)
  iblk = (ibits == fbits ? fblk : sblk);
  ablk = (abits == sbits ? sblk : fblk);

  // max 16*4 -> 64bits
  // no scattered writes on SSE
  // movd: 3clk 1/1, extrq: 2clk 1/2
  Col4 c0, c1, c2, c3,
       a0, a1, a2, a3;

  // m5:
  //  c00 = 0 -> [-][0] ->            , c08 = 1 -> [-][2] ->
  //  a00 = 0 -> [-][0] -> e7 d6 9c ff, a08 = 1 -> [-][1] -> ee d8 90 ff
  //  c01 = 1 -> [-][1] ->            , c09 = 3 -> [-][2] ->
  //  a01 = 1 -> [-][0] -> e7 d7 96 ff, a09 = 3 -> [-][1] -> ee d8 90 ff
  //  c02 = 3 -> [-][3] ->            , c10 = 3 -> [-][4] ->
  //  a02 = 3 -> [-][0] -> e7 d9 8b ff, a10 = 3 -> [-][1] -> ee db 84 ff
  //  c03 = 1 -> [-][4] ->            , c11 = 2 -> [-][4] ->
  //  a03 = 1 -> [-][1] -> ee db 84 ff, a11 = 2 -> [-][2] -> f4 db 84 ff
  //  c04 = 0 -> [-][2] ->            , c12 = 2 -> [-][3] ->
  //  a04 = 0 -> [-][0] -> e7 d8 90 ff, a12 = 2 -> [-][1] -> ee d9 8b ff
  //  c05 = 1 -> [-][3] ->            , c13 = 2 -> [-][3] ->
  //  a05 = 1 -> [-][0] -> e7 d9 8b ff, a13 = 2 -> [-][1] -> ee d9 8b ff
  //  c06 = 2 -> [-][3] ->            , c14 = 0 -> [-][4] ->
  //  a06 = 2 -> [-][0] -> e7 d9 8b ff, a14 = 0 -> [-][2] -> f4 db 84 ff
  //  c07 = 2 -> [-][2] ->            , c15 = 1 -> [-][4] ->
  //  a07 = 2 -> [-][2] -> f4 d8 90 ff, a15 = 1 -> [-][3] -> fb db 84 ff

  // m6:
  //  c00 = 0 -> [-][0] ->            , c08 = 1 -> [-][0] ->
  //  a00 = 0 -> [-][0] -> 00 00 00 00, a08 = 1 -> [-][0] -> 00 00 00 00
  //  c01 = 1 -> [-][0] ->            , c09 = 3 -> [-][0] ->
  //  a01 = 1 -> [-][0] -> 00 00 00 00, a09 = 3 -> [-][0] -> 00 00 00 00
  //  c02 = 3 -> [-][0] ->            , c10 = 3 -> [-][0] ->
  //  a02 = 3 -> [-][0] -> 00 00 00 00, a10 = 3 -> [-][0] -> 00 00 00 00
  //  c03 = 1 -> [-][0] ->            , c11 = 2 -> [-][0] ->
  //  a03 = 1 -> [-][0] -> 00 00 00 00, a11 = 2 -> [-][0] -> 00 00 00 00
  //  c04 = 0 -> [-][0] ->            , c12 = 2 -> [-][3] ->
  //  a04 = 0 -> [-][0] -> 00 00 00 00, a12 = 2 -> [-][3] -> 0b 00 00 0a
  //  c05 = 1 -> [-][0] ->            , c13 = 2 -> [-][1] ->
  //  a05 = 1 -> [-][0] -> 00 00 00 00, a13 = 2 -> [-][1] -> 04 00 00 03
  //  c06 = 2 -> [-][0] ->            , c14 = 0 -> [-][0] ->
  //  a06 = 2 -> [-][0] -> 00 00 00 00, a14 = 0 -> [-][0] -> 00 00 00 00
  //  c07 = 2 -> [-][0] ->            , c15 = 1 -> [-][0] ->
  //  a07 = 2 -> [-][0] -> 00 00 00 00, a15 = 1 -> [-][0] -> 00 00 00 00

  c0 = ExtrBits<ibits, ibits *  0>(iblk);
  a0 = ExtrBits<abits, abits *  0>(ablk);
  c1 = ExtrBits<ibits, ibits *  1>(iblk);
  a1 = ExtrBits<abits, abits *  1>(ablk);
  c2 = ExtrBits<ibits, ibits *  2>(iblk);
  a2 = ExtrBits<abits, abits *  2>(ablk);
  c3 = ExtrBits<ibits, ibits *  3>(iblk);
  a3 = ExtrBits<abits, abits *  3>(ablk);

  out[ 0] = icodes[0 + c0.GetLong()] +
            acodes[0 + a0.GetLong()];
  out[ 1] = icodes[0 + c1.GetLong()] +
            acodes[0 + a1.GetLong()];
  out[ 2] = icodes[0 + c2.GetLong()] +
            acodes[0 + a2.GetLong()];
  out[ 3] = icodes[0 + c3.GetLong()] +
            acodes[0 + a3.GetLong()];

  c0 = ExtrBits<ibits, ibits *  4>(iblk);
  a0 = ExtrBits<abits, abits *  4>(ablk);
  c1 = ExtrBits<ibits, ibits *  5>(iblk);
  a1 = ExtrBits<abits, abits *  5>(ablk);
  c2 = ExtrBits<ibits, ibits *  6>(iblk);
  a2 = ExtrBits<abits, abits *  6>(ablk);
  c3 = ExtrBits<ibits, ibits *  7>(iblk);
  a3 = ExtrBits<abits, abits *  7>(ablk);

  out[ 4] = icodes[0 + c0.GetLong()] +
            acodes[0 + a0.GetLong()];
  out[ 5] = icodes[0 + c1.GetLong()] +
            acodes[0 + a1.GetLong()];
  out[ 6] = icodes[0 + c2.GetLong()] +
            acodes[0 + a2.GetLong()];
  out[ 7] = icodes[0 + c3.GetLong()] +
            acodes[0 + a3.GetLong()];

  c0 = ExtrBits<ibits, ibits *  8>(iblk);
  a0 = ExtrBits<abits, abits *  8>(ablk);
  c1 = ExtrBits<ibits, ibits *  9>(iblk);
  a1 = ExtrBits<abits, abits *  9>(ablk);
  c2 = ExtrBits<ibits, ibits * 10>(iblk);
  a2 = ExtrBits<abits, abits * 10>(ablk);
  c3 = ExtrBits<ibits, ibits * 11>(iblk);
  a3 = ExtrBits<abits, abits * 11>(ablk);

  out[ 8] = icodes[0 + c0.GetLong()] +
            acodes[0 + a0.GetLong()];
  out[ 9] = icodes[0 + c1.GetLong()] +
            acodes[0 + a1.GetLong()];
  out[10] = icodes[0 + c2.GetLong()] +
            acodes[0 + a2.GetLong()];
  out[11] = icodes[0 + c3.GetLong()] +
            acodes[0 + a3.GetLong()];

  c0 = ExtrBits<ibits, ibits * 12>(iblk);
  a0 = ExtrBits<abits, abits * 12>(ablk);
  c1 = ExtrBits<ibits, ibits * 13>(iblk);
  a1 = ExtrBits<abits, abits * 13>(ablk);
  c2 = ExtrBits<ibits, ibits * 14>(iblk);
  a2 = ExtrBits<abits, abits * 14>(ablk);
  c3 = ExtrBits<ibits, ibits * 15>(iblk);
  a3 = ExtrBits<abits, abits * 15>(ablk);

  out[12] = icodes[0 + c0.GetLong()] +
            acodes[0 + a0.GetLong()];
  out[13] = icodes[0 + c1.GetLong()] +
            acodes[0 + a1.GetLong()];
  out[14] = icodes[0 + c2.GetLong()] +
            acodes[0 + a2.GetLong()];
  out[15] = icodes[0 + c3.GetLong()] +
            acodes[0 + a3.GetLong()];
#undef	fbits
#undef	sbits
}

/* -----------------------------------------------------------------------------
 */
void ReadPaletteBlock3_m1(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[3][FIELDN];
  Col4 b[3][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[3][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 1 mode bit, 4 partition bits
  partition = ExtrBits< 4,  1 -  0>(blkl).GetLong();

  // b2 -> set 3 stop  {0x0000000a, 0x0000000c, 0x0000000d, 0x00000000}
  // b1 -> set 2 stop  {0x00000008, 0x0000000c, 0x0000000e, 0x00000000}
  // b0 -> set 1 stop  {0x00000008, 0x0000000c, 0x0000000c, 0x00000000}
  // a2 -> set 3 start {0x00000008, 0x0000000c, 0x0000000c, 0x00000000}
  // a1 -> set 2 start {0x00000008, 0x0000000c, 0x0000000c, 0x00000000}
  // a0 -> set 1 start {0x00000008, 0x0000000d, 0x0000000e, 0x00000000}
  //
  // u  -> set 3 stop  {0x00000001, 0x00000001, 0x00000001, 0x00000000}
  // u  -> set 2 stop  {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 1 stop  {0x00000001, 0x00000001, 0x00000001, 0x00000000}
  // u  -> set 3 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 2 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 1 start {0x00000001, 0x00000001, 0x00000001, 0x00000000}
  //
  // b3 -> set 4 stop  {0x00000015, 0x00000019, 0x0000001b, 0x00000000}
  // b1 -> set 2 stop  {0x00000010, 0x00000018, 0x0000001c, 0x00000000}
  // b0 -> set 1 stop  {0x00000011, 0x00000019, 0x00000019, 0x00000000}
  // a2 -> set 3 start {0x00000010, 0x00000018, 0x00000018, 0x00000000}
  // a1 -> set 2 start {0x00000010, 0x00000018, 0x00000018, 0x00000000}
  // a0 -> set 1 start {0x00000011, 0x0000001b, 0x0000001d, 0x00000000}
  //
  // b2 -> set 3 stop  {0x000000ad, 0x000000ce, 0x000000de, 0x000000ff}
  // b1 -> set 2 stop  {0x00000084, 0x000000c6, 0x000000e7, 0x000000ff}
  // b0 -> set 1 stop  {0x0000008c, 0x000000ce, 0x000000ce, 0x000000ff}
  // a2 -> set 3 start {0x00000084, 0x000000c6, 0x000000c6, 0x000000ff}
  // a1 -> set 2 start {0x00000084, 0x000000c6, 0x000000c6, 0x000000ff}
  // a0 -> set 1 start {0x0000008c, 0x000000de, 0x000000ef, 0x000000ff}

  ExtrBits< 4, 73 - 64>(blkh, b[2][C]);	// 4 bits set 3 blue stop
  ExtrBits< 4, 65 - 64>(blkh, b[1][C]);	// 4 bits set 2 blue stop
  ExtrBits< 4, 57 -  0>(blkl, b[0][C]);	// 4 bits set 1 blue stop

  ExtrBits< 4, 69 - 64>(blkh, a[2][C]);	// 4 bits set 3 blue start
  ExtrBits< 3, 61 -  0>(blkl, a[1][C]);	// 4 bits set 2 blue start
  ExtrBits< 4, 53 -  0>(blkl, a[0][C]);	// 4 bits set 1 blue start
  a[1][C] =
  CopyBits< 1, 64 - 61>(a[1][C], blkh);	// 4 bits set 2 blue start

  ConcBits< 4, 49 -  0>(blkl, b[2][C]);	// 4 bits set 3 green stop
  ConcBits< 4, 41 -  0>(blkl, b[1][C]);	// 4 bits set 2 green stop
  ConcBits< 4, 33 -  0>(blkl, b[0][C]);	// 4 bits set 1 green stop

  ConcBits< 4, 45 -  0>(blkl, a[2][C]);	// 4 bits set 3 green start
  ConcBits< 4, 37 -  0>(blkl, a[1][C]);	// 4 bits set 2 green start
  ConcBits< 4, 29 -  0>(blkl, a[0][C]);	// 4 bits set 1 green start

  ConcBits< 4, 25 -  0>(blkl, b[2][C]);	// 4 bits set 3 red stop
  ConcBits< 4, 17 -  0>(blkl, b[1][C]);	// 4 bits set 2 red stop
  ConcBits< 4,  9 -  0>(blkl, b[0][C]);	// 4 bits set 1 red stop

  ConcBits< 4, 21 -  0>(blkl, a[2][C]);	// 4 bits set 3 red start
  ConcBits< 4, 13 -  0>(blkl, a[1][C]);	// 4 bits set 2 red start
  ConcBits< 4,  5 -  0>(blkl, a[0][C]);	// 4 bits set 1 red start

  ReplBits<-1, 82 - 64>(blkh, b[2][U]);	// 1 bits set 3 unique stop
  ReplBits<-1, 80 - 64>(blkh, b[1][U]);	// 1 bits set 2 unique stop
  ReplBits<-1, 78 - 64>(blkh, b[0][U]);	// 1 bits set 1 unique stop

  ReplBits<-1, 81 - 64>(blkh, a[2][U]);	// 1 bits set 3 unique start
  ReplBits<-1, 79 - 64>(blkh, a[1][U]);	// 1 bits set 2 unique start
  ReplBits<-1, 77 - 64>(blkh, a[0][U]);	// 1 bits set 1 unique start

  // insert the 1 unique bit & extend 4+1 bits to 8 bits
  UnpackFrom<4,4,4,0,1,0>(a);
  UnpackFrom<4,4,4,0,1,0>(b);

  // generate the midpoints
  CodebookP<3>(codes[2], a[2][C], b[2][C]);
  CodebookP<3>(codes[1], a[1][C], b[1][C]);
  CodebookP<3>(codes[0], a[0][C], b[0][C]);

  // 128 - 83 -> 45 index bits + 3 bit from 3 set start/end order -> 16 * 3bit
  ReadPaletteBlock<3, 3, 83>(partition, (unsigned int *)codes, blkl, blkh, (int *)rgba);
}

void ReadPaletteBlock3_m2(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[2][1 << 3];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 2 mode bit, 6 partition bits
  partition = ExtrBits< 6,  2 -  0>(blkl).GetLong();

  // b1 -> set 2 stop  {0x0000001b, 0x00000030, 0x00000032, 0x00000000}
  // b0 -> set 1 stop  {0x0000001f, 0x00000030, 0x00000033, 0x00000000}
  // a1 -> set 2 start {0x00000021, 0x00000030, 0x00000032, 0x00000000}
  // a0 -> set 1 start {0x00000025, 0x00000033, 0x00000038, 0x00000000}
  //
  // s -> set 12 stop  {0x00000001, 0x00000001, 0x00000001, 0x00000000}
  // s -> set 12 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  //
  // b1 -> set 2 stop  {0x00000037, 0x00000061, 0x00000065, 0x00000000}
  // b0 -> set 1 stop  {0x0000003e, 0x00000060, 0x00000066, 0x00000000}
  // a1 -> set 2 start {0x00000043, 0x00000061, 0x00000065, 0x00000000}
  // a0 -> set 1 start {0x0000004a, 0x00000066, 0x00000070, 0x00000000}
  //
  // b1 -> set 2 stop  {0x0000006e, 0x000000c3, 0x000000cb, 0x000000ff}
  // b0 -> set 1 stop  {0x0000007c, 0x000000c1, 0x000000cd, 0x000000ff}
  // a1 -> set 2 start {0x00000087, 0x000000c3, 0x000000cb, 0x000000ff}
  // a0 -> set 1 start {0x00000095, 0x000000cd, 0x000000e1, 0x000000ff}

  ExtrBits< 6, 74 - 64>(blkh, b[1][C]);	// 4 bits set 2 blue stop
  ExtrBits< 2, 62 -  0>(blkl, b[0][C]);	// 4 bits set 1 blue stop
  b[0][C] =
  CopyBits< 4, 64 - 62>(b[0][C], blkh);	// 4 bits set 1 blue stop

  ExtrBits< 6, 68 - 64>(blkh, a[1][C]);	// 4 bits set 2 blue start
  ExtrBits< 6, 56 -  0>(blkl, a[0][C]);	// 4 bits set 1 blue start

  ConcBits< 6, 50 -  0>(blkl, b[1][C]);	// 4 bits set 2 green stop
  ConcBits< 6, 38 -  0>(blkl, b[0][C]);	// 4 bits set 1 green stop

  ConcBits< 6, 44 -  0>(blkl, a[1][C]);	// 4 bits set 2 green start
  ConcBits< 6, 32 -  0>(blkl, a[0][C]);	// 4 bits set 1 green start

  ConcBits< 6, 26 -  0>(blkl, b[1][C]);	// 4 bits set 2 red stop
  ConcBits< 6, 14 -  0>(blkl, b[0][C]);	// 4 bits set 1 red stop

  ConcBits< 6, 20 -  0>(blkl, a[1][C]);	// 4 bits set 2 red start
  ConcBits< 6,  8 -  0>(blkl, a[0][C]);	// 4 bits set 1 red start

  ReplBits<-1, 81 - 64>(blkh, a[1][S]);	// 1 bits set 2 shared start/stop
  ReplBits<-1, 80 - 64>(blkh, a[0][S]);	// 1 bits set 1 shared start/stop

  b[1][S] = a[1][S];
  b[0][S] = a[0][S];

  // insert the 1 shared bit & extend 6+1 bits to 8 bits
  UnpackFrom<6,6,6,0,0,1>(a);
  UnpackFrom<6,6,6,0,0,1>(b);

  // generate the midpoints
  CodebookP<3>(codes[1], a[1][C], b[1][C]);
  CodebookP<3>(codes[0], a[0][C], b[0][C]);

  // 128 - 82 -> 46 index bits + 2 bit from 2 set start/end order -> 16 * 3bit
  ReadPaletteBlock<2, 3, 82>(partition, (unsigned int *)codes, blkl, blkh, (int *)rgba);
}

void ReadPaletteBlock3_m3(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[3][FIELDN];
  Col4 b[3][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[3][1 << 2];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 3 mode bit, 6 partition bits
  partition = ExtrBits< 6,  3 -  0>(blkl).GetLong();

  // b2 -> set 3 stop  {0x00000015, 0x00000016, 0x00000011, 0x00000000}
  // b1 -> set 2 stop  {0x00000016, 0x00000016, 0x00000012, 0x00000000}
  // b0 -> set 1 stop  {0x00000016, 0x00000017, 0x00000011, 0x00000000}
  // a2 -> set 3 start {0x00000016, 0x00000018, 0x00000014, 0x00000012}
  // a1 -> set 2 start {0x00000012, 0x00000015, 0x00000012, 0x00000000}
  // a0 -> set 1 start {0x00000013, 0x00000016, 0x00000012, 0x00000000}
  //
  // b2 -> set 3 stop  {0x000000ad, 0x000000b5, 0x0000008c, 0x000000ff}
  // b1 -> set 2 stop  {0x000000b5, 0x000000b5, 0x00000094, 0x000000ff}
  // b0 -> set 1 stop  {0x000000b5, 0x000000bd, 0x0000008c, 0x000000ff}
  // a2 -> set 3 start {0x000000b5, 0x000000c6, 0x000000a5, 0x000000ff}
  // a1 -> set 2 start {0x00000094, 0x000000ad, 0x00000094, 0x000000ff}
  // a0 -> set 1 start {0x0000009c, 0x000000b5, 0x00000094, 0x000000ff}

  ExtrBits< 5, 94 - 64>(blkh, b[2][C]);	// 4 bits set 3 blue stop
  ExtrBits< 5, 84 - 64>(blkh, b[1][C]);	// 4 bits set 2 blue stop
  ExtrBits< 5, 74 - 64>(blkh, b[0][C]);	// 4 bits set 1 blue stop

  ExtrBits< 5, 89 - 64>(blkh, a[2][C]);	// 4 bits set 3 blue start
  ExtrBits< 5, 79 - 64>(blkh, a[1][C]);	// 4 bits set 2 blue start
  ExtrBits< 5, 69 - 64>(blkh, a[0][C]);	// 4 bits set 1 blue start

  ConcBits< 5, 64 - 64>(blkh, b[2][C]);	// 4 bits set 3 green stop
  ConcBits< 5, 54 -  0>(blkl, b[1][C]);	// 4 bits set 2 green stop
  ConcBits< 5, 44 -  0>(blkl, b[0][C]);	// 4 bits set 1 green stop

  ConcBits< 5, 59 -  0>(blkl, a[2][C]);	// 4 bits set 3 green start
  ConcBits< 5, 49 -  0>(blkl, a[1][C]);	// 4 bits set 2 green start
  ConcBits< 5, 39 -  0>(blkl, a[0][C]);	// 4 bits set 1 green start

  ConcBits< 5, 34 -  0>(blkl, b[2][C]);	// 4 bits set 3 red stop
  ConcBits< 5, 24 -  0>(blkl, b[1][C]);	// 4 bits set 2 red stop
  ConcBits< 5, 14 -  0>(blkl, b[0][C]);	// 4 bits set 1 red stop

  ConcBits< 5, 29 -  0>(blkl, a[2][C]);	// 4 bits set 3 red start
  ConcBits< 5, 19 -  0>(blkl, a[1][C]);	// 4 bits set 2 red start
  ConcBits< 5,  9 -  0>(blkl, a[0][C]);	// 4 bits set 1 red start

  // extend 5 bits to 8 bits
  UnpackFrom<5,5,5,0,0,0>(a);
  UnpackFrom<5,5,5,0,0,0>(b);

  // generate the midpoints
  CodebookP<2>(codes[2], a[2][C], b[2][C]);
  CodebookP<2>(codes[1], a[1][C], b[1][C]);
  CodebookP<2>(codes[0], a[0][C], b[0][C]);

  // 128 - 99 -> 29 index bits + 3 bit from 3 set start/end order -> 16 * 2bit
  ReadPaletteBlock<3, 2, 99>(partition, (unsigned int *)codes, blkl, blkh, (int *)rgba);
}

void ReadPaletteBlock3_m4(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[2][1 << 2];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 4 mode bit, 6 partition bits
  partition = ExtrBits< 6,  4 -  0>(blkl).GetLong();

  // b1 -> set 2 stop  {0x00000046, 0x00000064, 0x0000006e, 0x00000000}
  // b0 -> set 1 stop  {0x00000041, 0x00000061, 0x00000068, 0x00000000}
  // a1 -> set 2 start {0x00000045, 0x00000061, 0x00000069, 0x00000000}
  // a0 -> set 1 start {0x00000049, 0x00000064, 0x0000006e, 0x00000000}
  //
  // u  -> set 2 stop  {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 1 stop  {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 2 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 1 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  //
  // b1 -> set 2 stop  {0x00000046, 0x00000064, 0x0000006e, 0x00000000}
  // b0 -> set 1 stop  {0x00000041, 0x00000061, 0x00000068, 0x00000000}
  // a1 -> set 2 start {0x00000045, 0x00000061, 0x00000069, 0x00000000}
  // a0 -> set 1 start {0x00000049, 0x00000064, 0x0000006e, 0x00000000}
  //
  // b1 -> set 2 stop  {0x0000008c, 0x000000c8, 0x000000dc, 0x000000ff}
  // b0 -> set 1 stop  {0x00000082, 0x000000c2, 0x000000d0, 0x000000ff}
  // a1 -> set 2 start {0x0000008a, 0x000000c2, 0x000000d2, 0x000000ff}
  // a0 -> set 1 start {0x00000092, 0x000000c8, 0x000000dc, 0x000000ff}

  ExtrBits< 7, 87 - 64>(blkh, b[1][C]);	// 4 bits set 2 blue stop
  ExtrBits< 7, 73 - 64>(blkh, b[0][C]);	// 4 bits set 1 blue stop

  ExtrBits< 7, 80 - 64>(blkh, a[1][C]);	// 4 bits set 2 blue start
  ExtrBits< 7, 66 - 64>(blkh, a[0][C]);	// 4 bits set 1 blue start

  ConcBits< 5, 59 -  0>(blkl, b[1][C]);	// 4 bits set 2 green stop
  ConcBits< 7, 45 -  0>(blkl, b[0][C]);	// 4 bits set 1 green stop
  b[1][C] =
  InjtBits< 2, 64 - 59>(b[1][C], blkh);	// 4 bits set 2 green stop

  ConcBits< 7, 52 -  0>(blkl, a[1][C]);	// 4 bits set 2 green start
  ConcBits< 7, 38 -  0>(blkl, a[0][C]);	// 4 bits set 1 green start

  ConcBits< 7, 31 -  0>(blkl, b[1][C]);	// 4 bits set 2 red stop
  ConcBits< 7, 17 -  0>(blkl, b[0][C]);	// 4 bits set 1 red stop

  ConcBits< 7, 24 -  0>(blkl, a[1][C]);	// 4 bits set 2 red start
  ConcBits< 7, 10 -  0>(blkl, a[0][C]);	// 4 bits set 1 red start

  ReplBits<-1, 97 - 64>(blkh, b[1][U]);	// 1 bits set 2 unique stop
  ReplBits<-1, 95 - 64>(blkh, b[0][U]);	// 1 bits set 1 unique stop

  ReplBits<-1, 96 - 64>(blkh, a[1][U]);	// 1 bits set 2 unique start
  ReplBits<-1, 94 - 64>(blkh, a[0][U]);	// 1 bits set 1 unique start

  // insert the 1 shared bit & extend 7+1 bits to 8 bits
  UnpackFrom<7,7,7,0,1,0>(a);
  UnpackFrom<7,7,7,0,1,0>(b);

  // generate the midpoints
  CodebookP<2>(codes[1], a[1][C], b[1][C]);
  CodebookP<2>(codes[0], a[0][C], b[0][C]);

  // 128 - 98 -> 30 index bits + 2 bit from 2 set start/end order -> 16 * 2bit
  ReadPaletteBlock<2, 2, 98>(partition, (unsigned int *)codes, blkl, blkh, (int *)rgba);
}

void ReadPaletteBlock4_m5(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[2][1 << 3];
  int rix;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 5 mode bit, 2 rotation bits, 1 index bit
  rix = ExtrBits< 3,  5 -  0>(blkl).GetLong();

  // b1 -> set 2 stop  {0x0000003e, 0x00000000, 0x00000000, 0x00000000}
  // b0 -> set 1 stop  {0x00000000, 0x0000001b, 0x0000000e, 0x0000001f}
  // a1 -> set 2 start {0x00000039, 0x00000000, 0x00000000, 0x00000000}
  // a0 -> set 1 start {0x00000000, 0x0000001a, 0x00000013, 0x0000001f}
  //
  // b1 -> set 2 stop  {0x000000fb, 0x00000000, 0x00000000, 0x00000000}
  // b0 -> set 1 stop  {0x00000000, 0x000000de, 0x00000073, 0x000000ff}
  // a1 -> set 2 start {0x000000e7, 0x00000000, 0x00000000, 0x00000000}
  // a0 -> set 1 start {0x00000000, 0x000000d6, 0x0000009c, 0x000000ff}

  // rotate colors
  assume(rix >= 0 && rix <= 7);
  switch (rix) {
    case 0: case 4:
  	ExtrBits< 6, 44 -  0>(blkl, b[1][C]);	// 6 bits set 1 alpha stop
  	ExtrBits< 6, 38 -  0>(blkl, a[1][C]);	// 6 bits set 1 alpha start

	b[1][C] = ShiftLeft<96>(b[1][C]);	// 6 bits set 1 alpha stop
	a[1][C] = ShiftLeft<96>(a[1][C]);	// 6 bits set 1 alpha start

	b[0][C] = ShiftLeft<0>(b[0][C]);	// 6 bits set 1 alpha stop
	a[0][C] = ShiftLeft<0>(a[0][C]);	// 6 bits set 1 alpha start

  	ExtrBits< 5, 33 -  0>(blkl, b[0][C]);	// 5 bits set 1 blue stop
  	ExtrBits< 5, 28 -  0>(blkl, a[0][C]);	// 5 bits set 1 blue start

  	ConcBits< 5, 23 -  0>(blkl, b[0][C]);	// 5 bits set 1 green stop
  	ConcBits< 5, 18 -  0>(blkl, a[0][C]);	// 5 bits set 1 green start

  	ConcBits< 5, 13 -  0>(blkl, b[0][C]);	// 5 bits set 1 red stop
  	ConcBits< 5,  8 -  0>(blkl, a[0][C]);	// 5 bits set 1 red start

	// extend 5/6 bits to 8 bits
	UnpackFrom<5,5,5,6,0,0>(a);
	UnpackFrom<5,5,5,6,0,0>(b);
    	break;
    case 1: case 5:
  	ExtrBits< 6, 44 -  0>(blkl, b[1][C]);	// 6 bits set 1 red stop
  	ExtrBits< 6, 38 -  0>(blkl, a[1][C]);	// 6 bits set 1 red start

	b[1][C] = ShiftLeft<0>(b[1][C]);	// 6 bits set 1 red stop
	a[1][C] = ShiftLeft<0>(a[1][C]);	// 6 bits set 1 red start

  	ExtrBits< 5, 13 -  0>(blkl, b[0][C]);	// 5 bits set 1 alpha stop
  	ExtrBits< 5,  8 -  0>(blkl, a[0][C]);	// 5 bits set 1 alpha start

  	ConcBits< 5, 33 -  0>(blkl, b[0][C]);	// 5 bits set 1 blue stop
  	ConcBits< 5, 28 -  0>(blkl, a[0][C]);	// 5 bits set 1 blue start

  	ConcBits< 5, 23 -  0>(blkl, b[0][C]);	// 5 bits set 1 green stop
  	ConcBits< 5, 18 -  0>(blkl, a[0][C]);	// 5 bits set 1 green start

	b[0][C] = ShiftLeft<32>(b[0][C]);	// 6 bits set 1 red stop
	a[0][C] = ShiftLeft<32>(a[0][C]);	// 6 bits set 1 red start

	// extend 5/6 bits to 8 bits
	UnpackFrom<6,5,5,5,0,0>(a);
	UnpackFrom<6,5,5,5,0,0>(b);
    	break;
    case 2: case 6:
  	ExtrBits< 6, 44 -  0>(blkl, b[1][C]);	// 6 bits set 1 green stop
  	ExtrBits< 6, 38 -  0>(blkl, a[1][C]);	// 6 bits set 1 green start

	b[1][C] = ShiftLeft<32>(b[1][C]);	// 6 bits set 1 green stop
	a[1][C] = ShiftLeft<32>(a[1][C]);	// 6 bits set 1 green start

  	ExtrBits< 5, 23 -  0>(blkl, b[0][C]);	// 5 bits set 1 alpha stop
  	ExtrBits< 5, 18 -  0>(blkl, a[0][C]);	// 5 bits set 1 alpha start

  	ConcBits< 5, 33 -  0>(blkl, b[0][C]);	// 5 bits set 1 blue stop
  	ConcBits< 5, 28 -  0>(blkl, a[0][C]);	// 5 bits set 1 blue start

	b[0][C] = ShiftLeft<32>(b[0][C]);	// 6 bits set 1 green stop
	a[0][C] = ShiftLeft<32>(a[0][C]);	// 6 bits set 1 green start

  	ConcBits< 5, 13 -  0>(blkl, b[0][C]);	// 5 bits set 1 red stop
  	ConcBits< 5,  8 -  0>(blkl, a[0][C]);	// 5 bits set 1 red start

	// extend 5/6 bits to 8 bits
	UnpackFrom<5,6,5,5,0,0>(a);
	UnpackFrom<5,6,5,5,0,0>(b);
    	break;
    case 3: case 7:
  	ExtrBits< 6, 44 -  0>(blkl, b[1][C]);	// 6 bits set 1 blue stop
  	ExtrBits< 6, 38 -  0>(blkl, a[1][C]);	// 6 bits set 1 blue start

	b[1][C] = ShiftLeft<64>(b[1][C]);	// 6 bits set 1 blue stop
	a[1][C] = ShiftLeft<64>(a[1][C]);	// 6 bits set 1 blue start

  	ExtrBits< 5, 33 -  0>(blkl, b[0][C]);	// 5 bits set 1 alpha stop
  	ExtrBits< 5, 28 -  0>(blkl, a[0][C]);	// 5 bits set 1 alpha start

	b[0][C] = ShiftLeft<32>(b[0][C]);	// 6 bits set 1 blue stop
	a[0][C] = ShiftLeft<32>(a[0][C]);	// 6 bits set 1 blue start

  	ConcBits< 5, 23 -  0>(blkl, b[0][C]);	// 5 bits set 1 green stop
  	ConcBits< 5, 18 -  0>(blkl, a[0][C]);	// 5 bits set 1 green start

  	ConcBits< 5, 13 -  0>(blkl, b[0][C]);	// 5 bits set 1 red stop
  	ConcBits< 5,  8 -  0>(blkl, a[0][C]);	// 5 bits set 1 red start

	// extend 5/6 bits to 8 bits
	UnpackFrom<5,5,6,5,0,0>(a);
	UnpackFrom<5,5,6,5,0,0>(b);
    	break;
  }

  rix >>= 2;

  // generate the midpoints
  assume(rix >= 0 && rix <= 1);
  if (!rix) {
    CodebookP<2>(codes[0], a[0][C], b[0][C]);
    CodebookP<3>(codes[1], a[1][C], b[1][C]);
  }
  else {
    CodebookP<3>(codes[0], a[0][C], b[0][C]);
    CodebookP<2>(codes[1], a[1][C], b[1][C]);
  }

  // 128 - 50 -> 78 index bits -> 31 + 47 index bits + 2 bit from 2 set start/end order -> 16 * 2bit + 16 * 3bit
  if (!rix)
    ReadPaletteBlock<1, 2,3, 50>(0, (unsigned int *)codes[0], (unsigned int *)codes[1], blkl, blkh, (int *)rgba);
  else
    ReadPaletteBlock<1, 3,2, 50>(0, (unsigned int *)codes[0], (unsigned int *)codes[1], blkl, blkh, (int *)rgba);
}

void ReadPaletteBlock4_m6(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[2][1 << 2];
  int rix;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 6 mode bit, 2 rotation bits
  rix = ExtrBits< 2,  6 -  0>(blkl).GetLong();

  // b1 -> set 2 stop  {0x0000000b, 0x00000000, 0x00000000, 0x00000000}
  // b0 -> set 1 stop  {0x00000000, 0x00000000, 0x00000000, 0x00000005}
  // a1 -> set 2 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // a0 -> set 1 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  //
  // b1 -> set 2 stop  {0x0000000b, 0x00000000, 0x00000000, 0x00000000}
  // b0 -> set 1 stop  {0x00000000, 0x00000000, 0x00000000, 0x0000000a}
  // a1 -> set 2 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // a0 -> set 1 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}

  // rotate colors
  assume(rix >= 0 && rix <= 3);
  switch (rix) {
    case 0:
  	ExtrBits< 6, 58 -  0>(blkl, b[1][C]);	// 8 bits set 1 alpha stop
  	ExtrBits< 8, 50 -  0>(blkl, a[1][C]);	// 8 bits set 1 alpha start
	b[1][C] =
	CopyBits< 2, 64 - 58>(b[1][C], blkh);	// 8 bits set 1 alpha stop

	b[1][C] = ShiftLeft<96>(b[1][C]);	// 8 bits set 1 alpha stop
	a[1][C] = ShiftLeft<96>(a[1][C]);	// 8 bits set 1 alpha start

	b[0][C] = ShiftLeft<0>(b[0][C]);	// 8 bits set 1 alpha stop
	a[0][C] = ShiftLeft<0>(a[0][C]);	// 8 bits set 1 alpha start

  	ExtrBits< 7, 43 -  0>(blkl, b[0][C]);	// 7 bits set 1 blue stop
  	ExtrBits< 7, 36 -  0>(blkl, a[0][C]);	// 7 bits set 1 blue start

  	ConcBits< 7, 29 -  0>(blkl, b[0][C]);	// 7 bits set 1 green stop
  	ConcBits< 7, 22 -  0>(blkl, a[0][C]);	// 7 bits set 1 green start

  	ConcBits< 7, 15 -  0>(blkl, b[0][C]);	// 7 bits set 1 red stop
  	ConcBits< 7,  8 -  0>(blkl, a[0][C]);	// 7 bits set 1 red start

	// extend 7/8 bits to 8 bits
	UnpackFrom<7,7,7,8,0,0>(a);
	UnpackFrom<7,7,7,8,0,0>(b);
    	break;
    case 1:
  	ExtrBits< 6, 58 -  0>(blkl, b[1][C]);	// 8 bits set 1 red stop
  	ExtrBits< 8, 50 -  0>(blkl, a[1][C]);	// 8 bits set 1 red start
	b[1][C] =
	CopyBits< 2, 64 - 58>(b[1][C], blkh);	// 8 bits set 1 red stop

	b[1][C] = ShiftLeft<0>(b[1][C]);	// 8 bits set 1 red stop
	a[1][C] = ShiftLeft<0>(a[1][C]);	// 8 bits set 1 red start

  	ExtrBits< 7, 15 -  0>(blkl, b[0][C]);	// 7 bits set 1 alpha stop
  	ExtrBits< 7,  8 -  0>(blkl, a[0][C]);	// 7 bits set 1 alpha start

  	ConcBits< 7, 43 -  0>(blkl, b[0][C]);	// 7 bits set 1 blue stop
  	ConcBits< 7, 36 -  0>(blkl, a[0][C]);	// 7 bits set 1 blue start

  	ConcBits< 7, 29 -  0>(blkl, b[0][C]);	// 7 bits set 1 green stop
  	ConcBits< 7, 22 -  0>(blkl, a[0][C]);	// 7 bits set 1 green start

	b[0][C] = ShiftLeft<32>(b[0][C]);	// 7 bits set 1 red stop
	a[0][C] = ShiftLeft<32>(a[0][C]);	// 7 bits set 1 red start

	// extend 7/8 bits to 8 bits
	UnpackFrom<8,7,7,7,0,0>(a);
	UnpackFrom<8,7,7,7,0,0>(b);
    	break;
    case 2:
  	ExtrBits< 6, 58 -  0>(blkl, b[1][C]);	// 8 bits set 1 green stop
  	ExtrBits< 8, 50 -  0>(blkl, a[1][C]);	// 8 bits set 1 green start
	b[1][C] =
	CopyBits< 2, 64 - 58>(b[1][C], blkh);	// 8 bits set 1 green stop

	b[1][C] = ShiftLeft<32>(b[1][C]);	// 8 bits set 1 red stop
	a[1][C] = ShiftLeft<32>(a[1][C]);	// 8 bits set 1 red start

  	ExtrBits< 7, 29 -  0>(blkl, b[0][C]);	// 7 bits set 1 alpha stop
  	ExtrBits< 7, 22 -  0>(blkl, a[0][C]);	// 7 bits set 1 alpha start

  	ConcBits< 7, 43 -  0>(blkl, b[0][C]);	// 7 bits set 1 blue stop
  	ConcBits< 7, 36 -  0>(blkl, a[0][C]);	// 7 bits set 1 blue start

	b[0][C] = ShiftLeft<32>(b[0][C]);	// 7 bits set 1 green stop
	a[0][C] = ShiftLeft<32>(a[0][C]);	// 7 bits set 1 green start

  	ConcBits< 7, 15 -  0>(blkl, b[0][C]);	// 7 bits set 1 red stop
  	ConcBits< 7,  8 -  0>(blkl, a[0][C]);	// 7 bits set 1 red start

	// extend 7/8 bits to 8 bits
	UnpackFrom<7,8,7,7,0,0>(a);
	UnpackFrom<7,8,7,7,0,0>(b);
    	break;
    case 3:
  	ExtrBits< 6, 58 -  0>(blkl, b[1][C]);	// 8 bits set 1 blue stop
  	ExtrBits< 8, 50 -  0>(blkl, a[1][C]);	// 8 bits set 1 blue start
	b[1][C] =
	CopyBits< 2, 64 - 58>(b[1][C], blkh);	// 8 bits set 1 blue stop

	b[1][C] = ShiftLeft<64>(b[1][C]);	// 8 bits set 1 blue stop
	a[1][C] = ShiftLeft<64>(a[1][C]);	// 8 bits set 1 blue start

  	ExtrBits< 7, 43 -  0>(blkl, b[0][C]);	// 7 bits set 1 alpha stop
  	ExtrBits< 7, 36 -  0>(blkl, a[0][C]);	// 7 bits set 1 alpha start

	b[0][C] = ShiftLeft<32>(b[0][C]);	// 7 bits set 1 blue stop
	a[0][C] = ShiftLeft<32>(a[0][C]);	// 7 bits set 1 blue start

  	ConcBits< 7, 29 -  0>(blkl, b[0][C]);	// 7 bits set 1 green stop
  	ConcBits< 7, 22 -  0>(blkl, a[0][C]);	// 7 bits set 1 green start

  	ConcBits< 7, 15 -  0>(blkl, b[0][C]);	// 7 bits set 1 red stop
  	ConcBits< 7,  8 -  0>(blkl, a[0][C]);	// 7 bits set 1 red start

	// extend 7/8 bits to 8 bits
	UnpackFrom<7,7,8,7,0,0>(a);
	UnpackFrom<7,7,8,7,0,0>(b);
    	break;
  }

  // generate the midpoints
  CodebookP<2>(codes[0], a[0][C], b[0][C]);
  CodebookP<2>(codes[1], a[1][C], b[1][C]);

  // 128 - 66 -> 62 index bits -> 31 + 31 index bits + 2 bit from 2 set start/end order -> 16 * 2bit + 16 * 3bit
  ReadPaletteBlock<1, 2,2, 66>(0, (unsigned int *)codes[0], (unsigned int *)codes[1], blkl, blkh, (int *)rgba);
}

void ReadPaletteBlock4_m7(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[1][FIELDN];
  Col4 b[1][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[1][1 << 4];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 7 mode bit, 0 partition bits
  partition = 0;//ExtrBits< 0,  7 -  0>(blkl).GetLong();

  // b0 -> set 1 stop  {0x0000006e, 0x0000007b, 0x0000007e, 0x0000007f}
  // a0 -> set 1 start {0x0000005f, 0x00000076, 0x0000007c, 0x0000007f}
  //
  // u  -> set 1 stop  {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 1 start {0x00000001, 0x00000001, 0x00000001, 0x00000001}
  //
  // b0 -> set 1 stop  {0x000000dc, 0x000000f6, 0x000000fc, 0x000000fe}
  // a0 -> set 1 start {0x000000bf, 0x000000ed, 0x000000f9, 0x000000ff}
  //
  // b0 -> set 1 stop  {0x000000dc, 0x000000f6, 0x000000fc, 0x000000fe}
  // a0 -> set 1 start {0x000000bf, 0x000000ed, 0x000000f9, 0x000000ff}

  ExtrBits< 7, 56 -  0>(blkl, b[0][C]);	// 7 bits set 1 alpha stop
  ExtrBits< 7, 49 -  0>(blkl, a[0][C]);	// 7 bits set 1 alpha start

  ConcBits< 7, 42 -  0>(blkl, b[0][C]);	// 7 bits set 1 blue stop
  ConcBits< 7, 35 -  0>(blkl, a[0][C]);	// 7 bits set 1 blue start

  ConcBits< 7, 28 -  0>(blkl, b[0][C]);	// 7 bits set 1 green stop
  ConcBits< 7, 21 -  0>(blkl, a[0][C]);	// 7 bits set 1 green start

  ConcBits< 7, 14 -  0>(blkl, b[0][C]);	// 7 bits set 1 red stop
  ConcBits< 7,  7 -  0>(blkl, a[0][C]);	// 7 bits set 1 red start

  ReplBits< 1, 64 - 64>(blkh, b[0][U]);	// 1 bits set 1 unique stop
  ReplBits< 1, 63 -  0>(blkl, a[0][U]);	// 1 bits set 1 unique start

  // insert the 1 shared bit & extend 7+1 bits to 8 bits
  UnpackFrom<7,7,7,7,1,0>(a);
  UnpackFrom<7,7,7,7,1,0>(b);

  // generate the midpoints
  CodebookP<4>(codes[0], a[0][C], b[0][C]);

  // 128 - 65 -> 63 index bits + 1 bit from 1 set start/end order -> 16 * 4bit
  ReadPaletteBlock<1, 4, 65>(partition, (unsigned int *)codes, blkl, blkh, (int *)rgba);
}

void ReadPaletteBlock4_m8(u8* rgba, void const* block) {
  // get the packed values
  Col4 a[2][FIELDN];
  Col4 b[2][FIELDN];
  Col4 blkl, blkh;

  // remap the indices
  unsigned int codes[2][1 << 2];
  int partition;

  /* read in */
  LoadUnaligned(blkl, blkh, block);

  // 8 mode bit, 6 partition bits
  partition = ExtrBits< 6,  8 -  0>(blkl).GetLong();

  // b1 -> set 2 stop  {0x00000008, 0x00000000, 0x00000000, 0x00000010}
  // b0 -> set 1 stop  {0x0000000a, 0x00000000, 0x00000000, 0x0000001a}
  // a1 -> set 2 start {0x0000001d, 0x0000000e, 0x00000004, 0x0000001e}
  // a0 -> set 1 start {0x0000001f, 0x0000000b, 0x00000003, 0x0000001f}
  //
  // u  -> set 2 stop  {0x00000001, 0x00000001, 0x00000001, 0x00000001}
  // u  -> set 1 stop  {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 2 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  // u  -> set 1 start {0x00000000, 0x00000000, 0x00000000, 0x00000000}
  //
  // b1 -> set 2 stop  {0x00000011, 0x00000001, 0x00000001, 0x00000021}
  // b0 -> set 1 stop  {0x00000014, 0x00000000, 0x00000000, 0x00000034}
  // a1 -> set 2 start {0x0000003a, 0x0000001c, 0x00000008, 0x0000003c}
  // a0 -> set 1 start {0x0000003e, 0x00000016, 0x00000006, 0x0000003e}
  //
  // b1 -> set 2 stop  {0x00000045, 0x00000004, 0x00000004, 0x00000086}
  // b0 -> set 1 stop  {0x00000051, 0x00000000, 0x00000000, 0x000000d3}
  // a1 -> set 2 start {0x000000eb, 0x00000071, 0x00000020, 0x000000f3}
  // a0 -> set 1 start {0x000000fb, 0x00000059, 0x00000018, 0x000000fb}

  ExtrBits< 5, 89 - 64>(blkh, b[1][C]);	// 7 bits set 2 alpha stop
  ExtrBits< 5, 79 - 64>(blkh, b[0][C]);	// 7 bits set 1 alpha stop

  ExtrBits< 5, 84 - 64>(blkh, a[1][C]);	// 7 bits set 2 alpha start
  ExtrBits< 5, 74 - 64>(blkh, a[0][C]);	// 7 bits set 1 alpha start

  ConcBits< 5, 69 - 64>(blkh, b[1][C]);	// 7 bits set 2 blue stop
  ConcBits< 5, 59 -  0>(blkl, b[0][C]);	// 7 bits set 1 blue stop

  ConcBits< 5, 64 - 64>(blkh, a[1][C]);	// 7 bits set 2 blue start
  ConcBits< 5, 54 -  0>(blkl, a[0][C]);	// 7 bits set 1 blue start

  ConcBits< 5, 49 -  0>(blkl, b[1][C]);	// 7 bits set 2 green stop
  ConcBits< 5, 39 -  0>(blkl, b[0][C]);	// 7 bits set 1 green stop

  ConcBits< 5, 44 -  0>(blkl, a[1][C]);	// 7 bits set 2 green start
  ConcBits< 5, 34 -  0>(blkl, a[0][C]);	// 7 bits set 1 green start

  ConcBits< 5, 29 -  0>(blkl, b[1][C]);	// 7 bits set 2 red stop
  ConcBits< 5, 19 -  0>(blkl, b[0][C]);	// 7 bits set 1 red stop

  ConcBits< 5, 24 -  0>(blkl, a[1][C]);	// 7 bits set 2 red start
  ConcBits< 5, 14 -  0>(blkl, a[0][C]);	// 7 bits set 1 red start

  ReplBits< 1, 97 - 64>(blkh, b[1][U]);	// 1 bits set 2 unique stop
  ReplBits< 1, 95 - 64>(blkh, b[0][U]);	// 1 bits set 1 unique stop

  ReplBits< 1, 96 - 64>(blkh, a[1][U]);	// 1 bits set 2 unique start
  ReplBits< 1, 94 - 64>(blkh, a[0][U]);	// 1 bits set 1 unique start

  // insert the 1 shared bit & extend 7+1 bits to 8 bits
  UnpackFrom<5,5,5,5,1,0>(a);
  UnpackFrom<5,5,5,5,1,0>(b);

  // generate the midpoints
  CodebookP<2>(codes[1], a[1][C], b[1][C]);
  CodebookP<2>(codes[0], a[0][C], b[0][C]);

  // 128 - 98 -> 30 index bits + 2 bit from 2 set start/end order -> 16 * 2bit
  ReadPaletteBlock<2, 2, 98>(partition, (unsigned int *)codes, blkl, blkh, (int *)rgba);
}

#undef	C
#undef	U
#undef	S

void DecompressColoursBtc7u(u8* rgba, void const* block)
{
  // get the block bytes
  u8 const* bytes = reinterpret_cast< u8 const* >(block);

  // get the mode
#ifdef __GNUC__
  unsigned long mode = __builtin_ctz(*((int *)bytes));
#else
  unsigned long mode; _BitScanForward(&mode, *((int *)bytes));
#endif


  assume(mode >= 0 && mode <= 7);
  switch (mode) {
    case 0: ReadPaletteBlock3_m1(rgba, block); break;
    case 1: ReadPaletteBlock3_m2(rgba, block); break;
    case 2: ReadPaletteBlock3_m3(rgba, block); break;
    case 3: ReadPaletteBlock3_m4(rgba, block); break;
    case 4: ReadPaletteBlock4_m5(rgba, block); break;
    case 5: ReadPaletteBlock4_m6(rgba, block); break;
    case 6: ReadPaletteBlock4_m7(rgba, block); break;
    case 7: ReadPaletteBlock4_m8(rgba, block); break;
  }
}

void DecompressColoursBtc7u(u16* rgba, void const* block)
{
  u8 bytes[4 * 16]; DecompressColoursBtc7u(bytes, block);

  for (int v = 0; v < (4 * 16); v++)
    rgba[v] = bytes[v] * (65535 / 255);
}

void DecompressColoursBtc7u(f23* rgba, void const* block)
{
  u8 bytes[4 * 16]; DecompressColoursBtc7u(bytes, block);

  for (int v = 0; v < (4 * 16); v++)
    rgba[v] = bytes[v] * (1.0f / 255.0f);
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
static void WritePaletteBlock(tile_barrier barrier, const int thread, lineI2 palette, index16 indices, out code64 block) amp_restricted
{
  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  threaded_cse(0) {
    // write the endpoints
    block[0] =
      ((palette[CSTRT] /*& 0xFFFF*/) <<  0) |
      ((palette[CSTOP] /*& 0xFFFF*/) << 16);

    block[1] =
      (((indices[ 0] << 0) + (indices[ 1] << 2) + (indices[ 2] << 4) + (indices[ 3] << 6)) <<  0) +
      (((indices[ 4] << 0) + (indices[ 5] << 2) + (indices[ 6] << 4) + (indices[ 7] << 6)) <<  8) +
      (((indices[ 8] << 0) + (indices[ 9] << 2) + (indices[10] << 4) + (indices[11] << 6)) << 16) +
      (((indices[12] << 0) + (indices[13] << 2) + (indices[14] << 4) + (indices[15] << 6)) << 24);
  }
}

void WritePaletteBlock3(tile_barrier barrier, const int thread, lineC2 cline, inout index16 indices, out code64 block,
		       IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 2, 3, 4, 5, 6, 7 };	// (a >  b)

  // degenerate case "start > stop" (AMP: local register, not enough for group-shared)
  int palette[CVALS];
  int sorted[CVALS];

  // get the packed values (AMP: cline is guaranteed to be valid/written)
  FloatTo565(cline, palette);

  // palette[CSTRT] > palette[CSTOP] := sorted[CSTRT] == palette[CSTRT]
  sorted[CSTRT] = min(palette[CSTRT], palette[CSTOP]);
  sorted[CSTOP] = max(palette[CSTRT], palette[CSTOP]);

  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  // remap the indices
  wavefrnt_for(i, 16) {
#if	defined(SQUISH_USE_COMPUTE)
    indices[i] = (sorted[CSTRT] == palette[CSTOP] ? indices[i] : yArr[IBL_COLOR3][indices[i]].mapped);
#else
    indices[i] = (sorted[CSTRT] == palette[CSTOP] ? indices[i] : yArr(IBL_COLOR3, indices[i]).mapped);
#endif
  }

  // write the block
  WritePaletteBlock(barrier, thread, sorted, indices, block);
}

void WritePaletteBlock4(tile_barrier barrier, const int thread, lineC2 cline, inout index16 indices, out code64 block,
		       IndexBlockLUT yArr) amp_restricted
{
//static ccr8 slut[8] = { 1, 0, 3, 2, 5, 4, 7, 6 };	// (a <  b)
//static ccr8 slut[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };	// (a == b)

  // degenerate case "start > stop"
  int palette[CVALS];
  int sorted[CVALS];

  // get the packed values (AMP: cline is guaranteed to be valid/written)
  FloatTo565(cline, palette);

  // palette[CSTRT] < palette[CSTOP] := sorted[CSTRT] == palette[CSTOP]
  sorted[CSTRT] = max(palette[CSTRT], palette[CSTOP]);
  sorted[CSTOP] = min(palette[CSTRT], palette[CSTOP]);

  // AMP: make "indices"-writes visible
  tile_static_memory_fence(barrier);

  // remap the indices
  wavefrnt_for(i, 16) {
    indices[i] = (indices[i] ^ (sorted[CSTRT] == palette[CSTOP] ? 1 : 0)) & (sorted[CSTRT] == sorted[CSTOP] ? 0 : 3);
  }

  // write the block
  WritePaletteBlock(barrier, thread, palette, indices, block);
}
#endif

} // namespace squish
