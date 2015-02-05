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

#include "alphanormalfit.h"
#include "maths.h"
#include "simd.h"

#if	!defined(SQUISH_USE_COMPUTE)
#include <cmath>
#include <algorithm>
#endif

#pragma warning ( disable : 4239 )

#include "inlineables.cpp"

//#define PYRAMID_PROJECTION
//#define Complement  ComplementPyramidal

// possibly defined in some std-headers
#undef FSCALE
#undef FOFFSET
#undef CSCALE
#undef COFFSET

#undef ISCALE
#undef IOFFSET
#undef OSCALE
#undef OOFFSET

// Codebook precision bits, up to 5 (+8)
#define CBLB	0	// low precision
#define CBHB	5	// high precision

// convert from floats to codes, codes are only u8/s8
#define FSCALE	( 0.5f * ((max - min)))
#define FOFFSET	(-0.5f * ((max - min) * (min == 0)))
// convert from codes to floats, codes are only u8/s8
#define CSCALE	( 2.0f / ((max - min) * (1 << prc)))			//  127.5, 127
#define COFFSET	(-0.5f * ((max - min) * (1 << prc) * (min == 0)))	// -127.5,   0

// convert from input to floats
#define ISCALE	( 2.0f / ((upr - lwr)))
#define IOFFSET	(-0.5f * ((upr - lwr) * (lwr == 0)))
// convert from floats to output
#define OSCALE	( 0.5f * ((upr - lwr)))
#define OOFFSET	(-0.5f * ((upr - lwr) * (lwr == 0)))

namespace squish {

/* *****************************************************************************
 */
#if	!defined(SQUISH_USE_PRE)
/*static void FixRange(int& min, int& max, int steps)
{
  if (max - min < steps)
    max = std::min<int>(min + steps, 0xFF);
  if (max - min < steps)
    min = std::max<int>(0x00, max - steps);
}*/

template<const int min, const int max, const int prc, typename otyp>
static Scr4 FitCodes(Vec4 const* xyz, int mask,
		     Col8 const &codesx, u8* indicesx,
		     Col8 const &codesy, u8* indicesy)
{
  // fit each coord value to the codebook
  Scr4 error = Scr4(DEVIANCE_BASE);

  for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
    // check this pixel is valid
    if ((imask & 1) == 0) {
      // use the first code
      indicesx[i] = 0;
      indicesy[i] = 0;
      continue;
    }

    // fetch floating point vector
    Vec4 valx = xyz[i].SplatX();
    Vec4 valy = xyz[i].SplatY();
    Vec4 valz = xyz[i].SplatZ();

    // find the closest code
    Scr4 dist = Vec4(DEVIANCE_MAX);

    int idxx = 0;
    int idxy = 0;

    // find the least error and corresponding index
    // prefer lower indices over higher ones
    
    // Col4(codesy[k    ], codesy[k    ], codesy[k    ], codesy[k    ]);
    int k = 7; Col8 _cy = codesy; do { Col4 _yy = RepeatUpper(_cy, otyp(0)); _cy = ShiftUp<1>(_cy);
    // Col4(codesx[j - 0], codesx[j - 1], codesx[j - 2], codesx[j - 3]);
    int j = 7; Col8 _cx = codesx; do { Col4 _xx = ExpandUpper(_cx, otyp(0)); _cx = ShiftUp<4>(_cx);
    
      Vec4 cxx = (_xx + Vec4(COFFSET)) * Vec4(CSCALE);
      Vec4 cyy = (_yy + Vec4(COFFSET)) * Vec4(CSCALE);
      Vec4 czz = Complement<DISARM>(cxx, cyy);

      // measure absolute angle-deviation (cosine)
      Vec4 dst = (valx * cxx) + (valy * cyy) + (valz * czz);

      // select the smallest deviation (NaN as first arg is ignored!)
      dist = HorizontalMax(Max(dst, dist));
      int match = CompareEqualTo(dst, dist);

      // will cause VS to make them all cmovs
      if (match & 0xF) { idxy = k; }
      if (match & 0x8) { idxx = j; } j--;
      if (match & 0x4) { idxx = j; } j--;
      if (match & 0x2) { idxx = j; } j--;
      if (match & 0x1) { idxx = j; } j--;
    } while (  j >= 0);
    } while (--k >= 0);

    // save the index
    indicesx[i] = (u8)idxx;
    indicesy[i] = (u8)idxy;

    // accumulate the error (sine)
    AddDeviance(dist, error);
  }

  // return the total error
  return error;
}

template<const int min, const int max, const int prc, typename otyp>
static Vec4 GetError(Vec4 const* xyz, int mask,
		     Col8 const &codes5x, Col8 const &codes7x,
		     Col8 const &codes5y, Col8 const &codes7y)
{
  // initial values
  Vec4 error = Vec4(DEVIANCE_BASE);

  for (int v = 0; v < 16; v++) {
    // find the closest code
    Vec4 dists = Vec4(DEVIANCE_MAX);

    // fetch floating point vector
    Vec4 valx = xyz[v].SplatX();
    Vec4 valy = xyz[v].SplatY();
    Vec4 valz = xyz[v].SplatZ();

    // for all possible codebook-entries
    
    // Col4(codes55y[g], codes77y[g], codes55y[g], codes77y[g]);
    int g = 7; Col8 _cy5 = codes5y, _cy7 = codes7y; do { Col4 _yy = InterleaveUpper(_cy5, _cy7, otyp(0)); _cy5 = ShiftUp<1>(_cy5); _cy7 = ShiftUp<1>(_cy7);
    // Col4(codes55x[f], codes55x[f], codes77x[f], codes77x[f]);
    int f = 7; Col8 _cx5 = codes5x, _cx7 = codes7x; do { Col4 _xx = ReplicateUpper (_cx5, _cx7, otyp(0)); _cx5 = ShiftUp<1>(_cx5); _cx7 = ShiftUp<1>(_cx7);

      Vec4 cxx = (_xx + Vec4(COFFSET)) * Vec4(CSCALE);
      Vec4 cyy = (_yy + Vec4(COFFSET)) * Vec4(CSCALE);
      Vec4 czz = Complement<DISARM>(cxx, cyy);

      // measure absolute angle-deviation (cosine)
      Vec4 dst = (valx * cxx) + (valy * cyy) + (valz * czz);

      // select the smallest deviation (NaN as first arg is ignored!)
      dists = Max(dst, dists);

    } while (--f >= 0);
    } while (--g >= 0);

    // accumulate the error (sine)
    AddDeviance(dists, error);
  }
  
  // return the total error
  return error;
}

#define FIT_THRESHOLD 1e-15f

template<const int min, const int max, const int prc, const int stepx, const int stepy>
static Scr4 FitError(Vec4 const* xyz, Col4 &minXY, Col4 &maxXY, Scr4 &errXY) {
  const Vec4 mprc = Vec4(1.0f * (1 << prc));

  // for 5-step
  const Vec4 fouf = Vec4(4.0f * (1 << prc) / 5.0f);
  const Vec4 thrf = Vec4(3.0f * (1 << prc) / 5.0f);
  const Vec4 twof = Vec4(2.0f * (1 << prc) / 5.0f);
  const Vec4 onef = Vec4(1.0f * (1 << prc) / 5.0f);
  
  // for 7-step
  const Vec4 sixs = Vec4(6.0f * (1 << prc) / 7.0f);
  const Vec4 fivs = Vec4(5.0f * (1 << prc) / 7.0f);
  const Vec4 fous = Vec4(4.0f * (1 << prc) / 7.0f);
  const Vec4 thrs = Vec4(3.0f * (1 << prc) / 7.0f);
  const Vec4 twos = Vec4(2.0f * (1 << prc) / 7.0f);
  const Vec4 ones = Vec4(1.0f * (1 << prc) / 7.0f);

  // for 5-step positions 1 & 6 become 0
  const Vec4 mix1 = Vec4(0.0f * (1 << prc) / 5.0f, 0.0f * (1 << prc) / 5.0f,   6.0f * (1 << prc) / 7.0f, 6.0f * (1 << prc) / 7.0f);
  const Vec4 mix2 = Vec4(4.0f * (1 << prc) / 5.0f, 4.0f * (1 << prc) / 5.0f,   5.0f * (1 << prc) / 7.0f, 5.0f * (1 << prc) / 7.0f);
  const Vec4 mix3 = Vec4(3.0f * (1 << prc) / 5.0f, 3.0f * (1 << prc) / 5.0f,   4.0f * (1 << prc) / 7.0f, 4.0f * (1 << prc) / 7.0f);
  const Vec4 mix4 = Vec4(2.0f * (1 << prc) / 5.0f, 2.0f * (1 << prc) / 5.0f,   3.0f * (1 << prc) / 7.0f, 3.0f * (1 << prc) / 7.0f);
  const Vec4 mix5 = Vec4(1.0f * (1 << prc) / 5.0f, 1.0f * (1 << prc) / 5.0f,   2.0f * (1 << prc) / 7.0f, 2.0f * (1 << prc) / 7.0f);
  const Vec4 mix6 = Vec4(0.0f * (1 << prc) / 5.0f, 0.0f * (1 << prc) / 5.0f,   1.0f * (1 << prc) / 7.0f, 1.0f * (1 << prc) / 7.0f);

  // for 5-step positions 1&6 become -1.0 & +1.0f
  const Vec4 scale  = Vec4(CSCALE);
  const Vec4 offset = Vec4(COFFSET);
  const Vec4 scmix  = Vec4( 1.0f,  1.0f, CSCALE , CSCALE );
  const Vec4 offmxn = Vec4( 1.0f,  1.0f, COFFSET, COFFSET);
  const Vec4 offmxp = Vec4(-1.0f, -1.0f, COFFSET, COFFSET);

#if	(SQUISH_USE_SIMD > 0)
  Vec4 minxy = minXY;
  Vec4 maxxy = maxXY;
  Vec4 rngxy = (maxXY - minXY) >> 1;	// max 127

  Scr4 errC = Scr4(errXY);
  Vec4 sx   = minxy.SplatX();
  Vec4 sy   = minxy.SplatY();
  Vec4 ex   = maxxy.SplatX();
  Vec4 ey   = maxxy.SplatY();
  Vec4 rx   = rngxy.SplatX() * Vec4(1.0f, -1.0f, 0.0f, 0.0f);
  Vec4 ry   = rngxy.SplatY() * Vec4(1.0f, -1.0f, 0.0f, 0.0f);
  Vec4 sxy  = Vec4(sx, sy, false, false);
  Vec4 exy  = Vec4(ex, ey, false, false);

  while ((rx + ry) > Vec4(0.0f)) {
    // s + os, s + os, s + os - r, s + os + r
    // e - oe - r, e - oe + r, e - oe, e - oe
    Vec4 cssmpx = sx + rx       ; cssmpx = Max(Min(cssmpx, Vec4(max)), Vec4(min));
    Vec4 cmpeex = ex + rx.Swap(); cmpeex = Max(Min(cmpeex, Vec4(max)), Vec4(min));

    Vec4 cssmpy = sy + ry       ; cssmpy = Max(Min(cssmpy, Vec4(max)), Vec4(min));
    Vec4 cmpeey = ey + ry.Swap(); cmpeey = Max(Min(cmpeey, Vec4(max)), Vec4(min));

    Vec4 cssss_ = sxy;
    Vec4 ceeee_ = exy;
    
    // construct unordered code-book
    Vec4 _cbs[8];
    Vec4 _cbx[8];
    Vec4 _cby[8];

    // for all possible codebook-entries of xy (XY==xy5, ZW==xy7)
    {
      _cbs[0] = (cssss_ * mprc)                  ; _cbs[0] = (Truncate(_cbs[0]) + offset) * scale;
      _cbs[1] = (cssss_ * mix1) + (ceeee_ * mix6); _cbs[1] = (Truncate(_cbs[1]) + offmxn) * scmix;
      _cbs[2] = (cssss_ * mix2) + (ceeee_ * mix5); _cbs[2] = (Truncate(_cbs[2]) + offset) * scale;
      _cbs[3] = (cssss_ * mix3) + (ceeee_ * mix4); _cbs[3] = (Truncate(_cbs[3]) + offset) * scale;
      _cbs[4] = (cssss_ * mix4) + (ceeee_ * mix3); _cbs[4] = (Truncate(_cbs[4]) + offset) * scale;
      _cbs[5] = (cssss_ * mix5) + (ceeee_ * mix2); _cbs[5] = (Truncate(_cbs[5]) + offset) * scale;
      _cbs[6] = (cssss_ * mix6) + (ceeee_ * mix1); _cbs[6] = (Truncate(_cbs[6]) + offmxp) * scmix;
      _cbs[7] =                   (ceeee_ * mprc); _cbs[7] = (Truncate(_cbs[7]) + offset) * scale;
    }

    // for all possible codebook-entries of x
    if (stepx == 5) {
      _cbx[0] = (cssmpx * mprc)                  ; _cbx[0] = (Truncate(_cbx[0]) + offset) * scale;
      _cbx[1] = (cssmpx * fouf) + (cmpeex * onef); _cbx[1] = (Truncate(_cbx[1]) + offset) * scale;
      _cbx[2] = (cssmpx * thrf) + (cmpeex * twof); _cbx[2] = (Truncate(_cbx[2]) + offset) * scale;
      _cbx[3] = (cssmpx * twof) + (cmpeex * thrf); _cbx[3] = (Truncate(_cbx[3]) + offset) * scale;
      _cbx[4] = (cssmpx * onef) + (cmpeex * fouf); _cbx[4] = (Truncate(_cbx[4]) + offset) * scale;
      _cbx[5] =                   (cmpeex * mprc); _cbx[5] = (Truncate(_cbx[5]) + offset) * scale;
      _cbx[6] = Vec4(-1.0f);
      _cbx[7] = Vec4( 1.0f);
    }
    else if (stepx == 7) {
      _cbx[0] = (cssmpx * mprc)                  ; _cbx[0] = (Truncate(_cbx[0]) + offset) * scale;
      _cbx[1] = (cssmpx * sixs) + (cmpeex * ones); _cbx[1] = (Truncate(_cbx[1]) + offset) * scale;
      _cbx[2] = (cssmpx * fivs) + (cmpeex * twos); _cbx[2] = (Truncate(_cbx[2]) + offset) * scale;
      _cbx[3] = (cssmpx * fous) + (cmpeex * thrs); _cbx[3] = (Truncate(_cbx[3]) + offset) * scale;
      _cbx[4] = (cssmpx * thrs) + (cmpeex * fous); _cbx[4] = (Truncate(_cbx[4]) + offset) * scale;
      _cbx[5] = (cssmpx * twos) + (cmpeex * fivs); _cbx[5] = (Truncate(_cbx[5]) + offset) * scale;
      _cbx[6] = (cssmpx * ones) + (cmpeex * sixs); _cbx[6] = (Truncate(_cbx[6]) + offset) * scale;
      _cbx[7] =                   (cmpeex * mprc); _cbx[7] = (Truncate(_cbx[7]) + offset) * scale;
    }

    // for all possible codebook-entries of y
    if (stepy == 5) {
      _cby[0] = (cssmpy * mprc)                  ; _cby[0] = (Truncate(_cby[0]) + offset) * scale;
      _cby[1] = (cssmpy * fouf) + (cmpeey * onef); _cby[1] = (Truncate(_cby[1]) + offset) * scale;
      _cby[2] = (cssmpy * thrf) + (cmpeey * twof); _cby[2] = (Truncate(_cby[2]) + offset) * scale;
      _cby[3] = (cssmpy * twof) + (cmpeey * thrf); _cby[3] = (Truncate(_cby[3]) + offset) * scale;
      _cby[4] = (cssmpy * onef) + (cmpeey * fouf); _cby[4] = (Truncate(_cby[4]) + offset) * scale;
      _cby[5] =                   (cmpeey * mprc); _cby[5] = (Truncate(_cby[5]) + offset) * scale;
      _cby[6] = Vec4(-1.0f);
      _cby[7] = Vec4( 1.0f);
    }
    else if (stepy == 7) {
      _cby[0] = (cssmpy * mprc)                  ; _cby[0] = (Truncate(_cby[0]) + offset) * scale;
      _cby[1] = (cssmpy * sixs) + (cmpeey * ones); _cby[1] = (Truncate(_cby[1]) + offset) * scale;
      _cby[2] = (cssmpy * fivs) + (cmpeey * twos); _cby[2] = (Truncate(_cby[2]) + offset) * scale;
      _cby[3] = (cssmpy * fous) + (cmpeey * thrs); _cby[3] = (Truncate(_cby[3]) + offset) * scale;
      _cby[4] = (cssmpy * thrs) + (cmpeey * fous); _cby[4] = (Truncate(_cby[4]) + offset) * scale;
      _cby[5] = (cssmpy * twos) + (cmpeey * fivs); _cby[5] = (Truncate(_cby[5]) + offset) * scale;
      _cby[6] = (cssmpy * ones) + (cmpeey * sixs); _cby[6] = (Truncate(_cby[6]) + offset) * scale;
      _cby[7] =                   (cmpeey * mprc); _cby[7] = (Truncate(_cby[7]) + offset) * scale;
    }

    Vec4 error0x  = Vec4(DEVIANCE_BASE);
    Vec4 error0y  = Vec4(DEVIANCE_BASE);
    Vec4 error1xy = Vec4(DEVIANCE_BASE);
    Vec4 error2xy = Vec4(DEVIANCE_BASE);
    Vec4 error3xy = Vec4(DEVIANCE_BASE);
    Vec4 error4xy = Vec4(DEVIANCE_BASE);
    for (int v = 0; v < 16; v++) {
      // find the closest code
      Vec4 valx = xyz[v].SplatX();
      Vec4 valy = xyz[v].SplatY();
      Vec4 valz = xyz[v].SplatZ();

      Vec4 dist0x  = Vec4(DEVIANCE_MAX);
      Vec4 dist0y  = Vec4(DEVIANCE_MAX);
      Vec4 dist1xy = Vec4(DEVIANCE_MAX);
      Vec4 dist2xy = Vec4(DEVIANCE_MAX);
      Vec4 dist3xy = Vec4(DEVIANCE_MAX);
      Vec4 dist4xy = Vec4(DEVIANCE_MAX);

      int j = 7; do {
      int i = 7; do {
	Vec4 _cxx = (stepx == 5 ? _cbs[i].SplatX() : _cbs[i].SplatZ());
	Vec4 _cyy = (stepy == 5 ? _cbs[j].SplatY() : _cbs[j].SplatW());

	//           ey -= ry	   =>   error0x1y -> _cb_[?].x _cby[?].w
	//           ey += ry	   =>   error0x2y -> _cb_[?].x _cby[?].z
	//           sy -= ry	   =>   error0x3y -> _cb_[?].x _cby[?].y
	//           sy += ry	   =>   error0x4y -> _cb_[?].x _cby[?].x
	Vec4 x0x  = _cxx;	        Vec4 y0x  = _cby[j];		Vec4 z0x  = Complement<DISARM>(x0x, y0x);

	// ex -= rx          	   =>   error1x0y -> _cbx[?].w _cb_[?].y
	// ex += rx          	   =>   error2x0y -> _cbx[?].z _cb_[?].y
	// sx -= rx          	   =>   error3x0y -> _cbx[?].y _cb_[?].y
	// sx += rx          	   =>   error4x0y -> _cbx[?].x _cb_[?].y
	Vec4 x0y  = _cbx[i];		Vec4 y0y  = _cyy;	        Vec4 z0y  = Complement<DISARM>(x0y, y0y);

	// ex -= rx, ey -= ry	   =>   error1x1y -> _cbx[?].w _cby[?].w
	// ex -= rx, ey += ry	   =>   error1x2y -> _cbx[?].w _cby[?].z
	// ex -= rx, sy -= ry	   =>   error1x3y -> _cbx[?].w _cby[?].y
	// ex -= rx, sy += ry	   =>   error1x4y -> _cbx[?].w _cby[?].x
	Vec4 x1xy = _cbx[i].SplatW();	Vec4 y1xy = _cby[j];		Vec4 z1xy = Complement<DISARM>(x1xy, y1xy);

	// ex += rx, ey -= ry	   =>   error2x1y -> _cbx[?].z _cby[?].w
	// ex += rx, ey += ry	   =>   error2x2y -> _cbx[?].z _cby[?].z
	// ex += rx, sy -= ry	   =>   error2x3y -> _cbx[?].z _cby[?].y
	// ex += rx, sy += ry	   =>   error2x4y -> _cbx[?].z _cby[?].x
	Vec4 x2xy = _cbx[i].SplatZ();	Vec4 y2xy = _cby[j];		Vec4 z2xy = Complement<DISARM>(x2xy, y2xy);

	// sx -= rx, ey -= ry	   =>   error3x1y -> _cbx[?].y _cby[?].w
	// sx -= rx, ey += ry	   =>   error3x2y -> _cbx[?].y _cby[?].z
	// sx -= rx, sy -= ry	   =>   error3x3y -> _cbx[?].y _cby[?].y
	// sx -= rx, sy += ry	   =>   error3x4y -> _cbx[?].y _cby[?].x
	Vec4 x3xy = _cbx[i].SplatY();	Vec4 y3xy = _cby[j];		Vec4 z3xy = Complement<DISARM>(x3xy, y3xy);

	// sx += rx, ey -= ry	   =>   error4x1y -> _cbx[?].x _cby[?].w
	// sx += rx, ey += ry	   =>   error4x2y -> _cbx[?].x _cby[?].z
	// sx += rx, sy -= ry	   =>   error4x3y -> _cbx[?].x _cby[?].y
	// sx += rx, sy += ry	   =>   error4x4y -> _cbx[?].x _cby[?].x
	Vec4 x4xy = _cbx[i].SplatX();	Vec4 y4xy = _cby[j];		Vec4 z4xy = Complement<DISARM>(x4xy, y4xy);

        // measure absolute angle-deviation (cosine)
	Vec4 d0x  = (valx * x0x ) + (valy * y0x ) + (valz * z0x );
	Vec4 d0y  = (valx * x0y ) + (valy * y0y ) + (valz * z0y );
	Vec4 d1xy = (valx * x1xy) + (valy * y1xy) + (valz * z1xy);
	Vec4 d2xy = (valx * x2xy) + (valy * y2xy) + (valz * z2xy);
	Vec4 d3xy = (valx * x3xy) + (valy * y3xy) + (valz * z3xy);
	Vec4 d4xy = (valx * x4xy) + (valy * y4xy) + (valz * z4xy);

	// select the smallest deviation (NaN as first arg is ignored!)
	dist0x  = Max(d0x , dist0x );
	dist0y  = Max(d0y , dist0y );
	dist1xy = Max(d1xy, dist1xy);
	dist2xy = Max(d2xy, dist2xy);
	dist3xy = Max(d3xy, dist3xy);
	dist4xy = Max(d4xy, dist4xy);
      } while(--i >= 0);
      } while(--j >= 0);

      // accumulate the error (sine)
      AddDeviance(dist0x , error0x );
      AddDeviance(dist0y , error0y );
      AddDeviance(dist1xy, error1xy);
      AddDeviance(dist2xy, error2xy);
      AddDeviance(dist3xy, error3xy);
      AddDeviance(dist4xy, error4xy);
    }

    // encourage OoO
    Vec4 error0  = Min(error0x , error0y );
    Vec4 error14 = Min(error1xy, error4xy);
    Vec4 error23 = Min(error2xy, error3xy);
    Vec4 errors  = Min(error0, Min(error14, error23));
    Scr4 merr    = HorizontalMin(errors);

    if (!(merr < errC)) {
      // half range
      if (rx > ry)
	rx = Truncate(rx * Vec4(0.5f));
      else
	ry = Truncate(ry * Vec4(0.5f));

      continue;
    }
    else if (CompareEqualTo(error0 , merr) != 0x00) {
      int value0x = CompareEqualTo(error0x, merr);
      int value0y = CompareEqualTo(error0y, merr);

      /**/ if (value0x & 0x8) ey = cmpeey.SplatW();	//           ey -= ry;// range up
      else if (value0x & 0x1) sy = cssmpy.SplatX();	//           sy += ry;// range dn
      else if (value0y & 0x8) ex = cmpeex.SplatW();	// ex -= rx          ;// range up
      else if (value0y & 0x1) sx = cssmpx.SplatX();	// sx += rx          ;// range dn
      else if (value0x & 0x4) ey = cmpeey.SplatZ();	//           ey += ry;// range up
      else if (value0x & 0x2) sy = cssmpy.SplatY();	//           sy -= ry;// range dn
      else if (value0y & 0x4) ex = cmpeex.SplatZ();	// ex += rx          ;// range up
      else /* (value0y & 2)*/ sx = cssmpx.SplatY();	// sx -= rx          ;// range dn
    }
    else if (CompareEqualTo(error14, merr) != 0x00) {
      int value1xy = CompareEqualTo(error1xy, merr);
      int value14  = CompareEqualTo(error14 , merr);

      /**/ if (value1xy != 0) ex = cmpeex.SplatW();	// ex -= rx, ........;// range up
      else /* (value4xy != */ sx = cssmpx.SplatX();	// sx += rx, ........;// range dn

      /**/ if (value14 & 0x8) ey = cmpeey.SplatW();	// ........, ey -= ry;// range +-
      else if (value14 & 0x1) sy = cssmpy.SplatX();	// ........, sy += ry;// range +-
      else if (value14 & 0x4) ey = cmpeey.SplatZ();	// ........, ey += ry;// range +-
      else /* (value14 & 2)*/ sy = cssmpy.SplatY();	// ........, sy -= ry;// range +-
    }
    else /*if (CompareEqualTo(error23, merr) != 0x00)*/ {
      int value2xy = CompareEqualTo(error2xy, merr);
      int value23  = CompareEqualTo(error23 , merr);

      /**/ if (value2xy != 0) ex = cmpeex.SplatZ();	// ex += rx, ........;// range up
      else /* (value3xy != */ sx = cssmpx.SplatY();	// sx -= rx, ........;// range dn

      /**/ if (value23 & 0x8) ey = cmpeey.SplatW();	// ........, ey -= ry;// range +-
      else if (value23 & 0x1) sy = cssmpy.SplatX();	// ........, sy += ry;// range +-
      else if (value23 & 0x4) ey = cmpeey.SplatZ();	// ........, ey += ry;// range +-
      else /* (value23 & 2)*/ sy = cssmpy.SplatY();	// ........, sy -= ry;// range +-
    }

    sxy  = Vec4(sx, sy, false, false);
    exy  = Vec4(ex, ey, false, false);
    errC = merr;
    
    /*
#if defined(TRACK_STATISTICS)
    static int counter = 0; counter++;
    static int matches[2][128] = {0};
    matches[0][(int)rx.X()]++;
    matches[1][(int)ry.X()]++;
#endif
     */

    // lossless
    if (!(errC > Scr4(FIT_THRESHOLD)))
      break;
  }

#if defined(TRACK_STATISTICS)
  /* there is a clear skew towards unweighted clusterfit (all weights == 1)
    *
    * C == 3, numset ==
    *  [0]	0x00f796e0 {124800, 15616}
    */
  if (steps == 5) {
    gstat.coord[0] += errXY;
    gstat.coord[2] += 1;
  }
  else if (steps == 7) {
    gstat.coord[3] += errXY;
    gstat.coord[5] += 1;
  }
#endif

  // final match
  minXY = FloatToInt<false>(sxy);
  maxXY = FloatToInt<false>(exy);
  errXY = errC;

#if defined(TRACK_STATISTICS)
  /* there is a clear skew towards unweighted clusterfit (all weights == 1)
    *
    * C == 3, numset ==
    *  [0]	0x00f796e0 {124800, 15616}
    */
  if (steps == 5)
    gstat.coord[1] += errXY;
  else if (steps == 7)
    gstat.coord[4] += errXY;
#endif
#else
  Scr4 error0x0y = Scr4(errXY);

  int sx = minXY.R();
  int ex = maxXY.R();
  int rx = (ex - sx) >> 1;	// max 127

  int sy = minXY.G();
  int ey = maxXY.G();
  int ry = (ey - sy) >> 1;	// max 127

  while ((rx != 0) || (ry != 0)) {
    int msx = std::max(std::min(sx - rx, max), min);  // os - r
    int csx = std::max(std::min(sx     , max), min);  // os
    int psx = std::max(std::min(sx + rx, max), min);  // os + r
    int mex = std::max(std::min(ex - rx, max), min);  // oe + r
    int cex = std::max(std::min(ex     , max), min);  // oe
    int pex = std::max(std::min(ex + rx, max), min);  // oe - r

    int msy = std::max(std::min(sy - ry, max), min);  // os - r
    int csy = std::max(std::min(sy     , max), min);  // os
    int psy = std::max(std::min(sy + ry, max), min);  // os + r
    int mey = std::max(std::min(ey - ry, max), min);  // oe + r
    int cey = std::max(std::min(ey     , max), min);  // oe
    int pey = std::max(std::min(ey + ry, max), min);  // oe - r

    // construct code-books
    Col8 cb0y, cb1y, cb2y, cb3y, cb4y;
    Col8 cb0x, cb1x, cb2x, cb3x, cb4x;

    if (stepy == 5) {
      Codebook6<min,max,prc>(cb0y, Col8(csy), Col8(cey));
      Codebook6<min,max,prc>(cb1y, Col8(csy), Col8(mey));
      Codebook6<min,max,prc>(cb2y, Col8(csy), Col8(pey));
      Codebook6<min,max,prc>(cb3y, Col8(msy), Col8(cey));
      Codebook6<min,max,prc>(cb4y, Col8(psy), Col8(cey));
    }
    else if (stepy == 7) {
      Codebook8<min,max,prc>(cb0y, Col8(csy), Col8(cey));
      Codebook8<min,max,prc>(cb1y, Col8(csy), Col8(mey));
      Codebook8<min,max,prc>(cb2y, Col8(csy), Col8(pey));
      Codebook8<min,max,prc>(cb3y, Col8(msy), Col8(cey));
      Codebook8<min,max,prc>(cb4y, Col8(psy), Col8(cey));
    }

    if (stepx == 5) {
      Codebook6<min,max,prc>(cb0x, Col8(csx), Col8(cex));
      Codebook6<min,max,prc>(cb1x, Col8(csx), Col8(mex));
      Codebook6<min,max,prc>(cb2x, Col8(csx), Col8(pex));
      Codebook6<min,max,prc>(cb3x, Col8(msx), Col8(cex));
      Codebook6<min,max,prc>(cb4x, Col8(psx), Col8(cex));
    }
    else if (stepx == 7) {
      Codebook8<min,max,prc>(cb0x, Col8(csx), Col8(cex));
      Codebook8<min,max,prc>(cb1x, Col8(csx), Col8(mex));
      Codebook8<min,max,prc>(cb2x, Col8(csx), Col8(pex));
      Codebook8<min,max,prc>(cb3x, Col8(msx), Col8(cex));
      Codebook8<min,max,prc>(cb4x, Col8(psx), Col8(cex));
    }

    Scr4 error0x1y, error0x2y, error0x3y, error0x4y;
    Scr4 error1x0y, error2x0y, error3x0y, error4x0y;
    Scr4 error1x1y, error1x2y, error1x3y, error1x4y;
    Scr4 error2x1y, error2x2y, error2x3y, error2x4y;
    Scr4 error3x1y, error3x2y, error3x3y, error3x4y;
    Scr4 error4x1y, error4x2y, error4x3y, error4x4y;

    error0x1y = error0x2y = error0x3y = error0x4y =
    error1x0y = error2x0y = error3x0y = error4x0y =
    error1x1y = error1x2y = error1x3y = error1x4y =
    error2x1y = error2x2y = error2x3y = error2x4y =
    error3x1y = error3x2y = error3x3y = error3x4y =
    error4x1y = error4x2y = error4x3y = error4x4y = Scr4(DEVIANCE_BASE);

    for (int v = 0; v < 16; v++) {
      // find the closest code
      Scr4 dist0x1y, dist0x2y, dist0x3y, dist0x4y;
      Scr4 dist1x0y, dist2x0y, dist3x0y, dist4x0y;
      Scr4 dist1x1y, dist1x2y, dist1x3y, dist1x4y;
      Scr4 dist2x1y, dist2x2y, dist2x3y, dist2x4y;
      Scr4 dist3x1y, dist3x2y, dist3x3y, dist3x4y;
      Scr4 dist4x1y, dist4x2y, dist4x3y, dist4x4y;

      dist0x1y = dist0x2y = dist0x3y = dist0x4y =
      dist1x0y = dist2x0y = dist3x0y = dist4x0y =
      dist1x1y = dist1x2y = dist1x3y = dist1x4y =
      dist2x1y = dist2x2y = dist2x3y = dist2x4y =
      dist3x1y = dist3x2y = dist3x3y = dist3x4y =
      dist4x1y = dist4x2y = dist4x3y = dist4x4y = Scr4(DEVIANCE_MAX);

      // fetch floating point vector
      Vec4 value = xyz[v];

      // for all possible codebook-entries
      for (int g = 0; g < 8; g++) {
      for (int f = 0; f < 8; f++) {
	// complement z
	Col4 _xyz0x1y = Col4(cb0x[f], cb1y[g]); Vec4 cxyz0x1y = scale * (offset + _xyz0x1y); cxyz0x1y = Complement<true>(cxyz0x1y);
	Col4 _xyz0x2y = Col4(cb0x[f], cb2y[g]); Vec4 cxyz0x2y = scale * (offset + _xyz0x2y); cxyz0x2y = Complement<true>(cxyz0x2y);
	Col4 _xyz0x3y = Col4(cb0x[f], cb3y[g]); Vec4 cxyz0x3y = scale * (offset + _xyz0x3y); cxyz0x3y = Complement<true>(cxyz0x3y);
	Col4 _xyz0x4y = Col4(cb0x[f], cb4y[g]); Vec4 cxyz0x4y = scale * (offset + _xyz0x4y); cxyz0x4y = Complement<true>(cxyz0x4y);
	Col4 _xyz1x1y = Col4(cb1x[f], cb1y[g]); Vec4 cxyz1x1y = scale * (offset + _xyz1x1y); cxyz1x1y = Complement<true>(cxyz1x1y);
	Col4 _xyz1x2y = Col4(cb1x[f], cb2y[g]); Vec4 cxyz1x2y = scale * (offset + _xyz1x2y); cxyz1x2y = Complement<true>(cxyz1x2y);
	Col4 _xyz1x3y = Col4(cb1x[f], cb3y[g]); Vec4 cxyz1x3y = scale * (offset + _xyz1x3y); cxyz1x3y = Complement<true>(cxyz1x3y);
	Col4 _xyz1x4y = Col4(cb1x[f], cb4y[g]); Vec4 cxyz1x4y = scale * (offset + _xyz1x4y); cxyz1x4y = Complement<true>(cxyz1x4y);
	Col4 _xyz2x1y = Col4(cb2x[f], cb1y[g]); Vec4 cxyz2x1y = scale * (offset + _xyz2x1y); cxyz2x1y = Complement<true>(cxyz2x1y);
	Col4 _xyz2x2y = Col4(cb2x[f], cb2y[g]); Vec4 cxyz2x2y = scale * (offset + _xyz2x2y); cxyz2x2y = Complement<true>(cxyz2x2y);
	Col4 _xyz2x3y = Col4(cb2x[f], cb3y[g]); Vec4 cxyz2x3y = scale * (offset + _xyz2x3y); cxyz2x3y = Complement<true>(cxyz2x3y);
	Col4 _xyz2x4y = Col4(cb2x[f], cb4y[g]); Vec4 cxyz2x4y = scale * (offset + _xyz2x4y); cxyz2x4y = Complement<true>(cxyz2x4y);
	Col4 _xyz3x1y = Col4(cb3x[f], cb1y[g]); Vec4 cxyz3x1y = scale * (offset + _xyz3x1y); cxyz3x1y = Complement<true>(cxyz3x1y);
	Col4 _xyz3x2y = Col4(cb3x[f], cb2y[g]); Vec4 cxyz3x2y = scale * (offset + _xyz3x2y); cxyz3x2y = Complement<true>(cxyz3x2y);
	Col4 _xyz3x3y = Col4(cb3x[f], cb3y[g]); Vec4 cxyz3x3y = scale * (offset + _xyz3x3y); cxyz3x3y = Complement<true>(cxyz3x3y);
	Col4 _xyz3x4y = Col4(cb3x[f], cb4y[g]); Vec4 cxyz3x4y = scale * (offset + _xyz3x4y); cxyz3x4y = Complement<true>(cxyz3x4y);
	Col4 _xyz4x1y = Col4(cb4x[f], cb1y[g]); Vec4 cxyz4x1y = scale * (offset + _xyz4x1y); cxyz4x1y = Complement<true>(cxyz4x1y);
	Col4 _xyz4x2y = Col4(cb4x[f], cb2y[g]); Vec4 cxyz4x2y = scale * (offset + _xyz4x2y); cxyz4x2y = Complement<true>(cxyz4x2y);
	Col4 _xyz4x3y = Col4(cb4x[f], cb3y[g]); Vec4 cxyz4x3y = scale * (offset + _xyz4x3y); cxyz4x3y = Complement<true>(cxyz4x3y);
	Col4 _xyz4x4y = Col4(cb4x[f], cb4y[g]); Vec4 cxyz4x4y = scale * (offset + _xyz4x4y); cxyz4x4y = Complement<true>(cxyz4x4y);
	Col4 _xyz1x0y = Col4(cb1x[f], cb0y[g]); Vec4 cxyz1x0y = scale * (offset + _xyz1x0y); cxyz1x0y = Complement<true>(cxyz1x0y);
	Col4 _xyz2x0y = Col4(cb2x[f], cb0y[g]); Vec4 cxyz2x0y = scale * (offset + _xyz2x0y); cxyz2x0y = Complement<true>(cxyz2x0y);
	Col4 _xyz3x0y = Col4(cb3x[f], cb0y[g]); Vec4 cxyz3x0y = scale * (offset + _xyz3x0y); cxyz3x0y = Complement<true>(cxyz3x0y);
	Col4 _xyz4x0y = Col4(cb4x[f], cb0y[g]); Vec4 cxyz4x0y = scale * (offset + _xyz4x0y); cxyz4x0y = Complement<true>(cxyz4x0y);

        // measure absolute angle-deviation (cosine)
	Scr4 d0x1y = Dot(value, cxyz0x1y);
	Scr4 d0x2y = Dot(value, cxyz0x2y);
	Scr4 d0x3y = Dot(value, cxyz0x3y);
	Scr4 d0x4y = Dot(value, cxyz0x4y);
	Scr4 d1x1y = Dot(value, cxyz1x1y);
	Scr4 d1x2y = Dot(value, cxyz1x2y);
	Scr4 d1x3y = Dot(value, cxyz1x3y);
	Scr4 d1x4y = Dot(value, cxyz1x4y);
	Scr4 d2x1y = Dot(value, cxyz2x1y);
	Scr4 d2x2y = Dot(value, cxyz2x2y);
	Scr4 d2x3y = Dot(value, cxyz2x3y);
	Scr4 d2x4y = Dot(value, cxyz2x4y);
	Scr4 d3x1y = Dot(value, cxyz3x1y);
	Scr4 d3x2y = Dot(value, cxyz3x2y);
	Scr4 d3x3y = Dot(value, cxyz3x3y);
	Scr4 d3x4y = Dot(value, cxyz3x4y);
	Scr4 d4x1y = Dot(value, cxyz4x1y);
	Scr4 d4x2y = Dot(value, cxyz4x2y);
	Scr4 d4x3y = Dot(value, cxyz4x3y);
	Scr4 d4x4y = Dot(value, cxyz4x4y);
	Scr4 d1x0y = Dot(value, cxyz1x0y);
	Scr4 d2x0y = Dot(value, cxyz2x0y);
	Scr4 d3x0y = Dot(value, cxyz3x0y);
	Scr4 d4x0y = Dot(value, cxyz4x0y);

	// select the smallest deviation
	if (dist0x1y < d0x1y) dist0x1y = d0x1y;
	if (dist0x2y < d0x2y) dist0x2y = d0x2y;
	if (dist0x3y < d0x3y) dist0x3y = d0x3y;
	if (dist0x4y < d0x4y) dist0x4y = d0x4y;
	if (dist1x1y < d1x1y) dist1x1y = d1x1y;
	if (dist1x2y < d1x2y) dist1x2y = d1x2y;
	if (dist1x3y < d1x3y) dist1x3y = d1x3y;
	if (dist1x4y < d1x4y) dist1x4y = d1x4y;
	if (dist2x1y < d2x1y) dist2x1y = d2x1y;
	if (dist2x2y < d2x2y) dist2x2y = d2x2y;
	if (dist2x3y < d2x3y) dist2x3y = d2x3y;
	if (dist2x4y < d2x4y) dist2x4y = d2x4y;
	if (dist3x1y < d3x1y) dist3x1y = d3x1y;
	if (dist3x2y < d3x2y) dist3x2y = d3x2y;
	if (dist3x3y < d3x3y) dist3x3y = d3x3y;
	if (dist3x4y < d3x4y) dist3x4y = d3x4y;
	if (dist4x1y < d4x1y) dist4x1y = d4x1y;
	if (dist4x2y < d4x2y) dist4x2y = d4x2y;
	if (dist4x3y < d4x3y) dist4x3y = d4x3y;
	if (dist4x4y < d4x4y) dist4x4y = d4x4y;
	if (dist1x0y < d1x0y) dist1x0y = d1x0y;
	if (dist2x0y < d2x0y) dist2x0y = d2x0y;
	if (dist3x0y < d3x0y) dist3x0y = d3x0y;
	if (dist4x0y < d4x0y) dist4x0y = d4x0y;
      }
      }

      // accumulate the error (sine)
      AddDeviance(dist0x1y, error0x1y);
      AddDeviance(dist0x2y, error0x2y);
      AddDeviance(dist0x3y, error0x3y);
      AddDeviance(dist0x4y, error0x4y);
      AddDeviance(dist1x1y, error1x1y);
      AddDeviance(dist1x2y, error1x2y);
      AddDeviance(dist1x3y, error1x3y);
      AddDeviance(dist1x4y, error1x4y);
      AddDeviance(dist2x1y, error2x1y);
      AddDeviance(dist2x2y, error2x2y);
      AddDeviance(dist2x3y, error2x3y);
      AddDeviance(dist2x4y, error2x4y);
      AddDeviance(dist3x1y, error3x1y);
      AddDeviance(dist3x2y, error3x2y);
      AddDeviance(dist3x3y, error3x3y);
      AddDeviance(dist3x4y, error3x4y);
      AddDeviance(dist4x1y, error4x1y);
      AddDeviance(dist4x2y, error4x2y);
      AddDeviance(dist4x3y, error4x3y);
      AddDeviance(dist4x4y, error4x4y);
      AddDeviance(dist1x0y, error1x0y);
      AddDeviance(dist2x0y, error2x0y);
      AddDeviance(dist3x0y, error3x0y);
      AddDeviance(dist4x0y, error4x0y);
    }

    Scr4                   merrx = error0x0y;
    if (merrx > error0x1y) merrx = error0x1y;
    if (merrx > error0x4y) merrx = error0x4y;
    if (merrx > error1x0y) merrx = error1x0y;
    if (merrx > error4x0y) merrx = error4x0y;
    if (merrx > error0x2y) merrx = error0x2y;
    if (merrx > error0x3y) merrx = error0x3y;
    if (merrx > error2x0y) merrx = error2x0y;
    if (merrx > error3x0y) merrx = error3x0y;
    if (merrx > error1x1y) merrx = error1x1y;
    if (merrx > error1x4y) merrx = error1x4y;
    if (merrx > error1x2y) merrx = error1x2y;
    if (merrx > error1x3y) merrx = error1x3y;
    if (merrx > error4x1y) merrx = error4x1y;
    if (merrx > error4x4y) merrx = error4x4y;
    if (merrx > error4x2y) merrx = error4x2y;
    if (merrx > error4x3y) merrx = error4x3y;
    if (merrx > error2x1y) merrx = error2x1y;
    if (merrx > error2x4y) merrx = error2x4y;
    if (merrx > error2x2y) merrx = error2x2y;
    if (merrx > error2x3y) merrx = error2x3y;
    if (merrx > error3x1y) merrx = error3x1y;
    if (merrx > error3x4y) merrx = error3x4y;
    if (merrx > error3x2y) merrx = error3x2y;
    if (merrx > error3x3y) merrx = error3x3y;

//  /**/ if (merrx == error0x0y) rx >>= 1, ry >>= 1;				// half range
    /**/ if (merrx == error0x0y) rx > ry ? rx >>= 1 : ry >>= 1;
    else if (merrx == error0x1y) 	   ey = mey;//          ey -= ry;	// range up
    else if (merrx == error0x4y) 	   sy = psy;//          sy += ry;	// range up
    else if (merrx == error1x0y) ex = mex;	    //ex -= rx          ;	// range up
    else if (merrx == error4x0y) sx = psx;	    //sx += rx          ;	// range up
    else if (merrx == error0x2y) 	   ey = pey;//          ey += ry;	// range up
    else if (merrx == error0x3y) 	   sy = msy;//          sy -= ry;	// range up
    else if (merrx == error2x0y) ex = pex;	    //ex += rx          ;	// range up
    else if (merrx == error3x0y) sx = msx;	    //sx -= rx          ;	// range up
    else if (merrx == error1x1y) ex = mex, ey = mey;//ex -= rx, ey -= ry;	// range up
    else if (merrx == error1x4y) ex = mex, sy = psy;//ex -= rx, sy += ry;	// range up
    else if (merrx == error1x2y) ex = mex, ey = pey;//ex -= rx, ey += ry;	// range up
    else if (merrx == error1x3y) ex = mex, sy = msy;//ex -= rx, sy -= ry;	// range up
    else if (merrx == error4x1y) sx = psx, ey = mey;//sx += rx, ey -= ry;	// range dn
    else if (merrx == error4x4y) sx = psx, sy = psy;//sx += rx, sy += ry;	// range dn
    else if (merrx == error4x2y) sx = psx, ey = pey;//sx += rx, ey += ry;	// range dn
    else if (merrx == error4x3y) sx = psx, sy = msy;//sx += rx, sy -= ry;	// range dn
    else if (merrx == error2x1y) ex = pex, ey = mey;//ex += rx, ey -= ry;	// range up
    else if (merrx == error2x4y) ex = pex, sy = psy;//ex += rx, sy += ry;	// range up
    else if (merrx == error2x2y) ex = pex, ey = pey;//ex += rx, ey += ry;	// range up
    else if (merrx == error2x3y) ex = pex, sy = msy;//ex += rx, sy -= ry;	// range up
    else if (merrx == error3x1y) sx = msx, ey = mey;//sx -= rx, ey -= ry;	// range dn
    else if (merrx == error3x4y) sx = msx, sy = psy;//sx -= rx, sy += ry;	// range dn
    else if (merrx == error3x2y) sx = msx, ey = pey;//sx -= rx, ey += ry;	// range dn
    else if (merrx == error3x3y) sx = msx, sy = msy;//sx -= rx, sy -= ry;	// range dn

    // lossless
    error0x0y = merrx;
    if (!(error0x0y > Scr4(FIT_THRESHOLD)))
      break;
  }

#if defined(TRACK_STATISTICS)
  /* way better!
    * [0]	17302780	int
    * [1]	3868483		int
    * [2]	69408		int
    */
  if (steps == 5) {
    gstat.coord[0] += errXY;
    gstat.coord[2] += 1;
  }
  else if (steps == 7) {
    gstat.coord[3] += errXY;
    gstat.coord[5] += 1;
  }
#endif

  // final match
  minXY = Col4(sx, sy);
  maxXY = Col4(ex, ey);
  errXY = error0x0y;

#if defined(TRACK_STATISTICS)
  if (steps == 5)
    gstat.coord[1] += errXY;
  else if (steps == 7)
    gstat.coord[4] += errXY;
#endif
#endif

  return errXY;
}

/* -----------------------------------------------------------------------------
 */
static void WriteNormalBlock(int coord0, int coord1, u8 const* indices, void* block)
{
  u8* bytes = reinterpret_cast< u8* >(block);

  // write the first two bytes
  bytes[0] = (u8)coord0;
  bytes[1] = (u8)coord1;

  // pack the indices with 3 bits each
  u8* dest = bytes + 2;
  u8 const* src = indices;
  for (int i = 0; i < 2; ++i) {
    // pack 8 3-bit values
    int value = 0;
    for (int j = 0; j < 8; ++j) {
      int index = *src++;
      value += (index << (3 * j));
    }

    // store in 3 bytes
    for (int j = 0; j < 3; ++j) {
      int byte = (value >> (8 * j)) & 0xFF;
      *dest++ = (u8)byte;
    }

    // 77766655.54443332.22111000, FFFEEEDD.DCCCBBBA.AA999888
    // 22111000.54443332.77766655, AA999888.DCCCBBBA.FFFEEEDD
  }
}

static doinline void WriteNormalBlock5(int coord0, int coord1, u8 const* indices, void* block)
{
  if (coord0 == coord1) {
    assert((indices[ 0] % 6) < 2); assert((indices[ 1] % 6) < 2);
    assert((indices[ 2] % 6) < 2); assert((indices[ 3] % 6) < 2);
    assert((indices[ 4] % 6) < 2); assert((indices[ 5] % 6) < 2);
    assert((indices[ 6] % 6) < 2); assert((indices[ 7] % 6) < 2);
    assert((indices[ 8] % 6) < 2); assert((indices[ 9] % 6) < 2);
    assert((indices[10] % 6) < 2); assert((indices[11] % 6) < 2);
    assert((indices[12] % 6) < 2); assert((indices[13] % 6) < 2);
    assert((indices[14] % 6) < 2); assert((indices[15] % 6) < 2);
  }

  // write the block
  WriteNormalBlock(coord0, coord1, indices, block);
}

static doinline void WriteNormalBlock7(int coord0, int coord1, u8 const* indices, void* block)
{
  if (coord0 == coord1) {
    assert(indices[ 0] < 2); assert(indices[ 1] < 2);
    assert(indices[ 2] < 2); assert(indices[ 3] < 2);
    assert(indices[ 4] < 2); assert(indices[ 5] < 2);
    assert(indices[ 6] < 2); assert(indices[ 7] < 2);
    assert(indices[ 8] < 2); assert(indices[ 9] < 2);
    assert(indices[10] < 2); assert(indices[11] < 2);
    assert(indices[12] < 2); assert(indices[13] < 2);
    assert(indices[14] < 2); assert(indices[15] < 2);
  }

  // write the block
  WriteNormalBlock(coord0, coord1, indices, block);
}

/* -----------------------------------------------------------------------------
 */

template<const int min, const int max, const int prc, typename otyp, typename dtyp>
static void CompressNormalBtc5v(Vec4 const* xyz, int mask, void* blockx, void* blocky, int flags)
{
  Col8 codes55x, codes57x, codes75x, codes77x;
  Col8 codes55y, codes57y, codes75y, codes77y;

  // get the range for 5-coord and 7-coord interpolation
  Col4 min55 = Col4(max), min57;
  Col4 max55 = Col4(min), max57;
  Col4 min77 = Col4(max), min75;
  Col4 max77 = Col4(min), max75;

  if (!((flags & kNormalIterativeFit) && prc)) {
    for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
      // check this pixel is valid
      if ((imask & 1) == 0)
        continue;
    
      Vec4 norm = xyz[i];
#ifdef PYRAMID_PROJECTION
      Vec4 nabs = Abs(norm);
      Vec4 len = HorizontalMaxXY(nabs) + nabs.SplatZ();

      norm /= len;
#endif
    
      // incorporate into the min/max
      {
        // create 8bit integer vector
        Col4 value = 
          FloatToInt<true>((norm * Vec4(FSCALE)) - Vec4(FOFFSET));
        
        Col4 vmask =
          IsValue<min>(value) |
          IsValue<max>(value);
      
        min77 = Min(min77, value);
        max77 = Max(max77, value);
        min55 = Min(min55, (vmask % value) + (vmask & min55));
        max55 = Max(max55, (vmask % value) + (vmask & max55));
      }
    }
  }
  else {
    for (int i = 0, imask = mask; i < 16; ++i, imask >>= 1) {
      // check this pixel is valid
      if ((imask & 1) == 0)
        continue;
    
      Vec4 norm = xyz[i];
#ifdef PYRAMID_PROJECTION
      Vec4 nabs = Abs(norm);
      Vec4 len = HorizontalMaxXY(nabs) + nabs.SplatZ();

      norm /= len;
#endif
    
      // calculate fudge value which guarantees "255.0f + fudge" round to 255 (tie-even makes 255.5f become 256)
      const float fudge = 0.5f - 0.5f / (1 << prc);

      // incorporate into the min/max
      {
        // create 8bit integer vector, do floor() with round-nearest-tie-even() fudged
        Col4 value = 
          FloatToInt<true>((norm * Vec4(FSCALE)) - Vec4(FOFFSET + fudge));
        
        Col4 vmask =
          IsValue<min>(value) |
          IsValue<max>(value);
      
        min77 = Min(min77, value);
        min55 = Min(min55, (vmask % value) + (vmask & min55));
      }
      
      {
        // create 8bit integer vector, do ceil() with round-nearest-tie-even() fudged
        Col4 value = 
          FloatToInt<true>((norm * Vec4(FSCALE)) - Vec4(FOFFSET - fudge));
        
        Col4 vmask =
          IsValue<min>(value) |
          IsValue<max>(value);
      
        max77 = Max(max77, value);
        max55 = Max(max55, (vmask % value) + (vmask & max55));
      }
    }
  }

  // handle the case that no valid range was found
  min55 = Min(min55, max55);
  min77 = Min(min77, max77);
  
  // expand empty range in 5-code
  Col4 eq55 =  CompareAllEqualTo_M4(min55, max55);
  min55 += CompareAllGreaterThan_M4(min55, Col4(min)) & eq55;
  max55 -=    CompareAllLessThan_M4(max55, Col4(max)) & eq55;
  
  // fix the range to be the minimum in each case
//FixRange(min55, max55, 5);
//FixRange(min77, max77, 7);

  // construct code-book at desired precision
  Codebook6<min,max,prc>(codes55x, Col8(min55.SplatR()), Col8(max55.SplatR()));
  Codebook6<min,max,prc>(codes55y, Col8(min55.SplatG()), Col8(max55.SplatG()));
  Codebook8<min,max,prc>(codes77x, Col8(max77.SplatR()), Col8(min77.SplatR()));
  Codebook8<min,max,prc>(codes77y, Col8(max77.SplatG()), Col8(min77.SplatG()));

  // construct cross-over ranges/codebooks
  // min57 = Col4(min55.R(), min77.G());
  // max57 = Col4(max55.R(), max77.G());
  // min75 = Col4(min77.R(), min55.G());
  // max75 = Col4(max77.R(), max55.G());
  Interleave(min57, min55, min77); codes57x = codes55x;
  Interleave(max57, max55, max77); codes57y = codes77y;
  Interleave(min75, min77, min55); codes75x = codes77x;
  Interleave(max75, max77, max55); codes75y = codes55y;

  // do the iterative tangent search
  if (flags & kNormalIterativeFit) {
    // get best error for the current min/max
    Vec4 error = GetError<min,max,prc,otyp>(xyz, mask, codes55x, codes77x, codes55y, codes77y);

    // binary search, tangent-fitting
    Scr4 err55 = error.SplatX();
    Scr4 err57 = error.SplatY();
    Scr4 err75 = error.SplatZ();
    Scr4 err77 = error.SplatW();
    Scr4 errM0x0y = Min(Min(err55, err77), Min(err57, err75));

    // !lossless
    if (errM0x0y > Scr4(FIT_THRESHOLD)) {
      errM0x0y = FitError<min,max,prc,7,7>(xyz, min77, max77, err77);

      // if better, reconstruct code-book
      Codebook8<min,max,prc>(codes77x, Col8(max77.SplatR()), Col8(min77.SplatR()));
      Codebook8<min,max,prc>(codes77y, Col8(max77.SplatG()), Col8(min77.SplatG()));

    // !lossless
    if (errM0x0y > Scr4(FIT_THRESHOLD)) {
      errM0x0y = FitError<min,max,prc,7,5>(xyz, min75, max75, err75);

      // if better, reconstruct code-book
      Codebook8<min,max,prc>(codes75x, Col8(max75.SplatR()), Col8(min75.SplatR()));
      Codebook6<min,max,prc>(codes75y, Col8(min75.SplatG()), Col8(max75.SplatG()));

    // !lossless
    if (errM0x0y > Scr4(FIT_THRESHOLD)) {
      errM0x0y = FitError<min,max,prc,5,7>(xyz, min57, max57, err57);

      // if better, reconstruct code-book
      Codebook6<min,max,prc>(codes57x, Col8(min57.SplatR()), Col8(max57.SplatR()));
      Codebook8<min,max,prc>(codes57y, Col8(max57.SplatG()), Col8(min57.SplatG()));

    // !lossless
    if (errM0x0y > Scr4(FIT_THRESHOLD)) {
      errM0x0y = FitError<min,max,prc,5,5>(xyz, min55, max55, err55);

      // if better, reconstruct code-book
      Codebook6<min,max,prc>(codes55x, Col8(min55.SplatR()), Col8(max55.SplatR()));
      Codebook6<min,max,prc>(codes55y, Col8(min55.SplatG()), Col8(max55.SplatG()));

    }}}}
  }

  // fit the data to both code books
  u8 indices55x[16], indices57x[16], indices75x[16], indices77x[16];
  u8 indices55y[16], indices57y[16], indices75y[16], indices77y[16];

  Scr4 err55 = FitCodes<min,max,prc,otyp>(xyz, mask, codes55x, indices55x, codes55y, indices55y);
  Scr4 err57 = FitCodes<min,max,prc,otyp>(xyz, mask, codes57x, indices57x, codes57y, indices57y);
  Scr4 err75 = FitCodes<min,max,prc,otyp>(xyz, mask, codes75x, indices75x, codes75y, indices75y);
  Scr4 err77 = FitCodes<min,max,prc,otyp>(xyz, mask, codes77x, indices77x, codes77y, indices77y);
  Scr4 err   = Min(Min(err55, err77), Min(err57, err75));
  
  /*
#if defined(TRACK_STATISTICS)
  static int matches[5] = {0}; matches[0]++;
  /  / if (err == err77) matches[1]++;
  else if (err == err75) matches[2]++;
  else if (err == err57) matches[3]++;
  else                   matches[4]++;
#endif
   */

#if 0
  /* save the block with least error (prefer 77, capture min==max)
   *     | Brickwork | TerrainRk |
   * all |    211288 |    312500 |
   *  77 |    198531 |    286150 |
   *  75 |      5923 |      3755 |
   *  57 |      6791 |     20994 |
   *  55 |        42 |      1601 |
   **/ if (err == err77) {
    WriteNormalBlock7(max77.R(), min77.R(), indices77x, blockx);
    WriteNormalBlock7(max77.G(), min77.G(), indices77y, blocky);
  }
  else if (err == err75) {
    WriteNormalBlock7(max75.R(), min75.R(), indices75x, blockx);
    WriteNormalBlock5(min75.G(), max75.G(), indices75y, blocky);
  }
  else if (err == err57) {
    WriteNormalBlock5(min57.R(), max57.R(), indices57x, blockx);
    WriteNormalBlock7(max57.G(), min57.G(), indices57y, blocky);
  }
  else {
    WriteNormalBlock5(min55.R(), max55.R(), indices55x, blockx);
    WriteNormalBlock5(min55.G(), max55.G(), indices55y, blocky);
  }
#else
  Col4 min__x, max__x; u8* indices__x;
  Col4 min__y, max__y; u8* indices__y;

  /**/ if (err == err77) {
    min__x = max77, max__x = min77; indices__x = indices77x;
    min__y = max77, max__y = min77; indices__y = indices77y;
  }
  else if (err == err75) {
    min__x = max75, max__x = min75; indices__x = indices75x;
    min__y = min75, max__y = max75; indices__y = indices75y;
  }
  else if (err == err57) {
    min__x = min57, max__x = max57; indices__x = indices57x;
    min__y = max57, max__y = min57; indices__y = indices57y;
  }
  else {
    min__x = min55, max__x = max55; indices__x = indices55x;
    min__y = min55, max__y = max55; indices__y = indices55y;
  }

  WriteNormalBlock(min__x.R(), max__x.R(), indices__x, blockx);
  WriteNormalBlock(min__y.G(), max__y.G(), indices__y, blocky);
#endif
}

#undef	FIT_THRESHOLD

template<const int min, const int max, const int prc, const int lwr, const int upr, typename otyp, typename dtyp>
static void CompressNormalsBtc5i(dtyp const* xyzd, int mask, void* blockx, void* blocky, int flags)
{
  Vec4 xyz[16];

  for (int i = 0; i < 16; ++i) {
    // create floating point vector
    Col4 value; LoadUnaligned(value, &xyzd[i * 4]);

    // create floating point vector
    xyz[i] = Normalize(KillW(value + Vec4(IOFFSET)) * Vec4(ISCALE));
  }
    
  CompressNormalBtc5v<min,max,prc,otyp,dtyp>(xyz, mask, blockx, blocky, flags);
}

template<const int min, const int max, const int prc, typename otyp, typename dtyp>
static void CompressNormalsBtc5f(dtyp const* xyzd, int mask, void* blockx, void* blocky, int flags)
{
  Vec4 xyz[16];

  for (int i = 0; i < 16; ++i) {
    // create floating point vector
    Vec4 value; LoadUnaligned(value, &xyzd[4 * i]);
    
    // create floating point vector
    xyz[i] = Normalize(KillW(value));
  }
    
  CompressNormalBtc5v<min,max,prc,otyp,dtyp>(xyz, mask, blockx, blocky, flags);
}

void CompressNormalsBtc5u(u8 const* xyzd, int mask, void* blockx, void* blocky, int flags) {
  CompressNormalsBtc5i<   0,255, CBLB,    0,255, unsigned>(xyzd, mask, blockx, blocky, flags); }
void CompressNormalsBtc5s(s8 const* xyzd, int mask, void* blockx, void* blocky, int flags) {
  CompressNormalsBtc5i<-127,127, CBLB, -127,127,   signed>(xyzd, mask, blockx, blocky, flags); }

void CompressNormalsBtc5u(u16 const* xyzd, int mask, void* blockx, void* blocky, int flags) {
  CompressNormalsBtc5i<   0,255, CBHB,      0,65535, unsigned>(xyzd, mask, blockx, blocky, flags); }
void CompressNormalsBtc5s(s16 const* xyzd, int mask, void* blockx, void* blocky, int flags) {
  CompressNormalsBtc5i<-127,127, CBHB, -32767,32767,   signed>(xyzd, mask, blockx, blocky, flags); }

void CompressNormalsBtc5u(f23 const* xyzd, int mask, void* blockx, void* blocky, int flags) {
  CompressNormalsBtc5f<   0,255, CBHB, unsigned>(xyzd, mask, blockx, blocky, flags); }
void CompressNormalsBtc5s(f23 const* xyzd, int mask, void* blockx, void* blocky, int flags) {
  CompressNormalsBtc5f<-127,127, CBHB,   signed>(xyzd, mask, blockx, blocky, flags); }

/* *****************************************************************************
 */
template<const int prc, typename ctyp, typename etyp>
static void ReadNormalBlock(
  ctyp (&codesx  )[ 8], ctyp (&codesy  )[ 8],
  u8   (&indicesx)[16], u8   (&indicesy)[16],
  void const* blockx, void const* blocky)
{
  // get the two coord values
  etyp const* bytesx = reinterpret_cast< etyp const* >(blockx);
  etyp const* bytesy = reinterpret_cast< etyp const* >(blocky);
  etyp coord0x = bytesx[0];
  etyp coord1x = bytesx[1];
  etyp coord0y = bytesy[0];
  etyp coord1y = bytesy[1];

  // compare the values to build the codebook
  codesx[0] = (ctyp)coord0x << prc;
  codesx[1] = (ctyp)coord1x << prc;
  codesy[0] = (ctyp)coord0y << prc;
  codesy[1] = (ctyp)coord1y << prc;

  // use 5-coord or 7-coord  codebook
  Codebook6or8<prc>(codesx, coord0x <= coord1x);
  Codebook6or8<prc>(codesy, coord0y <= coord1y);

  // decode the indices
  u8 const* srcx = (u8*)bytesx + 2;
  u8 const* srcy = (u8*)bytesy + 2;
  u8* destx = indicesx;
  u8* desty = indicesy;
  for (int i = 0; i < 2; ++i) {
    // grab 3 bytes
    int valuex = 0;
    int valuey = 0;
    for (int j = 0; j < 3; ++j) {
      int bytex = *srcx++;
      int bytey = *srcy++;
      valuex += (bytex << 8 * j);
      valuey += (bytey << 8 * j);
    }

    // unpack 8 3-bit values from it
    for (int j = 0; j < 8; ++j) {
      int indexx = (valuex >> 3 * j) & 0x7;
      int indexy = (valuey >> 3 * j) & 0x7;
      *destx++ = (u8)indexx;
      *desty++ = (u8)indexy;
    }
  }
}
  
template<const int min, const int max, const int prc, typename dtyp, typename ctyp, typename etyp, const int lwr, const int upr>
static void DecompressNormalsBtc5i(dtyp* xyzd, void const* blockx, void const* blocky)
{
  ctyp codesx[8]; u8 indicesx[16];
  ctyp codesy[8]; u8 indicesy[16];

  ReadNormalBlock<prc,ctyp,etyp>(codesx, codesy, indicesx, indicesy, blockx, blocky);

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i) {
    Col3 _xyz0  = Col3(codesx[indicesx[i]], codesy[indicesy[i]]);
    Vec3 cxyz0  = (_xyz0 + Vec3(COFFSET)) * Vec3(CSCALE);
         cxyz0  = Complement<DISARM>(cxyz0);
	 cxyz0  = (cxyz0 * Vec3(OSCALE)) - Vec3(OOFFSET);
	 _xyz0  = FloatToInt<true>(cxyz0);

    StoreUnaligned(_xyz0, &xyzd[4 * i]);
  }
}

template<const int min, const int max, const int prc, typename dtyp, typename ctyp, typename etyp>
static void DecompressNormalsBtc5f(dtyp* xyzd, void const* blockx, void const* blocky)
{
  ctyp codesx[8]; u8 indicesx[16];
  ctyp codesy[8]; u8 indicesy[16];

  ReadNormalBlock<prc,ctyp,etyp>(codesx, codesy, indicesx, indicesy, blockx, blocky);

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i) {
    Col3 _xyz0  = Col3(codesx[indicesx[i]], codesy[indicesy[i]]);
    Vec3 cxyz0  = (_xyz0 + Vec3(COFFSET)) * Vec3(CSCALE);
         cxyz0  = Complement<DISARM>(cxyz0);
//	 cxyz0  = (cxyz0 * Vec4(OSCALE)) - Vec4(OOFFSET);
//	 _xyz0  = FloatToFloat<false>(cxyz0);

    StoreUnaligned(cxyz0, &xyzd[4 * i]);
  }
}

void DecompressNormalsBtc5u(u8* xyzd, void const* blockx, void const* blocky) {
  DecompressNormalsBtc5i<   0,255, CBLB, u8 , u8 ,u8,      0,  255>(xyzd, blockx, blocky); }
void DecompressNormalsBtc5s(s8* xyzd, void const* blockx, void const* blocky) {
  DecompressNormalsBtc5i<-127,127, CBLB, s8 , s8 ,s8,   -127,  127>(xyzd, blockx, blocky); }

void DecompressNormalsBtc5u(u16* xyzd, void const* blockx, void const* blocky) {
  DecompressNormalsBtc5i<   0,255, CBHB, u16, u16,u8,      0,65535>(xyzd, blockx, blocky); }
void DecompressNormalsBtc5s(s16* xyzd, void const* blockx, void const* blocky) {
  DecompressNormalsBtc5i<-127,127, CBHB, s16, s16,s8, -32767,32767>(xyzd, blockx, blocky); }

void DecompressNormalsBtc5u(f23* xyzd, void const* blockx, void const* blocky) {
  DecompressNormalsBtc5f<   0,255, CBHB, f23, u16,u8>(xyzd, blockx, blocky); }
void DecompressNormalsBtc5s(f23* xyzd, void const* blockx, void const* blocky) {
  DecompressNormalsBtc5f<-127,127, CBHB, f23, s16,s8>(xyzd, blockx, blocky); }

#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#undef CBLB
#undef CBHB

#undef FSCALE
#undef FOFFSET
#undef CSCALE
#undef COFFSET

#undef ISCALE
#undef IOFFSET
#undef OSCALE
#undef OOFFSET
