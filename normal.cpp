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

#if 1
#include "normal.h"
#include "maths.h"
#include "simd.h"

#if	!defined(SQUISH_USE_COMPUTE)
#include <cmath>
#include <algorithm>
#endif

#include "inlineables.cpp"

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

static Scr4 FitCodes(Vec4 const* xyz, int mask,
		     u8 const* codesx, u8* indicesx,
		     u8 const* codesy, u8* indicesy)
{
  Vec4 scale  = Vec4( 1.0f / 127.5f);
  Vec4 offset = Vec4(-1.0f * 127.5f);

  // fit each coord value to the codebook
  Scr4 err = Scr4(16.0f);

  for (int i = 0; i < 16; ++i) {
    // check this pixel is valid
    int bit = 1 << i;
    if ((mask & bit) == 0) {
      // use the first code
      indicesx[i] = 0;
      indicesy[i] = 0;
      continue;
    }

    // fetch floating point vector
    Vec4 value = xyz[i];

    // find the closest code
    Scr4 dist = Scr4(-1.0f);

    int idxx = 0;
    int idxy = 0;

    for (int k = 0; k < 8; k += 1)
    for (int j = 0; j < 8; j += 0) {
      // complement z
      Col4 _xyz0 = Col4(codesx[j + 0], codesy[k]); Vec4 cxyz0 = scale * (offset + _xyz0); cxyz0 = Complement<true>(cxyz0);
      Col4 _xyz1 = Col4(codesx[j + 1], codesy[k]); Vec4 cxyz1 = scale * (offset + _xyz1); cxyz1 = Complement<true>(cxyz1);
      Col4 _xyz2 = Col4(codesx[j + 2], codesy[k]); Vec4 cxyz2 = scale * (offset + _xyz2); cxyz2 = Complement<true>(cxyz2);
      Col4 _xyz3 = Col4(codesx[j + 3], codesy[k]); Vec4 cxyz3 = scale * (offset + _xyz3); cxyz3 = Complement<true>(cxyz3);

      // measure absolute angle-deviation (cosine)
      Scr4 d0 = Dot(value, cxyz0);
      Scr4 d1 = Dot(value, cxyz1);
      Scr4 d2 = Dot(value, cxyz2);
      Scr4 d3 = Dot(value, cxyz3);

      // encourage OoO
      Scr4 da = Max(d0, d1);
      Scr4 db = Max(d2, d3);
      dist = Max(da, dist);
      dist = Max(db, dist);

      // will cause VS to make them all cmovs
      if (d0 == dist) { idxx = j; idxy = k; } j++;
      if (d1 == dist) { idxx = j; idxy = k; } j++;
      if (d2 == dist) { idxx = j; idxy = k; } j++;
      if (d3 == dist) { idxx = j; idxy = k; } j++;
    }

    // save the index
    indicesx[i] = (u8)idxx;
    indicesy[i] = (u8)idxy;

    // accumulate the error (sine)
    err -= dist;
  }

  // return the total error
  return err;
}

template<const int stepx, const int stepy>
static Scr4 FitError(Vec4 const* xyz, Col4 &minXY, Col4 &maxXY, Scr4 &errXY) {
  Vec4 scale  = Vec4( 1.0f / 127.5f);
  Vec4 offset = Vec4(-1.0f * 127.5f);

#if	(SQUISH_USE_SIMD > 0) && 0
  const Vec4 minxy = minXY;
  const Vec4 maxxy = maxXY;

  float rx = (float)((maxXY.R() - minXY.R()) >> 1);	// max 127
  float ry = (float)((maxXY.G() - minXY.G()) >> 1);	// max 127

  Scr4 errC   = Scr4(errXY);
  Vec4 sx     = minxy.SplatX();
  Vec4 sy     = minxy.SplatY();
  Vec4 ex     = maxxy.SplatX();
  Vec4 ey     = maxxy.SplatY();
  Vec4 rangex = Vec4(rx, -rx, 0.0f, 0.0f);
  Vec4 rangey = Vec4(ry, -ry, 0.0f, 0.0f);

  while ((rangex + rangey) > Vec4(0.0f)) {
    // s + os, s + os, s + os - r, s + os + r
    // e - oe - r, e - oe + r, e - oe, e - oe
    Vec4 rssmpx = sx + rangex;
    Vec4 rmpeex = ex + rangex.Swap();

    Vec4 cssmpx = Max(Min(rssmpx, Vec4(255.0f)), Vec4(0.0f));
    Vec4 cmpeex = Max(Min(rmpeex, Vec4(255.0f)), Vec4(0.0f));

    Vec4 rssmpy = sy + rangey;
    Vec4 rmpeey = ey + rangey.Swap();

    Vec4 cssmpy = Max(Min(rssmpy, Vec4(255.0f)), Vec4(0.0f));
    Vec4 cmpeey = Max(Min(rmpeey, Vec4(255.0f)), Vec4(0.0f));

    Vec4 cssss_ = Vec4(sx, sy);
    Vec4 ceeee_ = Vec4(ex, ey);

    Vec4 _cb5[8];
    Vec4 _cb7[8];
    Vec4 _cbx[8];
    Vec4 _cby[8];

    // for all possible codebook-entries of xy
    if ((stepx == 5) || (stepy == 5)) {
      _cb5[0] =  cssss_                                                    ; //b5[0] = scale * (offset + Truncate(_cb5[0]));
      _cb5[1] = (cssss_ * Vec4(4.0f / 5.0f)) + (ceeee_ * Vec4(1.0f / 5.0f)); _cb5[1] = scale * (offset + Truncate(_cb5[1]));
      _cb5[2] = (cssss_ * Vec4(3.0f / 5.0f)) + (ceeee_ * Vec4(2.0f / 5.0f)); _cb5[2] = scale * (offset + Truncate(_cb5[2]));
      _cb5[3] = (cssss_ * Vec4(2.0f / 5.0f)) + (ceeee_ * Vec4(3.0f / 5.0f)); _cb5[3] = scale * (offset + Truncate(_cb5[3]));
      _cb5[4] = (cssss_ * Vec4(1.0f / 5.0f)) + (ceeee_ * Vec4(4.0f / 5.0f)); _cb5[4] = scale * (offset + Truncate(_cb5[4]));
      _cb5[5] =                                 ceeee_                     ; //b5[5] = scale * (offset + Truncate(_cb5[5]));
      _cb5[6] = Vec4(-1.0f);
      _cb5[7] = Vec4( 1.0f);
    }

    // for all possible codebook-entries of xy
    if ((stepx == 7) || (stepy == 7)) {
      _cb7[0] =  cssss_                                                    ; //b7[0] = scale * (offset + Truncate(_cb7[0]));
      _cb7[1] = (cssss_ * Vec4(6.0f / 5.0f)) + (ceeee_ * Vec4(1.0f / 5.0f)); _cb7[1] = scale * (offset + Truncate(_cb7[1]));
      _cb7[2] = (cssss_ * Vec4(5.0f / 5.0f)) + (ceeee_ * Vec4(2.0f / 5.0f)); _cb7[2] = scale * (offset + Truncate(_cb7[2]));
      _cb7[3] = (cssss_ * Vec4(4.0f / 5.0f)) + (ceeee_ * Vec4(3.0f / 5.0f)); _cb7[3] = scale * (offset + Truncate(_cb7[3]));
      _cb7[4] = (cssss_ * Vec4(3.0f / 5.0f)) + (ceeee_ * Vec4(4.0f / 5.0f)); _cb7[4] = scale * (offset + Truncate(_cb7[4]));
      _cb7[5] = (cssss_ * Vec4(2.0f / 5.0f)) + (ceeee_ * Vec4(5.0f / 5.0f)); _cb7[5] = scale * (offset + Truncate(_cb7[5]));
      _cb7[6] = (cssss_ * Vec4(1.0f / 5.0f)) + (ceeee_ * Vec4(6.0f / 5.0f)); _cb7[6] = scale * (offset + Truncate(_cb7[6]));
      _cb7[7] =                                 ceeee_                     ; //b7[7] = scale * (offset + Truncate(_cb7[7]));
    }

    // for all possible codebook-entries of x
    if (stepx == 5) {
      _cbx[0] =  cssmpx                                                    ; //bx[0] = scale * (offset + Truncate(_cbx[0]));
      _cbx[1] = (cssmpx * Vec4(4.0f / 5.0f)) + (cmpeex * Vec4(1.0f / 5.0f)); _cbx[1] = scale * (offset + Truncate(_cbx[1]));
      _cbx[2] = (cssmpx * Vec4(3.0f / 5.0f)) + (cmpeex * Vec4(2.0f / 5.0f)); _cbx[2] = scale * (offset + Truncate(_cbx[2]));
      _cbx[3] = (cssmpx * Vec4(2.0f / 5.0f)) + (cmpeex * Vec4(3.0f / 5.0f)); _cbx[3] = scale * (offset + Truncate(_cbx[3]));
      _cbx[4] = (cssmpx * Vec4(1.0f / 5.0f)) + (cmpeex * Vec4(4.0f / 5.0f)); _cbx[4] = scale * (offset + Truncate(_cbx[4]));
      _cbx[5] =                                 cmpeex                     ; //bx[5] = scale * (offset + Truncate(_cbx[5]));
      _cbx[6] = Vec4(-1.0f);
      _cbx[7] = Vec4( 1.0f);
    }
    else if (stepx == 7) {
      _cbx[0] =  cssmpx                                                    ; //bx[0] = scale * (offset + Truncate(_cbx[0]));
      _cbx[1] = (cssmpx * Vec4(6.0f / 5.0f)) + (cmpeex * Vec4(1.0f / 5.0f)); _cbx[1] = scale * (offset + Truncate(_cbx[1]));
      _cbx[2] = (cssmpx * Vec4(5.0f / 5.0f)) + (cmpeex * Vec4(2.0f / 5.0f)); _cbx[2] = scale * (offset + Truncate(_cbx[2]));
      _cbx[3] = (cssmpx * Vec4(4.0f / 5.0f)) + (cmpeex * Vec4(3.0f / 5.0f)); _cbx[3] = scale * (offset + Truncate(_cbx[3]));
      _cbx[4] = (cssmpx * Vec4(3.0f / 5.0f)) + (cmpeex * Vec4(4.0f / 5.0f)); _cbx[4] = scale * (offset + Truncate(_cbx[4]));
      _cbx[5] = (cssmpx * Vec4(2.0f / 5.0f)) + (cmpeex * Vec4(5.0f / 5.0f)); _cbx[5] = scale * (offset + Truncate(_cbx[5]));
      _cbx[6] = (cssmpx * Vec4(1.0f / 5.0f)) + (cmpeex * Vec4(6.0f / 5.0f)); _cbx[6] = scale * (offset + Truncate(_cbx[6]));
      _cbx[7] =                                 cmpeex                     ; //bx[7] = scale * (offset + Truncate(_cbx[7]));
    }

    // for all possible codebook-entries of y
    if (stepy == 5) {
      _cby[0] =  cssmpy                                                    ; //by[0] = scale * (offset + Truncate(_cby[0]));
      _cby[1] = (cssmpy * Vec4(4.0f / 5.0f)) + (cmpeey * Vec4(1.0f / 5.0f)); _cby[1] = scale * (offset + Truncate(_cby[1]));
      _cby[2] = (cssmpy * Vec4(3.0f / 5.0f)) + (cmpeey * Vec4(2.0f / 5.0f)); _cby[2] = scale * (offset + Truncate(_cby[2]));
      _cby[3] = (cssmpy * Vec4(2.0f / 5.0f)) + (cmpeey * Vec4(3.0f / 5.0f)); _cby[3] = scale * (offset + Truncate(_cby[3]));
      _cby[4] = (cssmpy * Vec4(1.0f / 5.0f)) + (cmpeey * Vec4(4.0f / 5.0f)); _cby[4] = scale * (offset + Truncate(_cby[4]));
      _cby[5] =                                 cmpeey                     ; //by[5] = scale * (offset + Truncate(_cby[5]));
      _cby[6] = Vec4(-1.0f);
      _cby[7] = Vec4( 1.0f);
    }
    else if (stepy == 7) {
      _cby[0] =  cssmpy                                                    ; //by[0] = scale * (offset + Truncate(_cby[0]));
      _cby[1] = (cssmpy * Vec4(6.0f / 5.0f)) + (cmpeey * Vec4(1.0f / 5.0f)); _cby[1] = scale * (offset + Truncate(_cby[1]));
      _cby[2] = (cssmpy * Vec4(5.0f / 5.0f)) + (cmpeey * Vec4(2.0f / 5.0f)); _cby[2] = scale * (offset + Truncate(_cby[2]));
      _cby[3] = (cssmpy * Vec4(4.0f / 5.0f)) + (cmpeey * Vec4(3.0f / 5.0f)); _cby[3] = scale * (offset + Truncate(_cby[3]));
      _cby[4] = (cssmpy * Vec4(3.0f / 5.0f)) + (cmpeey * Vec4(4.0f / 5.0f)); _cby[4] = scale * (offset + Truncate(_cby[4]));
      _cby[5] = (cssmpy * Vec4(2.0f / 5.0f)) + (cmpeey * Vec4(5.0f / 5.0f)); _cby[5] = scale * (offset + Truncate(_cby[5]));
      _cby[6] = (cssmpy * Vec4(1.0f / 5.0f)) + (cmpeey * Vec4(6.0f / 5.0f)); _cby[6] = scale * (offset + Truncate(_cby[6]));
      _cby[7] =                                 cmpeey                     ; //by[7] = scale * (offset + Truncate(_cby[7]));
    }

    Vec4 error0x  = Vec4(16.0f);
    Vec4 error0y  = Vec4(16.0f);
    Vec4 error1xy = Vec4(16.0f);
    Vec4 error2xy = Vec4(16.0f);
    Vec4 error3xy = Vec4(16.0f);
    Vec4 error4xy = Vec4(16.0f);
    for (int v = 0; v < 16; v++) {
      // find the closest code
      Vec4 valx = xyz[v].SplatX();
      Vec4 valy = xyz[v].SplatY();
      Vec4 valz = xyz[v].SplatZ();

      Vec4 dist0x  = Vec4(-1.0f);
      Vec4 dist0y  = Vec4(-1.0f);
      Vec4 dist1xy = Vec4(-1.0f);
      Vec4 dist2xy = Vec4(-1.0f);
      Vec4 dist3xy = Vec4(-1.0f);
      Vec4 dist4xy = Vec4(-1.0f);

      for (int i = 0; i < 8; i++) {
      for (int j = 0; j < 8; j++) {
	Vec4 (&_cxx)[8] = (stepx == 5 ? _cb5 : _cb7);
	Vec4 (&_cyy)[8] = (stepy == 5 ? _cb5 : _cb7);

	//           ey -= ry	   =>   error0x1y -> _cb_[?].x _cby[?].w
	//           ey += ry	   =>   error0x2y -> _cb_[?].x _cby[?].z
	//           sy -= ry	   =>   error0x3y -> _cb_[?].x _cby[?].y
	//           sy += ry	   =>   error0x4y -> _cb_[?].x _cby[?].x
	Vec4 x0x  = _cxx[i].SplatX();	Vec4 y0x  = _cby[j];		Vec4 z0x  = Complement(x0x, y0x);

	// ex -= rx          	   =>   error1x0y -> _cbx[?].w _cb_[?].y
	// ex += rx          	   =>   error2x0y -> _cbx[?].z _cb_[?].y
	// sx -= rx          	   =>   error3x0y -> _cbx[?].y _cb_[?].y
	// sx += rx          	   =>   error4x0y -> _cbx[?].x _cb_[?].y
	Vec4 x0y  = _cbx[i];		Vec4 y0y  = _cyy[j].SplatY();	Vec4 z0y  = Complement(x0y, y0y);

	// ex -= rx, ey -= ry	   =>   error1x1y -> _cbx[?].w _cby[?].w
	// ex -= rx, ey += ry	   =>   error1x2y -> _cbx[?].w _cby[?].z
	// ex -= rx, sy -= ry	   =>   error1x3y -> _cbx[?].w _cby[?].y
	// ex -= rx, sy += ry	   =>   error1x4y -> _cbx[?].w _cby[?].x
	Vec4 x1xy = _cbx[i].SplatW();	Vec4 y1xy = _cby[j];		Vec4 z1xy = Complement(x1xy, y1xy);

	// ex += rx, ey -= ry	   =>   error2x1y -> _cbx[?].z _cby[?].w
	// ex += rx, ey += ry	   =>   error2x2y -> _cbx[?].z _cby[?].z
	// ex += rx, sy -= ry	   =>   error2x3y -> _cbx[?].z _cby[?].y
	// ex += rx, sy += ry	   =>   error2x4y -> _cbx[?].z _cby[?].x
	Vec4 x2xy = _cbx[i].SplatZ();	Vec4 y2xy = _cby[j];		Vec4 z2xy = Complement(x2xy, y2xy);

	// sx -= rx, ey -= ry	   =>   error3x1y -> _cbx[?].y _cby[?].w
	// sx -= rx, ey += ry	   =>   error3x2y -> _cbx[?].y _cby[?].z
	// sx -= rx, sy -= ry	   =>   error3x3y -> _cbx[?].y _cby[?].y
	// sx -= rx, sy += ry	   =>   error3x4y -> _cbx[?].y _cby[?].x
	Vec4 x3xy = _cbx[i].SplatY();	Vec4 y3xy = _cby[j];		Vec4 z3xy = Complement(x3xy, y3xy);

	// sx += rx, ey -= ry	   =>   error4x1y -> _cbx[?].x _cby[?].w
	// sx += rx, ey += ry	   =>   error4x2y -> _cbx[?].x _cby[?].z
	// sx += rx, sy -= ry	   =>   error4x3y -> _cbx[?].x _cby[?].y
	// sx += rx, sy += ry	   =>   error4x4y -> _cbx[?].x _cby[?].x
	Vec4 x4xy = _cbx[i].SplatX();	Vec4 y4xy = _cby[j];		Vec4 z4xy = Complement(x4xy, y4xy);

        // measure absolute angle-deviation (cosine)
	Vec4 d0x  = (valx * x0x ) + (valy * y0x ) + (valz * z0x );
	Vec4 d0y  = (valx * x0y ) + (valy * y0y ) + (valz * z0y );
	Vec4 d1xy = (valx * x1xy) + (valy * y1xy) + (valz * z1xy);
	Vec4 d2xy = (valx * x2xy) + (valy * y2xy) + (valz * z2xy);
	Vec4 d3xy = (valx * x3xy) + (valy * y3xy) + (valz * z3xy);
	Vec4 d4xy = (valx * x4xy) + (valy * y4xy) + (valz * z4xy);

	// select the smallest deviation
	dist0x  = Max(dist0x , d0x );
	dist0y  = Max(dist0y , d0y );
	dist1xy = Max(dist1xy, d1xy);
	dist2xy = Max(dist2xy, d2xy);
	dist3xy = Max(dist3xy, d3xy);
	dist4xy = Max(dist4xy, d4xy);
      }
      }

      // accumulate the error (sine)
      error0x  -= dist0x ;
      error0y  -= dist0y ;
      error1xy -= dist1xy;
      error2xy -= dist2xy;
      error3xy -= dist3xy;
      error4xy -= dist4xy;
    }

    // encourage OoO
    Vec4 error0  = Min(error0x , error0y );
    Vec4 error14 = Min(error1xy, error4xy);
    Vec4 error23 = Min(error2xy, error3xy);
    Vec4 errors  = Min(error0, Min(error14, error23));
    Scr4 merr    = HorizontalMin(errors);

    if (!(merr < errC)) {
      // half range
      if (rangex > rangey)
	rangex = Truncate(rangex * Vec4(0.5f));
      else
	rangey = Truncate(rangey * Vec4(0.5f));
    }
    else if (CompareEqualTo(error0 , merr) != 0x00) {
      int value0x = CompareEqualTo(error0x, merr);
      int value0y = CompareEqualTo(error0y, merr);

      /**/ if (value0x & 0x8) ey = rmpeey.SplatW();	//           ey -= ry;// range up
      else if (value0x & 0x1) sy = rssmpy.SplatX();	//           sy += ry;// range dn
      else if (value0y & 0x8) ex = rmpeex.SplatW();	// ex -= rx          ;// range up
      else if (value0y & 0x1) sx = rssmpx.SplatX();	// sx += rx          ;// range dn
      else if (value0x & 0x4) ey = rmpeey.SplatZ();	//           ey += ry;// range up
      else if (value0x & 0x2) sy = rssmpy.SplatY();	//           sy -= ry;// range dn
      else if (value0y & 0x4) ex = rmpeex.SplatZ();	// ex += rx          ;// range up
      else if (value0y & 0x2) sx = rssmpx.SplatY();	// sx -= rx          ;// range dn

      errC = merr;
    }
    else if (CompareEqualTo(error14, merr) != 0x00) {
      int value1xy = CompareEqualTo(error1xy, merr);
      int value14  = CompareEqualTo(error14 , merr);

      /**/ if (value1xy != 0) ex = rmpeex.SplatW();	// ex -= rx, ........;// range up
      else /* (value4xy != */ sx = rssmpx.SplatX();	// sx += rx, ........;// range dn

      /**/ if (value14 & 0x8) ey = rmpeey.SplatW();	// ........, ey -= ry;// range +-
      else if (value14 & 0x1) sy = rssmpy.SplatX();	// ........, sy += ry;// range +-
      else if (value14 & 0x4) ey = rmpeey.SplatZ();	// ........, ey += ry;// range +-
      else /* (value14 & 0x*/ sy = rssmpy.SplatY();	// ........, sy -= ry;// range +-

      errC = merr;
    }
    else /*if (CompareEqualTo(error23, merr) != 0x00)*/ {
      int value2xy = CompareEqualTo(error2xy, merr);
      int value23  = CompareEqualTo(error23 , merr);

      /**/ if (value2xy != 0) ex = rmpeex.SplatY();	// ex -= rx, ........;// range up
      else /* (value3xy != */ sx = rssmpx.SplatZ();	// sx += rx, ........;// range dn

      /**/ if (value23 & 0x8) ey = rmpeey.SplatW();	// ........, ey -= ry;// range +-
      else if (value23 & 0x1) sy = rssmpy.SplatX();	// ........, sy += ry;// range +-
      else if (value23 & 0x4) ey = rmpeey.SplatZ();	// ........, ey += ry;// range +-
      else /* (value23 & 0x*/ sy = rssmpy.SplatY();	// ........, sy -= ry;// range +-

      errC = merr;
    }

    // lossless
    if (!(errC > Vec4(0.0f)))
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
  minXY = FloatToInt<false>(Max(Min(Vec4(sx, sy), Vec4(255.0f)), Vec4(0.0f)));
  maxXY = FloatToInt<false>(Max(Min(Vec4(ex, ey), Vec4(255.0f)), Vec4(0.0f)));
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
    int msx = std::max(std::min(sx - rx, 0xFF), 0x00);  // os - r
    int csx = std::max(std::min(sx     , 0xFF), 0x00);  // os
    int psx = std::max(std::min(sx + rx, 0xFF), 0x00);  // os + r
    int mex = std::max(std::min(ex - rx, 0xFF), 0x00);  // oe + r
    int cex = std::max(std::min(ex     , 0xFF), 0x00);  // oe
    int pex = std::max(std::min(ex + rx, 0xFF), 0x00);  // oe - r

    int msy = std::max(std::min(sy - ry, 0xFF), 0x00);  // os - r
    int csy = std::max(std::min(sy     , 0xFF), 0x00);  // os
    int psy = std::max(std::min(sy + ry, 0xFF), 0x00);  // os + r
    int mey = std::max(std::min(ey - ry, 0xFF), 0x00);  // oe + r
    int cey = std::max(std::min(ey     , 0xFF), 0x00);  // oe
    int pey = std::max(std::min(ey + ry, 0xFF), 0x00);  // oe - r

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
    error4x1y = error4x2y = error4x3y = error4x4y = Scr4(16.0f);

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
      dist4x1y = dist4x2y = dist4x3y = dist4x4y = Scr4(-1.0f);

      // fetch floating point vector
      Vec4 value = xyz[v];

      // for all possible codebook-entries
      for (int g = 0; g < 8; g++) {
	int cb0y, cb1y, cb2y, cb3y, cb4y;
	
      	// http://embeddedgurus.com/stack-overflow/2009/06/division-of-integers-by-constants/
	// Divide by 5:  (((uint32_t)A * (uint32_t)0xCCCD) >> 16) >> 2
	// Divide by 7: ((((uint32_t)A * (uint32_t)0x2493) >> 16) + A) >> 1) >> 2
	if (stepy == 5) {
	  cb0y = (((5 - g) * csy + g * cey) / 5);
	  cb1y = (((5 - g) * csy + g * mey) / 5);
	  cb2y = (((5 - g) * csy + g * pey) / 5);
	  cb3y = (((5 - g) * msy + g * cey) / 5);
	  cb4y = (((5 - g) * psy + g * cey) / 5);

	  if (g >= 6) cb0y = cb1y = cb2y = cb3y = cb4y = 0x00;
	  if (g >= 7) cb0y = cb1y = cb2y = cb3y = cb4y = 0xFF;
	}
	else if (stepy == 7) {
	  cb0y = (((7 - g) * csy + g * cey) / 7);
	  cb1y = (((7 - g) * csy + g * mey) / 7);
	  cb2y = (((7 - g) * csy + g * pey) / 7);
	  cb3y = (((7 - g) * msy + g * cey) / 7);
	  cb4y = (((7 - g) * psy + g * cey) / 7);
	}
	else
	  abort();

      for (int f = 0; f < 8; f++) {
	int cb0x, cb1x, cb2x, cb3x, cb4x;

	if (stepx == 5) {
	  cb0x = (((5 - f) * csx + f * cex) / 5);
	  cb1x = (((5 - f) * csx + f * mex) / 5);
	  cb2x = (((5 - f) * csx + f * pex) / 5);
	  cb3x = (((5 - f) * msx + f * cex) / 5);
	  cb4x = (((5 - f) * psx + f * cex) / 5);

	  if (f >= 6) cb0x = cb1x = cb2x = cb3x = cb4x = 0x00;
	  if (f >= 7) cb0x = cb1x = cb2x = cb3x = cb4x = 0xFF;
	}
	else if (stepx == 7) {
	  cb0x = (((7 - f) * csx + f * cex) / 7);
	  cb1x = (((7 - f) * csx + f * mex) / 7);
	  cb2x = (((7 - f) * csx + f * pex) / 7);
	  cb3x = (((7 - f) * msx + f * cex) / 7);
	  cb4x = (((7 - f) * psx + f * cex) / 7);
	}
	else
	  abort();

	// complement z
	Col4 _xyz0x1y = Col4(cb0x, cb1y); Vec4 cxyz0x1y = scale * (offset + _xyz0x1y); cxyz0x1y = Complement<true>(cxyz0x1y);
	Col4 _xyz0x2y = Col4(cb0x, cb2y); Vec4 cxyz0x2y = scale * (offset + _xyz0x2y); cxyz0x2y = Complement<true>(cxyz0x2y);
	Col4 _xyz0x3y = Col4(cb0x, cb3y); Vec4 cxyz0x3y = scale * (offset + _xyz0x3y); cxyz0x3y = Complement<true>(cxyz0x3y);
	Col4 _xyz0x4y = Col4(cb0x, cb4y); Vec4 cxyz0x4y = scale * (offset + _xyz0x4y); cxyz0x4y = Complement<true>(cxyz0x4y);
	Col4 _xyz1x1y = Col4(cb1x, cb1y); Vec4 cxyz1x1y = scale * (offset + _xyz1x1y); cxyz1x1y = Complement<true>(cxyz1x1y);
	Col4 _xyz1x2y = Col4(cb1x, cb2y); Vec4 cxyz1x2y = scale * (offset + _xyz1x2y); cxyz1x2y = Complement<true>(cxyz1x2y);
	Col4 _xyz1x3y = Col4(cb1x, cb3y); Vec4 cxyz1x3y = scale * (offset + _xyz1x3y); cxyz1x3y = Complement<true>(cxyz1x3y);
	Col4 _xyz1x4y = Col4(cb1x, cb4y); Vec4 cxyz1x4y = scale * (offset + _xyz1x4y); cxyz1x4y = Complement<true>(cxyz1x4y);
	Col4 _xyz2x1y = Col4(cb2x, cb1y); Vec4 cxyz2x1y = scale * (offset + _xyz2x1y); cxyz2x1y = Complement<true>(cxyz2x1y);
	Col4 _xyz2x2y = Col4(cb2x, cb2y); Vec4 cxyz2x2y = scale * (offset + _xyz2x2y); cxyz2x2y = Complement<true>(cxyz2x2y);
	Col4 _xyz2x3y = Col4(cb2x, cb3y); Vec4 cxyz2x3y = scale * (offset + _xyz2x3y); cxyz2x3y = Complement<true>(cxyz2x3y);
	Col4 _xyz2x4y = Col4(cb2x, cb4y); Vec4 cxyz2x4y = scale * (offset + _xyz2x4y); cxyz2x4y = Complement<true>(cxyz2x4y);
	Col4 _xyz3x1y = Col4(cb3x, cb1y); Vec4 cxyz3x1y = scale * (offset + _xyz3x1y); cxyz3x1y = Complement<true>(cxyz3x1y);
	Col4 _xyz3x2y = Col4(cb3x, cb2y); Vec4 cxyz3x2y = scale * (offset + _xyz3x2y); cxyz3x2y = Complement<true>(cxyz3x2y);
	Col4 _xyz3x3y = Col4(cb3x, cb3y); Vec4 cxyz3x3y = scale * (offset + _xyz3x3y); cxyz3x3y = Complement<true>(cxyz3x3y);
	Col4 _xyz3x4y = Col4(cb3x, cb4y); Vec4 cxyz3x4y = scale * (offset + _xyz3x4y); cxyz3x4y = Complement<true>(cxyz3x4y);
	Col4 _xyz4x1y = Col4(cb4x, cb1y); Vec4 cxyz4x1y = scale * (offset + _xyz4x1y); cxyz4x1y = Complement<true>(cxyz4x1y);
	Col4 _xyz4x2y = Col4(cb4x, cb2y); Vec4 cxyz4x2y = scale * (offset + _xyz4x2y); cxyz4x2y = Complement<true>(cxyz4x2y);
	Col4 _xyz4x3y = Col4(cb4x, cb3y); Vec4 cxyz4x3y = scale * (offset + _xyz4x3y); cxyz4x3y = Complement<true>(cxyz4x3y);
	Col4 _xyz4x4y = Col4(cb4x, cb4y); Vec4 cxyz4x4y = scale * (offset + _xyz4x4y); cxyz4x4y = Complement<true>(cxyz4x4y);
	Col4 _xyz1x0y = Col4(cb1x, cb0y); Vec4 cxyz1x0y = scale * (offset + _xyz1x0y); cxyz1x0y = Complement<true>(cxyz1x0y);
	Col4 _xyz2x0y = Col4(cb2x, cb0y); Vec4 cxyz2x0y = scale * (offset + _xyz2x0y); cxyz2x0y = Complement<true>(cxyz2x0y);
	Col4 _xyz3x0y = Col4(cb3x, cb0y); Vec4 cxyz3x0y = scale * (offset + _xyz3x0y); cxyz3x0y = Complement<true>(cxyz3x0y);
	Col4 _xyz4x0y = Col4(cb4x, cb0y); Vec4 cxyz4x0y = scale * (offset + _xyz4x0y); cxyz4x0y = Complement<true>(cxyz4x0y);

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
      error0x1y -= dist0x1y;
      error0x2y -= dist0x2y;
      error0x3y -= dist0x3y;
      error0x4y -= dist0x4y;
      error1x1y -= dist1x1y;
      error1x2y -= dist1x2y;
      error1x3y -= dist1x3y;
      error1x4y -= dist1x4y;
      error2x1y -= dist2x1y;
      error2x2y -= dist2x2y;
      error2x3y -= dist2x3y;
      error2x4y -= dist2x4y;
      error3x1y -= dist3x1y;
      error3x2y -= dist3x2y;
      error3x3y -= dist3x3y;
      error3x4y -= dist3x4y;
      error4x1y -= dist4x1y;
      error4x2y -= dist4x2y;
      error4x3y -= dist4x3y;
      error4x4y -= dist4x4y;
      error1x0y -= dist1x0y;
      error2x0y -= dist2x0y;
      error3x0y -= dist3x0y;
      error4x0y -= dist4x0y;
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

//	/**/ if (merrx == error0x0y) rx >>= 1, ry >>= 1;	// half range
    /**/ if (merrx == error0x0y) rx > ry ? rx >>= 1 : ry >>= 1;
    else if (merrx == error0x1y)           ey -= ry;	// range up
    else if (merrx == error0x4y)           sy += ry;	// range up
    else if (merrx == error1x0y) ex -= rx          ;	// range up
    else if (merrx == error4x0y) sx += rx          ;	// range up
    else if (merrx == error0x2y)           ey += ry;	// range up
    else if (merrx == error0x3y)           sy -= ry;	// range up
    else if (merrx == error2x0y) ex += rx          ;	// range up
    else if (merrx == error3x0y) sx -= rx          ;	// range up
    else if (merrx == error1x1y) ex -= rx, ey -= ry;	// range up
    else if (merrx == error1x4y) ex -= rx, sy += ry;	// range up
    else if (merrx == error1x2y) ex -= rx, ey += ry;	// range up
    else if (merrx == error1x3y) ex -= rx, sy -= ry;	// range up
    else if (merrx == error4x1y) sx += rx, ey -= ry;	// range dn
    else if (merrx == error4x4y) sx += rx, sy += ry;	// range dn
    else if (merrx == error4x2y) sx += rx, ey += ry;	// range dn
    else if (merrx == error4x3y) sx += rx, sy -= ry;	// range dn
    else if (merrx == error2x1y) ex += rx, ey -= ry;	// range up
    else if (merrx == error2x4y) ex += rx, sy += ry;	// range up
    else if (merrx == error2x2y) ex += rx, ey += ry;	// range up
    else if (merrx == error2x3y) ex += rx, sy -= ry;	// range up
    else if (merrx == error3x1y) sx -= rx, ey -= ry;	// range dn
    else if (merrx == error3x4y) sx -= rx, sy += ry;	// range dn
    else if (merrx == error3x2y) sx -= rx, ey += ry;	// range dn
    else if (merrx == error3x3y) sx -= rx, sy -= ry;	// range dn

    // lossless
    error0x0y = merrx;
    if (!(error0x0y > Scr4(0.0f)))
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
  minXY = Max(Col4(sx, sy), Col4(0x00));
  maxXY = Min(Col4(ex, ey), Col4(0xFF));
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
      value += (index << 3 * j);
    }

    // store in 3 bytes
    for (int j = 0; j < 3; ++j) {
      int byte = (value >> 8 * j) & 0xFF;
      *dest++ = (u8)byte;
    }

    // 77766655.54443332.22111000, FFFEEEDD.DCCCBBBA.AA999888
    // 22111000.54443332.77766655, AA999888.DCCCBBBA.FFFEEEDD
  }
}

static void WriteNormalBlock5(int coord0, int coord1, u8 const* indices, void* block)
{
  // check the relative values of the endpoints
  if (coord0 > coord1) {
    // swap the indices
    u8 swapped[16];
    for (int i = 0; i < 16; ++i) {
      u8 index = indices[i];
      if (index == 0)
	swapped[i] = 1;
      else if (index == 1)
	swapped[i] = 0;
      else if (index <= 5)
	swapped[i] = 7 - index;
      else
	swapped[i] = index;
    }

    // write the block
    WriteNormalBlock(coord1, coord0, swapped, block);
  }
  else {
    // write the block
    WriteNormalBlock(coord0, coord1, indices, block);
  }
}

static void WriteNormalBlock7(int coord0, int coord1, u8 const* indices, void* block)
{
  // check the relative values of the endpoints
  if (coord0 < coord1) {
    // swap the indices
    u8 swapped[16];
    for (int i = 0; i < 16; ++i) {
      u8 index = indices[i];
      if (index == 0)
	swapped[i] = 1;
      else if (index == 1)
	swapped[i] = 0;
      else
	swapped[i] = 9 - index;
    }

    // write the block
    WriteNormalBlock(coord1, coord0, swapped, block);
  }
  else {
    // write the block
    WriteNormalBlock(coord0, coord1, indices, block);
  }
}

void CompressNormalBtc5(u8 const* xyzd, int mask, void* blockx, void* blocky, int flags)
{
  Vec4 scale  = Vec4( 1.0f / 127.5f);
  Vec4 offset = Vec4(-1.0f * 127.5f);
  Vec4 xyz[16];

  // get the range for 5-coord and 7-coord interpolation
  Col4 min5 = Col4(0xFF);
  Col4 max5 = Col4(0x00);
  Col4 min7 = Col4(0xFF);
  Col4 max7 = Col4(0x00);

  for (int i = 0; i < 16; ++i) {
    // check this pixel is valid
    int bit = 1 << i;
    if ((mask & bit) == 0)
      continue;

    // create integer vector
    Col4 value = Col4(
      xyzd[4 * i + 0],
      xyzd[4 * i + 1],
      xyzd[4 * i + 2]
    );

    // create floating point vector
    xyz[i] = Normalize(scale * KillW(offset + value));

    Col4 mask =
      IsZero(value) |
      IsOne (value);

    // incorporate into the min/max
    min5 = Min(min5, (mask % value) + (mask & min5));
    max5 = Max(max5, (mask % value) + (mask & max5));
    min7 = Min(min7, value);
    max7 = Max(max7, value);
  }

  // handle the case that no valid range was found
  min5 = Min(min5, max5);
  max5 = Max(min5, max5);
  min7 = Min(min7, max7);
  max7 = Max(min7, max7);

  // do the iterative tangent search
  if (flags & kNormalIterativeFit) {
    // initial values
    Scr4 err5 = Scr4(16.0f);
    Scr4 err7 = Scr4(16.0f);
    
#if defined(TRACK_STATISTICS)
    for (int v = 0; v < 16; v++) {
      // find the closest code
      Scr4 dist5, dist7;
      dist5 = dist7 = Scr4(-1.0f);

      // fetch floating point vector
      Vec4 value = xyz[v];

      // for all possible codebook-entries
      for (int g = 0; g < 8; g++) {
	int cb5y = ((min5.G() * (5 - g) + max5.G() * g) / 5);
	int cb7y = ((min7.G() * (7 - g) + max7.G() * g) / 7);

	if (g >= 6) cb5y = 0x00;
	if (g >= 7) cb5y = 0xFF;

	for (int f = 0; f < 8; f++) {
      	  // http://embeddedgurus.com/stack-overflow/2009/06/division-of-integers-by-constants/
	  // Divide by 5:  (((uint32_t)A * (uint32_t)0xCCCD) >> 16) >> 2
	  // Divide by 7: ((((uint32_t)A * (uint32_t)0x2493) >> 16) + A) >> 1) >> 2
	  int cb5x = ((min5.R() * (5 - f) + max5.R() * f) / 5);
	  int cb7x = ((min7.R() * (7 - f) + max7.R() * f) / 7);

	  if (f >= 6) cb5x = 0x00;
	  if (f >= 7) cb5x = 0xFF;

	  // complement z
	  Col4 _xyz0 = Col4(cb5x, cb5y); Vec4 cxyz0 = scale * (offset + _xyz0); cxyz0 = Complement<true>(cxyz0);
	  Col4 _xyz1 = Col4(cb7x, cb7y); Vec4 cxyz1 = scale * (offset + _xyz1); cxyz1 = Complement<true>(cxyz1);

	  // measure absolute angle-deviation (cosine)
	  Scr4 d5 = Dot(value, cxyz0);
	  Scr4 d7 = Dot(value, cxyz1);

	  // select the smallest deviation
	  if (dist5 < d5) dist5 = d5;
	  if (dist7 < d7) dist7 = d7;
	}
      }

      // accumulate the error (sine)
      err5 -= dist5;
      err7 -= dist7;
    }

    // binary search, tangent-fitting
    Scr4 errM0x0y = Min(err5, err7);

    // !lossless, !lossless
    if (errM0x0y > Scr4(0.0f))
      errM0x0y = FitError<5,5>(xyz, min5, max5, err5);
    if (errM0x0y > Scr4(0.0f))
      errM0x0y = FitError<7,7>(xyz, min7, max7, err7);
#else
    // binary search, tangent-fitting
    Scr4 errM0x0y = FitError<5,5>(xyz, min5, max5, err5);
    // !lossless
    if (errM0x0y > Scr4(0.0f))
      errM0x0y = FitError<7,7>(xyz, min7, max7, err7);
#endif
  }

  // fix the range to be the minimum in each case
//FixRange(min5, max5, 5);
//FixRange(min7, max7, 7);

  // set up the 5-coord code book
  u8 codes5x[8];
  u8 codes5y[8];
  codes5x[0] = (u8)min5.R();
  codes5x[1] = (u8)max5.R();
  codes5y[0] = (u8)min5.G();
  codes5y[1] = (u8)max5.G();
  for (int i = 1; i < 5; ++i) {
    codes5x[1 + i] = (u8)(((5 - i) * min5.R() + i * max5.R()) / 5);
    codes5y[1 + i] = (u8)(((5 - i) * min5.G() + i * max5.G()) / 5);
  }
  codes5x[6] = 0x00;
  codes5y[6] = 0x00;
  codes5x[7] = 0xFF;
  codes5y[7] = 0xFF;

  // set up the 7-coord code book
  u8 codes7x[8];
  u8 codes7y[8];
  codes7x[0] = (u8)min7.R();
  codes7x[1] = (u8)max7.R();
  codes7y[0] = (u8)min7.G();
  codes7y[1] = (u8)max7.G();
  for (int i = 1; i < 7; ++i) {
    codes7x[1 + i] = (u8)(((7 - i) * min7.R() + i * max7.R()) / 7);
    codes7y[1 + i] = (u8)(((7 - i) * min7.G() + i * max7.G()) / 7);
  }

  // fit the data to both code books
  u8 indices5x[16];
  u8 indices5y[16];
  u8 indices7x[16];
  u8 indices7y[16];

  Scr4 err5 = FitCodes(xyz, mask, codes5x, indices5x, codes5y, indices5y);
  Scr4 err7 = FitCodes(xyz, mask, codes7x, indices7x, codes7y, indices7y);

  // save the block with least error
  if (!(err5 > err7)) {
    WriteNormalBlock5(min5.R(), max5.R(), indices5x, blockx);
    WriteNormalBlock5(min5.G(), max5.G(), indices5y, blocky);
  }
  else {
    WriteNormalBlock7(min7.R(), max7.R(), indices7x, blockx);
    WriteNormalBlock7(min7.G(), max7.G(), indices7y, blocky);
  }
}

void DecompressNormalBtc5(u8* xyzd, void const* blockx, void const* blocky)
{
  // get the two coord values
  u8 const* bytesx = reinterpret_cast< u8 const* >(blockx);
  u8 const* bytesy = reinterpret_cast< u8 const* >(blocky);
  int coord0x = bytesx[0];
  int coord1x = bytesx[1];
  int coord0y = bytesy[0];
  int coord1y = bytesy[1];

  // compare the values to build the codebook
  u8 codesx[8];
  codesx[0] = (u8)coord0x;
  codesx[1] = (u8)coord1x;
  if (coord0x <= coord1x) {
    // use 5-coord codebook
    for (int i = 1; i < 5; ++i)
      codesx[1 + i] = (u8)(((5 - i) * coord0x + i * coord1x) / 5);

    codesx[6] = 0x00;
    codesx[7] = 0xFF;
  }
  else {
    // use 7-coord codebook
    for (int i = 1; i < 7; ++i)
      codesx[1 + i] = (u8)(((7 - i) * coord0x + i * coord1x) / 7);
  }

  u8 codesy[8];
  codesy[0] = (u8)coord0y;
  codesy[1] = (u8)coord1y;
  if (coord0y <= coord1y) {
    // use 5-coord codebook
    for (int i = 1; i < 5; ++i)
      codesy[1 + i] = (u8)(((5 - i) * coord0y + i * coord1y) / 5);

    codesy[6] = 0x00;
    codesy[7] = 0xFF;
  }
  else {
    // use 7-coord codebook
    for (int i = 1; i < 7; ++i)
      codesy[1 + i] = (u8)(((7 - i) * coord0y + i * coord1y) / 7);
  }

  // decode the indices
  u8 indicesx[16];
  u8 indicesy[16];
  u8 const* srcx = bytesx + 2;
  u8 const* srcy = bytesy + 2;
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

  // write out the indexed codebook values
  for (int i = 0; i < 16; ++i) {
    Vec3 scale  = Vec3( 1.0f / 127.5f);
    Vec3 offset = Vec3(-1.0f * 127.5f);
    Col3 _xyz0  = Col3(codesx[indicesx[i]], codesy[indicesy[i]]);
    Vec3 cxyz0  = scale * (offset + _xyz0); cxyz0 = Complement(cxyz0);
    Vec3 scalei = Vec3( 1.0f * 127.5f);
	 cxyz0  = (scalei * cxyz0) - offset;
	 _xyz0  = FloatToInt<true>(cxyz0);

    PackBytes(_xyz0, *((int *)xyzd[4 * i]));
  }
}
#endif

/* *****************************************************************************
 */
#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
#endif

} // namespace squish

#endif
