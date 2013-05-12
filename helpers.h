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

#ifndef SQUISH_HELPERS_H
#define SQUISH_HELPERS_H

#include <squish.h>

namespace squish {
  
/* *****************************************************************************
 */
template<typename dtyp>
class Weight;

template<>
class Weight<u8> {
  
private:
  u8   w;
  Scr4 W;

public:
  doinline Weight(const u8 *rgba, int pos, const u8 wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = rgba[(4 * pos) + 3] | wgt;
    W = Scr4(w + 1) * Scr4(1.0f / 256.0f);
  }
  
  doinline Weight(const u8 a, const u8 wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = a | wgt;
    W = Scr4(w + 1) * Scr4(1.0f / 256.0f);
  }

  doinline float GetWeight() const {
    return W.X(); }
  doinline Scr4 GetWeights() const {
    return W; }

  doinline bool IsOne() const {
    return !(u8)(~w); }
};

template<>
class Weight<u16> {
  
private:
  u16  w;
  Scr4 W;

public:
  doinline Weight(const u16 *rgba, int pos, const u16 wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = rgba[(4 * pos) + 3] | wgt;
    W = Scr4(w + 1) * Scr4(1.0f / 65536.0f);
  }
  
  doinline Weight(const u16 a, const u16 wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = a | wgt;
    W = Scr4(w + 1) * Scr4(1.0f / 65536.0f);
  }

  doinline float GetWeight() const {
    return W.X(); }
  doinline Scr4 GetWeights() const {
    return W; }

  doinline bool IsOne() const {
    return !(u16)(~w); }
};

template<>
class Weight<f23> {
  
private:
  Scr4 w;
  Scr4 W;

public:
  doinline Weight(const f23 *rgba, int pos, const Scr4 &wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = Max(Scr4(rgba[(4 * pos) + 3]), wgt);
    W = Scr4(w);
  }
  
  doinline Weight(const f23 a, const f23 wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = Max(Scr4(a), Scr4(wgt));
    W = w;
  }
  
  doinline Weight(const f23 a, const Scr4 &wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = Max(Scr4(a), wgt);
    W = w;
  }

  doinline float GetWeight() const {
    return W.X(); }
  doinline Scr4 GetWeights() const {
    return W; }

  doinline bool IsOne() const {
    return !CompareFirstLessThan(w, Scr4(1.0f)); }
};

template<>
class Weight<Scr4> {
  
private:
  Scr4 w;
  Scr4 W;

public:
  doinline Weight(const Scr4 (&weights)[4][16], int pos, const Scr4 &wgt) {
    // ensure there is always non-zero weight even for zero alpha
    w = Max(weights[3][pos], wgt);
    W = Scr4(w);
  }

  doinline float GetWeight() const {
    return W.X(); }
  doinline Scr4 GetWeights() const {
    return W; }

  doinline bool IsOne() const {
    return !CompareFirstLessThan(w, Scr4(1.0f)); }
};

} // namespace squish

#endif // ndef SQUISH_HELPERS_H
