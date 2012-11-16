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

namespace squish {

#if	defined(SQUISH_USE_AMP) || defined(SQUISH_USE_COMPUTE)
void ColourSingleFit_CCR::AssignSet(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, const int metric, const int fit) amp_restricted
{
  ColourFit_CCR::AssignSet(barrier, thread, m_colours, metric, fit);

  // grab the single colour
  point16 values = m_colours.GetPoints();

  threaded_cse(0) {
    // initialise the color
//  m_colour = FloatToInt(values[0] * 255.0f, 255);
    m_colour = QuantizeFloatToInt(values[0], 255);

    // initialise the index
//  m_index = 0;
  }
}

void ColourSingleFit_CCR::Compress(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block, const bool trans,
				   IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted
{
  /* all or nothing branches, OK, same for all threads */
  const bool isBtc1 = (trans);
  if (isBtc1 && m_colours.IsTransparent())
    Compress3 (barrier, thread, m_colours, block, yArr, lArr);
  else if (!isBtc1)
    Compress4 (barrier, thread, m_colours, block, yArr, lArr);
  else
    Compress34(barrier, thread, m_colours, block, yArr, lArr);
}

#define CASE3	0
#define CASE4	1

void ColourSingleFit_CCR::Compress3(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
				    IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted
{
  // find the best end-points and index
  ComputeEndPoints(barrier, thread, CASE3, lArr);

  // build the block if we win
  {
    // remap the one index (AMP: prevent writing to same location)
//  m_matches[1][   0  ] = m_index;
    m_matches[1][thread] = m_index;

    // save the block
    ColourFit_CCR::Compress3(barrier, thread, m_colours, block, yArr);
  }
}

void ColourSingleFit_CCR::Compress4(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
				    IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted
{
  // find the best end-points and index
  ComputeEndPoints(barrier, thread, CASE4, lArr);

  // build the block if we win
  {
    // remap the one index (AMP: prevent writing to same location)
//  m_matches[1][   0  ] = m_index;
    m_matches[1][thread] = m_index;

    // save the block
    ColourFit_CCR::Compress4(barrier, thread, m_colours, block, yArr);
  }
}

void ColourSingleFit_CCR::Compress34(tile_barrier barrier, const int thread, ColourSet_CCRr m_colours, out code64 block,
				     IndexBlockLUT yArr, ColourSingleLUT lArr) amp_restricted
{
  // find the best end-points and index
  int is4 = ComputeEndPoints(barrier, thread, lArr);

  // build the block if we win
  {
    // remap the one index (AMP: prevent writing to same location)
//  m_matches[1][   0  ] = m_index;
    m_matches[1][thread] = m_index;

    // save the block
    ColourFit_CCR::Compress34(barrier, thread, m_colours, block, is4, yArr);
  }
}

#if	defined(SQUISH_USE_COMPUTE)
  tile_static ::ColourSingleLookup_CCR sources[3 * 2];
  tile_static float3 lines[4][2];
  tile_static int /*error[4],*/ err[2];
#endif

#define SIDEL	0
#define SIDEH	1
#define SIDES	2
#define CBITS	1	// number of bits
#define CMASK	1	// mask

void ColourSingleFit_CCR::ComputeEndPoints(tile_barrier barrier, const int thread, const int is4,
					   ColourSingleLUT lArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static ::ColourSingleLookup_CCR sources[3];
  tile_static float3 lines[SIDES][CVALS];
  tile_static int error[SIDES];
#endif
  int e;

  /* 3 colors, by 2 sides (this or next) */
  threaded_cse(0) {
    // AMP: pull from constant to shared
    {
      // grab the lookup table and index for this channel
      // store a pointer to the source for this channel
#if	defined(SQUISH_USE_COMPUTE)
      sources[(0 << 0)] = lArr[is4][m_colour.r];
      sources[(1 << 0)] = lArr[is4][m_colour.g];
      sources[(2 << 0)] = lArr[is4][m_colour.b];
#else
      sources[(0 << 0)] = lArr(is4, m_colour.r);
      sources[(1 << 0)] = lArr(is4, m_colour.g);
      sources[(2 << 0)] = lArr(is4, m_colour.b);
#endif
    }
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // calculate the interpolation lines (2 choices)
  threaded_for(ljoint, CVALS * SIDES) {
    const int side = ljoint &  CMASK;
    const int cpos = ljoint >> CBITS;

    lines[side][cpos] = float3(
      (float)SB_CLINE(sources[(0 << 0)].sources[side + 0], cpos),
      (float)SB_CLINE(sources[(1 << 0)].sources[side + 2], cpos),
      (float)SB_CLINE(sources[(2 << 0)].sources[side + 0], cpos)
    ) /
      float3(31.0f, 63.0f, 31.0f);
  }

  // accumulate the errors
  threaded_for(ejoint, SIDES) {
    const int index = ejoint;
    const int side  = ejoint;

    error[index] =
      SB_ERRSQR(sources[(0 << 0)].sources[side + 0]) +
      SB_ERRSQR(sources[(1 << 0)].sources[side + 2]) +
      SB_ERRSQR(sources[(2 << 0)].sources[side + 0]);
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // AMP: final reduction (make the choice)
  threaded_cse(0) {
    // keep it if the error is lower
    e = (error[SIDEH] < error[SIDEL] ? SIDEH : SIDEL);

    m_line[CSTRT] = lines[e][CSTRT];
    m_line[CSTOP] = lines[e][CSTOP];
    m_index = (ccr8)(2 * e);
  }
}

#undef	CASE3
#undef	CASE4
#undef	CBITS
#undef	CMASK
#undef	SIDEL
#undef	SIDEH
#undef	SIDES

#define SIDEL	0
#define SIDEH	2
#define CASE3	0
#define CASE4	1
#define CASES	2
#define CNUMS	4
#define CBITS	2	// number of bits
#define CMASK	3	// mask

int ColourSingleFit_CCR::ComputeEndPoints(tile_barrier barrier, const int thread,
					  ColourSingleLUT lArr) amp_restricted
{
#if	!defined(SQUISH_USE_COMPUTE)
  tile_static ::ColourSingleLookup_CCR sources[3 * CASES];
  tile_static float3 lines[CNUMS][CVALS];
  tile_static int error[CNUMS], err[CASES];
#endif

  /* 3 colors, by 2 sides (this or next), by 2 types (3 or 4) */
  threaded_for(sjoint, CASES) {
    const int which = sjoint;

    {
      // store a pointer to the source for this channel
#if	defined(SQUISH_USE_COMPUTE)
      sources[(0 << 1) + which] = lArr[which][m_colour.r];
      sources[(1 << 1) + which] = lArr[which][m_colour.g];
      sources[(2 << 1) + which] = lArr[which][m_colour.b];
#else
      sources[(0 << 1) + which] = lArr(which, m_colour.r);
      sources[(1 << 1) + which] = lArr(which, m_colour.g);
      sources[(2 << 1) + which] = lArr(which, m_colour.b);
#endif
    }
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // calculate the interpolation lines (4 choices)
  threaded_for(ljoint, CVALS * CNUMS) {
    const int wsde  = (ljoint & 3);
    const int which = (ljoint & 1);
    const int side  = (ljoint >> 1) & 1;
    const int cpos  = (ljoint >> CBITS);

    lines[wsde][cpos] = float3(
      (float)SB_CLINE(sources[(0 << 1) + which].sources[side + 0], cpos),
      (float)SB_CLINE(sources[(1 << 1) + which].sources[side + 2], cpos),
      (float)SB_CLINE(sources[(2 << 1) + which].sources[side + 0], cpos)
    ) /
      float3(31.0f, 63.0f, 31.0f);
  }

  // accumulate the errors
  threaded_for(ejoint, CNUMS) {
    const int which = (ejoint & 1);
    const int side  = (ejoint >> 1);

    error[ejoint] =
      SB_ERRSQR(sources[(0 << 1) + which].sources[side + 0]) +
      SB_ERRSQR(sources[(1 << 1) + which].sources[side + 2]) +
      SB_ERRSQR(sources[(2 << 1) + which].sources[side + 0]);
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // choose each side for 3 and 4
  threaded_for(which, CASES) {
    err[which] = (error[which + SIDEL] < error[which + SIDEH] ? which + SIDEL : which + SIDEH);
  }

  // make it visible
  tile_static_memory_fence(barrier);

  // AMP: all thread end up here, and all make this compare
  int e = (error[err[CASE3]] < error[err[CASE4]] ? err[CASE3] : err[CASE4]);

  // AMP: final reduction (make the choice)
  threaded_cse(0) {
    // keep it if the error is lower
    m_line[CSTRT] = lines[e][CSTRT];
    m_line[CSTOP] = lines[e][CSTOP];
    m_index = (ccr8)(e & SIDEH);
  }

  // is4
  return (e & CASE4);
}

#undef	CASE3
#undef	CASE4
#undef	CBITS
#undef	CMASK
#undef	SIDEL
#undef	SIDEH
#undef	SIDES
#endif

} // namespace squish
