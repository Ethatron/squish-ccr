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

#define	USE_PRE		// means: we include all sources (a must for AMP & COMPUTE)
#define	USE_COMPUTE	// we request the compute-compatible code
#include "squish.h"

#define	ONLY_ARRAY
#include "singlecolourlookup_ccr.inl"
#include "degeneracy_ccr.inl"

#include "maths.cpp"
#include "alpha.cpp"
#include "colourblock.cpp"

#include "colourset.cpp"
#include "colourfit.cpp"
#include "singlecolourfit.cpp"
#include "rangefit.cpp"
#include "clusterfit.cpp"
#include "squish.cpp"

RWStructuredBuffer<uint4> g_Input;
RWStructuredBuffer<uint2> g_Output;

#define BLOCK_SIZE_Y		4
#define BLOCK_SIZE_X		4
#define BLOCK_SIZE		(BLOCK_SIZE_Y * BLOCK_SIZE_X)

groupshared int ins[16][4];
groupshared unsigned int ous[2][2];

cbuffer cbCS
{
    uint bwidth;
    uint bheight;
};

[numthreads(BLOCK_SIZE_X, BLOCK_SIZE_Y, 1)]
void main(
	uint  TR : SV_GroupIndex,
	uint3 DI : SV_DispatchThreadID,
	uint3 TI : SV_GroupThreadID,
	uint3 GI : SV_GroupID
)
{
  const int posloc = (DI.z                      ) * (BLOCK_SIZE_Y * bheight * BLOCK_SIZE_X * bwidth) +
                     (DI.y                      ) * (                         BLOCK_SIZE_X * bwidth) +
                     (DI.x                      ) * (                                    1         );
//const int posloc = (GI.z *            1 + TI.z) * (BLOCK_SIZE_Y * bheight * BLOCK_SIZE_X * bwidth) +
//                   (GI.y * BLOCK_SIZE_Y + TI.y) * (                         BLOCK_SIZE_X * bwidth) +
//                   (GI.x * BLOCK_SIZE_X + TI.x) * (                                    1         );
  const int tilloc = (GI.z *            1       ) * (BLOCK_SIZE_Y * bheight * BLOCK_SIZE_X * bwidth) +
                     (GI.y * BLOCK_SIZE_Y       ) * (                         BLOCK_SIZE_X * bwidth) +
                     (GI.x * BLOCK_SIZE_X       ) * (                                    1         );
//const int thread = (                      TI.z) * (BLOCK_SIZE_Y           * BLOCK_SIZE_X         ) +
//                   (                      TI.y) * (                         BLOCK_SIZE_X         ) +
//                   (                      TI.x) * (                                    1         );
  const int thread = TR;

  /* read the pixel into thread-shared memory */
  uint4 pel = g_Input[posloc];

  ins[thread][0] = pel.a;
  ins[thread][1] = pel.r;
  ins[thread][2] = pel.g;
  ins[thread][3] = pel.b;

  /* wait for completion */
  GroupMemoryBarrier();

  /* TODO: check why this leads to compile errors */
#if 0
  /* compress */
  squish::CompressColorDxt1(
    0,
    thread,
    ins,
    0xFFFF,
    ous[1],

    // flags and parameters
    SQUISH_METRIC_PERCEPTUAL,
    true,
    0,

    // pass LUT arrays
    lookup_c34a57_ccr,
    lookup_34_56_ccr
  );
#else
  int barrier = 0;
  int mask = 0xFFFF;
  int metric = SQUISH_METRIC_PERCEPTUAL;
  bool trans = true;
  int fit = SQUISH_FIT_RANGE;
  IndexBlockLUT yArr = lookup_c34a57_ccr;
  SingleColourLUT lArr = lookup_34_56_ccr;

  // create the minimal point set
  squish::colours.CountSet(barrier, thread, ins, mask, true, trans);

  // compress colour
  squish::CompressColorDxt(barrier, thread, ins, squish::colours, ous[1], metric, trans, fit, yArr, lArr);
#endif

  /* wait for completion */
  GroupMemoryBarrierWithGroupSync();

  /* throw it out */
  if (thread == 0) {
//  threaded_xch(g_Output[0], ous[1][0]);
//  threaded_xch(g_Output[1], ous[1][1]);

//  threaded_xch(g_Output[posloc][thread], ous[1][thread]);
    g_Output[posloc] = uint2(ous[1][0], ous[1][1]);
  }
}
