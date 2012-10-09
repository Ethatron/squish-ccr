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

#define	SQUISH_USE_PRE		// means: we include all sources (a must for AMP & COMPUTE)
#define	SQUISH_USE_COMPUTE	// we request the compute-compatible code
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

// -----------------------------------------------------------------------------
// can be configured:
//
//  #define READ_STRUCTURED	read from a StructuredBuffer<uint4>
//  #undef  READ_STRUCTURED	read from a Texture2D<uint4>
//  #define WRTE_STRUCTURED	write to a RWStructuredBuffer<uint2>
//  #undef  WRTE_STRUCTURED	write to a RWTexture2D<uint2>
//  #define SIZE_CBUFFER	you pass dimensions yourself
//  #undef  SIZE_CBUFFER	we figure the dimensions out ourselfs
//
// Notes:
// - READ_STRUCTURED is bandwidth intensive, read 4 ints, no datatype
//   conversion
// - not WRTE_STRUCTURED is a bit odd, and you have to setup a 2 component
//   DXT/BTC/BTC imposter
// - not WRTE_STRUCTURED + not READ_STRUCTURED will make all i/o dimensionality
//   implicit, no need to calculate anything at all

// -----------------------------------------------------------------------------
// size of the input block in pixels (4x4 chars)
#define BLOCK_SIZE_Y		4
#define BLOCK_SIZE_X		4
#define BLOCK_SIZE		(BLOCK_SIZE_Y * BLOCK_SIZE_X)

#ifdef	READ_STRUCTURED
StructuredBuffer<uint4> g_Input;
#else
Texture2D<uint4> g_Input;
#endif

// -----------------------------------------------------------------------------
// size of the output codes in uint2s (2 or 4 longs)
#define CODES_SIZE_Y		1
#define CODES_SIZE_X		1
#define CODES_SIZE		(CODES_SIZE_Y * CODES_SIZE_X)

#ifdef	WRTE_STRUCTURED
RWStructuredBuffer<uint2> g_Output;
#else
RWTexture2D<uint2> g_Output;
#endif

// -----------------------------------------------------------------------------
#ifdef	SIZE_CBUFFER
cbuffer cbCS
{
    uint bwidth;
    uint bheight;
};
#endif

// -----------------------------------------------------------------------------
groupshared int ins[16][4];
groupshared unsigned int ous[2][2];

[numthreads(BLOCK_SIZE_X, BLOCK_SIZE_Y, 1)]
void main(
	uint  TR : SV_GroupIndex,
	uint3 DI : SV_DispatchThreadID,
	uint3 TI : SV_GroupThreadID,
	uint3 GI : SV_GroupID
)
{
  /* get the flattened thread-id -------------------------------------------- */

//const int thread = (                      TI.z) * (BLOCK_SIZE_Y           * BLOCK_SIZE_X         ) +
//                   (                      TI.y) * (                         BLOCK_SIZE_X         ) +
//                   (                      TI.x) * (                                    1         );
  const int thread = TR;

  /* get the dimension of the used input type ------------------------------- */

#ifndef	SIZE_CBUFFER
  uint bwidth, bheight;

#ifdef	READ_STRUCTURED
  g_Input.GetDimensions(bheight, bwidth); bheight /= bwidth;
#else
  g_Input.GetDimensions(bwidth, bheight);
#endif
#endif

  /* read the pixel(s) into thread-shared memory----------------------------- */

#ifdef	READ_STRUCTURED
  // calculate the input "array" location
  const int posloc = (DI.z                      ) * (BLOCK_SIZE_Y * bheight * BLOCK_SIZE_X * bwidth) +
                     (DI.y                      ) * (                         BLOCK_SIZE_X * bwidth) +
                     (DI.x                      ) * (                                    1         );
//const int posloc = (GI.z *            1 + TI.z) * (BLOCK_SIZE_Y * bheight * BLOCK_SIZE_X * bwidth) +
//                   (GI.y * BLOCK_SIZE_Y + TI.y) * (                         BLOCK_SIZE_X * bwidth) +
//                   (GI.x * BLOCK_SIZE_X + TI.x) * (                                    1         );
//const int tilloc = (GI.z *            1       ) * (BLOCK_SIZE_Y * bheight * BLOCK_SIZE_X * bwidth) +
//                   (GI.y * BLOCK_SIZE_Y       ) * (                         BLOCK_SIZE_X * bwidth) +
//                   (GI.x * BLOCK_SIZE_X       ) * (                                    1         );

  uint4 pel = g_Input[posloc];
#else
  // calculate the input "texture" position
//const uint3 texloc = uint3((GI.x * BLOCK_SIZE_X + TI.x),
//                           (GI.y * BLOCK_SIZE_Y + TI.y),
//                           (GI.z *            1 + TI.z));
  const uint2 texloc = uint2((GI.x * BLOCK_SIZE_X + TI.x),
                             (GI.y * BLOCK_SIZE_Y + TI.y));

//uint4 pel = g_Input.Load(texloc);
  uint4 pel = g_Input[texloc];
#endif

  ins[thread][0] = pel.a;
  ins[thread][1] = pel.r;
  ins[thread][2] = pel.g;
  ins[thread][3] = pel.b;

  /* wait for completion */
  GroupMemoryBarrier();

  /* TODO: check why this leads to compile errors */
#if 0
  /* compress */
  squish::CompressColorBtc1(
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
  squish::code64 block = ous[1];

  // create the minimal point set
  squish::colours.CountSet(barrier, thread, ins, mask, true, trans);

  // compress colour
  squish::CompressColorBtc(barrier, thread, ins, squish::colours, block, metric, trans, fit, yArr, lArr);

  if (thread == 0) {
    ous[1] = block;
  }
#endif

  /* wait for completion */
  GroupMemoryBarrierWithGroupSync();

  /* write the code(s) into destination ------------------------------------- */

  /* throw it out */
  if (thread == 0) {
#ifdef	WRTE_STRUCTURED
    // calculate the output "array" location
    const int cdeloc = (DI.z                      ) * (CODES_SIZE_Y * bheight * CODES_SIZE_X * bwidth) +
                       (DI.y                      ) * (                         CODES_SIZE_X * bwidth) +
                       (DI.x                      ) * (                                    1         );

    g_Output[cdeloc] = uint2(ous[1][0], ous[1][1]);
#else
    // calculate the output "texture" position
//  const uint3 tcdloc = uint3((GI.x * CODES_SIZE_X + TI.x),
//                             (GI.y * CODES_SIZE_Y + TI.y),
//                             (GI.z *            1 + TI.z));
    const uint2 tcdloc = uint2((GI.x * CODES_SIZE_X + TI.x),
                               (GI.y * CODES_SIZE_Y + TI.y));

    g_Output[tcdloc] = uint2(ous[1][0], ous[1][1]);
#endif
  }
}
