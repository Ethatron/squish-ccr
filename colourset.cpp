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

#include "colourset.h"

namespace squish {

#if	!defined(USE_PRE)
ColourSet::ColourSet( u8 const* rgba, int mask, int flags )
  : m_count( 0 ),
	m_transparent( false )
{
  // check the compression mode for dxt1
  bool isDxt1 = ( ( flags & kDxt1 ) != 0 );
  bool weightByAlpha = ( ( flags & kWeightColourByAlpha ) != 0 );

  /*
  static const float dw[] = {
    sqrtf(1.5 * 1.5 + 1.5 * 1.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 1.5 * 1.5),

    sqrtf(0.5 * 0.5 + 1.5 * 1.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 1.5 * 1.5),

    sqrtf(0.5 * 0.5 + 1.5 * 1.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 0.5 * 0.5),
    sqrtf(0.5 * 0.5 + 1.5 * 1.5),

    sqrtf(1.5 * 1.5 + 1.5 * 1.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 0.5 * 0.5),
    sqrtf(1.5 * 1.5 + 1.5 * 1.5)
  };
  */

  // create the minimal set
  for( int i = 0; i < 16; ++i )
  {
    // check this pixel is enabled
    int bit = 1 << i;
    if( ( mask & bit ) == 0 )
    {
      m_remap[i] = -1;
      continue;
    }

    // check for transparent pixels when using dxt1
    if( isDxt1 && rgba[4*i + 3] < 128 )
    {
      m_remap[i] = -1;
      m_transparent = true;
      continue;
    }

    // loop over previous points for a match
    for( int j = 0;; ++j )
    {
      // allocate a new point
      if( j == i )
      {
	// normalise coordinates to [0,1]
	float x = ( float )rgba[4*i] / 255.0f;
	float y = ( float )rgba[4*i + 1] / 255.0f;
	float z = ( float )rgba[4*i + 2] / 255.0f;

	// ensure there is always non-zero weight even for zero alpha
	float w = ( float )( rgba[4*i + 3] + 1 ) / 256.0f;

	// add the point
	m_points[m_count] = Vec3( x, y, z );
	m_weights[m_count] = ( weightByAlpha ? w : 1.0f );
	m_remap[i] = m_count;

	// advance
	++m_count;
	break;
      }

      // check for a match
      int oldbit = 1 << j;
      bool match = ( ( mask & oldbit ) != 0 )
	&& ( rgba[4*i] == rgba[4*j] )
	&& ( rgba[4*i + 1] == rgba[4*j + 1] )
	&& ( rgba[4*i + 2] == rgba[4*j + 2] )
	&& ( rgba[4*j + 3] >= 128 || !isDxt1 );
      if( match )
      {
	// get the index of the match
	int index = m_remap[j];

	// ensure there is always non-zero weight even for zero alpha
	float w = ( float )( rgba[4*i + 3] + 1 ) / 256.0f;

	// map to this point and increase the weight
	m_weights[index] += ( weightByAlpha ? w : 1.0f );
	m_remap[i] = index;
	break;
      }
    }
  }

  // square root the weights
  for( int i = 0; i < m_count; ++i )
    m_weights[i] = std::sqrt( m_weights[i] );
}

void ColourSet::RemapIndices( u8 const* source, u8* target ) const
{
  for (int i = 0; i < 16; ++i) {
    u8 t = 3; t = ((m_remap[i] == -1) ? t : source[m_remap[i]]); target[i] = t;
  }
}
#endif

#if	defined(USE_AMP) || defined(USE_COMPUTE)
void ColourSet_CCR::CountSet(tile_barrier barrier, const int thread, pixel16 rgba, int mask, const bool tresh, const bool trans) amp_restricted
{
  // check the compression mode for dxt1
  const bool isDxt1 = (tresh);
  const bool weightByAlpha = (trans);

  // clear counters
  threaded_cse(0) {
    m_count = 0;
    m_transparent = 0;
  }

  wavefrnt_for(rscan, 16) {
    // clear out all values, so we may use all 16 of
    // them, instead of variable count at any time
    m_points [rscan] = 0.0f;
    m_weights[rscan] = 0;

    // use no code by default
    m_remap  [rscan] = -1;
  }

  // create the minimal set
  for (int i = 0; i < 16; ++i) {
    // check this pixel is enabled
    int bit = 1 << i;
    if ((mask & bit) != 0) {
      // check for transparent pixels when using dxt1
      if (isDxt1 && (rgba[i][0] < 128)) {
        // AMP: all thread end up here (CXCHG won't block concurrent reads)
        int zero = 0; threaded_set(m_transparent, zero, 1);
      }
      else {
        // make writes to "m_remap" visible to all
        tile_static_memory_fence(barrier);

        // loop over previous points for a match (AMP: parallel search)
        threaded_for(j, i) {
          // check for a match
          int oldbit = 1 << j;

          // get the index of the match
          bool match = ((mask & oldbit) != 0)
	    && (rgba[i][3] == rgba[j][3])
	    && (rgba[i][2] == rgba[j][2])
	    && (rgba[i][1] == rgba[j][1])
	    && (!isDxt1 || (rgba[j][0] >= 128));

          // AMP: there is only one thread max. which could find a match
          //      so it's only one thread which writes to index, no
          //      atomic needed
          if (match)
	    m_remap[i] = m_remap[j];
        }

        // make writes to "m_remap" visible to thread 0
        tile_static_memory_fence(barrier);

        // allocate a new point
        threaded_cse(0) {
          int index = m_remap[i];
          if (index < 0) {
	    // get the index of the new point
	    index = m_count++;

	    // add the point
	    m_weights[index] = 0;
	    m_points [index] = float3(
	      (float)rgba[i][3],
	      (float)rgba[i][2],
	      (float)rgba[i][1]
	    );

	    m_remap[i] = index;
          }

          // ensure there is always non-zero weight even for zero alpha
          float w = (float)(rgba[i][0] + 1) / 256.0f;

          // map to this point and increase the weight
          m_weights[index] += (weightByAlpha ? w : 1.0f);
        }
      }
    }
  }

  // square root the weights
#if	!defined(USE_COMPUTE)
  wavefrnt_for(wscan, 16) {
    // normalise coordinates to [0,1]
    m_points [wscan] /= 255.0f;
    m_weights[wscan]  = sqrtf(m_weights[wscan]);
  }
#else
  // error X3695: race condition writing to shared memory detected, consider making this write conditional.
  threaded_for(wscan, m_count) {
    // normalise coordinates to [0,1]
    m_points [wscan] /= 255.0f;
    m_weights[wscan]  = sqrtf(m_weights[wscan]);
  }
#endif
}

int ColourSet_CCR::GetCount() amp_restricted {
  return m_count;
}

point16 ColourSet_CCR::GetPoints() amp_restricted {
  return m_points;
}

weight16 ColourSet_CCR::GetWeights() amp_restricted {
  return m_weights;
}

bool ColourSet_CCR::IsTransparent() amp_restricted {
  return !!m_transparent;
}

void ColourSet_CCR::WriteSet(tile_barrier barrier, const int thread, lineC2 cline, inout index16x2 source, out code64 block, const int is4,
			     IndexBlockLUT yArr) amp_restricted
{
  // build the block if we win
  {
    // remap the indices
    RemapIndices(barrier, thread, source);

    // save the block
    {
      if (!is4)
	WriteColourBlock3(barrier, thread, cline, m_indices, block, yArr);
      else
	WriteColourBlock4(barrier, thread, cline, m_indices, block, yArr);
    }
  }
}

void ColourSet_CCR::RemapIndices(tile_barrier barrier, const int thread, inout index16x2 source) amp_restricted
{
  // make writes to "source"/"m_remap" visible to all
  tile_static_memory_fence(barrier);

#if	!defined(USE_COMPUTE)
  // remap to palette-indices
  wavefrnt_for(rmscan, 16) {
    m_indices[rmscan] = source[1][m_remap[rmscan]];
  }
#else
  // error X3500: array reference cannot be used as an l-value not natively addressable
  threaded_cse(0) {
    m_indices[ 0] = source[1][m_remap[ 0]];
    m_indices[ 1] = source[1][m_remap[ 1]];
    m_indices[ 2] = source[1][m_remap[ 2]];
    m_indices[ 3] = source[1][m_remap[ 3]];
    m_indices[ 4] = source[1][m_remap[ 4]];
    m_indices[ 5] = source[1][m_remap[ 5]];
    m_indices[ 6] = source[1][m_remap[ 6]];
    m_indices[ 7] = source[1][m_remap[ 7]];
    m_indices[ 8] = source[1][m_remap[ 8]];
    m_indices[ 9] = source[1][m_remap[ 9]];
    m_indices[10] = source[1][m_remap[10]];
    m_indices[11] = source[1][m_remap[11]];
    m_indices[12] = source[1][m_remap[12]];
    m_indices[13] = source[1][m_remap[13]];
    m_indices[14] = source[1][m_remap[14]];
    m_indices[15] = source[1][m_remap[15]];
  }
#endif
}
#endif

} // namespace squish
