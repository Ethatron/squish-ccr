/* -----------------------------------------------------------------------------

	Copyright (c) 2006 Simon Brown                          si@sjbrown.co.uk

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

/*! @file

	@brief	This program tests the error for 1 and 2-colour DXT/BTC compression.

	This tests the effectiveness of the DXT/BTC compression algorithm for all
	possible 1 and 2-colour blocks of pixels.
*/

#include <squish.h>
#include <iostream>
#include <cmath>
#include <cfloat>

using namespace squish;

double GetColourError( unsigned int (&a)[4][4], unsigned int (&b)[4][4] )
{
	double error = 0.0;

	for( int i = 0; i < 4; ++i )
	for( int j = 0; j < 4; ++j )
	{
		for( int c = 0; c < 3; ++c )
		{
			int diff = ( int )((a[i][j] >> (c * 8)) & 0xFF)
			         - ( int )((b[i][j] >> (c * 8)) & 0xFF);
			error += ( double )( diff*diff );
		}
	}
	return error / 16.0;
}

void TestOneColour( int flags )
{
	unsigned int input[4][4];
	unsigned int output[4][4];
	unsigned int block[2][2];

	double avg = 0.0, min = DBL_MAX, max = -DBL_MAX;
	int counter = 0;

#if	defined(SQUISH_USE_AMP) && !defined(USE_AMP_DEBUG)
	Concurrency::array_view<const ColourSingleLookup_CCR, 2> lArr(2, 256, (const ColourSingleLookup_CCR *)::lookup_34_56_ccr);
	Concurrency::array_view<const   IndexBlockLookup_CCR, 2> yArr(4, 8,   (const   IndexBlockLookup_CCR *)::lookup_c34a57_ccr);
#endif

	// test all single-channel colours
	for( int i = 0; i < 4; ++i )
	for( int j = 0; j < 4; ++j )
		input[i][j] = 255U << 24;

	for( int channel = 0; channel < 3; ++channel )
	{
		for( int value = 0; value < 255; ++value )
		{
			// set the channnel value
			for( int i = 0; i < 4; ++i )
			for( int j = 0; j < 4; ++j )
				input[i][j] = (255 << 24) + (value << 16) + (value << 8) + (value << 0);

#if	defined(SQUISH_USE_AMP) && !defined(USE_AMP_DEBUG)
			/* constant buffer array */
			/* get a two-dimensional extend over the whole output (without re-cast to LONG),
			 * then get a tile-extend over that one ()
			 */
			Concurrency::extent<2> ee(4, 4);
			Concurrency::tiled_extent<4, 4> te(ee);

			Concurrency::array_view<const unsigned int, 2> sArr(4, 4, (const unsigned int *)input);
			Concurrency::array_view<      unsigned int, 2> dArr(1, 2, (      unsigned int *)block);

			Concurrency::parallel_for_each(te, [=](tiled_index<4, 4> elm) restrict(amp) {
			  /* generate this level's 4x4-block from the original surface */
			  const int y = elm.tile[0] * 4;
			  const int x = elm.tile[1] * 4;
			  const int ly = elm.local[0];
			  const int lx = elm.local[1];
			  const int lxy = ly * 4 + lx;

			  /* round down */
			  int posx = (x + lx) >> 2;
			  int posy = (y + ly) >> 2;

			  {
			    const int yl = ((y + ly) << 0);
			    const int xl = ((x + lx) << 0);

			    int pxl[4*4][DIM];
			    unsigned int blk[2][2];

			    for (int oy = 0; oy < 1; oy += 1)
			    for (int ox = 0; ox < 1; ox += 1) {
			      const int posx = (xl + ox) % 4;
			      const int posy = (yl + oy) % 4;

			      const ULONG &t = sArr(posy, posx);

			      pxl[oy][ox] = (t >> 24) & 0xFF;
			      pxl[oy][ox] = (t >> 16) & 0xFF;
			      pxl[oy][ox] = (t >>  8) & 0xFF;
			      pxl[oy][ox] = (t >>  0) & 0xFF;
			    }

			    // compress ...
			    squish::CompressColourBtc1u(elm.barrier, lxy, pxl, 0xFFFF, blk[1], SQUISH_METRIC_PERCEPTUAL, false, 0, yArr, lArr);

			    dArr(posy, posx + 0) = blk[1][0];
			    dArr(posy, posx + 1) = blk[1][1];
			  }
			});

			dArr.synchronize();
#else
#endif

			// ... and decompress
			Decompress( output, block, flags );

			// test the results
			double rm = GetColourError( input, output );
			double rms = std::sqrt( rm );

			// accumulate stats
			min = std::min( min, rms );
			max = std::max( max, rms );
			avg += rm;
			++counter;
		}
	}

	// finish stats
	avg = std::sqrt( avg/counter );

	// show stats
	std::cout << "one colour error (min, max, avg): "
		<< min << ", " << max << ", " << avg << std::endl;
}

void TestOneColourRandom( int flags )
{
	u8 input[4*16];
	u8 output[4*16];
	u8 block[16];

	double avg = 0.0, min = DBL_MAX, max = -DBL_MAX;
	int counter = 0;

	// test all single-channel colours
	for( int test = 0; test < 1000; ++test )
	{
		// set a constant random colour
		for( int channel = 0; channel < 3; ++channel )
		{
			u8 value = ( u8 )( rand() & 0xff );
			for( int i = 0; i < 16; ++i )
				input[4*i + channel] = value;
		}
		for( int i = 0; i < 16; ++i )
			input[4*i + 3] = 255;

		// compress and decompress
		Compress( input, block, flags );
		Decompress( output, block, flags );

		// test the results
		double rm = GetColourError( input, output );
		double rms = std::sqrt( rm );

		// accumulate stats
		min = std::min( min, rms );
		max = std::max( max, rms );
		avg += rm;
		++counter;
	}

	// finish stats
	avg = std::sqrt( avg/counter );

	// show stats
	std::cout << "random one colour error (min, max, avg): "
		<< min << ", " << max << ", " << avg << std::endl;
}

void TestTwoColour( int flags )
{
	u8 input[4*16];
	u8 output[4*16];
	u8 block[16];

	double avg = 0.0, min = DBL_MAX, max = -DBL_MAX;
	int counter = 0;

	// test all single-channel colours
	for( int i = 0; i < 16*4; ++i )
		input[i] = ( ( i % 4 ) == 3 ) ? 255 : 0;
	for( int channel = 0; channel < 3; ++channel )
	{
		for( int value1 = 0; value1 < 255; ++value1 )
		{
			for( int value2 = value1 + 1; value2 < 255; ++value2 )
			{
				// set the channnel value
				for( int i = 0; i < 16; ++i )
					input[4*i + channel] = ( u8 )( ( i < 8 ) ? value1 : value2 );

				// compress and decompress
				Compress( input, block, flags );
				Decompress( output, block, flags );

				// test the results
				double rm = GetColourError( input, output );
				double rms = std::sqrt( rm );

				// accumulate stats
				min = std::min( min, rms );
				max = std::max( max, rms );
				avg += rm;
				++counter;
			}
		}

		// reset the channel value
		for( int i = 0; i < 16; ++i )
			input[4*i + channel] = 0;
	}

	// finish stats
	avg = std::sqrt( avg/counter );

	// show stats
	std::cout << "two colour error (min, max, avg): "
		<< min << ", " << max << ", " << avg << std::endl;
}

int main()
{
	TestOneColourRandom( kBtc1 + kColourRangeFit );
	TestOneColour( kBtc1 );
	TestTwoColour( kBtc1 );
}
