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

#include <iostream>
#include <string>

struct SourceBlock
{
	int start;
	int end;
	int error;
};

struct TargetValue
{
	SourceBlock sources[8];
};

static void GenerateData( std::string const& name, int bits, int colours, bool dxt )
{
	TargetValue values[256];
	int indices = 2;

	if (dxt) {
		// start + 1 intermediate
		if( colours == 3 )
			indices = 2;
		// start + 1 intermediate (codebook is symmetric)
		else
			indices = 2;
	}
	else {
		// start + 1 intermediate (codebook is symmetric)
		if( colours == 4 )
			indices = 2;
		// start + 3 intermediates (codebook is symmetric)
		else  if( colours == 8 )
			indices = 4;
		// start + 7 intermediates (codebook is symmetric)
		else
			indices = 8;
	}

	// initialise the data
	for( int target = 0; target < 256; ++target )
		for( int index = 0; index < indices; ++index )
			values[target].sources[index].error = 255;

	// loop over all possible source points
	int count = ( 1 << bits );
	for( int value1 = 0; value1 < count; ++value1 )
	{
		for( int value2 = 0; value2 < count; ++value2 )
		{
			// compute the 8-bit endpoints
			int a = ( value1 << ( 8 - bits ) ) | ( value1 >> ( 2*bits - 8 ) );
			int b = ( value2 << ( 8 - bits ) ) | ( value2 >> ( 2*bits - 8 ) );

			// fill in the codebook with the these and intermediates
			int codes[8];
			codes[0] = a;

			if (dxt) {
				// start + 1 intermediate
				if( colours == 3 )
					codes[1] = ( a + b )/2;
				// start + 1 intermediate (codebook is symmetric)
				else
					codes[1] = ( 2*a + b )/3;
			}
			else {
				// start + 1 intermediate (codebook is symmetric)
				if( colours == 4 )
					codes[1] = ( 43*a + 21*b )/64;
				// start + 3 intermediates (codebook is symmetric)
				else  if( colours == 8 ) {
					codes[1] = ( 55*a + 9*b )/64;
					codes[2] = ( 46*a + 18*b )/64;
					codes[3] = ( 37*a + 27*b )/64;
				}
				// start + 7 intermediates (codebook is symmetric)
				else {
					codes[1] = ( 60*a + 4*b )/64;
					codes[2] = ( 55*a + 9*b )/64;
					codes[3] = ( 51*a + 13*b )/64;
					codes[4] = ( 47*a + 17*b )/64;
					codes[5] = ( 43*a + 21*b )/64;
					codes[6] = ( 38*a + 26*b )/64;
					codes[7] = ( 34*a + 30*b )/64;
				}
			}

			// mark each target point with the endpoints and index needed for it
			for( int index = 0; index < indices; ++index )
			{
				int target = codes[index];

				SourceBlock& block = values[target].sources[index];
				if( block.error != 0 )
				{
					block.start = value1;
					block.end = value2;
					block.error = 0;
				}
			}
		}
	}

	// iteratively fill in the missing values
	for( ;; )
	{
		bool stable = true;
		for( int index = 0; index < indices; ++index )
		{
			for( int target = 0; target < 256; ++target )
			{
				if( target != 255 )
				{
					SourceBlock& current = values[target].sources[index];
					SourceBlock& next = values[target + 1].sources[index];
					if( current.error > next.error + 1 )
					{
						current.start = next.start;
						current.end = next.end;
						current.error = next.error + 1;
						stable = false;
					}
				}
				if( target != 0 )
				{
					SourceBlock& current = values[target].sources[index];
					SourceBlock& previous = values[target - 1].sources[index];
					if( current.error > previous.error + 1 )
					{
						current.start = previous.start;
						current.end = previous.end;
						current.error = previous.error + 1;
						stable = false;
					}
				}
			}
		}
		if( stable )
			break;
	}

	// debug
	if (dxt)
		std::cout << "\nstatic SingleColourLookup" << " const " << name << "[] = \n{\n";
	else
		std::cout << "\nstatic SinglePaletteLookup" << indices << " const " << name << "[] = \n{\n";

	for( int i = 0;; )
	{
		std::cout << "  {{";
		for( int j = 0;; )
		{
			SourceBlock const& block = values[i].sources[j];
			if( j < colours )
				std::cout << "{"
				<<     (block.start <= 9 ? " " : "") << block.start << ","
				<<     (block.end   <= 9 ? " " : "") << block.end   << ","
				<<     (block.error <= 9 ? " " : "") << block.error << "}";
			else
				std::cout << "{0,0,0}";
			if( ++j == indices )
				break;
			std::cout << ",";
		}
		std::cout << "}}";
		if( ++i == 256 )
			break;
		std::cout << ",";
		if(!(i & 3))
			std::cout << "\n";
	}
	std::cout << "\n};\n";
}

int main()
{
	// DXT types
	GenerateData( "lookup_5_3", 5, 3, true );
	GenerateData( "lookup_6_3", 6, 3, true );
	GenerateData( "lookup_5_4", 5, 4, true );
	GenerateData( "lookup_6_4", 6, 4, true );

	// BTC types
	GenerateData( "lookup_5_4", 5, 4, false );
	GenerateData( "lookup_6_4", 6, 4, false );	// 5u1
	GenerateData( "lookup_7_4", 7, 4, false );
	GenerateData( "lookup_8_4", 8, 4, false );	// 7u1

	GenerateData( "lookup_5_8", 5, 8, false );	// 4u1
	GenerateData( "lookup_6_8", 6, 8, false );	// 6s1
	GenerateData( "lookup_7_8", 7, 8, false );

	GenerateData( "lookup_8_16", 8, 16, false );	//7u1
}
