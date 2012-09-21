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
	SourceBlock sources[2];
};

static void GenerateDataDxt( std::string const& name, int bits, int colours )
{
	TargetValue values[256];
	int count = ( 1 << bits );
	int indices = 2;

	// start + 1 intermediate
	if( colours == 3 )
		indices = 2;
	// start + 1 intermediate (codebook is symmetric)
	else
		indices = 2;

	// initialise the data
	for( int target = 0; target < 256; ++target )
		for( int index = 0; index < indices; ++index )
			values[target].sources[index].error = 255;

	// loop over all possible source points
	for( int value1 = 0; value1 < count; ++value1 )
	{
		for( int value2 = 0; value2 < count; ++value2 )
		{
			// compute the 8-bit endpoints
			int a = ( value1 << ( 8 - bits ) ) |
				( value1 >> ( 2*bits - 8 ) );
			int b = ( value2 << ( 8 - bits ) ) |
				( value2 >> ( 2*bits - 8 ) );

			// fill in the codebook with the these and intermediates
			int codes[2];
			codes[0] = a;

			// start + 1 intermediate
			if( colours == 3 )
				codes[1] = ( a + b )/2;
			// start + 1 intermediate (codebook is symmetric)
			else
				codes[1] = ( 2*a + b )/3;

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

	std::cout << "\nstatic SingleColourLookup" << " const " << name << "[] = \n{\n";

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

#include <set>

struct SSourceBlock
{
	int error;
	int index;
	int start;
	int end;

	int interleaved;

	void seedValue(int _start, int _stop) {
		start = _start;
		end   = _stop;

#if 0
		// absolute values
		interleaved  = 0;
		interleaved |= ((start >> 7) & 1) << 23;
		interleaved |= ((end   >> 7) & 1) << 22;
		interleaved |= ((start >> 6) & 1) << 21;
		interleaved |= ((end   >> 6) & 1) << 20;
		interleaved |= ((start >> 5) & 1) << 19;
		interleaved |= ((end   >> 5) & 1) << 18;
		interleaved |= ((start >> 4) & 1) << 17;
		interleaved |= ((end   >> 4) & 1) << 16;
		interleaved |= ((start >> 3) & 1) << 15;
		interleaved |= ((end   >> 3) & 1) << 14;
		interleaved |= ((start >> 2) & 1) << 13;
		interleaved |= ((end   >> 2) & 1) << 12;
		interleaved |= ((start >> 1) & 1) << 11;
		interleaved |= ((end   >> 1) & 1) << 10;
		interleaved |= ((start >> 0) & 1) <<  9;
		interleaved |= ((end   >> 0) & 1) <<  8;
#endif
		// min deviance
		interleaved = ((      std::abs(start - end)) << 8) & 0x00FFFF00;
		interleaved = ((      std::abs(start - end)) << 0) & 0x0000FFFF;

		// max deviance
		// guaranteed to just produce indices below half
		// less error-prone to the introduction of random bits
		interleaved = ((255 - std::abs(start - end)) << 8) & 0x00FFFF00;
		interleaved = ((255 - std::abs(start - end)) << 0) & 0x0000FFFF;

		interleaved = ((start << 8) | (end << 0)) & 0x0000FFFF;
	}

	void setValue(int _index, int _error) {
		index = _index;
		error = std::abs(_error);

		interleaved &= 0x00FFFF00;
		interleaved |= error << 24;
		interleaved |= index <<  0;

		// gear the index towards the higher endpoint
		// random indices
		/*
		interleaved &= 0x00FFFF00;
		interleaved |= error << 24;
		interleaved |= (start < end ? index : 255 - index);
		 */

		interleaved &= 0x0000FFFF;
		interleaved |= error << 24;
		interleaved |= index << 16;
	}

	bool operator < (const SSourceBlock &other) const {
	  return this->interleaved < other.interleaved;
	}

};

struct S {
  bool operator()(const SSourceBlock *s1, const SSourceBlock *s2) const {
    return s1->interleaved < s2->interleaved;
  }
};

struct SSourceSet
{
	std::set<SSourceBlock> values[256][16];
	int maxgap;

	void setWrong(int _gap) {
		// four is half a bin of (1 << 3)
		maxgap = _gap;
	}

	void addValue(SSourceBlock &sb, int value, int index) {
		int lo = std::max(  0 - value, -maxgap);
		int hi = std::min(255 - value,  maxgap);

		for (int e = lo; e <= hi; e++) {
			sb.setValue(index, e);

			values[value + e][index].insert(sb);
		}
	}
};

SSourceSet sets[4];

static void GenerateDataBtc( std::string const& name, int bits, int sharedbits, int colours )
{
	int count = ( 1 << bits );
	int permutations = ( 1 << (sharedbits * 2) );
	int maxgap = 1 << (8 - bits - sharedbits - 0);
	int indices = 2;

	// start + 1 intermediate (codebook is symmetric)
	if( colours == 4 )
		indices = 2;
	// start + 3 intermediates (codebook is symmetric)
	else  if( colours == 8 )
		indices = 4;
	// start + 7 intermediates (codebook is symmetric)
	else
		indices = 8;
	
	for( int cases = 0; cases < 4; ++cases )
		for( int value  = 0; value < 256; ++value  )
			for( int index  = 0; index < 16; ++index  )
				sets[cases].values[value][index].clear();

	// loop over all possible source points
	for( int cases = 0; cases < permutations; ++cases )
	{
		int casev1 = ( ( cases >> 0 ) & 1 );
		int casev2 = ( ( cases >> 1 ) & 1 );

		SSourceSet &cset = sets[cases];

		cset.setWrong(maxgap);

		for( int value1 = 0; value1 < count; ++value1 )
		{
			for( int value2 = 0; value2 < count; ++value2 )
			{
				SSourceBlock sb;

				// compute the 8-bit endpoints
				int a = ( value1 << ( 8 - bits ) ) |
					( casev1 << ( 8 - bits - sharedbits ) );
				int b = ( value2 << ( 8 - bits ) ) |
					( casev2 << ( 8 - bits - sharedbits ) );

				a |= a >> ( bits + sharedbits );
				b |= b >> ( bits + sharedbits );

				// compute the real endpoints
				int s = ( value1 << ( sharedbits ) ) |
					( casev1 << ( 0 ) );
				int e = ( value2 << ( sharedbits ) ) |
					( casev2 << ( 0 ) );

				// fill in the codebook with the these and intermediates
				sb.seedValue(s, e);

				cset.addValue(sb, a, 0);

				// start + 1 intermediate (codebook is symmetric)
				if( colours == 4 ) {
					cset.addValue(sb, ( 43*a + 21*b )/64, 1);
				//	cset.addValue(sb, ( 21*a + 43*b )/64, 2);
				}
				// start + 3 intermediates (codebook is symmetric)
				else  if( colours == 8 ) {
					cset.addValue(sb, ( 55*a +  9*b )/64, 1);
					cset.addValue(sb, ( 46*a + 18*b )/64, 2);
					cset.addValue(sb, ( 37*a + 27*b )/64, 3);
				//	cset.addValue(sb, ( 27*a + 37*b )/64, 4);
				//	cset.addValue(sb, ( 18*a + 46*b )/64, 5);
				//	cset.addValue(sb, (  9*a + 55*b )/64, 6);
				}
				// start + 7 intermediates (codebook is symmetric)
				else {
					cset.addValue(sb, ( 60*a +  4*b )/64,  1);
					cset.addValue(sb, ( 55*a +  9*b )/64,  2);
					cset.addValue(sb, ( 51*a + 13*b )/64,  3);
					cset.addValue(sb, ( 47*a + 17*b )/64,  4);
					cset.addValue(sb, ( 43*a + 21*b )/64,  5);
					cset.addValue(sb, ( 38*a + 26*b )/64,  6);
					cset.addValue(sb, ( 34*a + 30*b )/64,  7);
				//	cset.addValue(sb, ( 30*a + 34*b )/64,  8);
				//	cset.addValue(sb, ( 26*a + 38*b )/64,  9);
				//	cset.addValue(sb, ( 21*a + 43*b )/64, 10);
				//	cset.addValue(sb, ( 17*a + 47*b )/64, 11);
				//	cset.addValue(sb, ( 13*a + 51*b )/64, 12);
				//	cset.addValue(sb, (  9*a + 55*b )/64, 13);
				//	cset.addValue(sb, (  4*a + 60*b )/64, 14);
				}

			//	cset.addValue(sb, b, indices - 1);
			}
		}
	}
	
	std::cout << "\n/* 256 * " << indices << " * 3 * " << permutations << " = " << (256 * indices * 3 * permutations) << " */";
	if ( permutations == 1 )
	  std::cout << "\nstatic SinglePaletteLookup" << indices << " const " << name << "[256] = \n{\n";
	else
	  std::cout << "\nstatic SinglePaletteLookup" << indices << " const " << name << "[4][256] = \n{\n";
	
	for( int k = 0;; )
	{
		if( permutations > 1 )
			std::cout << "{\n";

		for( int i = 0;; )
		{
			std::cout << "  {{";
		
			for( int j = 0;; )
			{
				// lowest error, lowest values, lowest index
				SSourceSet const& cset = sets[k];
				std::set<SSourceBlock, struct S>::iterator val = cset.values[i][j].begin();

				if( j < colours )
					std::cout << "{"
					  << (val->start <= 99 ? " " : "") << (val->start <= 9 ? " " : "") << val->start << ","
					  << (val->end   <= 99 ? " " : "") << (val->end   <= 9 ? " " : "") << val->end   << ","
					  <<                                                                  val->error << "}";

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
			
		if( permutations > 1 )
		  std::cout << "\n}";

		if( ++k == permutations )
			break;
		std::cout << ",";
	}

	std::cout << "\n};\n";
}

int main()
{
	// DXT types
	GenerateDataDxt( "lookup_5_3",    5,     3 );
	GenerateDataDxt( "lookup_6_3",    6,     3 );
	GenerateDataDxt( "lookup_5_4",    5,     4 );
	GenerateDataDxt( "lookup_6_4",    6,     4 );

	// BTC types
	GenerateDataBtc( "lookup_5_4",    5, 0,  4 );
	GenerateDataBtc( "lookup_6_4",    6, 0,  4 );
	GenerateDataBtc( "lookup_7_4",    7, 0,  4 );
	GenerateDataBtc( "lookup_8_4",    8, 0,  4 );
	GenerateDataBtc( "lookup_5_8",    5, 0,  8 );
	GenerateDataBtc( "lookup_6_8",    6, 0,  8 );
	GenerateDataBtc( "lookup_7_8",    7, 0,  8 );	// needed
	GenerateDataBtc( "lookup_8_16",   8, 0, 16 );	// needed

	GenerateDataBtc( "lookup_5u1_4",  5, 1,  4 );
	GenerateDataBtc( "lookup_7u1_4",  7, 1,  4 );
	GenerateDataBtc( "lookup_4u1_8",  4, 1,  8 );
	GenerateDataBtc( "lookup_6s1_8",  6, 1,  8 );
	GenerateDataBtc( "lookup_7u1_16", 7, 1, 16 );

	//4,0,  1,0,  3,0	4u1	8
	//5,0,  0,0,  2,0       5	4
	//5,5,  1,0,  2,0       5u1	4
	//5,6,  0,0,  2,3       5,6	4,8
	//6,0,  0,1,  3,0       6s1	8
	//7,0,  1,0,  2,0       7u1	4
	//7,7,  1,0,  4,0       7u1	16
	//7,8,  0,0,  2,2       7,8	4,4
}
