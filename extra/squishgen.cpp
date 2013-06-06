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

#pragma warning(disable: 4100)
#pragma warning(disable: 4189)
#pragma warning(disable: 4127)

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

static void GenerateDataDxt(std::string const& name, int bits, int colours, bool direct = false)
{
  TargetValue values[256];
  int count = ( 1 << bits );
  int indices = 2;

  // start + 1 intermediate
  if (colours == 3)
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
      if (colours == 3)
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

  if (!direct) {
    std::cout << "\nstatic ColourSingleLookup" << " const sc_" << name << "[] = \n{\n";

    for(int i = 0;;) {
      std::cout << "  {{";

      for (int j = 0;;) {
        SourceBlock const& block = values[i].sources[j];
        if( j < colours )
	  std::cout << "{"
	  <<     (block.start <= 9 ? " " : "") << block.start << ","
	  <<     (block.end   <= 9 ? " " : "") << block.end   << ","
	  <<     (block.error <= 9 ? " " : "") << block.error << "}";
        else
	  std::cout << "{0,0,0}";

        if (++j == indices)
	  break;
        std::cout << ",";
      }

      std::cout << "}}";
      if (++i == 256)
        break;

      std::cout << ",";
      if (!(i & 3))
        std::cout << "\n";
    }

    std::cout << "\n};\n";
  }
  else {
    std::cout << "\nstatic ColourSingleLookup" << " const sc_" << name << "[] = \n{\n";

    for(int i = 0;;) {
      std::cout << "  ";

      if (values[i].sources[0].error <
      	  values[i].sources[1].error) {
      	fprintf(stderr, "assumption (error0 >= error1) incorrect!\n");
      }

      int j = 1; {
        SourceBlock const& block = values[i].sources[j];

        int a = ( block.start << ( 8 - bits ) ) |
	        ( block.start >> ( 2*bits - 8 ) );
        int b = ( block.end   << ( 8 - bits ) ) |
	        ( block.end   >> ( 2*bits - 8 ) );

        if( j < colours )
	  std::cout << "{"
	  <<     (a <= 99 ? " " : "") << (a <= 9 ? " " : "") << a << ","
	  <<     (b <= 99 ? " " : "") << (b <= 9 ? " " : "") << b   << "}";
        else
	  std::cout << "{0,0}";
      }

      std::cout << "";
      if (++i == 256)
        break;

      std::cout << ",";
      if (!(i & 7))
        std::cout << "\n";
    }

    std::cout << "\n};\n";
  }
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

    interleaved = ((start << 8) + (end << 0)) & 0x0000FFFF;
    interleaved = ((start << 8) + (std::abs(start - end) << 0)) & 0x0000FFFF;
  }

  void setValue(int _index, int _error) {
    index = _index;
    error = std::abs(_error);

    /*
    interleaved &= 0x00FFFF00;
    interleaved |= error << 24;
    interleaved |= index <<  0;
    */

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

static void GenerateDataBtc( std::string const& name, int bits, int sharedbits, int colours, int exactcases = 0, bool direct = false )
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
    if ( ( exactcases == 2 ) && ( cases != 0 ) && ( cases != 3 ) )
      continue;

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
//	sb.seedValue(s, e);
	sb.seedValue(a, b);

	cset.addValue(sb, a, 0);

	// start + 1 intermediate (codebook is symmetric)
	if( colours == 4 ) {
	  cset.addValue(sb, ( 43*a + 21*b )/64, 1);
//	  cset.addValue(sb, ( 21*a + 43*b )/64, 2);
	}
	// start + 3 intermediates (codebook is symmetric)
	else  if( colours == 8 ) {
	  cset.addValue(sb, ( 55*a +  9*b )/64, 1);
	  cset.addValue(sb, ( 46*a + 18*b )/64, 2);
	  cset.addValue(sb, ( 37*a + 27*b )/64, 3);
//	  cset.addValue(sb, ( 27*a + 37*b )/64, 4);
//	  cset.addValue(sb, ( 18*a + 46*b )/64, 5);
//	  cset.addValue(sb, (  9*a + 55*b )/64, 6);
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
//	  cset.addValue(sb, ( 30*a + 34*b )/64,  8);
//	  cset.addValue(sb, ( 26*a + 38*b )/64,  9);
//	  cset.addValue(sb, ( 21*a + 43*b )/64, 10);
//	  cset.addValue(sb, ( 17*a + 47*b )/64, 11);
//	  cset.addValue(sb, ( 13*a + 51*b )/64, 12);
//	  cset.addValue(sb, (  9*a + 55*b )/64, 13);
//	  cset.addValue(sb, (  4*a + 60*b )/64, 14);
	}

//	cset.addValue(sb, b, indices - 1);
      }
    }
  }

  if (!direct) {
    std::cout << "\n/* 256 * " << indices << " * 3 * " << permutations << " = " << (256 * indices * 3 * permutations) << " */";
    if ( permutations == 1 )
      std::cout << "\nstatic PaletteSingleLookup" << indices << " const " << name << "[256] = \n{\n";
    else if ( exactcases == 2 )
      std::cout << "\nstatic PaletteSingleLookup" << indices << " const " << name << "[2][256] = \n{\n";
    else
      std::cout << "\nstatic PaletteSingleLookup" << indices << " const " << name << "[4][256] = \n{\n";

    for( int k = 0;; )
    {
      if ( ( exactcases == 2 ) && ( k != 0 ) && ( k != 3 ) ) {
        if( ++k == permutations )
	  break;
        continue;
      }

      if( permutations > 1 )
        std::cout << "{\n";

      for( int i = 0;; )
      {
        std::cout << "  {{";

        int min_error = 0xFFFF;
        int min_index = 0;

        for( int j = 0;; )
        {
	  // lowest error, lowest values, lowest index
	  SSourceSet const& cset = sets[k];
	  std::set<SSourceBlock, struct S>::iterator val = cset.values[i][j].begin();

	  if( j < colours )
	    if( min_error > val->error ) {
	      min_error = val->error;
	      min_index = j;
	    }
	    else if( min_error == val->error ) {
	      min_index = !min_index ? j : min_index;
	    }
	  if( ++j == indices )
	    break;
        }

        std::cout << "/*" << min_index << "*/";

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
  else {
    std::cout << "\n/* 256 * 2 * " << permutations << " = " << (256 * 2 * permutations) << " */";
    if ( permutations == 1 )
      std::cout << "\nstatic PaletteSingleLookup" << indices << " const " << name << "[256] = \n{\n";
    else if ( exactcases == 2 )
      std::cout << "\nstatic PaletteSingleLookup" << indices << " const " << name << "[2][256] = \n{\n";
    else
      std::cout << "\nstatic PaletteSingleLookup" << indices << " const " << name << "[4][256] = \n{\n";

    for( int k = 0;; )
    {
      if ( ( exactcases == 2 ) && ( k != 0 ) && ( k != 3 ) ) {
        if( ++k == permutations )
	  break;
        continue;
      }

      if( permutations > 1 )
        std::cout << "{\n";

      for( int i = 0;; )
      {
        std::cout << "  ";

        int j = 1; {
	  // lowest error, lowest values, lowest index
	  SSourceSet const& cset = sets[k];
	  std::set<SSourceBlock, struct S>::iterator val = cset.values[i][j].begin();

	  if( j < colours )
	    std::cout << "{"
	    << (val->start <= 99 ? " " : "") << (val->start <= 9 ? " " : "") << val->start << ","
	    << (val->end   <= 99 ? " " : "") << (val->end   <= 9 ? " " : "") << val->end   << "}";

	  else
	    std::cout << "{0,0}";
        }

        std::cout << "";
        if( ++i == 256 )
	  break;

        std::cout << ",";
        if(!(i & 7))
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
}

struct ASourceBlock
{
  int unique;
  int errorS, errorE;
  int start;
  int end;
  unsigned __int64 interleaved;

  void seedValue(int _unique, int _start, int _stop) {
    unique = _unique;
    start = _start;
    end   = _stop;

    interleaved  =    ((start << 8) + (end << 0))   & 0x0000FFFF;
    interleaved |= (std::max(errorS, errorE) << 24) & 0xFF000000;
    interleaved |= (std::min(errorS, errorE) << 16) & 0x00FF0000;
  }

  void setValue(int _delta, int _errorS, int _errorE) {
    errorS = std::abs(_errorS);
    errorE = std::abs(_errorE);

    interleaved &= 0xFFFFFFFF;
    interleaved |= ((unsigned __int64)((errorS * errorS) + (errorE * errorE))) << 32; // 17bits
    interleaved |= ((unsigned __int64)_delta) << (32 + 17);
  }

  bool operator < (const ASourceBlock &other) const {
    return this->interleaved < other.interleaved;
  }
};

struct A {
  bool operator()(const ASourceBlock *s1, const ASourceBlock *s2) const {
    return s1->interleaved < s2->interleaved;
  }
};

struct ASourceSet
{
  std::set<ASourceBlock> values[1 << 16];

  void addCase(ASourceBlock &sb, int delta, int unique, int value1, int value2, int start, int end, int errorS, int errorE) {
    std::set<ASourceBlock> &cset = values[(value2 << 8) | (value1 << 0)];

    {
      sb.seedValue(unique, start, end);
      sb.setValue(delta, errorS, errorE);

      std::set<ASourceBlock, struct A>::iterator val = cset.begin();

      if (val != cset.end())
	if ((val->interleaved >> 16) < (sb.interleaved >> 16))
	  return;

      cset.insert(sb);
    }
  }
};

ASourceSet asets[4];

static void GenerateDataRtg( std::string const& name, int bits, int sharedbits, int colours, int exactcases = 0 )
{
  int count = ( 1 << bits );
  int permutations = ( 1 << (sharedbits * 2) );
  int maxgap = 1 << (8 - bits - sharedbits - 0);
  int indices = 2;

  // start + 1 intermediate (codebook is symmetric)
  if (colours == 4)
    indices = 4;
  // start + 3 intermediates (codebook is symmetric)
  else if (colours == 8)
    indices = 8;
  // start + 7 intermediates (codebook is symmetric)
  else
    indices = 16;

  for (int cases = 0; cases < 4; ++cases)
    for (int value  = 0; value < 65536; ++value)
      asets[cases].values[value].clear();

  // loop over all possible source points
  for (int cases = 0; cases < permutations; ++cases) {
    if ((exactcases == 2) && (cases != 0) && (cases != 3))
      continue;

    int casev1 = ((cases >> 0) & 1);
    int casev2 = ((cases >> 1) & 1);

    ASourceSet &cset = asets[cases];

    for (int value1 = 0; value1 < count; ++value1) {
      for (int value2 = 0; value2 < count; ++value2) {
	ASourceBlock sb;

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
	int codes[16];

	codes[0] = a;

	// start + 1 intermediate (codebook is symmetric)
	if( colours == 4 ) {
	  codes[1] = (43*a + 21*b )/64;
	  codes[2] = (21*a + 43*b )/64;
	}
	// start + 3 intermediates (codebook is symmetric)
	else  if( colours == 8 ) {
	  codes[1] = ( 55*a +  9*b )/64;
	  codes[2] = ( 46*a + 18*b )/64;
	  codes[3] = ( 37*a + 27*b )/64;
	  codes[4] = ( 27*a + 37*b )/64;
	  codes[5] = ( 18*a + 46*b )/64;
	  codes[6] = (  9*a + 55*b )/64;
	}
	// start + 7 intermediates (codebook is symmetric)
	else {
	  codes[ 1] = ( 60*a +  4*b )/64;
	  codes[ 2] = ( 55*a +  9*b )/64;
	  codes[ 3] = ( 51*a + 13*b )/64;
	  codes[ 4] = ( 47*a + 17*b )/64;
	  codes[ 5] = ( 43*a + 21*b )/64;
	  codes[ 6] = ( 38*a + 26*b )/64;
	  codes[ 7] = ( 34*a + 30*b )/64;
	  codes[ 8] = ( 30*a + 34*b )/64;
	  codes[ 9] = ( 26*a + 38*b )/64;
	  codes[10] = ( 21*a + 43*b )/64;
	  codes[11] = ( 17*a + 47*b )/64;
	  codes[12] = ( 13*a + 51*b )/64;
	  codes[13] = (  9*a + 55*b )/64;
	  codes[14] = (  4*a + 60*b )/64;
	}

	codes[colours-1] = b;

	// count unique indices
	int unique = colours;
	for (int c1 = 0; c1 < colours; ++c1) {
	  for (int c2 = c1 + 1; c2 < colours; ++c2) {
	    if (codes[c1] == codes[c2]) {
	      unique--;
	      break;
	    }
	  }
	}

	/* we can parameterize a random distribution of points
	 * into three parameters:
	 *
	 * - range: start/end
	 * - median: deviation from the center
	 * - square magnitudes
	 *
	 * there are 16^255 possible incoming pixel combinations
	 */

	for (int valueS = 0; valueS <= 255; ++valueS) {
	  for (int valueE = valueS; valueE <= 255; ++valueE) {
	    int values[2];
	    int error[2];
	    int idxs[2];

	    values[0] = valueS;
	    values[1] = valueE;

	    int spread = std::min(colours, valueE - valueS);

	    // eliminate codebooks with less entries than the spread
//	    if (unique < spread)
//	      continue;

	    // find the closest codes
	    for (int i = 0; i < 2; ++i) {
	      int dist = 0x7FFFFFFF;
	      int idx = 0;

	      for (int j = 0; j < colours; ++j) {
		int d = std::abs(values[i] - codes[j]);
		if (d < dist) {
		  dist = d;
		  idx = j;
		}
	      }

	      idxs[i] = idx;
	      error[i] = dist;
	    }

	    // eliminate codebooks with less entries between the endpoints than the spread
	    int delta = 0;//-std::max(std::abs(idxs[0] - idxs[1]) - (spread - 2), 0);
//	    if (delta > 0)
//	      continue;

	    // requirement
	    int sp = std::max(spread - 2, 0);
	    // offer
	    int id = std::abs(idxs[0] - idxs[1]);

	    // less offered than requirement
	    if (id < sp)
	      delta = sp - id;

	    // penalize in the priorization
	    cset.addCase(sb, delta, unique, valueS, valueE, s, e, error[0], error[1]);
	  }
	}

	fprintf(stderr, "%f%%\r", 100.0f * ((value1 * count) + value2) / (count * count));
      }
    }
  }

  if ( permutations == 1 )
    std::cout << "\nstatic SingleChannelLookup" << indices << " const " << name << "[32768] = \n{\n";
  else if ( exactcases == 2 )
    std::cout << "\nstatic SingleChannelLookup" << indices << " const " << name << "[2][32768] = \n{\n";
  else
    std::cout << "\nstatic SingleChannelLookup" << indices << " const " << name << "[4][32768] = \n{\n";

  // [upper = 0] = 0 removed
  // [upper = 1] = 255 removed
  // [upper = 2] = 255 + 254 removed

  for( int k = 0;; )
  {
    if ( ( exactcases == 2 ) && ( k != 0 ) && ( k != 3 ) ) {
      if( ++k == permutations )
	break;
      continue;
    }

    if( permutations > 1 )
      std::cout << "{\n";

    int removed[256] = {0};

    for( int i = 0;; )
    {
      std::cout << "  /* max spread: [0," << i << "] */\n";
      std::cout << "  ";

      for( int j = 0;; )
      {
	if (i >= j) {
	  // lowest error, lowest values, lowest index
	  ASourceSet const& cset = asets[k];
	  std::set<ASourceBlock, struct A>::iterator val = cset.values[(i << 8) + (j)].begin();

	  if (val != cset.values[(i << 8) + (j)].end())
	    std::cout << "{"
	    << (val->start <= 99 ? " " : "") << (val->start <= 9 ? " " : "") << val->start  << ","
	    << (val->end   <= 99 ? " " : "") << (val->end   <= 9 ? " " : "") << val->end    << ","
	    <<                                                                  val->errorS << ","
	    <<                                                                  val->errorE << "/*"
	    <<                                                                  val->unique << "*/}";
	  else
	    std::cout << "{/*failed*/}";

	  std::cout << ",";
	  if((j & 15) == 15)
	    std::cout << "\n  ";
	}
	else
	  removed[i]++;

	if( ++j == 256 )
	  break;
      }

      std::cout << "\n";
      std::cout << "  /* removed: " << removed[i] << " */\n";

      if( ++i == 256 )
	break;

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
  if (0) {
    // DXT types
    GenerateDataDxt( "lookup_5_3",    5,     3 );
    GenerateDataDxt( "lookup_6_3",    6,     3 );
    GenerateDataDxt( "lookup_5_4",    5,     4 );
    GenerateDataDxt( "lookup_6_4",    6,     4 );
  }

  if (0) {
    // DXT types (with Ryg observation: error0 >= error1 always)
    GenerateDataDxt( "lookup_5_3",    5,     3 , true);
    GenerateDataDxt( "lookup_6_3",    6,     3 , true);
    GenerateDataDxt( "lookup_5_4",    5,     4 , true);
    GenerateDataDxt( "lookup_6_4",    6,     4 , true);
  }

  if (0) {
    // BTC types
    GenerateDataBtc( "lookup_5_4",    5, 0,  4, 0 );
    GenerateDataBtc( "lookup_6_4",    6, 0,  4, 0 );
    GenerateDataBtc( "lookup_7_4",    7, 0,  4, 0 );
    GenerateDataBtc( "lookup_8_4",    8, 0,  4, 0 );
    GenerateDataBtc( "lookup_5_8",    5, 0,  8, 0 );
    GenerateDataBtc( "lookup_6_8",    6, 0,  8, 0 );
    GenerateDataBtc( "lookup_7_8",    7, 0,  8, 0 );	// needed
    GenerateDataBtc( "lookup_8_16",   8, 0, 16, 0 );	// needed

    GenerateDataBtc( "lookup_5u1_4",  5, 1,  4, 4 );
    GenerateDataBtc( "lookup_7u1_4",  7, 1,  4, 4 );
    GenerateDataBtc( "lookup_4u1_8",  4, 1,  8, 4 );
    GenerateDataBtc( "lookup_6s1_8",  6, 1,  8, 2 );
    GenerateDataBtc( "lookup_7u1_16", 7, 1, 16, 4 );

    //4,0,  1,0,  3,0	    4u1	8
    //5,0,  0,0,  2,0       5	4
    //5,5,  1,0,  2,0       5u1	4
    //5,6,  0,0,  2,3       5,6	4,8
    //6,0,  0,1,  3,0       6s1	8
    //7,0,  1,0,  2,0       7u1	4
    //7,7,  1,0,  4,0       7u1	16
    //7,8,  0,0,  2,2       7,8	4,4
  }

  if (1) {
    // BTC types
    GenerateDataBtc( "lookup_5_4",    5, 0,  4, 0, true );
    GenerateDataBtc( "lookup_6_4",    6, 0,  4, 0, true );
    GenerateDataBtc( "lookup_7_4",    7, 0,  4, 0, true );
    GenerateDataBtc( "lookup_8_4",    8, 0,  4, 0, true );
    GenerateDataBtc( "lookup_5_8",    5, 0,  8, 0, true );
    GenerateDataBtc( "lookup_6_8",    6, 0,  8, 0, true );
    GenerateDataBtc( "lookup_7_8",    7, 0,  8, 0, true );	// needed
    GenerateDataBtc( "lookup_8_16",   8, 0, 16, 0, true );	// needed

    GenerateDataBtc( "lookup_5u1_4",  5, 1,  4, 4, true );
    GenerateDataBtc( "lookup_7u1_4",  7, 1,  4, 4, true );
    GenerateDataBtc( "lookup_4u1_8",  4, 1,  8, 4, true );
    GenerateDataBtc( "lookup_6s1_8",  6, 1,  8, 2, true );
    GenerateDataBtc( "lookup_7u1_16", 7, 1, 16, 4, true );
  }

  if (0) {
    // RTG types
    GenerateDataRtg( "lookup_5g_4",    5, 0,  4 );
//  GenerateDataRtg( "lookup_6g_4",    6, 0,  4 );
//  GenerateDataRtg( "lookup_7g_4",    6, 0,  4 );
  }
}