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

/*! @file

    @brief	Example program that converts between the PNG and DXT/BTC formats.

    This program requires libpng for PNG input and output, and is designed
    to show how to prepare data for the squish library when it is not simply
    a contiguous block of memory.
*/

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>
#include <assert.h>
#include <squish.h>
#include <png.h>
#include <intrin.h>

#if (_MSC_VER >= 1700) && defined(NDEBUG)
#include <chrono>
using namespace std::chrono;
#include <windows.h>
#include <concrt.h>
#include <ppl.h>
#undef min
#undef max
#endif

#ifdef _MSC_VER
#pragma warning( disable: 4511 4512 )
#endif // def _MSC_VER

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

  void PNGRead(png_structp png_ptr, png_bytep data, png_size_t length)
  {
    FILE *ar = (FILE *)png_ptr->io_ptr;

    fread(data, 1, length, ar);
  }

  void PNGWrite(png_structp png_ptr, png_bytep data, png_size_t length)
  {
    FILE *ar = (FILE *)png_ptr->io_ptr;

    fwrite(data, 1, length, ar);
  }

  void PNGFlush(png_structp png_ptr)
  {
    FILE *ar = (FILE *)png_ptr->io_ptr;

    fflush(ar);
  }

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

using namespace squish;

//! Simple exception class.
class Error : public std::exception
{
public:
  Error( std::string const& excuse ) : m_excuse( excuse ) {}
  ~Error() throw() {}

  virtual char const* what() const throw() { return m_excuse.c_str(); }

private:
  std::string m_excuse;
};

//! Base class to make derived classes non-copyable
class NonCopyable
{
public:
  NonCopyable() {}

private:
  NonCopyable( NonCopyable const& );
  NonCopyable& operator=( NonCopyable const& );
};

//! Memory object.
class Mem : NonCopyable
{
public:
  explicit Mem( int size ) : m_p( new unsigned char[size] ) {}
  ~Mem() { delete[] m_p; }

  unsigned char* Get() const { return m_p; }

private:
  unsigned char* m_p;
};

//! File object.
class File : NonCopyable
{
public:
  explicit File( FILE* fp ) : m_fp( fp ) {}
  ~File() { if( m_fp ) fclose( m_fp ); }

  bool IsValid() const { return m_fp != 0; }
  FILE* Get() const { return m_fp; }

private:
  FILE* m_fp;
};

//! PNG read object.
class PngReadStruct : NonCopyable
{
public:
  PngReadStruct()
    : m_png(0)
    , m_info(0)
    , m_end(0)
  {
    m_png = png_create_read_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0 );
    if (!m_png)
      throw Error( "failed to create png read struct" );

    m_info = png_create_info_struct(m_png);
    m_end = png_create_info_struct(m_png);
    if (!m_info || !m_end) {
      png_infopp info = m_info ? &m_info : 0;
      png_infopp end = m_end ? &m_end : 0;
      png_destroy_read_struct( &m_png, info, end );
      throw Error( "failed to create png info structs" );
    }
  }

  ~PngReadStruct()
  {
    png_destroy_read_struct( &m_png, &m_info, &m_end );
  }

  png_structp GetPng() const { return m_png; }
  png_infop GetInfo() const { return m_info; }

private:
  png_structp m_png;
  png_infop m_info, m_end;
};

//! PNG write object.
class PngWriteStruct : NonCopyable
{
public:
  PngWriteStruct()
    : m_png( 0 ),
    m_info( 0 )
  {
    m_png = png_create_write_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0 );
    if( !m_png )
      throw Error( "failed to create png read struct" );

    m_info = png_create_info_struct( m_png );
    if( !m_info )
    {
      png_infopp info = m_info ? &m_info : 0;
      png_destroy_write_struct( &m_png, info );
      throw Error( "failed to create png info structs" );
    }
  }

  ~PngWriteStruct()
  {
    png_destroy_write_struct( &m_png, &m_info );
  }

  png_structp GetPng() const { return m_png; }
  png_infop GetInfo() const { return m_info; }

private:
  png_structp m_png;
  png_infop m_info;
};

//! PNG rows object.
class PngRows : NonCopyable
{
public:
  PngRows( int width, int height, int stride ) : m_width( width ), m_height( height )
  {
    m_rows = (png_bytep *)malloc(m_height * sizeof(png_bytep));
    for (int i = 0; i < m_height; ++i)
      m_rows[i] = (png_bytep)malloc(m_width * stride);
  }

  ~PngRows()
  {
    for( int i = 0; i < m_height; ++i )
      free( m_rows[i] );
    free( m_rows );
  }

  png_bytep* Get() const {
    return m_rows;
  }
  png_bytep const GetRow( int row ) const {
    return (png_bytep)m_rows[row];
  }

  void Read( PngReadStruct *m_png ) {
    png_read_image( m_png->GetPng(), m_rows );
  }
  void Write( PngWriteStruct *m_png ) {
    png_write_image( m_png->GetPng(), m_rows );
  }

private:
  png_bytep* m_rows;
  int m_width, m_height;
};

class PngImage
{
public:
  explicit PngImage( std::string const& fileName );

  int GetWidth() const { return m_width; }
  int GetHeight() const { return m_height; }
  int GetStride() const { return m_stride; }
  bool IsColour() const { return m_colour; }
  bool IsAlpha() const { return m_alpha; }

#if (PNG_LIBPNG_VER_MINOR >= 3)
  png_bytep* Get() const { return (png_bytep *)m_rows->Get( ); }
  PngRows* GetRows( ) const { return m_rows; }
  u8 const* GetRow( int row ) const { return ( u8* )m_rows[row]; }
#else
  png_bytep* Get() const { return (png_bytep *)m_rows->Get( ); }
  PngRows* GetRows( ) const { return m_rows; }
  u8 const* GetRow( int row ) const { return ( u8* )m_rows->GetRow( row ); }
#endif

private:
  PngReadStruct m_png;

  int m_width;
  int m_height;
  int m_stride;
  bool m_colour;
  bool m_alpha;

#if (PNG_LIBPNG_VER_MINOR >= 3)
  png_bytep* m_rows;
#else
  PngRows* m_rows;
#endif
};

PngImage::PngImage( std::string const& fileName )
{
  // open the source file
  File file( fopen( fileName.c_str(), "rb" ) );
  if( !file.IsValid() )
  {
    std::ostringstream oss;
    oss << "failed to open \"" << fileName << "\" for reading";
    throw Error( oss.str() );
  }

  // check the signature bytes
  png_byte header[8];
  fread( header, 1, 8, file.Get() );
  if( png_sig_cmp( header, 0, 8 ) )
  {
    std::ostringstream oss;
    oss << "\"" << fileName << "\" does not look like a png file";
    throw Error( oss.str() );
  }

  // read the image into memory
  png_set_read_fn( m_png.GetPng(), (void *)file.Get(), PNGRead);
  //	png_init_io( m_png.GetPng(), file.Get() );
  png_set_sig_bytes( m_png.GetPng(), 8 );
#if (PNG_LIBPNG_VER_MINOR >= 3)
  png_read_png( m_png.GetPng(), m_png.GetInfo(), PNG_TRANSFORM_EXPAND, 0 );
#else
  png_read_info( m_png.GetPng(), m_png.GetInfo() );
#endif

  // get the image info
  png_uint_32 width;
  png_uint_32 height;
  int bitDepth;
  int colourType;

  png_get_IHDR( m_png.GetPng(), m_png.GetInfo(), &width, &height, &bitDepth, &colourType, 0, 0, 0 );

  // check the image is 8 bit
  if( bitDepth != 8 )
  {
    std::ostringstream oss;
    oss << "cannot process " << bitDepth << "-bit image (bit depth must be 8)";
    throw Error( oss.str() );
  }

  // save the info
  m_width = width;
  m_height = height;
  m_colour = ( ( colourType & PNG_COLOR_MASK_COLOR ) != 0 );
  m_alpha = ( ( colourType & PNG_COLOR_MASK_ALPHA ) != 0 );
  m_stride = ( ( m_colour ? 3 : 1 ) + ( m_alpha ? 1 : 0 ) ) * sizeof( u8 );

  // get the image rows
#if (PNG_LIBPNG_VER_MINOR >= 3)
  m_rows = png_get_rows( m_png.GetPng(), m_png.GetInfo() );
  if( !m_rows )
    throw Error( "failed to get image rows" );
#else
  m_rows = new PngRows( m_width, m_height, m_stride );
  m_rows->Read( &m_png );

  png_read_end( m_png.GetPng(), m_png.GetInfo() );
#endif
}

static void Compress(std::string const& sourceFileName, std::string const& targetFileName, int paint, int mapping, int flags)
{
  // load the source image
  PngImage sourceImage(sourceFileName);

  // get the image info
  int width = sourceImage.GetWidth();
  int height = sourceImage.GetHeight();
  int stride = sourceImage.GetStride();
  bool colour = sourceImage.IsColour();
  bool alpha = sourceImage.IsAlpha();

  /* kill alpha-channel */
  if (!alpha && (((flags & kBtcp) <= kBtc3) || ((flags & kBtcp) >= kBtc7)))
    flags = (flags | ((flags & kBtcp) > kBtc1 ? kExcludeAlphaFromPalette : 0)) & (~kWeightColourByAlpha) & (~kAlphaIterativeFit);
  if ((mapping == -1) && ((flags & kBtcp) == kBtc4))
    mapping = 0;
  if ((mapping == -1) && ((flags & kBtcp) == kBtc5))
    mapping = 3;

  // check the image dimensions
  if ((width % 4) != 0 || (height % 4) != 0) {
    std::ostringstream oss;
    oss << "cannot compress " << width << "x" << height
      << "image (dimensions must be multiples of 4)";

    throw Error(oss.str());
  }

  // create the target data
  int bytesPerBlock = 16;
  int targetDataSize = bytesPerBlock * width * height / 16;
  Mem targetData(targetDataSize);
  
  struct sqio s = GetSquishIO(width, height, sqio::dtp::DT_U8, flags);

  // loop over blocks and compress them
  clock_t start = std::clock();
  unsigned char* _targetBlock = targetData.Get();
  
#if (_MSC_VER >= 1700) && defined(NDEBUG)
  concurrency::parallel_for(0, height, 4, [&](int y) {
#else
  for (int y = 0; y < height; y += 4) {
#endif
    // process a row of blocks
    for (int x = 0; x < width; x += 4) {
      unsigned char* targetBlock = _targetBlock + ((y / 4) * (width / 4) + (x / 4)) * bytesPerBlock;

      // get the block data
      u8 sourceRgba[16 * 4];

      for (int py = 0, i = 0; py < 4; ++py) {
	u8 const* row = sourceImage.GetRow(y + py) + x * stride;

	for (int px = 0; px < 4; ++px, ++i) {
	  // get the pixel colour
	  if (colour) {
	    for (int j = 0; j < 3; ++j)
	      sourceRgba[4 * i + j] = *row++;
	  }
	  else {
	    for (int j = 0; j < 3; ++j)
	      sourceRgba[4 * i + j] = *row; ++row;
	  }

	  // skip alpha for now
	  if (alpha)
	    sourceRgba[4 * i + 3] = *row++;
	  else
	    sourceRgba[4 * i + 3] = 255;

	  if ((flags & kBtcp) == kBtc5) {
	    // duplicate alpha into g (required!)
	    // duplicate r into b (for symmetry)
	    if ((flags & kColourMetrics) == kColourMetricUnit)
	      ;
	    else {
	      sourceRgba[4 * i + 4 - mapping] = sourceRgba[4 * i + mapping];
	      sourceRgba[4 * i + 2          ] = sourceRgba[4 * i + 0      ];
	    }
	  }
	  else if ((flags & kBtcp) == kBtc4) {
	    sourceRgba[4 * i + 0] = sourceRgba[4 * i + mapping];
	    sourceRgba[4 * i + 1] = sourceRgba[4 * i + mapping];
	    sourceRgba[4 * i + 2] = sourceRgba[4 * i + mapping];
	  }
	}
      }

      if (flags & kSignedExternal) {
	for (int i = 0; i < 16 * 4; ++i)
	  sourceRgba[i] = sourceRgba[i] - 0x80 + (sourceRgba[i] < 0x80);
      }

      // compress this block
      s.encoder(sourceRgba, -1, targetBlock, flags);

#if (defined(VERIFY_QUANTIZER) || defined(VERIFY_ENCODER))
      // write the data into the target rows
      for (int py = 0, i = 0; py < 4; ++py) {
	u8 *row = (u8 *)sourceImage.Get()[y + py] + x * stride;
	for (int px = 0; px < 4; ++px, ++i) {
	  // get the pixel colour
	  if (colour) {
	    for (int j = 0; j < 3; ++j)
	      *row++ = sourceRgba[4 * i + j];
	  }

	  // skip alpha for now
	  if (alpha)
	    *row++ = sourceRgba[4 * i + 3];
	}
      }
#endif

      if (((flags & kBtcp) == kBtc7) && paint) {
	// draw the mode
	if (paint == 2 ) {
	  unsigned long btcvalue = *((int *)targetBlock);
	  unsigned long btcmode; _BitScanForward(&btcmode, btcvalue); btcmode += 1;

	  fprintf(stderr, "%1u", btcmode);
	}
	// draw the pattern
	else if (paint == 1) {
	  char *patchar = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzµß";
	  unsigned long btcvalue = *((int *)targetBlock);
	  unsigned long btcmode; _BitScanForward(&btcmode, btcvalue); btcmode += 1;
	  unsigned long btcpat = (btcmode == 5 || btcmode == 6 || btcmode == 7 ? 0 : (btcmode == 1 ? (1 << 4) - 1 : (1 << 6) - 1));

	  fprintf(stderr, "%c", patchar[(btcvalue >> btcmode) & btcpat]);
	}
      }

      // advance
      targetBlock += bytesPerBlock;
    }

    // draw the mode
    if (((flags & kBtcp) == kBtc7) && paint) {
      fprintf(stderr, "\n");
    }
  }
#if (_MSC_VER >= 1700) && defined(NDEBUG)
  );
#endif

  clock_t end = std::clock();
  double duration = (double)(end - start) / CLOCKS_PER_SEC;
  std::cout << "time taken: " << duration << " seconds" << std::endl;
  std::cout << "compression speed: " << ((double) (width * height) / (duration * 1024 * 1024)) << " MPixel/s" << std::endl;

  // open the target file
  File targetFile(fopen(targetFileName.c_str(), "wb"));
  if (!targetFile.IsValid()) {
    std::ostringstream oss;
    oss << "failed to open \"" << sourceFileName << "\" for writing";

    throw Error(oss.str());
  }

  // write the header
  fwrite(&width,  sizeof(int), 1, targetFile.Get());
  fwrite(&height, sizeof(int), 1, targetFile.Get());

  // write the data
  fwrite(targetData.Get(), 1, targetDataSize, targetFile.Get());

#if (defined(VERIFY_QUANTIZER) || defined(VERIFY_ENCODER))
  {
    // create the target PNG
    PngWriteStruct targetPng;

    // set up the image
    png_set_IHDR(
      targetPng.GetPng(), targetPng.GetInfo(), width, height,
      8*sizeof(u8), colour && alpha ? PNG_COLOR_TYPE_RGBA : (colour ? PNG_COLOR_TYPE_RGB : PNG_COLOR_TYPE_GRAY_ALPHA), PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT
      );

    // open the target file
    File targetFile( fopen( "verify.png", "wb" ) );
    if( !targetFile.IsValid() )
    {
      std::ostringstream oss;
      oss << "failed to open \"" << targetFileName << "\" for writing";
      throw Error( oss.str() );
    }

    // write the image
#if (PNG_LIBPNG_VER_MINOR >= 3)
    png_set_rows( targetPng.GetPng(), targetPng.GetInfo(), targetRows.Get() );
#endif

    png_set_write_fn( targetPng.GetPng(), (void *)targetFile.Get(), PNGWrite, PNGFlush);
//  png_init_io( targetPng.GetPng(), targetFile.Get() );

#if (PNG_LIBPNG_VER_MINOR >= 3)
    png_write_png( targetPng.GetPng(), targetPng.GetInfo(), PNG_TRANSFORM_IDENTITY, 0 );
#else
    png_write_info( targetPng.GetPng(), targetPng.GetInfo() );
    sourceImage.GetRows()->Write( &targetPng );
    png_write_end( targetPng.GetPng(), targetPng.GetInfo() );
#endif
  }
#endif
}

static void Decompress(std::string const& sourceFileName, std::string const& targetFileName, int paint, int mapping, int flags)
{
  // open the source file
  File sourceFile(fopen(sourceFileName.c_str(), "rb"));
  if (!sourceFile.IsValid()) {
    std::ostringstream oss;
    oss << "failed to open \"" << sourceFileName << "\" for reading";
    throw Error(oss.str());
  }

  // get the width and height
  int width, height;
  fread(&width , sizeof(int), 1, sourceFile.Get());
  fread(&height, sizeof(int), 1, sourceFile.Get());

  if ((mapping == -1) && ((flags & kBtcp) == kBtc4))
    mapping = 0;
  if ((mapping == -1) && ((flags & kBtcp) == kBtc5))
    mapping = 3;

  // work out the data size
  int bytesPerBlock = 16;
  int sourceDataSize = bytesPerBlock * width * height / 16;
  Mem sourceData(sourceDataSize);

  // read the source data
  fread(sourceData.Get(), 1, sourceDataSize, sourceFile.Get());

  // create the target rows
  PngRows targetRows(width, height, 4);
  
  struct sqio s = GetSquishIO(width, height, sqio::dtp::DT_U8, flags);

  // loop over blocks and compress them
  clock_t start = std::clock();
  unsigned char const* _sourceBlock = sourceData.Get();

#if (_MSC_VER >= 1700) && defined(NDEBUG)
  concurrency::parallel_for(0, height, 4, [&](int y) {
#else
  for (int y = 0; y < height; y += 4) {
#endif
    // process a row of blocks
    for (int x = 0; x < width; x += 4) {
      unsigned char const* sourceBlock = _sourceBlock + ((y / 4) * (width / 4) + (x / 4)) * bytesPerBlock;

      // decompress back
      u8 targetRgba[16 * 4];

      s.decoder(targetRgba, sourceBlock, flags);
      
      if (flags & kSignedExternal) {
	for (int i = 0; i < 16 * 4; ++i)
	  targetRgba[i] = targetRgba[i] + 0x80;
      }

      // write the data into the target rows
      for (int py = 0, i = 0; py < 4; ++py) {
	u8 *row = (u8 *)targetRows.Get()[y + py] + x * 4;

	for (int px = 0; px < 4; ++px, ++i) {
	  if ((flags & kBtcp) == kBtc5) {
	    // duplicate alpha from g (required!)
	    // duplicate g/b from r (for greyscale)
	    if ((flags & kColourMetrics) == kColourMetricUnit)
	      targetRgba[4 * i + 3] = 255;
	    else {
	      if (mapping == 3) {
		targetRgba[4 * i + 3] = targetRgba[4 * i + 1];
		targetRgba[4 * i + 2] = targetRgba[4 * i + 0];
		targetRgba[4 * i + 1] = targetRgba[4 * i + 0];
		targetRgba[4 * i + 0] = targetRgba[4 * i + 0];
	      }
	      else {
		targetRgba[4 * i + 3          ] = 255;
		targetRgba[4 * i + 0 + mapping] = targetRgba[4 * i + 1];
		targetRgba[4 * i + 3 - mapping] = targetRgba[4 * i + 0];
		targetRgba[4 * i + 0          ] = targetRgba[4 * i + 0];
	      }
	    }
	  }
	  else if ((flags & kBtcp) == kBtc4) {
	    if (mapping == 3) {
	      targetRgba[4 * i + 3] = targetRgba[4 * i + 0];
	      targetRgba[4 * i + 2] = 0;
	      targetRgba[4 * i + 1] = 0;
	      targetRgba[4 * i + 0] = 0;
	    }
	    else {
	      targetRgba[4 * i + 3] = 255;
	      targetRgba[4 * i + 2] = targetRgba[4 * i + 0];
	      targetRgba[4 * i + 1] = targetRgba[4 * i + 0];
	      targetRgba[4 * i + 0] = targetRgba[4 * i + 0];
	    }
	  }

	  if (flags & kExcludeAlphaFromPalette)
	    targetRgba[4 * i + 3] = 255;

	  for (int j = 0; j < 4; ++j)
	    *row++ = targetRgba[4 * i + j];
	}
      }

      if (((flags & kBtcp) == kBtc7) && paint) {
	// draw the mode
	if (paint == 2) {
	  unsigned long btcvalue = *((int *)sourceBlock);
	  unsigned long btcmode; _BitScanForward(&btcmode, btcvalue); btcmode += 1;

	  fprintf(stderr, "%1u", btcmode);
	}
	// draw the pattern
	else if (paint == 1) {
	  char *patchar = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzµß";
	  unsigned long btcvalue = *((int *)sourceBlock);
	  unsigned long btcmode; _BitScanForward(&btcmode, btcvalue); btcmode += 1;
	  unsigned long btcpat = (btcmode == 5 || btcmode == 6 || btcmode == 7 ? 0 : (btcmode == 1 ? (1 << 4) - 1 : (1 << 6) - 1));

	  fprintf(stderr, "%c", patchar[(btcvalue >> btcmode) & btcpat]);
	}
      }

      // advance
      sourceBlock += bytesPerBlock;
    }

    // draw the mode
    if (((flags & kBtcp) == kBtc7) && paint) {
      fprintf(stderr, "\n");
    }
  }
#if (_MSC_VER >= 1700) && defined(NDEBUG)
  );
#endif

  clock_t end = std::clock();
  double duration = (double)(end - start) / CLOCKS_PER_SEC;
  std::cout << "time taken: " << duration << " seconds" << std::endl;
  std::cout << "decompression speed: " << ((double)(width * height) / (duration * 1024 * 1024)) << " MPixel/s" << std::endl;

  // create the target PNG
  PngWriteStruct targetPng;

  // set up the image
  png_set_IHDR(
    targetPng.GetPng(), targetPng.GetInfo(), width, height,
    8 * sizeof(u8), PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT
  );

  // open the target file
  File targetFile( fopen( targetFileName.c_str(), "wb"));
  if (!targetFile.IsValid()) {
    std::ostringstream oss;
    oss << "failed to open \"" << targetFileName << "\" for writing";
    throw Error(oss.str());
  }

  // write the image
#if (PNG_LIBPNG_VER_MINOR >= 3)
  png_set_rows( targetPng.GetPng(), targetPng.GetInfo(), targetRows.Get() );
#endif

  png_set_write_fn(targetPng.GetPng(), (void *)targetFile.Get(), PNGWrite, PNGFlush);
//png_init_io( targetPng.GetPng(), targetFile.Get() );

#if (PNG_LIBPNG_VER_MINOR >= 3)
  png_write_png( targetPng.GetPng(), targetPng.GetInfo(), PNG_TRANSFORM_IDENTITY, 0 );
#else
  png_write_info(targetPng.GetPng(), targetPng.GetInfo());
  targetRows.Write(&targetPng);
  png_write_end(targetPng.GetPng(), targetPng.GetInfo());
#endif
}

static void Diff(std::string const& sourceFileName, std::string const& targetFileName)
{
  // load the images
  PngImage sourceImage(sourceFileName);
  PngImage targetImage(targetFileName);

  // get the image info
  int width = sourceImage.GetWidth();
  int height = sourceImage.GetHeight();
  int sourceStride = sourceImage.GetStride();
  int targetStride = targetImage.GetStride();
  int stride = std::min(sourceStride, targetStride);
  
  // check they match
  if (width != targetImage.GetWidth() || height != targetImage.GetHeight())
    throw Error("source and target dimensions do not match");

  // work out the error
  double error = 0.0;
  for (int y = 0; y < height; ++y) {
    u8 const* sourceRow = sourceImage.GetRow(y);
    u8 const* targetRow = targetImage.GetRow(y);

    for (int x = 0; x < width; ++x) {
      u8 const* sourcePixel = sourceRow + x * sourceStride;
      u8 const* targetPixel = targetRow + x * targetStride;

      for (int i = 0; i < stride; ++i) {
	int diff = (int)sourcePixel[i] - (int)targetPixel[i];
	error += (double)(diff * diff);
      }
    }
  }

  error = std::sqrt(error / (width * height));

  // print it out
  std::cout << "rms error: " << error << std::endl;
}

static void Benchmark(std::string const& sourceFileName, int mapping, int flags)
{
  // load the source image
  PngImage sourceImage(sourceFileName);

  // get the image info
  int width = sourceImage.GetWidth();
  int height = sourceImage.GetHeight();
  int stride = sourceImage.GetStride();
  bool colour = sourceImage.IsColour();
  bool alpha = sourceImage.IsAlpha();

  // create the target rows
  PngRows targetRows(width, height, 4);

  /* kill alpha-channel */
  if (!alpha && (((flags & kBtcp) <= kBtc3) || ((flags & kBtcp) >= kBtc7)))
    flags = (flags | kExcludeAlphaFromPalette) & (~kWeightColourByAlpha) & (~kAlphaIterativeFit);
  if ((mapping == -1) && ((flags & kBtcp) == kBtc4))
    mapping = 0;
  if ((mapping == -1) && ((flags & kBtcp) == kBtc5))
    mapping = 3;

  // check the image dimensions
  if ((width % 4) != 0 || (height % 4) != 0) {
    std::ostringstream oss;
    oss << "cannot compress " << width << "x" << height
      << "image (dimensions must be multiples of 4)";

    throw Error(oss.str());
  }

  // create the target data
  int bytesPerBlock = 16;
  int benchDataSize = bytesPerBlock * width * height / 16;
  Mem benchData(benchDataSize);

  int num = std::max(1, 8192 / width) * std::max(1, 8192 / width) * 2;
#if (_MSC_VER >= 1700) && defined(NDEBUG)
//high_resolution_clock::time_point start;
//high_resolution_clock::time_point end;
//duration<double> curduration;

  LARGE_INTEGER start;
  LARGE_INTEGER end;
  LARGE_INTEGER performanceFrequency; QueryPerformanceFrequency(&performanceFrequency);
  double frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
  double microsecondsToSeconds = 1000000.;

  double curduration;
  double minduration;
#else
  clock_t start;
  clock_t end;
  double curduration;
  double minduration;
#endif

  minduration = DBL_MAX;
  for (int l = 0; l < num; l += 1) {
    struct sqio s = GetSquishIO(width, height, sqio::dtp::DT_U8, flags);

    unsigned char* _targetBlock = benchData.Get();

    // loop over blocks and compress them
#if (_MSC_VER >= 1700) && defined(NDEBUG)
//  start = high_resolution_clock::now();
    QueryPerformanceCounter(&start);
#else
    start = std::clock();
#endif
    
#if (_MSC_VER >= 1700) && defined(NDEBUG)
    concurrency::parallel_for(0, height, 4, [&](int y) {
#else
    for (int y = 0; y < height; y += 4) {
#endif
      // process a row of blocks
      for (int x = 0; x < width; x += 4) {
	unsigned char* targetBlock = _targetBlock + ((y / 4) * (width / 4) + (x / 4)) * bytesPerBlock;

	// get the block data
	u8 sourceRgba[16 * 4];

	for (int py = 0, i = 0; py < 4; ++py) {
	  u8 const* row = sourceImage.GetRow(y + py) + x * stride;
	  for (int px = 0; px < 4; ++px, ++i) {
	    // get the pixel colour
	    if (colour) {
	      for (int j = 0; j < 3; ++j)
		sourceRgba[4 * i + j] = *row++;
	    }
	    else {
	      for (int j = 0; j < 3; ++j)
		sourceRgba[4 * i + j] = *row; ++row;
	    }

	    // skip alpha for now
	    if (alpha)
	      sourceRgba[4 * i + 3] = *row++;
	    else
	      sourceRgba[4 * i + 3] = 255;

	    if ((flags & kBtcp) == kBtc5) {
	      // duplicate alpha into g (required!)
	      // duplicate r into b (for symmetry)
	      if ((flags & kColourMetrics) == kColourMetricUnit)
		;
	      else {
		sourceRgba[4 * i + 4 - mapping] = sourceRgba[4 * i + mapping];
		sourceRgba[4 * i + 2          ] = sourceRgba[4 * i + 0      ];
	      }
	    }
	    else if ((flags & kBtcp) == kBtc4) {
	      sourceRgba[4 * i + 0] = sourceRgba[4 * i + mapping];
	      sourceRgba[4 * i + 1] = sourceRgba[4 * i + mapping];
	      sourceRgba[4 * i + 2] = sourceRgba[4 * i + mapping];
	    }
	  }
	}
	
	if (flags & kSignedExternal) {
	  for (int i = 0; i < 16 * 4; ++i)
	    sourceRgba[i] = sourceRgba[i] - 0x80;
	}

	// compress this block
	s.encoder(sourceRgba, -1, targetBlock, flags);

	// advance
	targetBlock += bytesPerBlock;
      }
    }
#if (_MSC_VER >= 1700) && defined(NDEBUG)
    );
#endif

#if (_MSC_VER >= 1700) && defined(NDEBUG)
//  end = high_resolution_clock::now();
//  curduration = duration_cast< duration<double> >(end - start);
//  minduration = std::min(minduration, curduration.count());

    QueryPerformanceCounter(&end);

    end.QuadPart -= start.QuadPart;
    curduration = (double)end.QuadPart / frequencyToMicroseconds;
    minduration = std::min(minduration, curduration / microsecondsToSeconds);
#else
    end = std::clock();
    curduration = (double)(end - start) / CLOCKS_PER_SEC;
    minduration = std::min(minduration, curduration);
#endif
  }

  std::cout << "time taken: " << minduration << " seconds" << std::endl;
  std::cout << "compression speed: " << ((double) (width * height) / (minduration * 1024 * 1024)) << " MPixel/s" << std::endl;

  minduration = DBL_MAX;
  for (int l = 0; l < num; l += 1) {
    struct sqio s = GetSquishIO(width, height, sqio::dtp::DT_U8, flags);

    unsigned char const* sourceBlock = benchData.Get();

    // loop over blocks and compress them
#if (_MSC_VER >= 1700) && defined(NDEBUG)
//  start = high_resolution_clock::now();
    QueryPerformanceCounter(&start);
#else
    start = std::clock();
#endif

    for (int y = 0; y < height; y += 4) {
      // process a row of blocks
      for (int x = 0; x < width; x += 4) {
	// decompress back
	u8 targetRgba[16 * 4];

	s.decoder(targetRgba, sourceBlock, flags);

	// write the data into the target rows
	for (int py = 0, i = 0; py < 4; ++py) {
	  u8* row = (u8*)targetRows.Get()[y + py] + x * 4;
	  for (int px = 0; px < 4; ++px, ++i) {
	    for (int j = 0; j < 4; ++j)
	      *row++ = targetRgba[4 * i + j];
	  }
	}

	// advance
	sourceBlock += bytesPerBlock;
      }
    }

#if (_MSC_VER >= 1700) && defined(NDEBUG)
//  end = high_resolution_clock::now();
//  curduration = duration_cast< duration<double> >(end - start);
//  minduration = std::min(minduration, curduration.count());

    QueryPerformanceCounter(&end);

    end.QuadPart -= start.QuadPart;
    curduration = (double)end.QuadPart / frequencyToMicroseconds;
    minduration = std::min(minduration, curduration / microsecondsToSeconds);
#else
    end = std::clock();
    curduration = (double)(end - start) / CLOCKS_PER_SEC;
    minduration = std::min(minduration, curduration);
#endif
  }

  std::cout << "time taken: " << minduration << " seconds" << std::endl;
  std::cout << "decompression speed: " << ((double)(width * height) / (minduration * 1024 * 1024)) << " MPixel/s" << std::endl;
}

enum Mode
{
  kCompress,
  kDecompress,
  kDiff,
  kBenchmark
};

int main(int argc, char* argv[]) {
  try {
    // parse the command-line
    std::string sourceFileName;
    std::string targetFileName;
    Mode mode = kBenchmark;//kCompress;
    int method = kBtc1;//kBtc7;
    int metric = kColourMetricPerceptual;
    int fit = kColourClusterFit;//kColourClusterFit;
    int alpha = 0;//kAlphaIterativeFit;
    int extra = 0;
    int paint = 0;
    int sign = 0;
    int mapping = -1;
    bool help = false;
    bool arguments = true;
    
    for (int i = 1; i < argc; ++i) {
      // check for options
      char const* word = argv[i];
      if (arguments && word[0] == '-') {
	for (int j = 1; word[j] != '\0'; ++j) {
	  switch (word[j]) {
	    case 'h': help = true; break;

	    case 'c': mode = kCompress; break;
	    case 'd': mode = kDecompress; break;
	    case 'e': mode = kDiff; break;
	    case 'b': mode = kBenchmark; break;

	    case 'p': paint = 1; break;
	    case 'm': paint = 2; break;

	    case '1': method = kBtc1; break;
	    case '2': method = kBtc2; break;
	    case '3': method = kBtc3; break;
	    case '4': method = kBtc4; break;
	    case '5': method = kBtc5; break;
	    case '7': method = kBtc7; break;
	    case '0': method = kCtx1; break;

	    case 'A': mapping = 3; break;
	    case 'B': mapping = 2; break;
	    case 'G': mapping = 1; break;
	    case 'R': mapping = 0; break;

	    case 'u': metric = kColourMetricUniform; break;
	    case 'l': metric = kColourMetricPerceptual; break;
	    case 'n': metric = kColourMetricUnit; break;

	    case 'a': alpha = kAlphaIterativeFit; break;
	    case 'r': fit = kColourRangeFit; break;
	    case 'i': fit = kColourIterativeClusterFit; break;
	    case 'x': fit = kColourClusterFit * 15; break;
	      
	    case 's': sign = kSignedExternal + kSignedInternal; break;

	    case 'w': extra = kWeightColourByAlpha; break;
	    case 'W': extra = kWeightColourByAlpha + kExcludeAlphaFromPalette; break;
	    case '-': arguments = false; break;
	    default:
	      std::cerr << "unknown option '" << word[j] << "'" << std::endl;
	      return -1;
	  }
	}
      }
      else {
	if (sourceFileName.empty())
	  sourceFileName.assign(word);
	else if (targetFileName.empty())
	  targetFileName.assign(word);
	else {
	  std::cerr << "unexpected argument \"" << word << "\"" << std::endl;
	}
      }
    }

    // check arguments
    if (help) {
      std::cout
	<< "SYNTAX" << std::endl
	<< "\tsquishldrpng [-bcde123457] <source> <target>" << std::endl
	<< "OPTIONS" << std::endl
	<< "\t-c\tCompress source png to target raw btc (default)" << std::endl
	<< "\t-d\tDecompress source raw btc to target png" << std::endl
	<< "\t-b\tBenchmark the chosen config" << std::endl
	<< "\t-e\tDiff source and target png" << std::endl
	<< "\t-123\tSpecifies whether to use DXT1/BC1, DXT3/BC2 or DXT5/BC3 compression" << std::endl
	<< "\t-45\tSpecifies whether to use ATI/BC4, ATI2/BC5 compression" << std::endl
	<< "\t-7\tSpecifies whether to use BC7 compression" << std::endl
	<< "\t-0\tSpecifies whether to use CTX1 compression" << std::endl
	<< "\t-s\tSpecifies whether to signed block compression" << std::endl
	<< "\t-a\tUse the slow iterative alpha/gray/normal compressor" << std::endl
	<< "\t-r\tUse the fast but inferior range-based colour compressor" << std::endl
	<< "\t-i\tUse the very slow but slightly better iterative colour compressor" << std::endl
	<< "\t-x\tUse the extreme slow but slightly better iterative colour compressor" << std::endl
	<< "\t-w\tWeight colour values by alpha in the cluster colour compressor" << std::endl
	<< "\t-W\tWeight colour values by alpha, but store everything opaque" << std::endl
	<< "\t-n\tUse normal map metrics and algorithms" << std::endl
	<< "\t-l\tUse CIE XYZ luminance metrics during colour compression" << std::endl
	<< "\t-u\tUse a uniform colour metric during colour compression" << std::endl
	<< "\t-p\tPaint the patterns used in BC7 to stderr" << std::endl
	<< "\t-m\tPaint the modes used in BC7 to stderr" << std::endl
	<< "\t-A\tStore RA into BC5 (default), store A into BC4" << std::endl
	<< "\t-G\tStore RG into BC5" << std::endl
	<< "\t-R\tStore R into BC4 (default)" << std::endl
//	<< "\t-g\tConvert png from gamma to linear and back to gamma png" << std::endl
//	<< "\t-s\tConvert unsigned png to signed and back to unsigned png" << std::endl
	;

      return 0;
    }

    if (sourceFileName.empty()) {
      std::cerr << "no source file given" << std::endl;
      return -1;
    }

    if (targetFileName.empty() && (mode != kBenchmark)) {
      std::cerr << "no target file given" << std::endl;
      return -1;
    }

    // do the work
    switch (mode) {
      case kCompress:
	Compress(sourceFileName, targetFileName, paint, mapping, method + sign + metric + fit + alpha + extra);
	break;

      case kDecompress:
	Decompress(sourceFileName, targetFileName, paint, mapping, method + sign + metric /*+ extra*/);
	break;

      case kDiff:
	Diff(sourceFileName, targetFileName);
	break;

      case kBenchmark:
	Benchmark(sourceFileName, mapping, method + metric + sign + fit + alpha + extra);
	break;

      default:
	std::cerr << "unknown mode" << std::endl;
	throw std::exception();
    }
  }
  catch (std::exception& excuse) {
    // complain
    std::cerr << "squishpng error: " << excuse.what() << std::endl;
    return -1;
  }

  // done
  return 0;
}