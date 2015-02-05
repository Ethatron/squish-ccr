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

#ifndef SQUISH_SIMD_SSE_H
#define SQUISH_SIMD_SSE_H

#include <xmmintrin.h>
#if ( SQUISH_USE_SSE > 1 )
#include <emmintrin.h>
#endif
#if ( SQUISH_USE_SSE >= 3 )
#include <pmmintrin.h>
#endif
#if ( SQUISH_USE_SSE >= 4 )
#include <smmintrin.h>
#endif
#if ( SQUISH_USE_XSSE == 3 )
#include <tmmintrin.h>
#endif
#if ( SQUISH_USE_XSSE == 4 )
#if defined(_MSC_VER) && (_MSC_VER > 1300)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>
#endif
#endif

#pragma warning(disable: 4127)

#define SQUISH_SSE_SPLAT( a )								\
	( ( a ) | ( ( a ) << 2 ) | ( ( a ) << 4 ) | ( ( a ) << 6 ) )

#define SQUISH_SSE_SHUF( x, y, z, w )							\
	( ( x ) | ( ( y ) << 2 ) | ( ( z ) << 4 ) | ( ( w ) << 6 ) )

#define SQUISH_SSE_SWAP64( )								\
	SQUISH_SSE_SHUF( 2, 3, 0, 1 )

#define SQUISH_SSE_SWAP32( )								\
	SQUISH_SSE_SHUF( 3, 2, 1, 0 )

#define SQUISH_SSE_SWAP16( )								\
	SQUISH_SSE_SHUF( 3, 2, 1, 0 )

#define _mm_shuffle_epi16(r,s)	_mm_shufflelo_epi16( _mm_shufflehi_epi16( res, s ), s )

namespace squish {

#define COL4_CONST( X ) Col4( X )


class Col3
{
public:
	typedef Col3 const& Arg;

	Col3() {}

	explicit Col3( __m128i v ) : m_v( v ) {}

	Col3( Arg arg ) : m_v( arg.m_v ) {}

	Col3& operator=( Arg arg )
	{
		m_v = arg.m_v;
		return *this;
	}

	explicit Col3(int s) : m_v( _mm_set1_epi32( s ) ) {}
	explicit Col3(float s) : m_v( _mm_cvttps_epi32( _mm_set1_ps ( s ) ) ) {}

	Col3( int r, int g, int b ) : m_v( _mm_setr_epi32( r, g, b, 0 ) ) {}
	Col3( int r, int g ) : m_v( _mm_setr_epi32( r, g, 0, 0 ) ) {}
	Col3( u16 r, u16 g, u16 b ) : m_v( _mm_setr_epi32( r, g, b, 0 ) ) {}
	Col3( s16 r, s16 g, s16 b ) : m_v( _mm_setr_epi32( r, g, b, 0 ) ) {}
	Col3( u16 r, u16 g ) : m_v( _mm_setr_epi32( r, g, 0, 0 ) ) {}
	Col3( s16 r, s16 g ) : m_v( _mm_setr_epi32( r, g, 0, 0 ) ) {}
	Col3( u8 r, u8 g, u8 b ) : m_v( _mm_setr_epi32( r, g, b, 0 ) ) {}
	Col3( s8 r, s8 g, s8 b ) : m_v( _mm_setr_epi32( r, g, b, 0 ) ) {}
	Col3( u8 r, u8 g ) : m_v( _mm_setr_epi32( r, g, 0, 0 ) ) {}
	Col3( s8 r, s8 g ) : m_v( _mm_setr_epi32( r, g, 0, 0 ) ) {}

	explicit Col3( unsigned int s ) : m_v( _mm_set1_epi32( s ) ) {}
	explicit Col3( const unsigned int (&_rgb)[3] ) : m_v( _mm_load_si128( (const __m128i *)&_rgb ) ) {}
	explicit Col3( u8 const *source ) : m_v( _mm_loadu_si128( (const __m128i *)source ) ) {}

	int GetM8() const
	{
		return _mm_movemask_epi8 ( m_v );
	}

	int GetLong() const
	{
		return _mm_cvtsi128_si32 ( m_v );
	}

	Col3 SetLong( int v ) const
	{
		return Col3 ( _mm_cvtsi32_si128 ( v ) );
	}

	int R() const { return _mm_extract_epi16( m_v, 0 ); }
	int G() const { return _mm_extract_epi16( m_v, 2 ); }
	int B() const { return _mm_extract_epi16( m_v, 4 ); }

	int &GetR() { return ((int *)&m_v)[0]; }
	int &GetG() { return ((int *)&m_v)[1]; }
	int &GetB() { return ((int *)&m_v)[2]; }
	// let the compiler figure this one out, probably spills to memory
	int &GetO(int o) { return ((int *)&m_v)[o]; }

	const int &GetR() const { return ((const int *)&m_v)[0]; }
	const int &GetG() const { return ((const int *)&m_v)[1]; }
	const int &GetB() const { return ((const int *)&m_v)[2]; }
	// let the compiler figure this one out, probably spills to memory
	const int &GetO(int o) const { return ((const int *)&m_v)[o]; }

	Col3 SplatR() const { return Col3( _mm_shuffle_epi32( m_v, SQUISH_SSE_SPLAT( 0 ) ) ); }
	Col3 SplatG() const { return Col3( _mm_shuffle_epi32( m_v, SQUISH_SSE_SPLAT( 1 ) ) ); }
	Col3 SplatB() const { return Col3( _mm_shuffle_epi32( m_v, SQUISH_SSE_SPLAT( 2 ) ) ); }

	template<const int inv>
	void SetRGB( int r, int g, int b )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, r, 0 );
		v = _mm_insert_epi16( v, g, 2 );
		v = _mm_insert_epi16( v, b, 4 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		m_v = v;
	}

	template<const int inv>
	void SetRGBpow2( int r, int g, int b )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, r, 0 );
		v = _mm_insert_epi16( v, g, 2 );
		v = _mm_insert_epi16( v, b, 4 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		v = _mm_slli_epi32( v, 23 );
		v = _mm_add_epi32( v, _mm_castps_si128( _mm_set1_ps(1.0f) ) );

		m_v = _mm_cvttps_epi32( _mm_castsi128_ps( v ) );
	}

	template<const int inv>
	void SetRGBpow2( int c )
	{
		__m128i v = _mm_shuffle_epi32( _mm_cvtsi32_si128( c ), SQUISH_SSE_SPLAT( 0 ) );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		v = _mm_slli_epi32( v, 23 );
		v = _mm_add_epi32( v, _mm_castps_si128( _mm_set1_ps(1.0f) ) );

		m_v = _mm_cvttps_epi32( _mm_castsi128_ps( v ) );
	}

	Col3& operator&=( Arg v )
	{
		m_v = _mm_and_si128( m_v, v.m_v );
		return *this;
	}

	Col3& operator%=( Arg v )
	{
		m_v = _mm_andnot_si128( m_v, v.m_v );
		return *this;
	}

	Col3& operator^=( Arg v )
	{
		m_v = _mm_xor_si128( m_v, v.m_v );
		return *this;
	}

	Col3& operator|=( Arg v )
	{
		m_v = _mm_or_si128( m_v, v.m_v );
		return *this;
	}

	Col3& operator>>=( const int n )
	{
		m_v = _mm_srli_epi32( m_v, n );
		return *this;
	}

	Col3& operator<<=( const int n )
	{
		m_v = _mm_slli_epi32( m_v, n );
		return *this;
	}

	Col3& operator+=( Arg v )
	{
		m_v = _mm_add_epi32( m_v, v.m_v );
		return *this;
	}

	Col3& operator-=( Arg v )
	{
		m_v = _mm_sub_epi32( m_v, v.m_v );
		return *this;
	}

	Col3& operator*=( Arg v )
	{
#if ( SQUISH_USE_SSE >= 4 )
		m_v = _mm_mullo_epi32( m_v, v.m_v );
#else
		m_v = _mm_mullo_epi16( m_v, v.m_v );
#endif
		return *this;
	}

	Col3& operator*=( short v )
	{
#if ( SQUISH_USE_SSE >= 4 )
		m_v = _mm_mullo_epi32( m_v, _mm_set1_epi32(v) );
#else
		m_v = _mm_mulhi_epu16( m_v, _mm_set1_epi32(v) );
#endif
		return *this;
	}

	Col3& operator/=( short v )
	{
		__m128
			
		fp = _mm_cvtepi32_ps(m_v);
		fp = _mm_div_ps(fp, _mm_set1_ps(v));
		m_v = _mm_cvttps_epi32(fp);

		return *this;
	}

	friend int operator<( Arg left, Arg right  )
	{
		return CompareFirstLessThan(left, right);
	}

	friend int operator>( Arg left, Arg right  )
	{
		return CompareFirstGreaterThan(left, right);
	}

	friend int operator==( Arg left, Arg right  )
	{
		return CompareFirstEqualTo(left, right);
	}

	friend Col3 operator~( Arg left )
	{
		return left ^ Col3(~0);
	}

	friend Col3 operator&( Arg left, Arg right  )
	{
		return Col3( _mm_and_si128( left.m_v, right.m_v ) );
	}

	friend Col3 operator%( Arg left, Arg right )
	{
		return Col3( _mm_andnot_si128( left.m_v, right.m_v ) );
	}

	friend Col3 operator^( Arg left, Arg right )
	{
		return Col3( _mm_xor_si128( left.m_v, right.m_v ) );
	}

	friend Col3 operator|( Arg left, Arg right )
	{
		return Col3( _mm_or_si128( left.m_v, right.m_v ) );
	}

	friend Col3 operator>>( Arg left, int right )
	{
		return Col3( _mm_srli_epi32( left.m_v, right ) );
	}

	friend Col3 operator<<( Arg left, int right )
	{
		return Col3( _mm_slli_epi32( left.m_v, right ) );
	}

	friend Col3 operator+( Arg left, Arg right )
	{
		return Col3( _mm_add_epi32( left.m_v, right.m_v ) );
	}

	friend Col3 operator-( Arg left, Arg right )
	{
		return Col3( _mm_sub_epi32( left.m_v, right.m_v ) );
	}

	friend Col3 operator*( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col3( _mm_mullo_epi32( left.m_v, right.m_v ) );
#else
		return Col3( _mm_mullo_epi16( left.m_v, right.m_v ) );
#endif
	}

	friend Col3 operator*( Arg left, int right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col3( _mm_mullo_epi32( left.m_v, _mm_set1_epi32( right ) ) );
#else
		return Col3( _mm_mullo_epi16( left.m_v, _mm_set1_epi32( right ) ) );
#endif
	}

	template<const int n>
	friend Col3 ShiftLeft( Arg a );
	template<const int n>
	friend Col3 ShiftLeft( Arg a )
	{
		if ((n) <= 0)
			return Col3( a.m_v );
		if ((n) <= 7)
			return Col3( _mm_slli_epi32( a.m_v, (n) & 7 ) );
		if ((n) & 7)
			return Col3( _mm_slli_epi32( _mm_slli_si128( a.m_v, (n) >> 3 ), (n) & 7 ) );

			return Col3( _mm_slli_si128( a.m_v, (n) >> 3 ) );
	}

	template<const int n>
	friend Col3 ShiftRight( Arg a );
	template<const int n>
	friend Col3 ShiftRight( Arg a )
	{
		if ((n) <= 0)
			return Col3( a.m_v );
		if ((n) <= 7)
			return Col3( _mm_srli_epi32( a.m_v, (n) & 7 ) );
		if ((n) & 7)
			return Col3( _mm_srli_epi32( _mm_srli_si128( a.m_v, (n) >> 3 ), (n) & 7 ) );

			return Col3( _mm_srli_si128( a.m_v, (n) >> 3 ) );
	}

	template<const int n>
	friend Col3 ShiftRightHalf( Arg a );
	template<const int n>
	friend Col3 ShiftRightHalf( Arg a )
	{
		return Col3( (n) > 0 ? _mm_srli_epi64( a.m_v, (n) ) : a.m_v );
	}

	friend Col3 ShiftRightHalf( Arg a, const int n )
	{
		return Col3( _mm_srl_epi64( a.m_v, _mm_cvtsi32_si128( n ) ) );
	}

	friend Col3 ShiftRightHalf( Arg a, Arg b )
	{
		return Col3( _mm_srl_epi64( a.m_v, b.m_v ) );
	}

	template<const int n>
	friend Col3 ShiftLeftHalf( Arg a );
	template<const int n>
	friend Col3 ShiftLeftHalf( Arg a )
	{
		return Col3( (n) > 0 ? _mm_slli_epi64( a.m_v, (n) ) : a.m_v );
	}

	friend Col3 ShiftLeftHalf( Arg a, const int n )
	{
		return Col3( _mm_sll_epi64( a.m_v, _mm_cvtsi32_si128( n ) ) );
	}

	template<const int r, const int g, const int b>
	friend Col3 ShiftLeftLo( Arg v )
	{
		// (1 << r, 1 << g, 1 << b);
		Col3 p2; p2.SetRGBpow2<0>(r, g, b);

	//	return Col3( _mm_mullo_epi32( v.m_v, p2.m_v ) );
		return Col3( _mm_mullo_epi16( v.m_v, p2.m_v ) );
	}

	template<const int n, const int p>
	friend Col3 MaskBits( Arg a );
	template<const int n, const int p>
	friend Col3 MaskBits( Arg a )
	{
		if ((p + n) <= 0)
			return Col3(0);
		if ((p + n) >= 64)
			return a;

		// compile time
		__int64 base = ~(0xFFFFFFFFFFFFFFFFULL << (     (p + n) & 63));
	//	__int64 base =  (0xFFFFFFFFFFFFFFFFULL >> (64 - (p + n) & 63));
		__m128i mask = _mm_setr_epi32(
		  (int)(base >>  0),
		  (int)(base >> 32), 0, 0
		);

		return Col3( _mm_and_si128( a.m_v, mask ) );
	}

	friend Col3 MaskBits( Arg a, const int n, const int p )
	{
		const int val = 64 - (p + n);

		__m128i shift = _mm_max_epi16( _mm_cvtsi32_si128( val ), _mm_set1_epi32( 0 ) );
		__m128i mask = _mm_setr_epi32(
		  0xFFFFFFFF,
		  0xFFFFFFFF, 0, 0
		);

		mask = _mm_srl_epi64( mask, shift );

		// (0xFFFFFFFFFFFFFFFFULL >> (64 - (p + n) & 63))
		return Col3( _mm_and_si128( a.m_v, mask ) );
	}

	template<const int n, const int p>
	friend Col3 CopyBits( Arg left, Arg right );
	template<const int n, const int p>
	friend Col3 CopyBits( Arg left, Arg right )
	{
		if (!(n))
			return left;
		if (!(p))
			return MaskBits<n, 0>(right);
		if (((p) + (n)) >= 64)
			return (left) + ShiftLeftHalf<p>(right);

#if ( SQUISH_USE_XSSE == 4 )
		return Col3( _mm_inserti_si64( left.m_v, right.m_v, n, p ) );
#else
		return MaskBits<p, 0>(left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
	//	return               (left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
#endif
	}

	friend Col3 CopyBits( Arg left, Col3 &right, const int n, const int p )
	{
#if ( SQUISH_USE_XSSE == 4 )
		/* ---- ---bl xxxx xxxx */
		const int val = (p << 8) + (n << 0);

		right.m_v = _mm_unpacklo_epi64( right.m_v, _mm_cvtsi32_si128( val ) );
		return Col3( _mm_insert_si64( left.m_v, right.m_v ) );
#else
		return MaskBits(left, p, 0) + MaskBits(ShiftLeftHalf(right, p), n, p);
	//	return         (left      ) + MaskBits(ShiftLeftHalf(right, p), n, p);
#endif
	}

	template<const int n, const int p>
	friend Col3 ExtrBits( Arg a );
	template<const int n, const int p>
	friend Col3 ExtrBits( Arg a )
	{
		if (!(n))
			return Col3(0);
		if (!(p))
			return MaskBits<n, 0>(a);
		if (((n) + (p)) >= 64)
			return ShiftRightHalf<p>(a);

#if ( SQUISH_USE_XSSE == 4 )
		return Col3( _mm_extracti_si64( a.m_v, n, p ) );
#else
		return MaskBits<n, 0>(ShiftRightHalf<p>(a));
#endif
	}

	friend Col3 ExtrBits( Arg a, const int n, const int p )
	{
#if ( SQUISH_USE_XSSE == 4 )
		/* ---- ----- ---- ---bl */
		const int val = (p << 8) + (n << 0);

		return Col3( _mm_extract_si64( a.m_v, _mm_cvtsi32_si128( val ) ) );
#else
		return MaskBits(ShiftRightHalf(a, p), n, 0);
#endif
	}

	template<const int n, const int p>
	friend void ExtrBits( Arg left, Col3 &right );
	template<const int n, const int p>
	friend void ExtrBits( Arg left, Col3 &right )
	{
		right  = ExtrBits<n, p>( left );
	}

	template<const int n, const int p>
	friend void ConcBits( Arg left, Col3 &right );
	template<const int n, const int p>
	friend void ConcBits( Arg left, Col3 &right )
	{
		right  = ShiftLeft<32>( right );
		if (n > 0)
			right += ExtrBits<n, p>( left );
	}

	template<const int n, const int p>
	friend void ReplBits( Arg left, Col3 &right );
	template<const int n, const int p>
	friend void ReplBits( Arg left, Col3 &right )
	{
		if (!n)
			return;
		if ((n < 0)) {
			right  = ExtrBits<-n, p>( left );
			right.m_v = _mm_shuffle_epi32( right.m_v, SQUISH_SSE_SHUF( 0, 0, 0, 3 ) );
		}
		else {
			right  = ExtrBits< n, p>( left );
			right.m_v = _mm_shuffle_epi32( right.m_v, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );
		}
	}

	friend Col3 Mul16x16u( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col3( _mm_mullo_epi32( a.m_v, b.m_v ) );
#else
		__m128i lo = _mm_mullo_epi16( a.m_v, b.m_v );
		__m128i hi = _mm_mulhi_epu16( a.m_v, b.m_v );

		return Col3( _mm_or_si128( lo, _mm_slli_si128( hi, 2 ) ) );
#endif
	}

	friend Col3 Mul16x16s( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col3( _mm_mullo_epi32( a.m_v, b.m_v ) );
#else
		__m128i lo = _mm_mullo_epi16( a.m_v, b.m_v );
		__m128i hi = _mm_mulhi_epi16( a.m_v, b.m_v );

		lo = _mm_and_si128( lo, _mm_set1_epi32(0x0000FFFF) );
		hi = _mm_slli_epi32( hi, 16 );

		return Col3( _mm_or_si128( lo, hi ) );
#endif
	}

	friend Col3 Div32x16u( Arg a, Arg b )
	{
		return Col3(
			(int)((unsigned int)a.GetR() / (unsigned int)b.GetR()),
			(int)((unsigned int)a.GetG() / (unsigned int)b.GetG()),
			(int)((unsigned int)a.GetB() / (unsigned int)b.GetB())
		);
	}

	friend Col3 Div32x16s( Arg a, Arg b )
	{
		return Col3(
			(int)((int)a.GetR() / (int)b.GetR()),
			(int)((int)a.GetG() / (int)b.GetG()),
			(int)((int)a.GetB() / (int)b.GetB())
		);
	}

	friend Col3 Div32x16sr( Arg a, Arg b )
	{
		Col3 r = a + (b >> 1);

		return Col3(
			(int)((int)r.GetR() / (int)b.GetR()),
			(int)((int)r.GetG() / (int)b.GetG()),
			(int)((int)r.GetB() / (int)b.GetB())
		);
	}

	//! Returns a*b + c
	friend Col3 MultiplyAdd( Arg a, Arg b, Arg c )
	{
#if ( SQUISH_USE_SSE >= 4 )
	  	return Col3( _mm_add_epi32( _mm_mullo_epi32( a.m_v, b.m_v ), c.m_v ) );
#else
		return Col3( _mm_add_epi32( _mm_mullo_epi16( a.m_v, b.m_v ), c.m_v ) );
#endif
	}

	//! Returns -( a*b - c )
	friend Col3 NegativeMultiplySubtract( Arg a, Arg b, Arg c )
	{
#if ( SQUISH_USE_SSE >= 4 )
	  	return Col3( _mm_sub_epi32( c.m_v, _mm_mullo_epi32( a.m_v, b.m_v ) ) );
#else
		return Col3( _mm_sub_epi32( c.m_v, _mm_mullo_epi16( a.m_v, b.m_v ) ) );
#endif
	}

	template<const int f, const int t>
	friend Col3 Shuffle( Arg a );
	template<const int f, const int t>
	friend Col3 Shuffle( Arg a )
	{
		if (f == t)
			return a;

		return Col3( _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SHUF(
			(t == 0 ? f : 0),
			(t == 1 ? f : 1),
			(t == 2 ? f : 2),
			(t == 3 ? f : 3)
		) ) );
	}

	template<const int f, const int t>
	friend Col3 Exchange( Arg a );
	template<const int f, const int t>
	friend Col3 Exchange( Arg a )
	{
		if (f == t)
			return a;

		return Col3( _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SHUF(
			(t == 0 ? f : (f == 0 ? t : 0)),
			(t == 1 ? f : (f == 1 ? t : 1)),
			(t == 2 ? f : (f == 2 ? t : 2)),
			(t == 3 ? f : (f == 3 ? t : 3))
		) ) );
	}

	friend Col3 HorizontalAdd( Arg a )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128i res = a.m_v;

		res = _mm_and_si128( res , _mm_setr_epi32( ~0, ~0, ~0, 0 ) );
		res = _mm_hadd_epi32( res, res );
		res = _mm_hadd_epi32( res, res );

		return Col3( res );
#else
		__m128i res = a.m_v;

		res = _mm_and_si128( res , _mm_setr_epi32( ~0, ~0, ~0, 0 ) );
		res = _mm_add_epi32( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP64() ) );
		res = _mm_add_epi32( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP32() ) );

		return Col3( res );
#endif
	}

	friend Col3 HorizontalAdd( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128i resc;

		resc = _mm_hadd_epi32( a.m_v, b.m_v );
		resc = _mm_and_si128( resc , _mm_setr_epi32( ~0, ~0, ~0, 0 ) );
		resc = _mm_hadd_epi32( resc, resc );
		resc = _mm_hadd_epi32( resc, resc );

		return Col3( resc );
#else
		__m128i resa = a.m_v;
		__m128i resb = b.m_v;
		__m128i resc;

		resc = _mm_add_epi32( resa, resb );
		resc = _mm_and_si128( resc , _mm_setr_epi32( ~0, ~0, ~0, 0 ) );
		resc = _mm_add_epi32( resc, _mm_shuffle_epi32( resc, SQUISH_SSE_SWAP64() ) );
		resc = _mm_add_epi32( resc, _mm_shuffle_epi32( resc, SQUISH_SSE_SWAP32() ) );

		return Col3( resc );
#endif
	}

	friend Col3 HorizontalAddTiny( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 ) && 0
		__m128 res = _mm_castsi128_ps ( a.m_v );

		res = _mm_dp_ps( res, _mm_set1_ps ( 1.0f ), 0x77 );

		return Col3( _mm_castps_si128 ( res ) );
#elif ( SQUISH_USE_SSE >= 3 )
		__m128 res = _mm_castsi128_ps ( a.m_v );

		// relies on correct de-normal floating-point treatment
		res = _mm_and_ps( res , _mm_castsi128_ps ( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
		res = _mm_hadd_ps( res, res );
		res = _mm_hadd_ps( res, res );

		return Col3( _mm_castps_si128 ( res ) );
#else
		return HorizontalAdd( a );
#endif
	}

	friend Col3 HorizontalAddTiny( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128 resa = _mm_castsi128_ps ( a.m_v );
		__m128 resb = _mm_castsi128_ps ( b.m_v );
		__m128 resc;

		// relies on correct de-normal floating-point treatment
		resc = _mm_hadd_ps( resa, resb );
		resc = _mm_and_ps( resc , _mm_castsi128_ps ( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
		resc = _mm_hadd_ps( resc, resc );
		resc = _mm_hadd_ps( resc, resc );

		return Col3( _mm_castps_si128 ( resc ) );
#else
		return HorizontalAdd( a, b );
#endif
	}
	
	friend Col3 HorizontalMaxTiny( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 ) && 0
		__m128i inv;
		__m128i res;

		inv = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SHUF( 0, 1, 2, 2 ) );
		inv = _mm_xor_si128( inv, _mm_set1_epi32 ( ~0 ) );
		inv = _mm_minpos_epu16 ( inv );
		res = _mm_xor_si128( inv, _mm_set1_epi32 ( ~0 ) );
		res = _mm_and_si128( inv, _mm_set1_epi32 ( 0x0000FFFF ) );

		return Col3( res ).SplatR();
#else
		__m128 res = _mm_castsi128_ps ( a.m_v );

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 2, 0, 1, 3 ) ) );
		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 2, 0, 3 ) ) );

		return Col3( _mm_castps_si128 ( res ) );
#endif
	}

	friend Col3 HorizontalMinTiny( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 )
		__m128i res;

		res = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SHUF( 0, 1, 2, 2 ) );
		res = _mm_minpos_epu16 ( res );
		res = _mm_and_si128( res , _mm_set1_epi32 ( 0x0000FFFF ) );

		return Col3( res ).SplatR();
#else
		__m128 res = _mm_castsi128_ps ( a.m_v );

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 2, 0, 1, 3 ) ) );
		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 2, 0, 3 ) ) );

		return Col3( _mm_castps_si128 ( res ) );
#endif
	}

	friend Col3 Dot( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return HorizontalAdd( Col3( _mm_mullo_epi32( left.m_v, right.m_v ) ) );
#else
		return HorizontalAdd( Col3( _mm_mullo_epi16( left.m_v, right.m_v ) ) );
#endif
	}

	friend Col3 DotTiny( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return HorizontalAdd    ( Col3( _mm_mullo_epi32( left.m_v, right.m_v ) ) );
#else
		return HorizontalAddTiny( Col3( _mm_mullo_epi16( left.m_v, right.m_v ) ) );
#endif
	}

	friend Col3 Abs( Arg a )
	{
		__m128i sign = _mm_srai_epi32( a.m_v, 31 );

		return Col3( _mm_sub_epi32( _mm_xor_si128( a.m_v, sign ), sign ) );
	}

	friend Col3 Min( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col3( _mm_min_epi32( left.m_v, right.m_v ) );
#else
		return Col3( _mm_min_epi16( left.m_v, right.m_v ) );
#endif
	}

	friend Col3 MinTiny( Arg left, Arg right )
	{
		__m128 resa = _mm_castsi128_ps( left.m_v );
		__m128 resb = _mm_castsi128_ps( right.m_v );
		__m128 resc;

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		resc = _mm_min_ps( resa, resb );

		return Col3( _mm_castps_si128 ( resc ) );
	}

	friend Col3 Max( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col3( _mm_max_epi32( left.m_v, right.m_v ) );
#else
		return Col3( _mm_max_epi16( left.m_v, right.m_v ) );
#endif
	}

	friend Col3 MaxTiny( Arg left, Arg right )
	{
	      __m128 resa = _mm_castsi128_ps( left.m_v );
	      __m128 resb = _mm_castsi128_ps( right.m_v );
	      __m128 resc;

	      // relies on correct de-normal floating-point treatment
	      // doesn't support negative values
	      resc = _mm_max_ps( resa, resb );

	      return Col3( _mm_castps_si128 ( resc ) );
	}
	
	friend bool CompareFirstLessThan( Arg left, Arg right )
	{
		__m128i bits = _mm_cmplt_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return !!(value & 0x000F);
	}

	friend bool CompareFirstGreaterThan( Arg left, Arg right )
	{
		__m128i bits = _mm_cmpgt_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return !!(value & 0x000F);
	}

	friend bool CompareFirstEqualTo( Arg left, Arg right )
	{
		__m128i bits = _mm_cmpeq_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return !!(value & 0x000F);
	}

	friend bool CompareAnyLessThan( Arg left, Arg right )
	{
		__m128i bits = _mm_cmpeq_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return (value & 0x0FFF) != 0x0000;
	}

	friend bool CompareAllEqualTo( Arg left, Arg right )
	{
		__m128i bits = _mm_cmpeq_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return (value & 0x0FFF) == 0x0FFF;
	}

	friend Col3 CompareAllEqualTo_M4( Arg left, Arg right )
	{
		return Col3( _mm_cmpeq_epi32( left.m_v, right.m_v ) );
	}

	friend int CompareGreaterThan( Arg left, Arg right )
	{
		return _mm_movemask_epi8( _mm_cmpgt_epi32( left.m_v, right.m_v ) );
	}

	friend Col3 IsOne( Arg v )
	{
		return Col3( _mm_cmpeq_epi32( v.m_v, _mm_set1_epi32( 0x000000FF ) ) );
	}

	friend Col3 IsZero( Arg v )
	{
		return Col3( _mm_cmpeq_epi32( v.m_v, _mm_setzero_si128( ) ) );
	}

	friend Col3 IsNotZero( Arg v )
	{
		return Col3( _mm_cmpgt_epi32( v.m_v, _mm_setzero_si128( ) ) );
	}

	friend void PackBytes( Arg a, unsigned int &loc )
	{
		__m128i

		r = _mm_packs_epi32( a.m_v, a.m_v );
		r = _mm_packus_epi16( r, r );

		loc = _mm_cvtsi128_si32( r );
	}
	
	friend void PackBytes( Arg a, int &loc )
	{
		__m128i

		r = _mm_packs_epi32( a.m_v, a.m_v );
		r = _mm_packs_epi16( r, r );

		loc = _mm_cvtsi128_si32( r );
	}
	
	friend void PackWords( Arg a, unsigned__int64 &loc )
	{
		__m128i

#if ( SQUISH_USE_SSE >= 4 )
		r = _mm_packus_epi32( a.m_v, a.m_v );
#else
		// fix-up/down
		r = _mm_sub_epi32( a.m_v, _mm_set1_epi32( 32768 ) );
		r = _mm_packs_epi32( r, r );
		r = _mm_add_epi16( r, _mm_set1_epi16( (short)-32768 ) );
#endif

//		loc = _mm_cvtsi128_si64( r );
		_mm_storel_epi64( (__m128i *)&loc, r );
	}
	
	friend void PackWords( Arg a, __int64 &loc )
	{
		__m128i
		  
		r = _mm_packs_epi32( a.m_v, a.m_v );

//		loc = _mm_cvtsi128_si64( r );
		_mm_storel_epi64( (__m128i *)&loc, r );
	}
	
	// clamp the output to [0, 1]
	Col3 Clamp() const {
		Col3 const one (0xFF);
		Col3 const zero(0x00);

		return Min(one, Max(zero, *this));
	}

	friend void LoadAligned( Col3 &a, Col3 &b, Arg c )
	{
	        a.m_v = c.m_v;
		b.m_v = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void LoadAligned( Col3 &a, void const *source )
	{
		a.m_v = _mm_load_si128( (__m128i const *)source );
	}

	friend void LoadAligned( Col3 &a, Col3 &b, void const *source )
	{
		a.m_v = _mm_load_si128( (__m128i const *)source );
		b.m_v = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void LoadUnaligned( Col3 &a, Col3 &b, void const *source )
	{
		a.m_v = _mm_loadu_si128( (__m128i const *)source );
		b.m_v = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void StoreAligned( Arg a, Arg b, Col3 &c )
	{
		c.m_v = _mm_unpacklo_epi64( a.m_v, b.m_v );
	}

	friend void StoreAligned( Arg a, void *destination )
	{
		_mm_store_si128( (__m128i *)destination, a.m_v );
	}

	friend void StoreAligned( Arg a, Arg b, void *destination )
	{
		_mm_store_si128( (__m128i *)destination, _mm_unpacklo_epi64( a.m_v, b.m_v ) );
	}
	
	friend void StoreUnaligned( Arg a, void *destination )
	{
		_mm_storeu_si128( (__m128i *)destination, a.m_v );
	}
	
	friend void StoreUnaligned( Arg a, Arg b, void *destination )
	{
		_mm_storeu_si128( (__m128i *)destination, _mm_unpacklo_epi64( a.m_v, b.m_v ) );
	}
	
	friend void StoreUnaligned( Arg a, u8* loc ) {
	  PackBytes( a, (unsigned int&) (*((unsigned int *)loc)) ); }
	friend void StoreUnaligned( Arg a, u16* loc ) {
	  PackWords( a, (unsigned__int64&) (*((unsigned__int64 *)loc)) ); }
	friend void StoreUnaligned( Arg a, s8* loc ) {
	  PackBytes( a, (int&) (*((int *)loc)) ); }
	friend void StoreUnaligned( Arg a, s16* loc ) {
	  PackWords( a, (__int64&) (*((__int64 *)loc)) ); }

private:
	__m128i m_v;

	friend class Col4;
	friend class Vec3;
};

class Col4
{
public:
	typedef Col4 const& Arg;

	Col4() {}

	explicit Col4( __m128i v ) : m_v( v ) {}
	explicit Col4( Col3::Arg v ) : m_v( v.m_v ) {}

	Col4( Arg arg ) : m_v( arg.m_v ) {}

	Col4& operator=( Arg arg )
	{
		m_v = arg.m_v;
		return *this;
	}

	explicit Col4(int s) : m_v( _mm_set1_epi32( s ) ) {}
	explicit Col4(float s) : m_v( _mm_cvttps_epi32( _mm_set1_ps ( s ) ) ) {}

	Col4( int r, int g, int b, int a ) : m_v( _mm_setr_epi32( r, g, b, a ) ) {}
	Col4( int r, int g, int b ) : m_v( _mm_setr_epi32( r, g, b, 0 ) ) {}
	Col4( int r, int g ) : m_v( _mm_setr_epi32( r, g, 0, 0 ) ) {}
	Col4( Col3::Arg v, int w ) : m_v( v.m_v ) { ((int *)&m_v)[3] = w; }

	explicit Col4( unsigned int s ) : m_v( _mm_set1_epi32( s ) ) {}
	explicit Col4( const unsigned int (&_rgba)[4] ) : m_v( _mm_load_si128( (const __m128i *)&_rgba ) ) {}
	explicit Col4( u8 const *source ) : m_v( _mm_loadu_si128( (const __m128i *)source ) ) {}

	Col3 GetCol3() const
	{
		return Col3( m_v );
	}

	int GetM8() const
	{
		return _mm_movemask_epi8 ( m_v );
	}

	int GetLong() const
	{
		return _mm_cvtsi128_si32 ( m_v );
	}

	Col4 SetLong( int v ) const
	{
		return Col4 ( _mm_cvtsi32_si128 ( v ) );
	}

	int R() const { return _mm_extract_epi16( m_v, 0 ); }
	int G() const { return _mm_extract_epi16( m_v, 2 ); }
	int B() const { return _mm_extract_epi16( m_v, 4 ); }
	int A() const { return _mm_extract_epi16( m_v, 6 ); }

	int &GetR() { return ((int *)&m_v)[0]; }
	int &GetG() { return ((int *)&m_v)[1]; }
	int &GetB() { return ((int *)&m_v)[2]; }
	int &GetA() { return ((int *)&m_v)[3]; }
	// let the compiler figure this one out, probably spills to memory
	int &GetO(int o) { return ((int *)&m_v)[o]; }

	const int &GetR() const { return ((const int *)&m_v)[0]; }
	const int &GetG() const { return ((const int *)&m_v)[1]; }
	const int &GetB() const { return ((const int *)&m_v)[2]; }
	const int &GetA() const { return ((const int *)&m_v)[3]; }
	// let the compiler figure this one out, probably spills to memory
	const int &GetO(int o) const { return ((const int *)&m_v)[o]; }

	Col4 SplatR() const { return Col4( _mm_shuffle_epi32( m_v, SQUISH_SSE_SPLAT( 0 ) ) ); }
	Col4 SplatG() const { return Col4( _mm_shuffle_epi32( m_v, SQUISH_SSE_SPLAT( 1 ) ) ); }
	Col4 SplatB() const { return Col4( _mm_shuffle_epi32( m_v, SQUISH_SSE_SPLAT( 2 ) ) ); }
	Col4 SplatA() const { return Col4( _mm_shuffle_epi32( m_v, SQUISH_SSE_SPLAT( 3 ) ) ); }

	template<const int inv>
	void SetRGBA( int r, int g, int b, int a )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, r, 0 );
		v = _mm_insert_epi16( v, g, 2 );
		v = _mm_insert_epi16( v, b, 4 );
		v = _mm_insert_epi16( v, a, 6 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		m_v = v;
	}

	template<const int inv>
	void SetRGBApow2( int r, int g, int b, int a )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, r, 0 );
		v = _mm_insert_epi16( v, g, 2 );
		v = _mm_insert_epi16( v, b, 4 );
		v = _mm_insert_epi16( v, a, 6 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		v = _mm_slli_epi32( v, 23 );
		v = _mm_add_epi32( v, _mm_castps_si128( _mm_set1_ps(1.0f) ) );

		m_v = _mm_cvttps_epi32( _mm_castsi128_ps( v ) );
	}

	template<const int inv>
	void SetRGBApow2( int c )
	{
		__m128i v = _mm_shuffle_epi32( _mm_cvtsi32_si128( c ), SQUISH_SSE_SPLAT( 0 ) );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		v = _mm_slli_epi32( v, 23 );
		v = _mm_add_epi32( v, _mm_castps_si128( _mm_set1_ps(1.0f) ) );

		m_v = _mm_cvttps_epi32( _mm_castsi128_ps( v ) );
	}

	Col4& operator&=( Arg v )
	{
		m_v = _mm_and_si128( m_v, v.m_v );
		return *this;
	}

	Col4& operator^=( Arg v )
	{
		m_v = _mm_xor_si128( m_v, v.m_v );
		return *this;
	}

	Col4& operator|=( Arg v )
	{
		m_v = _mm_or_si128( m_v, v.m_v );
		return *this;
	}

	Col4& operator>>=( const int n )
	{
		m_v = _mm_srli_epi32( m_v, n );
		return *this;
	}

	Col4& operator<<=( const int n )
	{
		m_v = _mm_slli_epi32( m_v, n );
		return *this;
	}

	Col4& operator+=( Arg v )
	{
		m_v = _mm_add_epi32( m_v, v.m_v );
		return *this;
	}

	Col4& operator-=( Arg v )
	{
		m_v = _mm_sub_epi32( m_v, v.m_v );
		return *this;
	}

	Col4& operator*=( Arg v )
	{
#if ( SQUISH_USE_SSE >= 4 )
		m_v = _mm_mullo_epi32( m_v, v.m_v );
#else
		m_v = _mm_mullo_epi16( m_v, v.m_v );
#endif
		return *this;
	}

	friend int operator<( Arg left, Arg right  )
	{
		return CompareFirstLessThan(left, right);
	}

	friend int operator>( Arg left, Arg right  )
	{
		return CompareFirstGreaterThan(left, right);
	}

	friend int operator==( Arg left, Arg right  )
	{
		return CompareFirstEqualTo(left, right);
	}

	friend Col4 operator~( Arg left )
	{
		return left ^ Col4(~0);
	}

	friend Col4 operator&( Arg left, Arg right  )
	{
		return Col4( _mm_and_si128( left.m_v, right.m_v ) );
	}

	friend Col4 operator%( Arg left, Arg right )
	{
		return Col4( _mm_andnot_si128( left.m_v, right.m_v ) );
	}

	friend Col4 operator^( Arg left, Arg right )
	{
		return Col4( _mm_xor_si128( left.m_v, right.m_v ) );
	}

	friend Col4 operator|( Arg left, Arg right )
	{
		return Col4( _mm_or_si128( left.m_v, right.m_v ) );
	}

	friend Col4 operator>>( Arg left, int right )
	{
		return Col4( _mm_srli_epi32( left.m_v, right ) );
	}

	friend Col4 operator<<( Arg left, int right )
	{
		return Col4( _mm_slli_epi32( left.m_v, right ) );
	}

	friend Col4 operator+( Arg left, Arg right )
	{
		return Col4( _mm_add_epi32( left.m_v, right.m_v ) );
	}

	friend Col4 operator-( Arg left, Arg right )
	{
		return Col4( _mm_sub_epi32( left.m_v, right.m_v ) );
	}

	friend Col4 operator*( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_mullo_epi32( left.m_v, right.m_v ) );
#else
		return Col4( _mm_mullo_epi16( left.m_v, right.m_v ) );
#endif
	}

	friend Col4 operator*( Arg left, int right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_mullo_epi32( left.m_v, _mm_set1_epi32( right ) ) );
#else
		return Col4( _mm_mullo_epi16( left.m_v, _mm_set1_epi32( right ) ) );
#endif
	}

	template<const int n>
	friend Col4 FillSign( Arg a );
	template<const int n>
	friend Col4 FillSign( Arg a )
	{
		return Col4( _mm_srai_epi32( _mm_slli_epi32( a.m_v, n ), n ) );
	}

	template<const int n>
	friend Col4 ExtendSign( Arg a );
	template<const int n>
	friend Col4 ExtendSign( Arg a )
	{
		return Col4( _mm_srai_epi32( a.m_v, n ) );
	}
	
	template<const int n>
	friend Col4 ShiftLeft( Arg a );
	template<const int n>
	friend Col4 ShiftLeft( Arg a )
	{
		if ((n) <= 0)
			return Col4( a.m_v );
		if ((n) <= 7)
			return Col4( _mm_slli_epi32( a.m_v, (n) & 7 ) );
		if ((n) & 7)
			return Col4( _mm_slli_epi32( _mm_slli_si128( a.m_v, (n) >> 3 ), (n) & 7 ) );

			return Col4( _mm_slli_si128( a.m_v, (n) >> 3 ) );
	}

	template<const int n>
	friend Col4 ShiftRight( Arg a );
	template<const int n>
	friend Col4 ShiftRight( Arg a )
	{
		if ((n) <= 0)
			return Col4( a.m_v );
		if ((n) <= 7)
			return Col4( _mm_srli_epi32( a.m_v, (n) & 7 ) );
		if ((n) & 7)
			return Col4( _mm_srli_epi32( _mm_srli_si128( a.m_v, (n) >> 3 ), (n) & 7 ) );

			return Col4( _mm_srli_si128( a.m_v, (n) >> 3 ) );
	}

	template<const int n>
	friend Col4 ShiftRightHalf( Arg a );
	template<const int n>
	friend Col4 ShiftRightHalf( Arg a )
	{
		return Col4( (n) > 0 ? _mm_srli_epi64( a.m_v, (n) ) : a.m_v );
	}

	friend Col4 ShiftRightHalf( Arg a, const int n )
	{
		return Col4( _mm_srl_epi64( a.m_v, _mm_cvtsi32_si128( n ) ) );
	}

	friend Col4 ShiftRightHalf( Arg a, Arg b )
	{
		return Col4( _mm_srl_epi64( a.m_v, b.m_v ) );
	}

	template<const int n>
	friend Col4 ShiftLeftHalf( Arg a );
	template<const int n>
	friend Col4 ShiftLeftHalf( Arg a )
	{
		return Col4( (n) > 0 ? _mm_slli_epi64( a.m_v, (n) ) : a.m_v );
	}

	friend Col4 ShiftLeftHalf( Arg a, const int n )
	{
		return Col4( _mm_sll_epi64( a.m_v, _mm_cvtsi32_si128( n ) ) );
	}

	template<const int r, const int g, const int b, const int a>
	friend Col4 ShiftLeftLo( Arg v );
	template<const int r, const int g, const int b, const int a>
	friend Col4 ShiftLeftLo( Arg v )
	{
		// (1 << r, 1 << g, 1 << b, 1 << a);
		Col4 p2; p2.SetRGBApow2<0>(r, g, b, a);

#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_mullo_epi32( v.m_v, p2.m_v ) );
#else
		return Col4( _mm_mullo_epi16( v.m_v, p2.m_v ) );
#endif
	}

	template<const int n, const int p>
	friend Col4 MaskBits( Arg a );
	template<const int n, const int p>
	friend Col4 MaskBits( Arg a )
	{
		if (((p) + (n)) <= 0)
			return Col4(0);
		if (((p) + (n)) >= 64)
			return a;

		// compile time
		__int64 base = ~(0xFFFFFFFFFFFFFFFFULL << (     ((p) + (n)) & 63));
	//	__int64 base =  (0xFFFFFFFFFFFFFFFFULL >> (64 - ((p) + (n)) & 63));
		__m128i mask = _mm_setr_epi32(
		  (int)(base >>  0),
		  (int)(base >> 32), 0, 0
		);

		return Col4( _mm_and_si128( a.m_v, mask ) );
	}

	friend Col4 MaskBits( Arg a, const int n, const int p )
	{
		const int val = 64 - ((p) + (n));

		__m128i shift = _mm_max_epi16( _mm_cvtsi32_si128( val ), _mm_set1_epi32( 0 ) );
		__m128i mask = _mm_setr_epi32(
		  0xFFFFFFFF,
		  0xFFFFFFFF, 0, 0
		);

		mask = _mm_srl_epi64( mask, shift );

		// (0xFFFFFFFFFFFFFFFFULL >> (64 - (p + n) & 63))
		return Col4( _mm_and_si128( a.m_v, mask ) );
	}

	template<const int n, const int p>
	friend Col4 CopyBits( Arg left, Arg right );
	template<const int n, const int p>
	friend Col4 CopyBits( Arg left, Arg right )
	{
		if (!(n))
			return left;
		if (!(p))
			return MaskBits<n, 0>(right);
		if (((p) + (n)) >= 64)
			return (left) + ShiftLeftHalf<p>(right);

#if ( SQUISH_USE_XSSE == 4 )
		return Col4( _mm_inserti_si64( left.m_v, right.m_v, n, p ) );
#else
		return MaskBits<p, 0>(left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
	//	return               (left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
#endif
	}

	friend Col4 CopyBits( Arg left, Col4& right, const int n, const int p )
	{
#if ( SQUISH_USE_XSSE == 4 )
		/* ---- ---bl xxxx xxxx */
		const int val = (p << 8) + (n << 0);

		right.m_v = _mm_unpacklo_epi64( right.m_v, _mm_cvtsi32_si128( val ) );
		return Col4( _mm_insert_si64( left.m_v, right.m_v ) );
#else
		return MaskBits(left, p, 0) + MaskBits(ShiftLeftHalf(right, p), n, p);
	//	return         (left      ) + MaskBits(ShiftLeftHalf(right, p), n, p);
#endif
	}

	template<const int n, const int p>
	friend Col4 KillBits( Arg a );
	template<const int n, const int p>
	friend Col4 KillBits( Arg a )
	{
		if (!n || (p >= 64))
			return a;
		if (!p && (n >= 64))
			return Col4(0);

		// compile time
		__int64 base1 =  (0xFFFFFFFFFFFFFFFFULL << (     (p + 0) & 63));
		__int64 base2 =  (0xFFFFFFFFFFFFFFFFULL >> (64 - (p + n) & 63));
	//	__int64 base1 = ~(0xFFFFFFFFFFFFFFFFULL >> (64 - (p + 0) & 63));
	//	__int64 base2 = ~(0xFFFFFFFFFFFFFFFFULL << (64 - (p + n) & 63));

		__m128i mask;

		if ((p + n) >= 64)
		  base2 = 0xFFFFFFFFFFFFFFFFULL;

		mask = _mm_setr_epi32(
		  (int)((base1 ^ base2) >>  0),
		  (int)((base1 ^ base2) >> 32), 0, 0
		);

		return Col4( _mm_and_si128( a.m_v, mask ) );
	}

	friend Col4 KillBits( Arg a, const int n, const int p )
	{
		const int val1 =      (p + 0);
		const int val2 = 64 - (p + n);

		__m128i shift1 = _mm_max_epi16( _mm_cvtsi32_si128( val1 ), _mm_set1_epi32( 0 ) );
		__m128i shift2 = _mm_max_epi16( _mm_cvtsi32_si128( val2 ), _mm_set1_epi32( 0 ) );
		__m128i mask1 = _mm_setr_epi32(
		  0xFFFFFFFF,
		  0xFFFFFFFF, 0, 0
		);
		__m128i mask2 = _mm_setr_epi32(
		  0xFFFFFFFF,
		  0xFFFFFFFF, 0, 0
		);

		mask1 = _mm_sll_epi64( mask1, shift1 );
		mask2 = _mm_srl_epi64( mask2, shift2 );

		return Col4( _mm_and_si128( a.m_v, _mm_xor_si128( mask1, mask2 ) ) );
	}

	template<const int n, const int p>
	friend Col4 InjtBits( Arg left, Arg right );
	template<const int n, const int p>
	friend Col4 InjtBits( Arg left, Arg right )
	{
		if (!n || (p >= 64))
			return right;
		if ((p + n) >= 64)
			return KillBits<n, p>(left) + ShiftLeftHalf<p>(right);
	//		return               (left) + ShiftLeftHalf<p>(right);


#if ( SQUISH_USE_XSSE == 4 )
		return Col4( _mm_inserti_si64( left.m_v, right.m_v, n, p ) );
#else
		return KillBits<n, p>(left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
	//	return               (left) + MaskBits<n, p>(ShiftLeftHalf<p>(right));
#endif
	}

	friend Col4 InjtBits( Arg left, Col4& right, const int n, const int p )
	{
#if ( SQUISH_USE_XSSE == 4 )
		/* ---- ---bl xxxx xxxx */
		const int val = (p << 8) + (n << 0);

		right.m_v = _mm_unpacklo_epi64( right.m_v, _mm_cvtsi32_si128( val ) );
		return Col4( _mm_insert_si64( left.m_v, right.m_v ) );
#else
		return KillBits(left, n, p) + MaskBits(ShiftLeftHalf(right, p), n, p);
	//	return         (left      ) + MaskBits(ShiftLeftHalf(right, p), n, p);
#endif
	}

	template<const int n, const int p>
	friend Col4 ExtrBits( Arg a );
	template<const int n, const int p>
	friend Col4 ExtrBits( Arg a )
	{
		if (!n)
			return Col4(0);
		if (!p)
			return MaskBits<n, 0>(a);
		if ((n + p) >= 64)
			return ShiftRightHalf<p>(a);

#if ( SQUISH_USE_XSSE == 4 )
		return Col4( _mm_extracti_si64( a.m_v, n, p ) );
#else
		return MaskBits<n, 0>(ShiftRightHalf<p>(a));
#endif
	}

	friend Col4 ExtrBits( Arg a, const int n, const int p )
	{
#if ( SQUISH_USE_XSSE == 4 )
		/* ---- ----- ---- ---bl */
		const int val = (p << 8) + (n << 0);

		return Col4( _mm_extract_si64( a.m_v, _mm_cvtsi32_si128( val ) ) );
#else
		return MaskBits(ShiftRightHalf(a, p), n, 0);
#endif
	}

	template<const int n, const int p>
	friend void ExtrBits( Arg left, Col4 &right );
	template<const int n, const int p>
	friend void ExtrBits( Arg left, Col4 &right )
	{
		right  = ExtrBits<n, p>( left );
	}

	template<const int n, const int p>
	friend void ConcBits( Arg left, Col4 &right );
	template<const int n, const int p>
	friend void ConcBits( Arg left, Col4 &right )
	{
		right  = ShiftLeft<32>( right );
		if (n > 0)
			right += ExtrBits<n, p>( left );
	}

	template<const int n, const int p>
	friend void ReplBits( Arg left, Col4 &right );
	template<const int n, const int p>
	friend void ReplBits( Arg left, Col4 &right )
	{
		if (!n)
			return;
		if ((n < 0)) {
			right  = ExtrBits<-n, p>( left );
			right.m_v = _mm_shuffle_epi32( right.m_v, SQUISH_SSE_SHUF( 0, 0, 0, 3 ) );
		}
		else {
			right  = ExtrBits< n, p>( left );
			right.m_v = _mm_shuffle_epi32( right.m_v, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );
		}
	}

	friend Col4 RevsBits( Col4::Arg v )
	{
		// reverse bits:
		// http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith32Bits
		// b = ((b * 0x0802LU & 0x22110LU) | (b * 0x8020LU & 0x88440LU)) * 0x10101LU >> 16;
		__m128i res = v.m_v;

		res = _mm_and_si128( res, _mm_set1_epi32( 0x000000FF ) );
	//	res = _mm_unpacklo_epi64( res, res );
		res = _mm_mul_epu32( res, _mm_setr_epi32( 0x00802, 0x0, 0x08020, 0x0 ) );
		res = _mm_and_si128( res, _mm_setr_epi32( 0x22110, 0x0, 0x88440, 0x0 ) );
		res = _mm_or_si128( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP64() ) );
		res = _mm_mul_epu32( res, _mm_setr_epi32( 0x10101, 0x0, 0x10101, 0x0 ) );
		res = _mm_srli_epi32( res, 16 );

		return Col4( res );
	}

	friend Col4 Mul16x16( Col4::Arg a, Col4::Arg b )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_mullo_epi32( a.m_v, b.m_v ) );
#else
		__m128i lo = _mm_mullo_epi16( a.m_v, b.m_v );
		__m128i hi = _mm_mulhi_epu16( a.m_v, b.m_v );

		return Col4( _mm_or_si128( lo, _mm_slli_si128( hi, 2 ) ) );
#endif
	}

	friend Col4 Div32x16( Col4::Arg a, Col4::Arg b )
	{
		return Col4(
			(int)a.GetR() / (int)b.GetR(),
			(int)a.GetG() / (int)b.GetG(),
			(int)a.GetB() / (int)b.GetB(),
			(int)a.GetA() / (int)b.GetA()
		);
	}

	//! Returns a*b + c
	friend Col4 MultiplyAdd( Arg a, Arg b, Arg c )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_add_epi32( _mm_mullo_epi32( a.m_v, b.m_v ), c.m_v ) );
#else
		return Col4( _mm_add_epi32( _mm_mullo_epi16( a.m_v, b.m_v ), c.m_v ) );
#endif
	}

	//! Returns -( a*b - c )
	friend Col4 NegativeMultiplySubtract( Arg a, Arg b, Arg c )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_sub_epi32( c.m_v, _mm_mullo_epi32( a.m_v, b.m_v ) ) );
#else
		return Col4( _mm_sub_epi32( c.m_v, _mm_mullo_epi16( a.m_v, b.m_v ) ) );
#endif
	}

	template<const int f, const int t>
	friend Col4 Shuffle( Arg a );
	template<const int f, const int t>
	friend Col4 Shuffle( Arg a )
	{
		if (f == t)
			return a;

		return Col4( _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SHUF(
			(t == 0 ? f : 0),
			(t == 1 ? f : 1),
			(t == 2 ? f : 2),
			(t == 3 ? f : 3)
		) ) );
	}

	template<const int f, const int t>
	friend Col4 Exchange( Arg a );
	template<const int f, const int t>
	friend Col4 Exchange( Arg a )
	{
		if (f == t)
			return a;

		return Col4( _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SHUF(
			(t == 0 ? f : (f == 0 ? t : 0)),
			(t == 1 ? f : (f == 1 ? t : 1)),
			(t == 2 ? f : (f == 2 ? t : 2)),
			(t == 3 ? f : (f == 3 ? t : 3))
		) ) );
	}

	friend Col4 HorizontalAdd( Arg a )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128i res = _mm_hadd_epi32( a.m_v, a.m_v );
		return Col4( _mm_hadd_epi32( res, res ) );
#else
		__m128i res = a.m_v;

		res = _mm_add_epi32( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP64() ) );
		res = _mm_add_epi32( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP32() ) );

		return Col4( res );
#endif
	}

	friend Col4 HorizontalAdd( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128i resc;

		resc = _mm_hadd_epi32( a.m_v, b.m_v );
		resc = _mm_hadd_epi32( resc, resc );
		resc = _mm_hadd_epi32( resc, resc );

		return Col4( resc );
#else
		__m128i resa = a.m_v;
		__m128i resb = b.m_v;
		__m128i resc;

		resc = _mm_add_epi32( resa, resb );
		resc = _mm_add_epi32( resc, _mm_shuffle_epi32( resc, SQUISH_SSE_SWAP64() ) );
		resc = _mm_add_epi32( resc, _mm_shuffle_epi32( resc, SQUISH_SSE_SWAP32() ) );

		return Col4( resc );
#endif
	}

	friend Col4 HorizontalAddTiny( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 ) && 0
		__m128 res = _mm_castsi128_ps ( a.m_v );

		res = _mm_dp_ps( res, _mm_set1_ps ( 1.0f ), 0xFF );

		return Col4( _mm_castps_si128 ( res ) );
#elif ( SQUISH_USE_SSE >= 3 )
		__m128 res = _mm_castsi128_ps ( a.m_v );

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		res = _mm_hadd_ps( res, res );
		res = _mm_hadd_ps( res, res );

		return Col4( _mm_castps_si128 ( res ) );
#else
		return HorizontalAdd( a );
#endif
	}

	friend Col4 HorizontalAddTiny( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128 resa = _mm_castsi128_ps ( a.m_v );
		__m128 resb = _mm_castsi128_ps ( b.m_v );
		__m128 resc;

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		resc = _mm_hadd_ps( resa, resb );
		resc = _mm_hadd_ps( resc, resc );
		resc = _mm_hadd_ps( resc, resc );

		return Col4( _mm_castps_si128 ( resc ) );
#else
		return HorizontalAdd( a, b );
#endif
	}

	friend Col4 HorizontalMaxTiny( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 )
		__m128i inv;
		__m128i res;

		inv = _mm_xor_si128( a.m_v, _mm_set1_epi32 ( ~0 ) );
		inv = _mm_minpos_epu16 ( inv );
		res = _mm_xor_si128(  inv , _mm_set1_epi32 ( ~0 ) );
		res = _mm_and_si128(  inv , _mm_set1_epi32 ( 0x0000FFFF ) );

		return Col4( res ).SplatR();
#else
		__m128 res = _mm_castsi128_ps ( a.m_v );

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP64() ) );
		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP32() ) );

		return Col4( _mm_castps_si128 ( res ) );
#endif
	}

	friend Col4 HorizontalMinTiny( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 )
		__m128i res;

		res = _mm_minpos_epu16 ( a.m_v );
		res = _mm_and_si128( res , _mm_set1_epi32 ( 0x0000FFFF ) );

		return Col4( res ).SplatR();
#else
		__m128 res = _mm_castsi128_ps ( a.m_v );

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP64() ) );
		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP32() ) );

		return Col4( _mm_castps_si128 ( res ) );
#endif
	}

	friend Col4 Dot( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
	  	return HorizontalAdd( Col4( _mm_mullo_epi32( left.m_v, right.m_v ) ) );
#else
		return HorizontalAdd( Col4( _mm_mullo_epi16( left.m_v, right.m_v ) ) );
#endif
	}

	friend Col4 DotTiny( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return HorizontalAddTiny( Col4( _mm_mullo_epi32( left.m_v, right.m_v ) ) );
#else
		return HorizontalAddTiny( Col4( _mm_mullo_epi16( left.m_v, right.m_v ) ) );
#endif
	}

	friend Col4 Abs( Col4::Arg a )
	{
		__m128i sign = _mm_srai_epi32( a.m_v, 31 );

		return Col4( _mm_sub_epi32( _mm_xor_si128( a.m_v, sign ), sign ) );
	}

	friend Col4 Min( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_min_epi32( left.m_v, right.m_v ) );
#else
		return Col4( _mm_min_epi16( left.m_v, right.m_v ) );
#endif
	}

	friend Col4 MinTiny( Arg left, Arg right )
	{
		__m128 resa = _mm_castsi128_ps( left.m_v );
		__m128 resb = _mm_castsi128_ps( right.m_v );
		__m128 resc;

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		resc = _mm_min_ps( resa, resb );

		return Col4( _mm_castps_si128 ( resc ) );
	}

	friend Col4 Max( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
	  	return Col4( _mm_max_epi32( left.m_v, right.m_v ) );
#else
		return Col4( _mm_max_epi16( left.m_v, right.m_v ) );
#endif
	}
	
	friend Col4 MaxTiny( Arg left, Arg right )
	{
		__m128 resa = _mm_castsi128_ps( left.m_v );
		__m128 resb = _mm_castsi128_ps( right.m_v );
		__m128 resc;

		// relies on correct de-normal floating-point treatment
		// doesn't support negative values
		resc = _mm_max_ps( resa, resb );

		return Col4( _mm_castps_si128 ( resc ) );
	}

#if 0
	friend Col4 AbsoluteDifference( Arg left, Arg right )
	{
		return Max( left, right ) - Min( left, right );
	}

	friend Col4 SummedAbsoluteDifference( Arg left, Arg right )
	{
		return HorizontalAdd( AbsoluteDifference( left, right ) );
	}

	friend Col4 MaximumAbsoluteDifference( Arg left, Arg right )
	{
		return HorizontalMax( AbsoluteDifference( left, right ) );
	}
#endif

	friend bool CompareFirstLessThan( Arg left, Arg right )
	{
		__m128i bits = _mm_cmplt_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return !!(value & 0x000F);
	}

	friend bool CompareFirstGreaterThan( Arg left, Arg right )
	{
		__m128i bits = _mm_cmpgt_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return !!(value & 0x000F);
	}

	friend bool CompareFirstEqualTo( Arg left, Arg right )
	{
		__m128i bits = _mm_cmpeq_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return !!(value & 0x000F);
	}

	friend bool CompareAnyLessThan( Arg left, Arg right )
	{
		__m128i bits = _mm_cmplt_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return value != 0x0000;
	}

	friend bool CompareAllEqualTo( Arg left, Arg right )
	{
		__m128i bits = _mm_cmpeq_epi32( left.m_v, right.m_v );
		int value = _mm_movemask_epi8( bits );
		return value == 0xFFFF;
	}

	friend Col4 CompareAllGreaterThan_M4( Arg left, Arg right )
	{
		return Col4( _mm_cmpgt_epi32( left.m_v, right.m_v ) );
	}

	friend Col4 CompareAllLessThan_M4( Arg left, Arg right )
	{
		return Col4( _mm_cmplt_epi32( left.m_v, right.m_v ) );
	}

	friend Col4 CompareAllEqualTo_M4( Arg left, Arg right )
	{
		return Col4( _mm_cmpeq_epi32( left.m_v, right.m_v ) );
	}

	friend Col4 CompareAllLessThan_M8( Arg left, Arg right )
	{
		return Col4( _mm_cmplt_epi8( left.m_v, right.m_v ) );
	}
	
	friend Col4 CompareAllEqualTo_M8( Arg left, Arg right )
	{
		return Col4( _mm_cmpeq_epi8( left.m_v, right.m_v ) );
	}

	friend Col4 IsNotZero( Arg v )
	{
		return Col4( _mm_cmpgt_epi32( v.m_v, _mm_setzero_si128( ) ) );
	}

	friend Col4 IsZero( Arg v )
	{
		return Col4( _mm_cmpeq_epi32( v.m_v, _mm_setzero_si128( ) ) );
	}

	friend Col4 IsOne( Arg v )
	{
		return Col4( _mm_cmpeq_epi32( v.m_v, _mm_set1_epi32( 0x000000FF ) ) );
	}

	template<const int value>
	friend Col4 IsValue( Arg v );
	template<const int value>
	friend Col4 IsValue( Arg v )
	{
		return Col4( _mm_cmpeq_epi32( v.m_v, _mm_set1_epi32( value ) ) );
	}

	friend Col4 TransferA( Arg left, Arg right )
	{
		__m128i l = _mm_and_si128( left.m_v , _mm_setr_epi32( ~0, ~0, ~0,  0 ) );
		__m128i r = _mm_and_si128( right.m_v, _mm_setr_epi32(  0,  0,  0, ~0 ) );

		return Col4( _mm_or_si128( l, r ) );
	}

	friend Col4 KillA( Arg left )
	{
		return Col4( _mm_or_si128( left.m_v, _mm_setr_epi32( 0x00, 0x00, 0x00, 0xFF ) ) );
	}
	
	friend Col4 CollapseA( Arg r, Arg g, Arg b, Arg a )
	{
		return Col4( _mm_packus_epi16(
			_mm_packs_epi32( _mm_srli_epi32( r.m_v, 24 ), _mm_srli_epi32( g.m_v, 24 ) ),
			_mm_packs_epi32( _mm_srli_epi32( b.m_v, 24 ), _mm_srli_epi32( a.m_v, 24 ) )
		) );
	}

	friend void PackBytes( Arg a, unsigned int &loc )
	{
		__m128i

		r = _mm_packs_epi32( a.m_v, a.m_v );
		r = _mm_packus_epi16( r, r );

		loc = _mm_cvtsi128_si32 ( r );
	}
	
	friend void PackBytes( Arg a, int &loc )
	{
		__m128i

		r = _mm_packs_epi32( a.m_v, a.m_v );
		r = _mm_packs_epi16( r, r );

		loc = _mm_cvtsi128_si32 ( r );
	}
	
	friend void PackWords( Arg a, unsigned__int64 &loc )
	{
		__m128i

#if ( SQUISH_USE_SSE >= 4 )
		r = _mm_packus_epi32( a.m_v, a.m_v );
#else
		// fix-up/down
		r = _mm_sub_epi32( a.m_v, _mm_set1_epi32( 32768 ) );
		r = _mm_packs_epi32( r, r );
		r = _mm_add_epi16( r, _mm_set1_epi16( (short)-32768 ) );
#endif

//		loc = _mm_cvtsi128_si64( r );
		_mm_storel_epi64( (__m128i *)&loc, r );
	}
	
	friend void PackWords( Arg a, __int64 &loc )
	{
		__m128i
		  
		r = _mm_packs_epi32( a.m_v, a.m_v );

//		loc = _mm_cvtsi128_si64( r );
		_mm_storel_epi64( (__m128i *)&loc, r );
	}

	friend Col4 PackSShorts( Arg a )
	{
		return Col4( _mm_packs_epi32( a.m_v, a.m_v ) );
	}

	friend Col4 PackUShorts( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Col4( _mm_packus_epi32( a.m_v, a.m_v ) );
#else
		// fix-up/down
		__m128i
		r = _mm_sub_epi32( a.m_v, _mm_set1_epi32( 32768 ) );
		r = _mm_packs_epi32( r, r );
		r = _mm_add_epi16( r, _mm_set1_epi16( (short)-32768 ) );

		return Col4(r);
#endif
	}

	friend void UnpackBytes( Col4 &a, const unsigned int &loc )
	{
		__m128i

		r = _mm_cvtsi32_si128 ( loc );
		r = _mm_unpacklo_epi8( r, _mm_setzero_si128() );
		r = _mm_unpacklo_epi16( r, _mm_setzero_si128() );

		a = Col4( r );
	}
	
	friend void UnpackBytes( Col4 &a, const int &loc )
	{
		__m128i

		r = _mm_cvtsi32_si128 ( loc );
		r = _mm_unpacklo_epi8( r, r );
		r = _mm_unpacklo_epi16( r, r );
		
		a = ExtendSign<24>( Col4( r ) );
	}
	
	friend void UnpackWords( Col4 &a, const unsigned__int64 &loc )
	{
		__m128i

		r = _mm_loadl_epi64( (__m128i *)&loc );
		r = _mm_unpacklo_epi16( r, _mm_setzero_si128() );

		a = Col4( r );
	}
	
	friend void UnpackWords( Col4 &a, const __int64 &loc )
	{
		__m128i

		r = _mm_loadl_epi64( (__m128i *)&loc );
		r = _mm_unpacklo_epi16( r, r );
		
		a = ExtendSign<16>( Col4( r ) );
	}
	
	// clamp the output to [0, 1]
	Col4 Clamp() const {
		Col4 const one (0xFF);
		Col4 const zero(0x00);

		return Min(one, Max(zero, *this));
	}

	friend void Interleave( Col4 &a, Arg b, Arg c )
	{
		a = Col4( _mm_shuffle_epi32( _mm_unpacklo_epi32( b.m_v , c.m_v ), SQUISH_SSE_SHUF(0, 3, 0, 3) ) );
	}

	friend void LoadAligned( Col4 &a, Col4 &b, Arg c )
	{
	        a.m_v = c.m_v;
		b.m_v = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void LoadAligned( Col4 &a, void const *source )
	{
		a.m_v = _mm_load_si128( (__m128i const *)source );
	}

	friend void LoadAligned( Col4 &a, Col4 &b, void const *source )
	{
		a.m_v = _mm_load_si128( (__m128i const *)source );
		b.m_v = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void LoadUnaligned( Col4 &a, Col4 &b, void const *source )
	{
		a.m_v = _mm_loadu_si128( (__m128i const *)source );
		b.m_v = _mm_shuffle_epi32( a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void StoreAligned( Arg a, Arg b, Col4 &c )
	{
		c.m_v = _mm_unpacklo_epi64( a.m_v, b.m_v );
	}

	friend void StoreAligned( Arg a, void *destination )
	{
		_mm_store_si128( (__m128i *)destination, a.m_v );
	}

	friend void StoreAligned( Arg a, Arg b, void *destination )
	{
		_mm_store_si128( (__m128i *)destination, _mm_unpacklo_epi64( a.m_v, b.m_v ) );
	}

	friend void StoreUnaligned( Arg a, void *destination )
	{
		_mm_storeu_si128( (__m128i *)destination, a.m_v );
	}

	friend void StoreUnaligned( Arg a, Arg b, void *destination )
	{
		_mm_storeu_si128( (__m128i *)destination, _mm_unpacklo_epi64( a.m_v, b.m_v ) );
	}
	
	friend void StoreUnaligned( Arg a, u8* loc ) {
	  PackBytes( a, (unsigned int&) (*((unsigned int *)loc)) ); }
	friend void StoreUnaligned( Arg a, u16* loc ) {
	  PackWords( a, (unsigned__int64&) (*((unsigned__int64 *)loc)) ); }
	friend void StoreUnaligned( Arg a, s8* loc ) {
	  PackBytes( a, (int&) (*((int *)loc)) ); }
	friend void StoreUnaligned( Arg a, s16* loc ) {
	  PackWords( a, (__int64&) (*((__int64 *)loc)) ); }
	
	friend void LoadUnaligned( Col4 &a, const u8* loc ) {
	  UnpackBytes( a, (const unsigned int&) (*((const unsigned int *)loc)) ); }
	friend void LoadUnaligned( Col4 &a, const u16* loc ) {
	  UnpackWords( a, (const unsigned__int64&) (*((const unsigned__int64 *)loc)) ); }
	friend void LoadUnaligned( Col4 &a, const s8* loc ) {
	  UnpackBytes( a, (const int&) (*((const int *)loc)) ); }
	friend void LoadUnaligned( Col4 &a, const s16* loc ) {
	  UnpackWords( a, (const __int64&) (*((const __int64 *)loc)) ); }

	void SwapRGBA( Col4 &with )
	{
	  /* inplace swap based on xors */
	       m_v = _mm_xor_si128( m_v, with.m_v );
	  with.m_v = _mm_xor_si128( with.m_v, m_v );
	       m_v = _mm_xor_si128( m_v, with.m_v );
	}

private:
	__m128i m_v;

	friend class Vec4;
	friend class Col8;
};

#if	!defined(SQUISH_USE_PRE)
inline Col3 LengthSquared( Col3::Arg v )
{
  return Dot( v, v );
}

inline Col3 LengthSquaredTiny( Col3::Arg v )
{
  return DotTiny( v, v );
}

inline Col4 LengthSquared( Col4::Arg v )
{
  return Dot( v, v );
}

inline Col4 LengthSquaredTiny( Col4::Arg v )
{
  return DotTiny( v, v );
}
#endif

class Col8
{
public:
	typedef Col8 const& Arg;

	Col8() {}

	explicit Col8( __m128i v ) : m_v( v ) {}

	Col8( Arg arg ) : m_v( arg.m_v ) {}

	Col8& operator=( Arg arg )
	{
		m_v = arg.m_v;
		return *this;
	}

	explicit Col8(Col4::Arg s) : m_v( s.m_v ) {
		m_v = _mm_packs_epi32( m_v, m_v );
	}

	explicit Col8(int s) : m_v( _mm_set1_epi16( (short)s ) ) {}
	explicit Col8(u16 s) : m_v( _mm_set1_epi16( s ) ) {}
	explicit Col8(s16 s) : m_v( _mm_set1_epi16( s ) ) {}
	explicit Col8(u8 s) : m_v( _mm_set1_epi16( s ) ) {}
	explicit Col8(s8 s) : m_v( _mm_set1_epi16( s ) ) {}

	Col8( int a, int b, int c, int d, int e, int f, int g, int h )
	  : m_v( _mm_setr_epi16( (short)a, (short)b, (short)c, (short)d,
				 (short)e, (short)f, (short)g, (short)h ) ) {}
	Col8( u16 a, u16 b, u16 c, u16 d, u16 e, u16 f, u16 g, u16 h )
	  : m_v( _mm_setr_epi16( a, b, c, d, e, f, g, h ) ) {}
	Col8( s16 a, s16 b, s16 c, s16 d, s16 e, s16 f, s16 g, s16 h )
	  : m_v( _mm_setr_epi16( a, b, c, d, e, f, g, h ) ) {}
	Col8( u8 a, u8 b, u8 c, u8 d, u8 e, u8 f, u8 g, u8 h )
	  : m_v( _mm_setr_epi16( a, b, c, d, e, f, g, h ) ) {}
	Col8( s8 a, s8 b, s8 c, s8 d, s8 e, s8 f, s8 g, s8 h )
	  : m_v( _mm_setr_epi16( a, b, c, d, e, f, g, h ) ) {}

	int Get0() const
	{
		return _mm_extract_epi16( m_v, 0 );
	}
	
#pragma warning ( push )
#pragma warning ( disable : 4100 )
	friend Col4 LoCol4(Arg v, const unsigned dummy)
	{
		return Col4( _mm_unpacklo_epi16( v.m_v, _mm_setzero_si128() ) );
	}
	
	friend Col4 HiCol4(Arg v, const unsigned dummy)
	{
		return Col4( _mm_unpackhi_epi16( v.m_v, _mm_setzero_si128() ) );
	}
	
	friend Col4 LoCol4(Arg v, const signed dummy)
	{
		return Col4( _mm_srai_epi32( _mm_unpacklo_epi16( _mm_setzero_si128(), v.m_v ), 16 ) );
	}
	
	friend Col4 HiCol4(Arg v, const signed dummy)
	{
		return Col4( _mm_srai_epi32( _mm_unpackhi_epi16( _mm_setzero_si128(), v.m_v ), 16 ) );
	}
#pragma warning ( pop )
	
	const u16 &operator[]( int pos ) const
	{
		return ((u16 *)&m_v)[pos];
	//	return m_v.m128i_u16[pos];
	}

	Col8& operator*=( Arg v )
	{
		m_v = _mm_mullo_epi16( m_v, v.m_v );
		return *this;
	}

	friend Col8 operator>>( Arg left, unsigned int right )
	{
		return Col8( _mm_srli_epi16( left.m_v, right ) );
	}
	
	friend Col8 operator>>( Arg left, int right )
	{
		return Col8( _mm_srai_epi16( left.m_v, right ) );
	}

	friend Col8 operator<<( Arg left, unsigned int right )
	{
		return Col8( _mm_slli_epi16( left.m_v, right ) );
	}
	
	friend Col8 operator<<( Arg left, int right )
	{
		return Col8( _mm_slli_epi16( left.m_v, right ) );
	}

	friend Col8 operator+( Arg left, Arg right )
	{
		return Col8( _mm_add_epi16( left.m_v, right.m_v ) );
	}

	friend Col8 operator-( Arg left, Arg right )
	{
		return Col8( _mm_sub_epi16( left.m_v, right.m_v ) );
	}

	friend Col8 operator*( Arg left, Arg right )
	{
		return Col8( _mm_mullo_epi16( left.m_v, right.m_v ) );
	}

	friend Col8 operator*( Arg left, unsigned int right )
	{
		return Col8( _mm_mulhi_epu16( left.m_v, _mm_set1_epi16( (unsigned short)right ) ) );
	}
	
	friend Col8 operator*( Arg left, int right )
	{
		return Col8( _mm_mulhi_epi16( left.m_v, _mm_set1_epi16( (short)right ) ) );
	}

	template<const int n>
	friend Col8 ExtendSign(Arg a);
	template<const int n>
	friend Col8 ExtendSign(Arg a)
	{
		return Col8( _mm_srai_epi16( a.m_v, n ) );
	}
	
	friend Col8 HorizontalMin(Arg a)
	{
		__m128i res = a.m_v;

#if ( SQUISH_USE_SSE >= 4 )
		res = _mm_min_epu16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP64() ) );
		res = _mm_min_epu16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP32() ) );
		res = _mm_min_epu16( res, _mm_shuffle_epi16( res, SQUISH_SSE_SWAP16() ) );
#else
		res = _mm_sub_epi16( res, _mm_set1_epi16( (short)-32768 ) );
		res = _mm_min_epi16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP64() ) );
		res = _mm_min_epi16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP32() ) );
		res = _mm_min_epi16( res, _mm_shuffle_epi16( res, SQUISH_SSE_SWAP16() ) );
		res = _mm_add_epi16( res, _mm_set1_epi16( (short)-32768 ) );
#endif

		return Col8( res );
	}

	friend Col8 HorizontalMax(Arg a)
	{
		__m128i res = a.m_v;

#if ( SQUISH_USE_SSE >= 4 )
		res = _mm_max_epu16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP64() ) );
		res = _mm_max_epu16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP32() ) );
		res = _mm_max_epu16( res, _mm_shuffle_epi16( res, SQUISH_SSE_SWAP16() ) );
#else
		res = _mm_sub_epi16( res, _mm_set1_epi16( (short)-32768 ) );
		res = _mm_max_epi16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP64() ) );
		res = _mm_max_epi16( res, _mm_shuffle_epi32( res, SQUISH_SSE_SWAP32() ) );
		res = _mm_max_epi16( res, _mm_shuffle_epi16( res, SQUISH_SSE_SWAP16() ) );
		res = _mm_add_epi16( res, _mm_set1_epi16( (short)-32768 ) );
#endif

		return Col8( res );
	}

	template<const int n>
	friend Col8 ShiftUp(Arg a);
	template<const int n>
	friend Col8 ShiftUp(Arg a)
	{
		return Col8( _mm_slli_si128( a.m_v, n << 1 ) );
	}
	
#pragma warning ( push )
#pragma warning ( disable : 4100 )
	friend Col4 ExpandUpper(Arg a, const unsigned dummy) {
		__m128i res = a.m_v;
		
		res = _mm_unpackhi_epi16( res, _mm_setzero_si128() );

#ifdef _MSV_VER
		assert(res.m128i_u32[0] == a.m_v.m128i_u16[4]);
		assert(res.m128i_u32[1] == a.m_v.m128i_u16[5]);
		assert(res.m128i_u32[2] == a.m_v.m128i_u16[6]);
		assert(res.m128i_u32[3] == a.m_v.m128i_u16[7]);
#endif

		return Col4( res );
	}

	friend Col4 RepeatUpper(Arg a, const unsigned dummy) {
		__m128i res = a.m_v;
		
		res = _mm_unpackhi_epi16( res, _mm_setzero_si128() );
		res = _mm_shuffle_epi32( res, SQUISH_SSE_SPLAT(3) );

#ifdef _MSV_VER
		assert(res.m128i_u32[0] == a.m_v.m128i_u16[7]);
		assert(res.m128i_u32[1] == a.m_v.m128i_u16[7]);
		assert(res.m128i_u32[2] == a.m_v.m128i_u16[7]);
		assert(res.m128i_u32[3] == a.m_v.m128i_u16[7]);
#endif

		return Col4( res );
	}
	
	friend Col4 InterleaveUpper(Arg a, Arg b, const unsigned dummy) {
		__m128i res;
		
		res = _mm_unpackhi_epi16( a.m_v, b.m_v );
		res = _mm_unpackhi_epi16( res, _mm_setzero_si128() );
		res = _mm_unpackhi_epi64( res, res );

#ifdef _MSV_VER
		assert(res.m128i_u32[0] == a.m_v.m128i_u16[7]);
		assert(res.m128i_u32[1] == b.m_v.m128i_u16[7]);
		assert(res.m128i_u32[2] == a.m_v.m128i_u16[7]);
		assert(res.m128i_u32[3] == b.m_v.m128i_u16[7]);
#endif

		return Col4( res );
	}

	friend Col4 ReplicateUpper(Arg a, Arg b, const unsigned dummy) {
		__m128i res;
		
		res = _mm_unpackhi_epi16( a.m_v, b.m_v );
		res = _mm_unpackhi_epi16( res, _mm_setzero_si128() );
		res = _mm_unpackhi_epi32( res, res );

#ifdef _MSV_VER
		assert(res.m128i_u32[0] == a.m_v.m128i_u16[7]);
		assert(res.m128i_u32[1] == a.m_v.m128i_u16[7]);
		assert(res.m128i_u32[2] == b.m_v.m128i_u16[7]);
		assert(res.m128i_u32[3] == b.m_v.m128i_u16[7]);
#endif

		return Col4( res );
	}

	friend Col4 ExpandUpper(Arg a, const signed dummy) {
		__m128i res = a.m_v;
		
		res = _mm_unpackhi_epi16( res, res );
		res = _mm_srai_epi32( res, 16 );

#ifdef _MSV_VER
		assert(res.m128i_i32[0] == a.m_v.m128i_i16[4]);
		assert(res.m128i_i32[1] == a.m_v.m128i_i16[5]);
		assert(res.m128i_i32[2] == a.m_v.m128i_i16[6]);
		assert(res.m128i_i32[3] == a.m_v.m128i_i16[7]);
#endif

		return Col4( res );
	}

	friend Col4 RepeatUpper(Arg a, const signed dummy) {
		__m128i res = a.m_v;

		res = _mm_srai_epi32( res, 16 );
		res = _mm_shuffle_epi32( res, SQUISH_SSE_SPLAT(3) );

#ifdef _MSV_VER
		assert(res.m128i_i32[0] == a.m_v.m128i_i16[7]);
		assert(res.m128i_i32[1] == a.m_v.m128i_i16[7]);
		assert(res.m128i_i32[2] == a.m_v.m128i_i16[7]);
		assert(res.m128i_i32[3] == a.m_v.m128i_i16[7]);
#endif

		return Col4( res );
	}
	
	friend Col4 InterleaveUpper(Arg a, Arg b, const signed dummy) {
		__m128i res;
		
		res = _mm_unpackhi_epi32( a.m_v, b.m_v );
		res = _mm_srai_epi32( res, 16 );
		res = _mm_unpackhi_epi64( res, res );

#ifdef _MSV_VER
		assert(res.m128i_i32[0] == a.m_v.m128i_i16[7]);
		assert(res.m128i_i32[1] == b.m_v.m128i_i16[7]);
		assert(res.m128i_i32[2] == a.m_v.m128i_i16[7]);
		assert(res.m128i_i32[3] == b.m_v.m128i_i16[7]);
#endif

		return Col4( res );
	}

	friend Col4 ReplicateUpper(Arg a, Arg b, const signed dummy) {
		__m128i res;
		
		res = _mm_unpackhi_epi32( a.m_v, b.m_v );
		res = _mm_srai_epi32( res, 16 );
		res = _mm_unpackhi_epi32( res, res );
		
#ifdef _MSV_VER
		assert(res.m128i_i32[0] == a.m_v.m128i_i16[7]);
		assert(res.m128i_i32[1] == a.m_v.m128i_i16[7]);
		assert(res.m128i_i32[2] == b.m_v.m128i_i16[7]);
		assert(res.m128i_i32[3] == b.m_v.m128i_i16[7]);
#endif

		return Col4( res );
	}
#pragma warning ( pop )
	
	/*
	friend Col4 Expand(Arg a, int ia) {
		__m128i res = _mm_setzero_si128();

		res = _mm_insert_epi16( res, a.m_v.m128i_u16[ia - 0], 0 );
		res = _mm_insert_epi16( res, a.m_v.m128i_u16[ia - 1], 2 );
		res = _mm_insert_epi16( res, a.m_v.m128i_u16[ia - 2], 4 );
		res = _mm_insert_epi16( res, a.m_v.m128i_u16[ia - 3], 6 );

		return Col4( res );
	}

	friend Col4 Repeat(Arg a, int ia) {
		__m128i res = _mm_setzero_si128();

		res = _mm_insert_epi16( res, a.m_v.m128i_u16[ia], 0 );
		res = _mm_shuffle_epi32( res, SQUISH_SSE_SPLAT(0) );

		return Col4( res );
	}

	friend Col4 Interleave(Arg a, Arg b, int ia, int ib) {
		__m128i res = _mm_setzero_si128();

		res = _mm_insert_epi16( res, a.m_v.m128i_u16[ia], 0 );
		res = _mm_insert_epi16( res, b.m_v.m128i_u16[ib], 2 );
		res = _mm_unpacklo_epi64( res, res );

		return Col4( res );
	}

	friend Col4 Replicate(Arg a, Arg b, int ia, int ib) {
		__m128i res = _mm_setzero_si128();

		res = _mm_insert_epi16( res, a.m_v.m128i_u16[ia], 0 );
		res = _mm_insert_epi16( res, b.m_v.m128i_u16[ib], 2 );
		res = _mm_unpacklo_epi32( res, res );

		return Col4( res );
	}
	*/
	
	friend int CompareEqualTo( Arg left, Arg right )
	{
		return _mm_movemask_epi8( _mm_cmpeq_epi16( left.m_v, right.m_v ) );
	}
	
	friend Col8 CompareAllEqualTo( Arg left, Arg right )
	{
		return Col8( _mm_cmpeq_epi16( left.m_v, right.m_v ) );
	}
	
	friend Col8 CompareAllLessThan( Arg left, Arg right )
	{
		return Col8( _mm_cmplt_epi16( left.m_v, right.m_v ) );
	}

private:
	__m128i m_v;

	friend class Vec4;
};

#define VEC4_CONST( X ) Vec4( X )

class Vec3
{
public:
	typedef Vec3 const& Arg;

	Vec3() {}

	explicit Vec3( __m128 v ) : m_v( v ) {}

	Vec3( Arg arg ) : m_v( arg.m_v ) {}

	Vec3& operator=( Arg arg )
	{
		m_v = arg.m_v;
		return *this;
	}

	explicit Vec3(float s) : m_v( _mm_set1_ps( s ) ) {}
	explicit Vec3(int s) : m_v( _mm_cvtepi32_ps( _mm_set1_epi32 ( s ) ) ) {}

	Vec3( const float* x, const float* y, const float* z ) {
		m_v = _mm_unpacklo_ps(_mm_load_ss(x), _mm_load_ss(y));
		m_v = _mm_movelh_ps(m_v, _mm_load_ss(z));
	}
	
	Vec3( bool x, bool y, bool z ) : m_v( _mm_castsi128_ps( _mm_setr_epi32( x ? ~0 : 0, y ? ~0 : 0, z ? ~0 : 0, 0 ) ) ) {}

	Vec3( float x, float y, float z ) : m_v( _mm_setr_ps( x, y, z, 0.0f ) ) {}
	Vec3( float x, float y ) : m_v( _mm_setr_ps( x, y, 0.0f, 0.0f ) ) {}
	Vec3( Arg x, Arg y, Arg z ) : m_v( _mm_unpacklo_ps( _mm_unpacklo_ps( x.m_v, z.m_v ), y.m_v ) ) {}
	Vec3( Arg x, Arg y ) : m_v( _mm_unpacklo_ps( _mm_unpacklo_ps( x.m_v, y.m_v ), _mm_set1_ps( 0.0f ) ) ) {}

	Vec3( Col3::Arg c ) : m_v( _mm_cvtepi32_ps( c.m_v ) ) {}

	void StoreX(float *x) const { _mm_store_ss(x, m_v); }
	void StoreY(float *y) const { _mm_store_ss(y, _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 1 ) )); }
	void StoreZ(float *z) const { _mm_store_ss(z, _mm_movehl_ps( m_v, m_v ) ); }
	
	float X() const { return ((float *)&m_v)[0]; }
	float Y() const { return ((float *)&m_v)[1]; }
	float Z() const { return ((float *)&m_v)[2]; }

	float &GetX() { return ((float *)&m_v)[0]; }
	float &GetY() { return ((float *)&m_v)[1]; }
	float &GetZ() { return ((float *)&m_v)[2]; }
	// let the compiler figure this one out, probably spills to memory
	float &GetO(int o) { return ((float *)&m_v)[o]; }

	Vec3 SplatX() const { return Vec3( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 0 ) ) ); }
	Vec3 SplatY() const { return Vec3( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 1 ) ) ); }
	Vec3 SplatZ() const { return Vec3( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 2 ) ) ); }

	template<const int inv>
	void SetXYZ( int x, int y, int z )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, x, 0 );
		v = _mm_insert_epi16( v, y, 2 );
		v = _mm_insert_epi16( v, z, 4 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		m_v = _mm_cvtepi32_ps( v );
	}

	template<const int inv>
	void SetXYZpow2( int x, int y, int z )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, x, 0 );
		v = _mm_insert_epi16( v, y, 2 );
		v = _mm_insert_epi16( v, z, 4 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		v = _mm_slli_epi32( v, 23 );
		v = _mm_add_epi32( v, _mm_castps_si128( _mm_set1_ps(1.0f) ) );

		m_v = _mm_castsi128_ps( v );
	}

	Vec3& operator+=( Arg v )
	{
		m_v = _mm_add_ps( m_v, v.m_v );
		return *this;
	}

	Vec3& operator-=( Arg v )
	{
		m_v = _mm_sub_ps( m_v, v.m_v );
		return *this;
	}

	Vec3& operator*=( Arg v )
	{
		m_v = _mm_mul_ps( m_v, v.m_v );
		return *this;
	}
	
	Vec3& operator/=( Arg v )
	{
		*this *= Reciprocal( v );
		return *this;
	}

	Vec3& operator/=( float v )
	{
		*this *= Reciprocal( Vec3( v ) );
		return *this;
	}

	Vec3& operator/=( int v )
	{
		*this *= Reciprocal( Vec3( v ) );
		return *this;
	}

	friend int operator!( Arg left )
	{
		return CompareFirstEqualTo(left, Vec3(0.0f));
	}

	friend int operator<( Arg left, Arg right )
	{
		return CompareFirstLessThan(left, right);
	}

	friend int operator>( Arg left, Arg right )
	{
		return CompareFirstGreaterThan(left, right);
	}

	friend int operator==( Arg left, Arg right )
	{
		return CompareFirstEqualTo(left, right);
	}

	friend Vec3 operator&( Arg left, Arg right )
	{
		return Vec3( _mm_and_ps( left.m_v, right.m_v ) );
	}

	friend Vec3 operator%( Arg left, Arg right )
	{
		return Vec3( _mm_andnot_ps( left.m_v, right.m_v ) );
	}

	friend Vec3 operator+( Arg left, Arg right )
	{
		return Vec3( _mm_add_ps( left.m_v, right.m_v ) );
	}

	friend Vec3 operator-( Arg left, Arg right )
	{
		return Vec3( _mm_sub_ps( left.m_v, right.m_v ) );
	}

	friend Vec3 operator*( Arg left, Arg right )
	{
		return Vec3( _mm_mul_ps( left.m_v, right.m_v ) );
	}

	friend Vec3 operator*( Arg left, float right )
	{
		return Vec3( _mm_mul_ps( left.m_v, _mm_set1_ps( right ) ) );
	}

	friend Vec3 operator*( float left, Arg right )
	{
		return Vec3( _mm_mul_ps( _mm_set1_ps( left ), right.m_v ) );
	}

	friend Vec3 operator/( Arg left, float right )
	{
		return left * Reciprocal( Vec3( right ) );
	}

	friend Vec3 operator*( Arg left, int right )
	{
#if ( SQUISH_USE_SSE == 1 )
		...
#else
		return Vec3( _mm_mul_ps( left.m_v, _mm_cvtepi32_ps( _mm_set1_epi32( right ) ) ) );
#endif
	}

	//! Returns a*b + c
	friend Vec3 MultiplyAdd( Arg a, Arg b, Arg c )
	{
		return Vec3( _mm_add_ps( _mm_mul_ps( a.m_v, b.m_v ), c.m_v ) );
	}

	//! Returns -( a*b - c )
	friend Vec3 NegativeMultiplySubtract( Arg a, Arg b, Arg c )
	{
		return Vec3( _mm_sub_ps( c.m_v, _mm_mul_ps( a.m_v, b.m_v ) ) );
	}

	template<const int f, const int t>
	friend Vec3 Shuffle( Arg a );
	template<const int f, const int t>
	friend Vec3 Shuffle( Arg a )
	{
		if (f == t)
			return a;

		return Vec3( _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( a.m_v ), SQUISH_SSE_SHUF(
			(t == 0 ? f : 0),
			(t == 1 ? f : 1),
			(t == 2 ? f : 2),
			(t == 3 ? f : 3)
		) ) ) );
	}

	template<const int f, const int t>
	friend Vec3 Exchange( Arg a );
	template<const int f, const int t>
	friend Vec3 Exchange( Arg a )
	{
		if (f == t)
			return a;

		return Vec3( _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( a.m_v ), SQUISH_SSE_SHUF(
			(t == 0 ? f : (f == 0 ? t : 0)),
			(t == 1 ? f : (f == 1 ? t : 1)),
			(t == 2 ? f : (f == 2 ? t : 2)),
			(t == 3 ? f : (f == 3 ? t : 3))
		) ) ) );
	}

	template<const int n>
	friend Vec3 RotateLeft( Arg a );
	template<const int n>
	friend Vec3 RotateLeft( Arg a )
	{
		return Vec3( _mm_shuffle_ps( a.m_v , a.m_v , SQUISH_SSE_SHUF(
			(n + 0) % 3,
			(n + 1) % 3,
			(n + 2) % 3,
			3
		) ) );
	}

	friend Vec3 HorizontalAdd( Arg a )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128 res = a.m_v;

		res = _mm_and_ps( res , _mm_castsi128_ps( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
		res = _mm_hadd_ps( res, res );
		res = _mm_hadd_ps( res, res );

		return Vec3( res );
#else
		__m128 res = a.m_v;

		res = _mm_and_ps( res , _mm_castsi128_ps( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
		res = _mm_add_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP64() ) );
		res = _mm_add_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP32() ) );

		return Vec3( res );
#endif
	}

	friend Vec3 HorizontalAdd( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128 resc;

		resc = _mm_hadd_ps( a.m_v, b.m_v );
		resc = _mm_and_ps( resc , _mm_castsi128_ps( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
		resc = _mm_hadd_ps( resc, resc );
		resc = _mm_hadd_ps( resc, resc );

		return Vec3( resc );
#else
		__m128 resc;

		resc = _mm_add_ps( a.m_v, b.m_v );
		resc = _mm_and_ps( resc , _mm_castsi128_ps( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
		resc = _mm_add_ps( resc, _mm_shuffle_ps( resc, resc, SQUISH_SSE_SWAP64() ) );
		resc = _mm_add_ps( resc, _mm_shuffle_ps( resc, resc, SQUISH_SSE_SWAP32() ) );

		return Vec3( resc );
#endif
	}

	friend Vec3 Select( Arg a, Arg b, Arg c )
	{
#if 0
		__m128 res;
		__m128 bits = _mm_cmpeq_ps( b.m_v, c.m_v );
		int mask = _mm_movemask_ps( bits );

		/* (1 >> 1) = 0
		 * (2 >> 1) = 1
		 * (4 >> 1) = 2
		mask = (mask & 7) >> 1;
		mask = (mask) * ((1 << 0) + (1 << 2) + (1 << 4) + (1 << 6));
		 */

		/**/ if (mask & 1)
		  res = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );
		else if (mask & 2)
		  res = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SHUF( 1, 1, 1, 1 ) );
		else
		  res = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SHUF( 2, 2, 2, 2 ) );

		return Vec3( res );
#else
		// branch free, and no CPU<->SSEunit transfer
		__m128 mask = _mm_cmpeq_ps( b.m_v, c.m_v );
		__m128 res = _mm_and_ps( a.m_v, mask );

		__m128 r0 = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );
		__m128 r1 = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 1, 1, 1 ) );
		__m128 r2 = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 2, 2, 2, 2 ) );

		res = _mm_or_ps( _mm_or_ps( r0, r1 ), r2 );

		return Vec3(res);
#endif
	}

	friend Vec3 HorizontalMin( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 2, 0, 1, 3 ) ) );
		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 2, 0, 3 ) ) );

		return Vec3( res );
	}

	friend Vec3 HorizontalMax( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 2, 0, 1, 3 ) ) );
		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 2, 0, 3 ) ) );

		return Vec3( res );
	}
	
	friend Vec3 HorizontalMaxXY( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_max_ps(
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 1, 0, 1 ) ),
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 0, 1, 0 ) )
		);

		return Vec3( res );
	}
	
	friend Vec3 HorizontalMinXY( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_min_ps(
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 1, 0, 1 ) ),
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 0, 1, 0 ) )
		);

		return Vec3( res );
	}

	friend Vec3 Reciprocal( Arg v )
	{
		// get the reciprocal estimate
		__m128 estimate = _mm_rcp_ps( v.m_v );

		// one round of Newton-Rhaphson refinement
		__m128 diff = _mm_sub_ps( _mm_set1_ps( 1.0f ), _mm_mul_ps( estimate, v.m_v ) );
		return Vec3( _mm_add_ps( _mm_mul_ps( diff, estimate ), estimate ) );
	}

	friend Vec3 ReciprocalSqrt( Arg v )
	{
		// get the reciprocal estimate
		__m128 estimate = _mm_rsqrt_ps( v.m_v );

		// one round of Newton-Rhaphson refinement
		__m128 diff = _mm_sub_ps( _mm_set1_ps( 3.0f ), _mm_mul_ps( estimate, _mm_mul_ps( estimate, v.m_v ) ) );
		return Vec3( _mm_mul_ps( _mm_mul_ps( diff, _mm_set1_ps( 0.5f ) ), estimate ) );
	}

	friend Vec3 Sqrt( Arg v )
	{
		return Vec3( _mm_sqrt_ps( v.m_v ) );
	}

	friend Vec3 Length( Arg left )
	{
		Vec3 sum = HorizontalAdd( Vec3( _mm_mul_ps( left.m_v, left.m_v ) ) );
		Vec3 sqt = Vec3( _mm_sqrt_ps( sum.m_v ) );

		return sqt;
	}

	friend Vec3 ReciprocalLength( Arg left )
	{
		Vec3 sum = HorizontalAdd( Vec3( _mm_mul_ps( left.m_v, left.m_v ) ) );
		Vec3 rsq = ReciprocalSqrt(sum);

		return rsq;
	}

	friend Vec3 Normalize( Arg left )
	{
		return left * ReciprocalLength(left);
	}

	friend Vec3 Normalize( Vec3 &x, Vec3 &y, Vec3 &z )
	{
		Vec3 xx = x * x;
		Vec3 yy = y * y;
		Vec3 zz = z * z;

		Vec3 sum = xx + yy + zz;
		Vec3 rsq = ReciprocalSqrt(sum);

		x = x * rsq;
		y = y * rsq;
		z = z * rsq;

		return rsq;
	}

	template<const bool disarm>
	friend Vec3 Complement( Arg left );
	template<const bool disarm>
	friend Vec3 Complement( Arg left )
	{
		__m128 ren, res, rez;

		ren = left.m_v;
		rez = _mm_set1_ps( 1.0f );
		res = _mm_mul_ps( left.m_v, left.m_v );
#if ( SQUISH_USE_SSE >= 3 )
		res = _mm_hadd_ps( res, res );
#else
		res = _mm_add_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 0, 1, 0 ) ) );
#endif
		if (!disarm) {
			// correct x² + y² > 1.0f by renormalization
			if ( _mm_comigt_ss( res, rez ) ) {
				res = ReciprocalSqrt( Vec3(res) ).m_v;
				res = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );

				ren = _mm_mul_ps( ren, res );
				res = rez;
			}
		}
		
		rez = _mm_sub_ps( rez, _mm_min_ps( rez, res ) );
		rez = _mm_sqrt_ps( rez );
		res = _mm_movelh_ps( left.m_v, rez );

		// sqrt(1.0f - (x*x + y*y))
		return Vec3( res );
	}

	template<const bool disarm>
	friend Vec3 Complement( Vec3 &left, Vec3 &right );
	template<const bool disarm>
	friend Vec3 Complement( Vec3 &left, Vec3 &right )
	{
		if (!disarm) {
			Vec3 len = (left * left) + (right * right);
			Vec3 adj = ReciprocalSqrt(Max(Vec3(1.0f), len));

			// correct x² + y² > 1.0f by renormalization
			left  *= adj;
			right *= adj;

			// sqrt(1.0f - (x² + y²))
			return Sqrt(Vec3(1.0f) - Min(Vec3(1.0f), len));
		}
		else {
			Vec3 len = (left * left) + (right * right);

			// disarm x² + y² > 1.0f by letting NaN happen
			// ...

			// sqrt(1.0f - (x² + y²))
			return Sqrt(Vec3(1.0f) - len);
		}
	}

	template<const bool disarm>
	friend Vec3 ComplementPyramidal( Vec3 &left );
	template<const bool disarm>
	friend Vec3 ComplementPyramidal( Vec3 &left )
	{
		// 1 - max(abs(l, r)) == 1 + min(-abs(l, r))
		Vec3 res = HorizontalMinXY(Neg(left)) + Vec3(1.0f);

		res = Normalize(TransferZ(left, res));

		return res;
	}

	template<const bool disarm>
	friend Vec3 ComplementPyramidal( Vec3 &left, Vec3 &right );
	template<const bool disarm>
	friend Vec3 ComplementPyramidal( Vec3 &left, Vec3 &right )
	{
		// 1 - max(abs(l, r)) == 1 + min(-abs(l, r))
		Vec3 res = Min(Neg(left), Neg(right)) + Vec3(1.0f);

		Normalize(left, right, res);

		return res;
	}

	friend Vec3 Dot( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Vec3( _mm_dp_ps ( left.m_v, right.m_v, 0x77 ) );
#else
		return HorizontalAdd( Vec3( _mm_mul_ps( left.m_v, right.m_v ) ) );
#endif
	}

	friend void Dot( Arg left, Arg right, float *r )
	{
		Vec3 res = Dot( left, right );

		_mm_store_ss( r, res.m_v );
	}

	friend Vec3 Abs( Arg a )
	{
		return Vec3( _mm_and_ps( a.m_v, _mm_castsi128_ps( _mm_set1_epi32( 0x7FFFFFFF ) ) ) );
	}
	
	friend Vec3 Neg( Arg a )
	{
		return Vec3( _mm_or_ps( a.m_v, _mm_castsi128_ps( _mm_set1_epi32( 0x80000000 ) ) ) );
	}

	friend Vec3 Min( Arg left, Arg right )
	{
		return Vec3( _mm_min_ps( left.m_v, right.m_v ) );
	}

	friend Vec3 Max( Arg left, Arg right )
	{
		return Vec3( _mm_max_ps( left.m_v, right.m_v ) );
	}

	// clamp the output to [0, 1]
	Vec3 Clamp() const {
		Vec3 const one (1.0f);
		Vec3 const zero(0.0f);

		return Min(one, Max(zero, *this));
	}

	template<const bool round>
	friend Col3 FloatToInt( Arg v );
	template<const bool round>
	friend Col3 FloatToInt( Arg v )
	{
#if ( SQUISH_USE_SSE == 1 )
		...
#else
		// use SSE2 instructions
		if (round)
		      return Col3( _mm_cvtps_epi32( v.m_v ) );
		else
		      return Col3( _mm_cvttps_epi32( v.m_v ) );
#endif
	}

	friend Vec3 Truncate( Arg v )
	{
#if ( SQUISH_USE_SSE == 1 )
		// convert to ints
		__m128 input = v.m_v;
		__m64 lo = _mm_cvttps_pi32( input );
		__m64 hi = _mm_cvttps_pi32( _mm_movehl_ps( input, input ) );

		// convert to floats
		__m128 part = _mm_movelh_ps( input, _mm_cvtpi32_ps( input, hi ) );
		__m128 truncated = _mm_cvtpi32_ps( part, lo );

		// clear out the MMX multimedia state to allow FP calls later
		_mm_empty();
		return Vec3( truncated );
#else
		// use SSE2 instructions
		return Vec3( _mm_cvtepi32_ps( _mm_cvttps_epi32( v.m_v ) ) );
#endif
	}

	friend Vec3 AbsoluteDifference( Arg left, Arg right )
	{
		__m128 diff = _mm_sub_ps( left.m_v, right.m_v );
		diff = _mm_and_ps( diff, _mm_castsi128_ps( _mm_set1_epi32( 0x7FFFFFFF ) ) );
		return Vec3( diff );
	}

	friend Vec3 SummedAbsoluteDifference( Arg left, Arg right )
	{
		return HorizontalAdd( AbsoluteDifference( left, right ) );
	}

	friend Vec3 MaximumAbsoluteDifference( Arg left, Arg right )
	{
		return HorizontalMax( AbsoluteDifference( left, right ) );
	}

	friend bool CompareAnyLessThan( Arg left, Arg right )
	{
		__m128 bits = _mm_cmplt_ps( left.m_v, right.m_v );
		int value = _mm_movemask_ps( bits );
		return (value & 0x7) != 0x0;
	}

	friend bool CompareAnyGreaterThan( Arg left, Arg right )
	{
		__m128 bits = _mm_cmpgt_ps( left.m_v, right.m_v );
		int value = _mm_movemask_ps( bits );
		return (value & 0x7) != 0x0;
	}

	friend bool CompareAllEqualTo( Arg left, Arg right )
	{
		__m128 bits = _mm_cmpeq_ps( left.m_v, right.m_v );
		int value = _mm_movemask_ps( bits );
		return (value & 0x7) == 0x7;
	}

	friend Col3 CompareAllEqualTo_M8( Arg left, Arg right )
	{
		return Col3( _mm_cmpeq_epi8( _mm_castps_si128 ( left.m_v ), _mm_castps_si128 ( right.m_v ) ) );
	}

	friend int CompareFirstLessThan( Arg left, Arg right )
	{
		return _mm_comilt_ss( left.m_v, right.m_v );
	}

	friend int CompareFirstGreaterThan( Arg left, Arg right )
	{
		return _mm_comigt_ss( left.m_v, right.m_v );
	}

	friend int CompareFirstEqualTo( Arg left, Arg right )
	{
		return _mm_comieq_ss( left.m_v, right.m_v );
	}

	Vec3 IsOne( ) const
	{
		return Vec3( _mm_cmpeq_ps( m_v, _mm_set1_ps( 1.0f ) ) );
	}

	Vec3 IsNotOne( ) const
	{
		return Vec3( _mm_cmpneq_ps( m_v, _mm_set1_ps( 1.0f ) ) );
	}
	
	friend Vec3 TransferZ( Arg left, Arg right )
	{
		return Vec3( _mm_shuffle_ps( left.m_v, right.m_v, SQUISH_SSE_SHUF( 0, 1, 2, 3 ) ) );
	}

	void SwapXYZ( Vec3 &with )
	{
		/* inplace swap based on xors */
		     m_v = _mm_xor_ps( m_v, with.m_v );
		with.m_v = _mm_xor_ps( with.m_v, m_v );
		     m_v = _mm_xor_ps( m_v, with.m_v );
	}

	friend void LoadAligned( Vec3 &a, Vec3 &b, Arg c )
	{
	        a.m_v = c.m_v;
		b.m_v = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SWAP64() );
		b.m_v = _mm_and_ps( b.m_v , _mm_castsi128_ps( _mm_setr_epi32( ~0, 0, 0, 0 ) ) );
	}

	friend void LoadAligned( Vec3 &a, void const *source )
	{
		a.m_v = _mm_load_ps( (float const *)source );
		a.m_v = _mm_and_ps( a.m_v , _mm_castsi128_ps( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
	}

	friend void LoadUnaligned( Vec3 &a, void const *source )
	{
		a.m_v = _mm_loadu_ps( (float const *)source );
		a.m_v = _mm_and_ps( a.m_v , _mm_castsi128_ps( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) );
	}

	friend void LoadAligned( Vec3 &a, Vec3 &b, void const *source )
	{
		a.m_v = _mm_load_ps( (float const *)source );
		b.m_v = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SWAP64() );
		b.m_v = _mm_and_ps( b.m_v , _mm_castsi128_ps( _mm_setr_epi32( ~0, 0, 0, 0 ) ) );
	}

	friend void LoadUnaligned( Vec3 &a, Vec3 &b, void const *source )
	{
		a.m_v = _mm_loadu_ps( (float const *)source );
		b.m_v = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SWAP64() );
		b.m_v = _mm_and_ps( b.m_v , _mm_castsi128_ps( _mm_setr_epi32( ~0, 0, 0, 0 ) ) );
	}

	friend void StoreUnaligned( Arg a, void *destination )
	{
		_mm_storeu_ps( (float *)destination, a.m_v );
	}

private:
	__m128 m_v;

	friend class Vec4;
};

template<const bool round>
Col3 FloatToUHalf( Vec3::Arg v );
template<const bool round>
Col3 FloatToUHalf( Vec3::Arg v )
{
	Col3 h;

	h.GetR() = FloatToUHalf(v.X());
	h.GetG() = FloatToUHalf(v.Y());
	h.GetB() = FloatToUHalf(v.Z());

	return h;
}

template<const bool round>
Col3 FloatToSHalf( Vec3::Arg v );
template<const bool round>
Col3 FloatToSHalf( Vec3::Arg v )
{
	Col3 h;

	h.GetR() = FloatToSHalf(v.X());
	h.GetG() = FloatToSHalf(v.Y());
	h.GetB() = FloatToSHalf(v.Z());

	return h;
}

Vec3 UHalfToFloat( Col3::Arg v )
{
	Vec3 f;

	f.GetX() = UHalfToFloat((u16)v.R());
	f.GetY() = UHalfToFloat((u16)v.G());
	f.GetZ() = UHalfToFloat((u16)v.B());

	return f;
}

Vec3 SHalfToFloat( Col3::Arg v )
{
	Vec3 f;

	f.GetX() = SHalfToFloat((u16)v.R());
	f.GetY() = SHalfToFloat((u16)v.G());
	f.GetZ() = SHalfToFloat((u16)v.B());

	return f;
}

class Vec4
{
public:
	typedef Vec4 const& Arg;

	Vec4() {}

	explicit Vec4( __m128 v ) : m_v( v ) {}

	Vec4( Vec4 const& arg ) : m_v( arg.m_v ) {}
	Vec4( Vec3 const& arg ) : m_v( arg.m_v ) {}

	Vec4& operator=( Vec4 const& arg )
	{
		m_v = arg.m_v;
		return *this;
	}

	Vec4& operator=( Vec3 const& arg )
	{
		m_v = arg.m_v;
		return *this;
	}
	
	operator Vec3()
	{
		return Vec3(m_v);
	}

	explicit Vec4(float s) : m_v( _mm_set1_ps( s ) ) {}
	explicit Vec4(int   s) : m_v( _mm_cvtepi32_ps( _mm_set1_epi32 ( s ) ) ) {}

	Vec4( const float* x, const float* y, const float* z, const float* w ) {
		__m128 m_w;

		m_v = _mm_unpacklo_ps(_mm_load_ss(x), _mm_load_ss(y));
		m_w = _mm_unpacklo_ps(_mm_load_ss(z), _mm_load_ss(w));
		m_v = _mm_movelh_ps(m_v, m_w);
	}

	Vec4( const float* x, const float* y, const float* z ) {
		m_v = _mm_unpacklo_ps(_mm_load_ss(x), _mm_load_ss(y));
		m_v = _mm_movelh_ps(m_v, _mm_load_ss(z));
	}

	Vec4( const float* x, const float* y ) {
		m_v = _mm_unpacklo_ps(_mm_load_ss(x), _mm_load_ss(y));
		m_v = _mm_movelh_ps(m_v, _mm_set1_ps( 0.0f ));
	}

	Vec4( const float* x ) {
		m_v = _mm_load_ss(x);
		m_v = _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 0 ) );
	}
	
	Vec4( const unsigned short* x ) {
		__m128i v = _mm_setzero_si128();

		m_v = _mm_cvtepi32_ps( _mm_insert_epi16( v, *x, 0 ) );
		m_v = _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 0 ) );
	}
	
	Vec4( const signed short* x ) {
		__m128i v = _mm_setzero_si128();

		m_v = _mm_cvtepi32_ps( _mm_insert_epi16( v, *x, 0 ) );
		m_v = _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 0 ) );
	}
	
	Vec4( bool x, bool y, bool z, bool w ) : m_v( _mm_castsi128_ps( _mm_setr_epi32( x ? ~0 : 0, y ? ~0 : 0, z ? ~0 : 0, w ? ~0 : 0 ) ) ) {}

	Vec4( int x, int y, int z, int w ) : m_v( _mm_cvtepi32_ps( _mm_setr_epi32( x, y, z, w ) ) ) {}
	Vec4( int x, int y, int z ) : m_v( _mm_cvtepi32_ps( _mm_setr_epi32( x, y, z, 0 ) ) ) {}
	Vec4( int x, int y ) : m_v( _mm_cvtepi32_ps( _mm_setr_epi32( x, y, 0, 0 ) ) ) {}

	Vec4( float x, float y, float z, float w ) : m_v( _mm_setr_ps( x, y, z, w ) ) {}
	Vec4( float x, float y, float z ) : m_v( _mm_setr_ps( x, y, z, 0.0f ) ) {}
	Vec4( float x, float y ) : m_v( _mm_setr_ps( x, y, 0.0f, 0.0f ) ) {}

	Vec4( Arg x, Arg y, Arg z, Arg w ) : m_v( _mm_unpacklo_ps( _mm_unpacklo_ps( x.m_v, z.m_v ), _mm_unpacklo_ps( y.m_v, w.m_v ) ) ) {}
	Vec4( Arg x, Arg y, Arg z ) : m_v( _mm_unpacklo_ps( _mm_unpacklo_ps( x.m_v, z.m_v ), _mm_unpacklo_ps( y.m_v,  _mm_set1_ps( 0.0f ) ) ) ) {}
	Vec4( Arg x, Arg y ) : m_v( _mm_movelh_ps( _mm_unpacklo_ps( x.m_v, y.m_v ), _mm_set1_ps( 0.0f ) ) ) {}

	Vec4( Arg x, Arg y, Arg z, bool w ) : m_v( _mm_unpacklo_ps( _mm_unpacklo_ps( x.m_v, z.m_v ), y.m_v ) ) {w=w;}
	Vec4( Arg x, Arg y, bool z, bool w ) : m_v( _mm_unpacklo_ps( x.m_v, y.m_v ) ) {z=z;w=w;}

	Vec4( Vec3::Arg v, float w ) : m_v( v.m_v ) { m_v = _mm_or_ps( _mm_and_ps( m_v, _mm_castsi128_ps( _mm_setr_epi32( ~0, ~0, ~0, 0 ) ) ),  _mm_setr_ps( 0.0f, 0.0f, 0.0f, w ) ); }

	Vec4( Col4::Arg c ) : m_v( _mm_cvtepi32_ps( c.m_v ) ) {}

	Vec3 GetVec3() const
	{
		return Vec3( m_v );
	}
	
	int GetM4() const
	{
		return _mm_movemask_ps( m_v );
	}

	template<class dtyp> friend Vec4 LoVec4(Col8 const&v, const dtyp& dummy);
	template<class dtyp> friend Vec4 LoVec4(Col8 const&v, const dtyp& dummy)
	{
		return Vec4( LoCol4( v, dummy ) );
	}

	template<class dtyp> friend Vec4 HiVec4(Col8 const&v, const dtyp& dummy);
	template<class dtyp> friend Vec4 HiVec4(Col8 const&v, const dtyp& dummy)
	{
		return Vec4( HiCol4( v, dummy ) );
	}

	void StoreX(float *x) const { _mm_store_ss(x, m_v); }
	void StoreY(float *y) const { _mm_store_ss(y, _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 1 ) )); }
	void StoreZ(float *z) const { _mm_store_ss(z, _mm_movehl_ps( m_v, m_v ) ); }
	void StoreW(float *w) const { _mm_store_ss(w, _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 3 ) )); }

	float X() const { return ((float *)&m_v)[0]; }
	float Y() const { return ((float *)&m_v)[1]; }
	float Z() const { return ((float *)&m_v)[2]; }
	float W() const { return ((float *)&m_v)[3]; }

	float &GetX() { return ((float *)&m_v)[0]; }
	float &GetY() { return ((float *)&m_v)[1]; }
	float &GetZ() { return ((float *)&m_v)[2]; }
	float &GetW() { return ((float *)&m_v)[3]; }
	// let the compiler figure this one out, probably spills to memory
	float &GetO(int o) { return ((float *)&m_v)[o]; }

	Vec4 Swap  () const { return Vec4( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SWAP64() ) ); }
	Vec4 SplatX() const { return Vec4( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 0 ) ) ); }
	Vec4 SplatY() const { return Vec4( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 1 ) ) ); }
	Vec4 SplatZ() const { return Vec4( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 2 ) ) ); }
	Vec4 SplatW() const { return Vec4( _mm_shuffle_ps( m_v, m_v, SQUISH_SSE_SPLAT( 3 ) ) ); }

	template<const int inv>
	void SetXYZW( int x, int y, int z, int w )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, x, 0 );
		v = _mm_insert_epi16( v, y, 2 );
		v = _mm_insert_epi16( v, z, 4 );
		v = _mm_insert_epi16( v, w, 6 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		m_v = _mm_cvtepi32_ps( v );
	}

	template<const int inv>
	void SetXYZWpow2( int x, int y, int z, int w )
	{
		__m128i v = _mm_setzero_si128();

		v = _mm_insert_epi16( v, x, 0 );
		v = _mm_insert_epi16( v, y, 2 );
		v = _mm_insert_epi16( v, z, 4 );
		v = _mm_insert_epi16( v, w, 6 );

		if (inv) {
			v = _mm_sub_epi32( _mm_set1_epi32( inv ), v );
		}

		v = _mm_slli_epi32( v, 23 );
		v = _mm_add_epi32( v, _mm_castps_si128( _mm_set1_ps(1.0f) ) );

		m_v = _mm_castsi128_ps( v );
	}

	template<const int p>
	void Set( const float val )
	{
		__m128 mask = _mm_castsi128_ps ( _mm_setr_epi32(
		  p != 0 ? ~0 : 0,
		  p != 1 ? ~0 : 0,
		  p != 2 ? ~0 : 0,
		  p != 3 ? ~0 : 0
		) );

		__m128 fill = _mm_setr_ps(
		  p != 0 ? 0.0f : val,
		  p != 1 ? 0.0f : val,
		  p != 2 ? 0.0f : val,
		  p != 3 ? 0.0f : val
		);

		m_v = _mm_or_ps( _mm_and_ps( m_v, mask ), fill );
	}

	Vec4& operator&=( Arg v )
	{
		m_v = _mm_and_ps( m_v, v.m_v );
		return *this;
	}

	Vec4& operator+=( Arg v )
	{
		m_v = _mm_add_ps( m_v, v.m_v );
		return *this;
	}

	Vec4& operator-=( Arg v )
	{
		m_v = _mm_sub_ps( m_v, v.m_v );
		return *this;
	}

	Vec4& operator*=( Arg v )
	{
		m_v = _mm_mul_ps( m_v, v.m_v );
		return *this;
	}
	
	Vec4& operator*=( float v )
	{
		m_v = _mm_mul_ps( m_v, Vec4( v ).m_v );
		return *this;
	}

	Vec4& operator/=( Arg v )
	{
		*this *= Reciprocal( v );
		return *this;
	}
	
	Vec4& operator/=( float v )
	{
		*this *= Reciprocal( Vec4( v ) );
		return *this;
	}

	Vec4& operator/=( int v )
	{
		*this *= Reciprocal( Vec4( v ) );
		return *this;
	}

	friend int operator!( Arg left )
	{
		return CompareFirstEqualTo(left, Vec4(0.0f));
	}

	friend int operator<( Arg left, Arg right )
	{
		return CompareFirstLessThan(left, right);
	}

	friend int operator>( Arg left, Arg right )
	{
		return CompareFirstGreaterThan(left, right);
	}

	friend int operator>=( Arg left, Arg right )
	{
		return CompareFirstGreaterEqualTo(left, right);
	}

	friend int operator==( Arg left, Arg right )
	{
		return CompareFirstEqualTo(left, right);
	}

	friend Vec4 operator&( Arg left, Arg right )
	{
		return Vec4( _mm_and_ps( left.m_v, right.m_v ) );
	}

	friend Vec4 operator%( Arg left, Arg right )
	{
		return Vec4( _mm_andnot_ps( left.m_v, right.m_v ) );
	}

	friend Vec4 operator+( Arg left, Arg right )
	{
		return Vec4( _mm_add_ps( left.m_v, right.m_v ) );
	}

	friend Vec4 operator-( Arg left, Arg right )
	{
		return Vec4( _mm_sub_ps( left.m_v, right.m_v ) );
	}

	friend Vec4 operator*( Arg left, Arg right )
	{
		return Vec4( _mm_mul_ps( left.m_v, right.m_v ) );
	}

	friend Vec4 operator*( Arg left, float right )
	{
		return Vec4( _mm_mul_ps( left.m_v, _mm_set1_ps( right ) ) );
	}

	friend Vec4 operator*( float left, Arg right )
	{
		return Vec4( _mm_mul_ps( _mm_set1_ps( left ), right.m_v ) );
	}

	friend Vec4 operator/( Arg left, float right )
	{
		return left * Reciprocal( Vec4( right ) );
	}

	friend Vec4 operator*( Arg left, int right )
	{
#if ( SQUISH_USE_SSE == 1 )
		...
#else
		return Vec4( _mm_mul_ps( left.m_v, _mm_cvtepi32_ps( _mm_set1_epi32( right ) ) ) );
#endif
	}

	//! Returns a*b + c
	friend Vec4 MultiplyAdd( Arg a, Arg b, Arg c )
	{
		return Vec4( _mm_add_ps( _mm_mul_ps( a.m_v, b.m_v ), c.m_v ) );
	}

	//! Returns -( a*b - c )
	friend Vec4 NegativeMultiplySubtract( Arg a, Arg b, Arg c )
	{
		return Vec4( _mm_sub_ps( c.m_v, _mm_mul_ps( a.m_v, b.m_v ) ) );
	}

	template<const int a, const int b, const int c, const int d>
	friend Vec4 Merge( Arg lo, Arg hi );
	template<const int a, const int b, const int c, const int d>
	friend Vec4 Merge( Arg lo, Arg hi )
	{
		return Vec4( _mm_shuffle_ps( lo.m_v , hi.m_v , SQUISH_SSE_SHUF(
			a % 4,
			b % 4,
			c % 4,
			d % 4
		) ) );
	}

	template<const int f, const int t>
	friend Vec4 Shuffle( Arg a );
	template<const int f, const int t>
	friend Vec4 Shuffle( Arg a )
	{
		if (f == t)
			return a;

		return Vec4( _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( a.m_v ), SQUISH_SSE_SHUF(
			(t == 0 ? f : 0),
			(t == 1 ? f : 1),
			(t == 2 ? f : 2),
			(t == 3 ? f : 3)
		) ) ) );
	}

	template<const int f, const int t>
	friend Vec4 Exchange( Arg a );
	template<const int f, const int t>
	friend Vec4 Exchange( Arg a )
	{
		if (f == t)
			return a;

		return Vec4( _mm_castsi128_ps( _mm_shuffle_epi32( _mm_castps_si128( a.m_v ), SQUISH_SSE_SHUF(
			(t == 0 ? f : (f == 0 ? t : 0)),
			(t == 1 ? f : (f == 1 ? t : 1)),
			(t == 2 ? f : (f == 2 ? t : 2)),
			(t == 3 ? f : (f == 3 ? t : 3))
		) ) ) );
	}

	template<const int n>
	friend Vec4 RotateLeft( Arg a );
	template<const int n>
	friend Vec4 RotateLeft( Arg a )
	{
		return Vec4( _mm_shuffle_ps( a.m_v , a.m_v , SQUISH_SSE_SHUF(
			(n + 0) % 4,
			(n + 1) % 4,
			(n + 2) % 4,
			(n + 3) % 4
		) ) );
	}

	friend Vec4 Threshold( Arg a, Arg b ) {
		__m128 mask = _mm_cmpge_ps( a.m_v, b.m_v );
		__m128 res = _mm_and_ps( _mm_set1_ps(1.0f), mask );

		return Vec4( res );
	}

	friend Vec4 Select( Arg a, Arg b, Arg c )
	{
#if 0
		__m128 res;
		__m128 bits = _mm_cmpeq_ps( b.m_v, c.m_v );
		int mask = _mm_movemask_ps( bits );

		/* (1 >> 1) = 0
		 * (2 >> 1) = 1
		 * (4 >> 1) = 2
		mask = (mask & 7) >> 1;
		mask = (mask) * ((1 << 0) + (1 << 2) + (1 << 4) + (1 << 6));
		 */

		/**/ if (mask & 1)
		  res = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );
		else if (mask & 2)
		  res = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SHUF( 1, 1, 1, 1 ) );
		else if (mask & 4)
		  res = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SHUF( 2, 2, 2, 2 ) );
		else
		  res = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SHUF( 3, 3, 3, 3 ) );

		return Vec3( res );
#else
		// branch free, and no CPU<->SSEunit transfer
		__m128 mask = _mm_cmpeq_ps( b.m_v, c.m_v );
		__m128 res = _mm_and_ps( a.m_v, mask );

		__m128 r0 = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );
		__m128 r1 = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 1, 1, 1 ) );
		__m128 r2 = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 2, 2, 2, 2 ) );
		__m128 r3 = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 3, 3, 3, 3 ) );

		res = _mm_or_ps( _mm_or_ps( r0, r1 ), _mm_or_ps( r2, r3 ) );

		return Vec4(res);
#endif
	}

	friend Vec4 HorizontalAdd( Arg a )
	{
#if ( SQUISH_USE_SSE >= 4 ) && 0
		__m128 res = a.m_v;

		res = _mm_dp_ps( res, _mm_set1_ps ( 1.0f ), 0xFF );

		return Vec4( res );
#elif ( SQUISH_USE_SSE >= 3 )
		__m128 res = a.m_v;

		res = _mm_hadd_ps( res, res );
		res = _mm_hadd_ps( res, res );

		return Vec4( res );
#else
		__m128 res = a.m_v;

		res = _mm_add_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP64() ) );
		res = _mm_add_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP32() ) );

		return Vec4( res );
#endif
	}

	friend Vec4 HorizontalAdd( Arg a, Arg b )
	{
#if ( SQUISH_USE_SSE >= 3 )
		__m128 resc;

		resc = _mm_hadd_ps( a.m_v, b.m_v );
		resc = _mm_hadd_ps( resc, resc );
		resc = _mm_hadd_ps( resc, resc );

		return Vec4( resc );
#else
		__m128 resc;

		resc = _mm_add_ps( a.m_v, b.m_v );
		resc = _mm_add_ps( resc, _mm_shuffle_ps( resc, resc, SQUISH_SSE_SWAP64() ) );
		resc = _mm_add_ps( resc, _mm_shuffle_ps( resc, resc, SQUISH_SSE_SWAP32() ) );

		return Vec4( resc );
#endif
	}

	friend Vec4 HorizontalMin( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP64() ) );
		res = _mm_min_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP32() ) );

		return Vec4( res );
	}

	friend Vec4 HorizontalMax( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP64() ) );
		res = _mm_max_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SWAP32() ) );

		return Vec4( res );
	}
	
	friend Vec4 HorizontalMaxXY( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_max_ps(
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 1, 0, 1 ) ),
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 0, 1, 0 ) )
		);

		return Vec4( res );
	}
	
	friend Vec4 HorizontalMinXY( Arg a )
	{
		__m128 res = a.m_v;

		res = _mm_min_ps(
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 1, 0, 1 ) ),
		  _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 0, 1, 0 ) )
		);

		return Vec4( res );
	}

	friend Vec4 Reciprocal( Arg v )
	{
		// get the reciprocal estimate
		__m128 estimate = _mm_rcp_ps( v.m_v );

		// one round of Newton-Rhaphson refinement
		__m128 diff = _mm_sub_ps( _mm_set1_ps( 1.0f ), _mm_mul_ps( estimate, v.m_v ) );
		return Vec4( _mm_add_ps( _mm_mul_ps( diff, estimate ), estimate ) );
	}

	friend Vec4 ReciprocalSqrt( Arg v )
	{
		// get the reciprocal estimate
		__m128 estimate = _mm_rsqrt_ps( v.m_v );

		// one round of Newton-Rhaphson refinement
		__m128 diff = _mm_sub_ps( _mm_set1_ps( 3.0f ), _mm_mul_ps( estimate, _mm_mul_ps( estimate, v.m_v ) ) );
		return Vec4( _mm_mul_ps( _mm_mul_ps( diff, _mm_set1_ps( 0.5f ) ), estimate ) );
	}

	friend Vec4 Sqrt( Arg v )
	{
		return Vec4( _mm_sqrt_ps( v.m_v ) );
	}

	friend Vec4 Length( Arg left )
	{
		Vec4 sum = HorizontalAdd( Vec4( _mm_mul_ps( left.m_v, left.m_v ) ) );
		Vec4 sqt = Vec4( _mm_sqrt_ps( sum.m_v ) );

		return sqt;
	}

	friend Vec4 ReciprocalLength( Arg left )
	{
		Vec4 sum = HorizontalAdd( Vec4( _mm_mul_ps( left.m_v, left.m_v ) ) );
		Vec4 rsq = ReciprocalSqrt(sum);

		return rsq;
	}
	
	friend Vec4 Normalize( Arg left )
	{
		Vec4 sum = HorizontalAdd( Vec4( _mm_mul_ps( left.m_v, left.m_v ) ) );
		Vec4 rsq = ReciprocalSqrt(sum);

		return left * rsq;
	}
	
	friend Vec4 Normalize( Vec4& x, Vec4& y, Vec4& z )
	{
		Vec4 xx = x * x;
		Vec4 yy = y * y;
		Vec4 zz = z * z;

		Vec4 sum = xx + yy + zz;
		Vec4 rsq = ReciprocalSqrt(sum);

		x = x * rsq;
		y = y * rsq;
		z = z * rsq;

		return rsq;
	}

	template<const bool disarm, const bool killw>
	friend Vec4 Complement( Arg left );
	template<const bool disarm, const bool killw>
	friend Vec4 Complement( Arg left )
	{
		__m128 ren, res, rez;

		ren = left.m_v;
		rez = _mm_set1_ps( 1.0f );
		res = _mm_mul_ps( left.m_v, left.m_v );
#if ( SQUISH_USE_SSE >= 3 )
		res = _mm_hadd_ps( res, res );
#else
		res = _mm_add_ps( res, _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 1, 0, 1, 0 ) ) );
#endif
		if (!disarm) {
			// correct x² + y² > 1.0f by renormalization
			if ( _mm_comigt_ss( res, rez ) ) {
				res = ReciprocalSqrt( Vec4(res) ).m_v;
				res = _mm_shuffle_ps( res, res, SQUISH_SSE_SHUF( 0, 0, 0, 0 ) );

				ren = _mm_mul_ps( ren, res );
				res = rez;
			}
		}

		rez = _mm_sub_ps( rez, res );
		rez = _mm_sqrt_ps( rez );

		if (!killw) {
			res = _mm_shuffle_ps( ren, rez, SQUISH_SSE_SHUF( 3, 3, 0, 0 ) );
			res = _mm_shuffle_ps( ren, res, SQUISH_SSE_SHUF( 0, 1, 2, 0 ) );
		}
		else {
			res = _mm_movelh_ps( ren, rez );
			res = _mm_and_ps( res, _mm_castsi128_ps ( _mm_setr_epi32( ~0, ~0, ~0,  0 ) ) );
		}

		// sqrt(1.0f - (x² + y²))
		return Vec4( res );
	}

	template<const bool disarm>
	friend Vec4 Complement( Vec4 &left, Vec4 &right );
	template<const bool disarm>
	friend Vec4 Complement( Vec4 &left, Vec4 &right )
	{
		if (!disarm) {
			Vec4 len = left * left + right * right;
			Vec4 adj = ReciprocalSqrt(Max(Vec4(1.0f), len));

			// correct x² + y² > 1.0f by renormalization
			left  *= adj;
			right *= adj;

			// sqrt(1.0f - (x² + y²))
			return Sqrt(Vec4(1.0f) - Min(Vec4(1.0f), len));
		}
		else {
			Vec4 len = (left * left) + (right * right);

			// disarm x² + y² > 1.0f by letting NaN happen
			// ...

			// sqrt(1.0f - (x² + y²))
			return Sqrt(Vec4(1.0f) - len);
		}
	}

	template<const bool disarm>
	friend Vec4 ComplementPyramidal( Vec4 &left );
	template<const bool disarm>
	friend Vec4 ComplementPyramidal( Vec4 &left )
	{
		// 1 - max(abs(l, r)) == 1 + min(-abs(l, r))
		Vec4 res = HorizontalMinXY(Neg(left)) + Vec4(1.0f);

		res = TransferW(Normalize(TransferZW(left, KillW(res))), left);

		return res;
	}

	template<const bool disarm>
	friend Vec4 ComplementPyramidal( Vec4 &left, Vec4 &right );
	template<const bool disarm>
	friend Vec4 ComplementPyramidal( Vec4 &left, Vec4 &right )
	{
		// 1 - max(abs(l, r)) == 1 + min(-abs(l, r))
		Vec4 res = Min(Neg(left), Neg(right)) + Vec4(1.0f);

		Normalize(left, right, res);

		return res;
	}

	friend Vec4 Dot( Arg left, Arg right )
	{
#if ( SQUISH_USE_SSE >= 4 )
		return Vec4( _mm_dp_ps ( left.m_v, right.m_v, 0xFF ) );
#else
		return HorizontalAdd( Vec4( _mm_mul_ps( left.m_v, right.m_v ) ) );
#endif
	}

	friend void Dot( Arg left, Arg right, float *r )
	{
		Vec4 res = Dot( left, right );

		_mm_store_ss( r, res.m_v );
	}

	friend Vec4 Abs( Arg a )
	{
		return Vec4( _mm_and_ps( a.m_v, _mm_castsi128_ps( _mm_set1_epi32( 0x7FFFFFFF ) ) ) );
	}
	
	friend Vec4 Neg( Arg a )
	{
		return Vec4( _mm_or_ps( a.m_v, _mm_castsi128_ps( _mm_set1_epi32( 0x80000000 ) ) ) );
	}

	friend Vec4 Min( Arg left, Arg right )
	{
		return Vec4( _mm_min_ps( left.m_v, right.m_v ) );
	}

	friend Vec4 Max( Arg left, Arg right )
	{
		return Vec4( _mm_max_ps( left.m_v, right.m_v ) );
	}

	// clamp the output to [0, 1]
	Vec4 Clamp() const {
		Vec4 const one (1.0f);
		Vec4 const zero(0.0f);

		return Min(one, Max(zero, *this));
	}

	template<const bool round>
	friend Col4 FloatToInt( Vec4::Arg v );
	template<const bool round>
	friend Col4 FloatToInt( Vec4::Arg v )
	{
#if ( SQUISH_USE_SSE == 1 )
		...
#else
		// use SSE2 instructions
		if (round)
		      return Col4( _mm_cvtps_epi32( v.m_v ) );
		else
		      return Col4( _mm_cvttps_epi32( v.m_v ) );
#endif
	}

	friend Vec4 Truncate( Arg v )
	{
#if ( SQUISH_USE_SSE == 1 )
		// convert to ints
		__m128 input = v.m_v;
		__m64 lo = _mm_cvttps_pi32( input );
		__m64 hi = _mm_cvttps_pi32( _mm_movehl_ps( input, input ) );

		// convert to floats
		__m128 part = _mm_movelh_ps( input, _mm_cvtpi32_ps( input, hi ) );
		__m128 truncated = _mm_cvtpi32_ps( part, lo );

		// clear out the MMX multimedia state to allow FP calls later
		_mm_empty();
		
		return Vec4( truncated );
#else
		// use SSE2 instructions
		return Vec4( _mm_cvtepi32_ps( _mm_cvttps_epi32( v.m_v ) ) );
#endif
	}

	friend Vec4 AbsoluteDifference( Arg left, Arg right )
	{
		__m128 diff = _mm_sub_ps( left.m_v, right.m_v );
		diff = _mm_and_ps( diff, _mm_castsi128_ps( _mm_set1_epi32( 0x7FFFFFFF ) ) );
		return Vec4( diff );
	}

	friend Vec4 SummedAbsoluteDifference( Arg left, Arg right )
	{
		return HorizontalAdd( AbsoluteDifference( left, right ) );
	}

	friend Vec4 MaximumAbsoluteDifference( Arg left, Arg right )
	{
		return HorizontalMax( AbsoluteDifference( left, right ) );
	}

	friend int CompareEqualTo( Arg left, Arg right )
	{
		return _mm_movemask_ps( _mm_cmpeq_ps( left.m_v, right.m_v ) );
	}
	
	friend int CompareNotEqualTo( Arg left, Arg right )
	{
		return _mm_movemask_ps( _mm_cmpneq_ps( left.m_v, right.m_v ) );
	}

	friend int CompareLessThan( Arg left, Arg right )
	{
		return _mm_movemask_ps( _mm_cmplt_ps( left.m_v, right.m_v ) );
	}
	
	friend int CompareGreaterThan( Arg left, Arg right )
	{
		return _mm_movemask_ps( _mm_cmpgt_ps( left.m_v, right.m_v ) );
	}

	friend int CompareGreaterEqual( Arg left, Arg right )
	{
		return _mm_movemask_ps( _mm_cmpge_ps( left.m_v, right.m_v ) );
	}

	friend bool CompareAnyLessThan( Arg left, Arg right )
	{
		__m128 bits = _mm_cmplt_ps( left.m_v, right.m_v );
		int value = _mm_movemask_ps( bits );
		return value != 0x0;
	}

	friend bool CompareAnyGreaterThan( Arg left, Arg right )
	{
		__m128 bits = _mm_cmpgt_ps( left.m_v, right.m_v );
		int value = _mm_movemask_ps( bits );
		return value != 0x0;
	}

	friend bool CompareAllEqualTo( Arg left, Arg right )
	{
		__m128 bits = _mm_cmpeq_ps( left.m_v, right.m_v );
		int value = _mm_movemask_ps( bits );
		return value == 0xF;
	}

	friend Col4 CompareAllEqualTo_M4( Arg left, Arg right )
	{
		return Col4( _mm_cmpeq_epi32( _mm_castps_si128 ( left.m_v ), _mm_castps_si128 ( right.m_v ) ) );
	}
	
	friend Col4 CompareAllEqualTo_M8( Arg left, Arg right )
	{
		return Col4( _mm_cmpeq_epi8( _mm_castps_si128 ( left.m_v ), _mm_castps_si128 ( right.m_v ) ) );
	}
	
	friend int CompareFirstLessThan( Arg left, Arg right )
	{
		return _mm_comilt_ss( left.m_v, right.m_v );
	}
	
	friend int CompareFirstLessEqualTo( Arg left, Arg right )
	{
		return _mm_comile_ss( left.m_v, right.m_v );
	}

	friend int CompareFirstGreaterThan( Arg left, Arg right )
	{
		return _mm_comigt_ss( left.m_v, right.m_v );
	}

	friend int CompareFirstGreaterEqualTo( Arg left, Arg right )
	{
		return _mm_comige_ss( left.m_v, right.m_v );
	}

	friend int CompareFirstEqualTo( Arg left, Arg right )
	{
		return _mm_comieq_ss( left.m_v, right.m_v );
	}
	
	friend Vec4 IsGreaterThan( Arg left, Arg right )
	{
		return Vec4( _mm_cmpgt_ps( left.m_v, right.m_v ) );
	}
	
	friend Vec4 IsGreaterEqual( Arg left, Arg right )
	{
		return Vec4( _mm_cmpge_ps( left.m_v, right.m_v ) );
	}
	
	friend Vec4 IsNotEqualTo( Arg left, Arg right )
	{
		return Vec4( _mm_cmpneq_ps( left.m_v, right.m_v ) );
	}

	Vec4 IsOne( ) const
	{
		return Vec4( _mm_cmpeq_ps( m_v, _mm_set1_ps( 1.0f ) ) );
	}

	Vec4 IsNotOne( ) const
	{
		return Vec4( _mm_cmpneq_ps( m_v, _mm_set1_ps( 1.0f ) ) );
	}

	Vec4 IsZero( ) const
	{
		return Vec4( _mm_cmpeq_ps( m_v, _mm_set1_ps( 0.0f ) ) );
	}

	Vec4 IsNotZero( ) const
	{
		return Vec4( _mm_cmpneq_ps( m_v, _mm_set1_ps( 0.0f ) ) );
	}

	friend Vec4 TransferW( Arg left, Arg right )
	{
		/* [new W, ....., ....., old Z] */
		//m128 u = _mm_unpackhi_ps( left.m_v, right.m_v );
		/* [new W, new W, old Z, old Z] */
		__m128 u = _mm_shuffle_ps( left.m_v, right.m_v, SQUISH_SSE_SHUF( 2, 2, 3, 3 ) );
		/* [new W, old Z, old Y, old X] */
		       u = _mm_shuffle_ps( left.m_v, u, SQUISH_SSE_SHUF( 0, 1, 0, 2 ) );

		return Vec4( u );
	}

	friend Vec4 TransferZW( Arg left, Arg right )
	{
		return Vec4( _mm_shuffle_ps( left.m_v, right.m_v, SQUISH_SSE_SHUF( 0, 1, 2, 3 ) ) );
	}

	friend Vec4 KillW( Arg left )
	{
		return Vec4( _mm_and_ps( left.m_v, _mm_castsi128_ps ( _mm_setr_epi32( ~0, ~0, ~0,  0 ) ) ) );
	}

	friend Vec4 OnlyW( Arg left )
	{
		return Vec4( _mm_and_ps( left.m_v, _mm_castsi128_ps ( _mm_setr_epi32(  0,  0,  0, ~0 ) ) ) );
	}
	
	friend Vec4 CollapseW( Arg x, Arg y, Arg z, Arg w )
	{
		return Vec4( _mm_unpackhi_ps( _mm_unpackhi_ps( x.m_v, z.m_v ), _mm_unpackhi_ps( y.m_v, w.m_v ) ) );
	}

	void SwapXYZW( Vec4 &with )
	{
		/* inplace swap based on xors */
		     m_v = _mm_xor_ps( m_v, with.m_v );
		with.m_v = _mm_xor_ps( with.m_v, m_v );
		     m_v = _mm_xor_ps( m_v, with.m_v );
	}

	void SwapXYZ ( Vec4 &with )
	{
		/* [old W, old W, new Z, new Z] */
		__m128 u = _mm_shuffle_ps( m_v, with.m_v, SQUISH_SSE_SHUF( 3, 3, 2, 2 ) );
		__m128 v = _mm_shuffle_ps( with.m_v, m_v, SQUISH_SSE_SHUF( 3, 3, 2, 2 ) );
		__m128 w = m_v;

		/* [new X, new Y, new Z, old W] */
		     m_v = _mm_shuffle_ps( with.m_v, u, SQUISH_SSE_SHUF( 0, 1, 2, 0 ) );
		with.m_v = _mm_shuffle_ps(        w, v, SQUISH_SSE_SHUF( 0, 1, 2, 0 ) );
	}

	void SwapW   ( Vec4 &with )
	{
		/* [old Z, old Z, new W, new W] */
		__m128 u = _mm_shuffle_ps( m_v, with.m_v, SQUISH_SSE_SHUF( 2, 2, 3, 3 ) );
		__m128 v = _mm_shuffle_ps( with.m_v, m_v, SQUISH_SSE_SHUF( 2, 2, 3, 3 ) );

		/* [old X, old Y, old Z, new W] */
		     m_v = _mm_shuffle_ps(      m_v, u, SQUISH_SSE_SHUF( 0, 1, 0, 2 ) );
		with.m_v = _mm_shuffle_ps( with.m_v, v, SQUISH_SSE_SHUF( 0, 1, 0, 2 ) );
	}

	friend void LoadAligned( Vec4 &a, Vec4 &b, Arg c )
	{
	        a.m_v = c.m_v;
		b.m_v = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void LoadAligned( Vec4 &a, void const *source )
	{
		a.m_v = _mm_load_ps( (float const *)source );
	}

	friend void LoadUnaligned( Vec4 &a, void const *source )
	{
		a.m_v = _mm_loadu_ps( (float const *)source );
	}

	friend void LoadAligned( Vec4 &a, Vec4 &b, void const *source )
	{
		a.m_v = _mm_load_ps( (float const *)source );
		b.m_v = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void LoadUnaligned( Vec4 &a, Vec4 &b, void const *source )
	{
		a.m_v = _mm_loadu_ps( (float const *)source );
		b.m_v = _mm_shuffle_ps( a.m_v, a.m_v, SQUISH_SSE_SWAP64() );
	}

	friend void StoreAligned( Arg a, Arg b, Vec4 &c )
	{
		c.m_v = _mm_unpacklo_ps( a.m_v, b.m_v );
	}

	friend void StoreAligned( Arg a, void *destination )
	{
		_mm_store_ps( (float *)destination, a.m_v );
	}

	friend void StoreAligned( Arg a, Arg b, void *destination )
	{
		_mm_store_ps( (float *)destination, _mm_unpacklo_ps( a.m_v, b.m_v ) );
	}

	friend void StoreUnaligned( Arg a, void *destination )
	{
		_mm_storeu_ps( (float *)destination, a.m_v );
	}

	friend void StoreUnaligned( Arg a, Arg b, void *destination )
	{
		_mm_storeu_ps( (float *)destination, _mm_unpacklo_ps( a.m_v, b.m_v ) );
	}

private:
	__m128 m_v;
};

template<const bool round>
Col4 FloatToUHalf( Vec4::Arg v );
template<const bool round>
Col4 FloatToUHalf( Vec4::Arg v )
{
	Col4 h;

	h.GetR() = FloatToUHalf(v.X());
	h.GetG() = FloatToUHalf(v.Y());
	h.GetB() = FloatToUHalf(v.Z());
	h.GetA() = FloatToUHalf(v.W());

	return h;
}

template<const bool round>
Col4 FloatToSHalf( Vec4::Arg v );
template<const bool round>
Col4 FloatToSHalf( Vec4::Arg v )
{
	Col4 h;

	h.GetR() = FloatToSHalf(v.X());
	h.GetG() = FloatToSHalf(v.Y());
	h.GetB() = FloatToSHalf(v.Z());
	h.GetA() = FloatToSHalf(v.W());

	return h;
}

Vec4 UHalfToFloat( Col4::Arg v )
{
	Vec4 f;

	f.GetX() = UHalfToFloat((u16)v.R());
	f.GetY() = UHalfToFloat((u16)v.G());
	f.GetZ() = UHalfToFloat((u16)v.B());
	f.GetW() = UHalfToFloat((u16)v.A());

	return f;
}

Vec4 SHalfToFloat( Col4::Arg v )
{
	Vec4 f;

	f.GetX() = SHalfToFloat((u16)v.R());
	f.GetY() = SHalfToFloat((u16)v.G());
	f.GetZ() = SHalfToFloat((u16)v.B());
	f.GetW() = SHalfToFloat((u16)v.A());

	return f;
}

// TODO: figure out how to put static const instances into an incomplete class body
namespace Vec4C {

  const Vec4 zero = Vec4(0.0f);
  const Vec4 one = Vec4(1.0f);
  const Vec4 half = Vec4(0.5f);

}

// scalar types
typedef	Vec3  Scr3;
typedef	Vec4  Scr4;

#if	!defined(SQUISH_USE_PRE)
inline Vec3 LengthSquared( Vec3::Arg v )
{
  return Dot( v, v );
}

inline void LengthSquared( Vec3::Arg v , float *r )
{
  Dot( v, v, r );
}

inline Vec4 LengthSquared( Vec4::Arg v )
{
  return Dot( v, v );
}

inline void LengthSquared( Vec4::Arg v , float *r )
{
  Dot( v, v, r );
}
#endif

} // namespace squish

#endif // ndef SQUISH_SIMD_SSE_H
