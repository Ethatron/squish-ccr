#ifndef SQUISH_CONFIG_H
#include "config.h"
#endif

#ifndef IBL_ONLY_STRUCT
#define IBL_ONLY_STRUCT
struct IndexBlockLookup_CCR
{
  int mapped;
};

#define	IBL_ALPHA5  0
#define	IBL_ALPHA7  1
#define	IBL_COLOR3  2
#define	IBL_COLOR4  3

#if	defined(USE_AMP) || defined(USE_COMPUTE)
#if	!defined(USE_COMPUTE)
typedef const ::Concurrency::array_view<const   IndexBlockLookup_CCR, 2>   &IndexBlockLUT;
#else
typedef const                                   IndexBlockLookup_CCR        IndexBlockLUT[4][8];
#endif
#endif
#endif

#if	defined(ONLY_ARRAY)
static_hlsl const IndexBlockLookup_CCR lookup_c34a57_ccr[4][8] =
{
  // all   : 0 is start-value
  // all   : 1 is stop-value
  // alpha5: 6 is white (can't be changed, fixed definition)
  // alpha5: 7 is black (can't be changed, fixed definition)
  // color3: 3 is trans (can't be changed, fixed definition)

  {  1, 0, 5, 4, 3, 2, 6, 7  },	// alpha5: (a >  b), 8 indices valid, 1 << 3
  {  1, 0, 7, 6, 5, 4, 3, 2  },	// alpha7: (a <  b), 8 indices valid, 1 << 3
  { 1, 0, 2, 3,   4, 5, 6, 7 },	// color3: (a >  b), 4 indices valid, 1 << 2
  { 1, 0, 3, 2,   5, 4, 7, 6 },	// color4: (a <  b), 4 indices valid, 1 << 2
//{ 0, 0, 0, 0,   0, 0, 0, 0 };	// color4: (a == b)
};
#elif	!defined(USE_COMPUTE)
extern const IndexBlockLookup_CCR lookup_c34a57_ccr[4][8];
#endif // LUT
