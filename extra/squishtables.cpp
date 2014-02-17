#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <set>
#include <map>

struct qerror {

  /* calculate the largest left-aligned multiplier which repeats
   * values at most "codebits" long:
   *  codebits = 8 -> 0x01010101 (repeated 4 times)
   *  codebits = 6 -> 0x04104104 (repeated 5 times)
   *  codebits = 4 -> 0x11111111 (repeated 8 times)
   */
  static unsigned int reconstructor(const int numbits, const int codebits) {
    if (numbits >= codebits)
      return (1U << (numbits - codebits)) + reconstructor(numbits - codebits, codebits);
    return 0U;
  }

  /* calculate a reconstruction of value p with "codebits" length
   * with "precisionbits" length at left-aligned "position":
   *  position = 16, precision 8: 0b0000000000000000RRRRRRRR00000000
   *  position = 12, precision 6: 0b00000000000000000000RRRRRR000000
   *  position =  6, precision 8: 0b00000000000000000000000000RRRRRRrr
   */
  template<const int position, const int codebits, const int precisionbits>
  static unsigned int reconstruct(unsigned int p) {
    const unsigned int rec = reconstructor(32, codebits);
    const unsigned int msk = (~0U) << (32 - precisionbits);
    const unsigned int val = (rec * p) >> (32 - position);
    const unsigned int cot = (msk    ) >> (32 - position);

    if (cot > ((1U << precisionbits) - 1))
      return val & cot;
    return val;
  }

  /* shovel the highest possible number off the given error-value
   * which fits into "codebits" at left-aligned "position",
   * the value is precisely reconstructed/evaluated as defined by the
   * block-coding spec.
   */
  template<const int position, const int codebits, const int precisionbits>
  static unsigned int code(unsigned int &e) {
    assert(codebits      <= 8);
    assert(precisionbits <= 8);
    assert(codebits      <= precisionbits);

    // dropbits
    const int CB =            codebits;
    const int DB = position - codebits;
    const int NB = position           ;
    const int NR = (1 << NB);

    // ...ggggg?........ 14 - 6
    // ........bbbb?.... 14 - 5 - 5
    // ............rrrrr 14 - 0
    // ...gggggbbbbrrrrr

    // top bits
    assert(e < NR);
    unsigned int top = (e >> DB);

    // reconstruction
    unsigned int rec = reconstruct<NB,CB,precisionbits>(  top);
    // ensure truncation to next lower reconstruction value
    if (rec > e) rec = reconstruct<NB,CB,precisionbits>(--top);

    // reduce remaining error
    assert(e >= rec);
    e -= rec;

    // return raw unreconstructed top value
    return top;
  }
};

// Associated to partition -1, 16 * 0 bit
static const unsigned int partitionmasks_1[1] =
{
  0xFFFF
};

// Associated to partition 0-63, 16 * 1 bit
static const unsigned int partitionmasks_2[64] =
{
  0xCCCC, 0x8888, 0xEEEE, 0xECC8,
  0xC880, 0xFEEC, 0xFEC8, 0xEC80,
  0xC800, 0xFFEC, 0xFE80, 0xE800,
  0xFFE8, 0xFF00, 0xFFF0, 0xF000,
  0xF710, 0x008E, 0x7100, 0x08CE,
  0x008C, 0x7310, 0x3100, 0x8CCE,
  0x088C, 0x3110, 0x6666, 0x366C,
  0x17E8, 0x0FF0, 0x718E, 0x399C,

  0xaaaa, 0xf0f0, 0x5a5a, 0x33cc,
  0x3c3c, 0x55aa, 0x9696, 0xa55a,
  0x73ce, 0x13c8, 0x324c, 0x3bdc,
  0x6996, 0xc33c, 0x9966, 0x0660,
  0x0272, 0x04e4, 0x4e40, 0x2720,
  0xc936, 0x936c, 0x39c6, 0x639c,
  0x9336, 0x9cc6, 0x817e, 0xe718,
  0xccf0, 0x0fcc, 0x7744, 0xee22,
};

// Associated to partition 64-127, 16 * 2 bit
static const unsigned int partitionmasks_3[64] =
{
  0xf60008cc, 0x73008cc8, 0x3310cc80, 0x00ceec00,
  0xcc003300, 0xcc0000cc, 0x00ccff00, 0x3300cccc,
  0xf0000f00, 0xf0000ff0, 0xff0000f0, 0x88884444,
  0x88886666, 0xcccc2222, 0xec80136c, 0x7310008c,
  0xc80036c8, 0x310008ce, 0xccc03330, 0x0cccf000,
  0xee0000ee, 0x77008888, 0xcc0022c0, 0x33004430,
  0x00cc0c22, 0xfc880344, 0x06606996, 0x66009960,
  0xc88c0330, 0xf9000066, 0x0cc0c22c, 0x73108c00,

  0xec801300, 0x08cec400, 0xec80004c, 0x44442222,
  0x0f0000f0, 0x49242492, 0x42942942, 0x0c30c30c,
  0x03c0c03c, 0xff0000aa, 0x5500aa00, 0xcccc3030,
  0x0c0cc0c0, 0x66669090, 0x0ff0a00a, 0x5550aaa0,
  0xf0000aaa, 0x0e0ee0e0, 0x88887070, 0x99906660,
  0xe00e0ee0, 0x88880770, 0xf0000666, 0x99006600,
  0xff000066, 0xc00c0cc0, 0xcccc0330, 0x90006000,
  0x08088080, 0xeeee1010, 0xfff0000a, 0x731008ce,
};

static const int shorterindex[64][/*6*/3] = {
  { /*0,*/  /*0,*/15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3, 8},
  { /*0,*/  /*0,*/15,  /*0,*/15, 8},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3},
  { /*0,*/  /*0,*/15,  /*0,*/15, 8},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3, 8},
  { /*0,*/  /*0,*/15,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 3, 8},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 3},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 3, 8},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 6,15},
  { /*0,*/  /*0,*/15,  /*0,*/10, 8},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 5, 3},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 8, 6},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 6,10},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15,10},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 8},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3},
  { /*0,*/  /*0,*/ 6,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 5,10},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 6,10},
  { /*0,*/  /*0,*/ 8,  /*0,*/10, 8},
  { /*0,*/  /*0,*/15,  /*0,*/ 8, 9},
  { /*0,*/  /*0,*/15,  /*0,*/15,10},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 6},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 8},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 3},
  { /*0,*/  /*0,*/15,  /*0,*/15, 6},
  { /*0,*/  /*0,*/15,  /*0,*/15, 6},
  { /*0,*/  /*0,*/ 6,  /*0,*/15, 8},
  { /*0,*/  /*0,*/ 6,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/15, 3},
  { /*0,*/  /*0,*/ 6,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 8,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/10,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 5,15},
  { /*0,*/  /*0,*/15,  /*0,*/10,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 8,15},
  { /*0,*/  /*0,*/15,  /*0,*/13,15},
  { /*0,*/  /*0,*/15,  /*0,*/15, 3},
  { /*0,*/  /*0,*/ 2,  /*0,*/12,15},
  { /*0,*/  /*0,*/ 2,  /*0,*/ 3,15},
  { /*0,*/  /*0,*/15,  /*0,*/ 3, 8}
};

int main(int argc, char **argv) {
#if 0
  printf("static const int firstofsetinpartition[64][6] = {\n");
  for (int p = 0; p < 64; p++) {
    printf("{ ");

    for (int s = 1; s < 4; s++) {
      for (int o = 0; o < s; o++) {
      	int masks = 0, mask[4];
      	int first = 0;

	// partition_1 mask is: bit cleared -> set 1, bit set -> set 2
	// partition_2 mask is: bit cleared -> set 1, bit set -> set 2, hi bit set -> set 3
      	if (s == 1) {
      		masks = partitionmasks_1[0];
      		first = 0;

		mask[0] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);	// color-set
		mask[1] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);	// alpha-set

		masks = mask[o];
		while (!(masks & 1))
			masks >>= 1, first++;
      	}
      	else if (s == 2) {
      		masks = partitionmasks_2[p];
      		first = 0;

		mask[0] = (~masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);
		mask[1] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);

		masks = mask[o];
		while (!(masks & 1))
			masks >>= 1, first++;
      	}
      	else if (s == 3) {
      		masks = partitionmasks_3[p];
      		first = 0;

		mask[0] = (~masks & 0xFFFF) & (~masks >> 16);
		mask[1] = ( masks & 0xFFFF) & (~masks >> 16);
		mask[2] = (0xFFFF & 0xFFFF) & ( masks >> 16);

		masks = mask[o];
		while (!(masks & 1))
			masks >>= 1, first++;
      	}

        printf("%2d%s", first, (s == 3 && o == 2 ? "" : ","));
      }

      printf("%s", (s == 3 ? "" : " "));
    }

    printf("}%s\n", (p == 63 ? "" : ","));
  }
  printf("};\n");
#endif

#if 0
//printf("static const unsigned int blockxor[2][64][/*6*/3][4] = {\n");
//for (int b = 2; b < 4; b++) {
//  printf("/* %d bit index */\n", b);
//  printf("{\n");
    printf("static const unsigned int blockxor[64][/*6*/5][4] = {\n");
    for (int p = 0; p < 64; p++) {
      printf("{ ");

      for (int s = 1; s < 4; s++) {
        for (int o = 0; o < s; o++) {
//        const char hibit = 1 << b;
//        const char hixor = hibit - 1;
      	  const unsigned int hixor = 0xFF;
      	  unsigned int masks = 0, mask[4];
          unsigned int bxor[4] = {0};

	  // partition_1 mask is: bit cleared -> set 1, bit set -> set 2
	  // partition_2 mask is: bit cleared -> set 1, bit set -> set 2, hi bit set -> set 3
      	  if (s == 1) {
      		masks = partitionmasks_1[0];

		mask[0] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);	// color-set
		mask[1] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);	// alpha-set

		masks = mask[o];
//fprintf(stderr, "%d 0x%04x %d\n", o, masks, (masks >> 0) & 1);
		for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j) {
		  int pos = ((i * 4) + j);
		  if (((masks >> pos) & 1))
		    bxor[i] |= (hixor << (j * 8));
		}
      	  }
      	  else if (s == 2) {
      		masks = partitionmasks_2[p];

		mask[0] = (~masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);
		mask[1] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);

		masks = mask[o];
//fprintf(stderr, "%d 0x%04x %d\n", o, masks, (masks >> 0) & 1);
		for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j) {
		  int pos = ((i * 4) + j);
		  if (((masks >> pos) & 1))
		    bxor[i] |= (hixor << (j * 8));
		}
      	  }
      	  else if (s == 3) {
      		masks = partitionmasks_3[p];

		mask[0] = (~masks & 0xFFFF) & (~masks >> 16);
		mask[1] = ( masks & 0xFFFF) & (~masks >> 16);
		mask[2] = (0xFFFF & 0xFFFF) & ( masks >> 16);

		masks = mask[o];
//fprintf(stderr, "%d 0x%04x %d\n", o, masks, (masks >> 0) & 1);
		for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j) {
		  int pos = ((i * 4) + j);
		  if (((masks >> pos) & 1))
		    bxor[i] |= (hixor << (j * 8));
		}
      	  }

          printf("%s{0x%08x,0x%08x,0x%08x,0x%08x}%s%s",
          	(s == 1 ? "/*" : ""),
          	bxor[0], bxor[1], bxor[2], bxor[3],
          	(s == 3 && o == 2 ? "" : ","),
          	(s == 1 ? "*/ " : " ")
          );
        }
      }

      printf("}%s\n", (p == 63 ? "" : ","));
    }
    printf("};\n");
//  printf("}%s\n", (b == 3 ? "" : ","));
//}
//printf("};\n");
#endif

#if 0
  printf("static const u8 whichsetinpartition[64][/*3*/2][16] = {\n");
  for (int p = 0; p < 64; p++) {
    printf("{ ");

    for (int s = 1; s < 4; s++) {
      int masks = 0, mask[4];
      int first = 0;

      // partition_1 mask is: bit cleared -> set 1, bit set -> set 2
      // partition_2 mask is: bit cleared -> set 1, bit set -> set 2, hi bit set -> set 3
      if (s == 1) {
      	masks = partitionmasks_1[0];

      	mask[0] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);	// color-set
      	mask[1] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);	// alpha-set

        printf("/*{0},*/ ");
      }
      else if (s == 2) {
      	masks = partitionmasks_2[p];

      	mask[0] = (~masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);
      	mask[1] = ( masks & 0xFFFF) & ( 0xFFFFFFFF >> 16);

        printf("{");
        for (int i = 0; i < 16; i++) {
          if (mask[0] & (1 << i))
            printf("%d", 0);
          else if (mask[1] & (1 << i))
            printf("%d", 1);

          printf("%s", (i == 15 ? "" : ","));
        }
        printf("}");
      }
      else if (s == 3) {
      	masks = partitionmasks_3[p];
      	first = 0;

      	mask[0] = (~masks & 0xFFFF) & (~masks >> 16);
      	mask[1] = ( masks & 0xFFFF) & (~masks >> 16);
      	mask[2] = (0xFFFF & 0xFFFF) & ( masks >> 16);

        printf("{");
        for (int i = 0; i < 16; i++) {
          if (mask[0] & (1 << i))
            printf("%d", 0);
          else if (mask[1] & (1 << i))
            printf("%d", 1);
          else if (mask[2] & (1 << i))
            printf("%d", 2);

          printf("%s", (i == 15 ? "" : ","));
        }
        printf("}");
      }

      printf("%s", (s != 2 ? "" : ","));
    }

    printf("}%s\n", (p == 63 ? "" : ","));
  }
  printf("};\n");
#endif

#if 0
  // getting some sRGB/gamma-tables
  float basefpartition = 0.0031308f;
  float baseipartition = 0.004045f;
  float basefslope = 12.92f / 1.0f;
  float baseislope = 1.0f / 12.92f;
  float basefgamma = 2.4f / 1.0f;
  float baseigamma = 1.0f / 2.4f;
  float baseoffset = 0.055f;
  float baseLUT[3][256];

  printf("float baseLUT[256] = {\n");
//for (int c = 0; c < 3; c++) {
    printf("   ");
//  printf("{");
    for (int v = 0; v < 256; v++) {
      double fv = (double)v / 255.0;

      if (fv <= baseipartition)
      	printf("%gf%s", (float)(fv / 12.92), (v == 255 ? "" : ","));
      else
      	printf("%gf%s", (float)pow((fv + 0.055) / (1.0 + 0.055), 2.4), (v == 255 ? "" : ","));

      if (((v & 7) == 7) && (v != 255))
        printf("\n   ");
    }
//  printf("}%s\n", (c == 2 ? "" : ","));
//}
  printf("};\n");

  printf("float baseLUT[256] = {\n");
//for (int c = 0; c < 3; c++) {
    printf("   ");
//  printf("{");
    for (int v = 0; v < 256; v++) {
      double fv = (double)v / 255.0;

      printf("%gf%s", (float)(fv), (v == 255 ? "" : ","));

      if (((v & 7) == 7) && (v != 255))
        printf("\n   ");
    }
//  printf("}%s\n", (c == 2 ? "" : ","));
//}
  printf("};\n");
#endif

#if 0
  // checking out the quantizer
  for (int b = 6; b <= 6; b++) {
    float grid    = 1.0f * ((1 << b) - 1);
    float gridrcp = 1.0f / ((1 << b) - 1);
//  float grid2    = 1.0f * ((1 << b) - 1) * ( 1.0f + gridrcp);
    float grid2    = 1.0f + ((1 << b) - 1);
//  float grid2off = ((128.0f / 255.0f) / (1.0f + gridrcp) - (128.0f / 255.0f));
//  float grid2off = (-(127.5f / 255.0f) * gridrcp) / (1.0f + gridrcp);
//  float grid2off = (-0.5f * gridrcp) / (1.0f + gridrcp);
//  float grid3off = 1.0f * ((1 << b) - 1) * (-0.5f * gridrcp);
    float grid3off = -0.5f;
    float grid2rcp = 1.0f / ((1 << b) - 1);

    for (int i = 0; i <= 255; i++) {
      int j = (i >> (8 - b));
      int k = (j << (8 - b)); k |= (k >> (b)); k |= (k >> (2 * b)); k |= (k >> (3 * b));

      float x = 1.0f / 255.0f * i;		// rgba[] / 255.0f
      float l = floor(grid * x + 0.5f) * gridrcp;
      int m = (int)floor(l * 255.0f);		// FloatToInt(x * 255.0f)

      printf("%3d %3d %3d %.8f %3d -> %3d", i, j, k, l, m, m - k);

      printf("    |    ");

//    float y = 1.0f / 255.0f * (128 - ((1.0f + gridrcp) * (128 - i)));
//    float y = (1.0f + gridrcp) * (1.0f / 255.0f * 128 / (1.0f + gridrcp) - (1.0f / 255.0f * 128 - 1.0f / 255.0f * i));
//    float n = floor(grid2 * (y + grid2off) + 0.5f) * gridrcp;

      float y = 1.0f / 255.0f * i;		// rgba[] / 255.0f
      float n = floor(grid2 * y + grid3off + 0.5f) * gridrcp;
      int o = (int)floor(n * 255.0f + 0.5f);	// FloatToInt(x * 255.0f)

      printf("%3d %3d %3d %.8f %3d -> %3d\n", i, j, k, n, o, o - k);
    }
  }
#endif

#if 0
  float qLUT[3][9][256];

  // checking out the quantizer
  for (int b = 1; b <= 8; b++) {
#if 0
    for (int i = 0; i < (1 << b); i++) {
      int j = ((i + 0) << (8 - b)); j |= (j >> (b)); j |= (j >> (2 * b)); j |= (j >> (3 * b));
      int k = ((i + 1) << (8 - b)); k |= (k >> (b)); k |= (k >> (2 * b)); k |= (k >> (3 * b));

      printf(" %3d [%2x, %3d]\n", i, j, k - j);
    }
#endif

    for (int i = 0; i <= 255; i++) {
      // half of a bin
      int   h =  (1 << (8 - b)) >> 1;
      float f = ((1 << (8 - b)) >> 1) / 255.0f;
      // range
      float r = (1 << b) - 1;
      float w =   1.0f / 255.0f;
      float fw = (((1 << (8 - b)) >> 1) - 1) / 255.0f;
      float fw_ = (((1 << (8 - b)) >> 1) - 1);
      float u =   1.0f / (1 <<      b );
      float v = 255.0f / (1 <<      b );
      float s =   1.0f / (1 << (8 - b));
      float z = 255.0f / (1 << (8 - b));
      float o = 127.0f / (255.0f * 255.0f);
      float x = f * s + o * s;
      float y = s - s * s;
      float mul  = ((float)(1 << (8 - b)) / (1 << (0 * b)));
            mul += ((float)(1 << (8 - b)) / (1 << (1 * b)));
            mul += ((float)(1 << (8 - b)) / (1 << (2 * b)));
            mul += ((float)(1 << (8 - b)) / (1 << (3 * b)));

      // truncating quantizer
      int j = (i >> (8 - b));
      int k = (j << (8 - b)); k |= (k >> (b)); k |= (k >> (2 * b)); k |= (k >> (3 * b));
      // nearest rounding quantizer
//    int l = ((i + (h - j) - 1 + (i >> 7)) >> (8 - b));
//    int l = ((i + (h - (i >> b)) - 1 + (i >> 7)) >> (8 - b));
      int l = ((i + (h - (i >> b))) >> (8 - b));
//    int l = ((i + (((h << b) - i - 1 + (1 << (b - 1))) >> b) - 0 + (0 >> 7)) >> (8 - b));
      int m = (l << (8 - b)); m |= (m >> (b)); m |= (m >> (2 * b)); m |= (m >> (3 * b));

      float val = i / 255.0f;

//    float p = (i + (h - (i >> b)) - 1 + (i >> 7)) >> (8 - b);
//    float p = (i + (h - (i >> b)) - 1 + (i >> 7)) * s;
//    float p = (i + (h - (i / (1 << b))) - 1 + (i >> 7)) * s;
//    float p = (i + (h - floor(i * u)) - 1 + (i >> 7)) * s;
//    float p = (i + (h - floor(val * v)) - 1 + (i / 127.5f)) * s;
//    float p = (i + (i / 128.0f) + (h - floor(val * v)) - 1) * s;
//    float q = floor(p) / r;
//    float q = (127.5f - (127.5f - floor(p)) * 1.001f) / r;	???
//    float q = floor(p * mul) / 255.0f;

//    float p = (i/255.0f + (i/255.0f / 128.0f) + (h/255.0f - floor(val * v)/255.0f) - 1/255.0f) * s;
//    float p = (val + (val / 128.0f) + (f - floor(val * v) / 255.0f) - w) * s;
//    float p = ((val * 129.0f) / 128.0f + (f - floor(val * v) / 255.0f) - w) * s;
//    float p = (val * (129.0f / 128.0f) + (f - floor(val * v) / 255.0f) - w) * s;
//    float p = (val * (129.0f / 128.0f) + fw - floor(val * v) / 255.0f) * s;
//    float p = (val + fw / (128.0f / 127.0f) - floor(val * v) / 255.0f / (128.0f / 127.0f)) * s * (128.0f / 127.0f);
//    float p = (val * (129.0f / 128.0f) * 255.0f + fw_ - floor(val * v)) * s;
//    float p = (val * 257.0f + fw_ - floor(val * v)) * s;
//    float q = floor(p) / r;
//    float q = floor(floor(p) * mul) / 255.0f;

      float newrange = (1 << b) - 1;
      float generalhalf = (0.5f * (1 << (8 - b))) * newrange / 255.0f;
      float p = val * newrange + generalhalf;
      float q = floor(p) / newrange;
//    float r = floor(floor(p * newrange) * (255.0f / newrange));

//    float p = (i/255.0f + (h/255.0f - (val * z)/255.0f) + (127.0f / 255.0f)/255.0f) * s;
//    float p = (val + (f - (val * s)) + o) * s;
//    float p = (val * (1.0f - s) + f + o) * s;
//    float p = val * y + x;
//    float q = floor(p * 255.0f) / r;

//    int t = (int)floor(p * 255.0f);
      int t = (int)floor(q * 255.0f) >> (8 - b);
      int d = (t << (8 - b)); d |= (d >> (1 * b)); d |= (d >> (2 * b)); d |= (d >> (3 * b));

      printf("%3d %.8f [%2x, %3d] [%2x, %3d] %.8f %3d %.8f\n", i, val, k, i - k, m, i - m, q, t, d / 255.0f);

      int lval = (int)floor(p);

      int _b = b - 1;
      float _newrange = (1 << _b) - 1;
      float _generalhalf = (0.5f * (1 << (8 - _b))) * _newrange / 255.0f;
      float _p = val * _newrange + _generalhalf;
      float _q = floor(_p) / _newrange;

      int _t = (int)floor(_q * 255.0f) >> (8 - _b);
      int dc = (_t << (8 - _b)); dc |= 0 << (8 - b); dc |= (dc >> (1 * b)); dc |= (dc >> (2 * b)); dc |= (dc >> (3 * b));
      int ds = (_t << (8 - _b)); ds |= 1 << (8 - b); ds |= (ds >> (1 * b)); ds |= (ds >> (2 * b)); ds |= (ds >> (3 * b));
      int distc = (i - dc); if (distc < 0) distc = -distc;
      int dists = (i - ds); if (dists < 0) dists = -dists;

      for (int j = (lval & (~1)) - 1; j <= (lval | ( 1)); j++) {
        int __t = lval >> 1;//j >> (8 - _b);
        int _dc = (__t << (8 - _b)); _dc |= 0 << (8 - b); _dc |= (_dc >> (1 * b)); _dc |= (_dc >> (2 * b)); _dc |= (_dc >> (3 * b));
        int _ds = (__t << (8 - _b)); _ds |= 1 << (8 - b); _ds |= (_ds >> (1 * b)); _ds |= (_ds >> (2 * b)); _ds |= (_ds >> (3 * b));

        int _distc = (i - _dc); if (_distc < 0) _distc = -_distc;
        int _dists = (i - _ds); if (_dists < 0) _dists = -_dists;

        if (distc > _distc) { distc = _distc; dc = _dc; }
        if (dists > _dists) { dists = _dists; ds = _ds; }
      }

      qLUT[0][b][lval] = d;
      qLUT[1][b][lval] = dc;
      qLUT[2][b][lval] = ds;
    }

//  printf("\n");
  }

  // checking out the quantizer
  for (int b = 1; b <= 8; b++) {
  for (int r = 0; r <= 2; r++) {
    printf("const a16 float qLUT_%d%s[%d] = {\n", b, r == 0 ? "all" : (r == 1 ? "clr" : "set"), 1 << b);
    for (int i = 0, j = 0; i < (1 << b); i++, j++) {
      if (j == 0)
        printf(" ");
      printf(" %3d.0f / 255.0f,", (int)qLUT[r][b][i]);
      if (j == 7)
        printf("\n"), j = -1;
    }
    printf("};\n\n", b);
  }
    printf("/* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */\n", b);
  }
#endif

#if 0
  // checking out the odd-test for fp
  float check = 1.0f;
  while (check < 256.0f) {
    // bit-flip approach
    //  when we have the integer value at 1 << (b+1) check
    //  if the number is even or odd
    //   clear bit: N - (odd + 0)
    //   fill  bit: N + (1 - odd)
    // we have to find a fast function for fp which tells us
    // if the number is odd after truncation
    //  after truncation means all the mantissa is free of
    //  fp-bits and only full of number-bits
    //  if the mantissa is zero it's a 2^x number
    //   (mantissa <<= (exp - 128)) & 1
    //   (mantissa <<= ((fp >> 23) & 0x7F)) & (1 << 22)
    //   ((mantissa <<= (1 + ((fp >> 23) & 0x7F)))) * 0x7F) & 0x3F800000
    //   ((mantissa <<= (1 + ((fp >> 23) & 0x7F)))) * 0x27F) for any != 1
    //   ((mantissa <<= (1 + (fp >> 23))) * 0x27F) for positive > 1
    //   ((mantissa <<= (1 + ((fp + 2) >> 23))) * 0x27F) for positive (is shifted to 3 doesn't matter for oddness)
    //  works for all except 1
    //  alternating between 0 and 1 is by xoring => (fp ^ 0x3F800000)
    //  other approach
    //   multiply by ???
    //    max(256, 255 * 256), test if <> 256
    //    max(256, 127 * 256), test if <> 256
    float adj = check + 2.0f;
    unsigned long ifp = *((unsigned long int *)&adj);

    unsigned long pow2 = (ifp >> 23) /*& 0x7F*/ /*shift modulo?*/;
    unsigned long mant = (ifp << pow2) & 0x00400000;
    unsigned long conv = (mant * 0xFE);

    float odd = *((float *)&conv);

    printf("%3g -> 0x%08x %3d 0x%08x 0x%08x -> %3g\n", adj, ifp, pow2, mant, conv, odd);

    check += 1.0f;
  }
#endif

#if 0
  // how much different codebooks has alpha5-mode?
  std::set<unsigned long long> trial5;
  std::map<unsigned long long, std::set<unsigned short> > match5;

  for (int s = 1; s <= 254; s++) {
    for (int e = s; e <= 254; e++) {
      int m1 = (((5 - 1) * s + 1 * e) / 5);
      int m2 = (((5 - 2) * s + 2 * e) / 5);
      int m3 = (((5 - 3) * s + 3 * e) / 5);
      int m4 = (((5 - 4) * s + 4 * e) / 5);

      unsigned long long cb;
      unsigned long long cba =
        (m4 << 24) | (m3 << 16) | (m2 << 8) | (m1 << 0);
      unsigned long long cbb =
        (m1 << 24) | (m2 << 16) | (m3 << 8) | (m4 << 0);

      if ((254 - e) < (s - 1))
      	continue;
      else
      	cb = cba;

      trial5.insert(cb);
      match5[cb].insert((s << 0) | (e << 8));
    }
  }

  int lookat = 0;

  for (int s = 1; s <= 254; s++) {
    for (int e = s; e <= 254; e++) {

    }
  }

  // well, that are a 19532 distinct ones!
  printf("alpha5: %d codebooks\n", trial5.size());

  // how much different codebooks has alpha7-mode?
  std::set<unsigned long long> trial7;
  std::map<unsigned long long, std::set<unsigned short> > match7;

  for (int s = 0; s <= 255; s++) {
    for (int e = s; e <= 255; e++) {
      int m1 = (((7 - 1) * s + 1 * e) / 7);
      int m2 = (((7 - 2) * s + 2 * e) / 7);
      int m3 = (((7 - 3) * s + 3 * e) / 7);
      int m4 = (((7 - 4) * s + 4 * e) / 7);
      int m5 = (((7 - 5) * s + 5 * e) / 7);
      int m6 = (((7 - 6) * s + 6 * e) / 7);
      unsigned long long cb1 =
        (m4 << 24) | (m3 << 16) | (m2 << 8) | (m1 << 0);
      unsigned long long cb2 =
        (m5 <<  8) | (m6 <<  0);
      unsigned long long cb =
        (cb2 << 32) | (cb1 <<  0);

      trial7.insert(cb);
      match7[cb].insert((s << 0) | (e << 8));
    }
  }

  // well, that are a 23569 distinct ones!
  printf("alpha7: %d codebooks\n", trial7.size());
#endif

#if 0
  int vmin = 0xFF, vmax = 0x00;
  int vals[16];
  for (int v = 0; v < 16; v++) {
    vals[v] = std::rand() & 0xFF;
    if (vmin > vals[v])
      vmin = vals[v];
    if (vmax < vals[v])
      vmax = vals[v];
  }

  // binary search, tangent-fitting
  int s = vmin;
  int e = vmax;
  int p = 1;
  int l = 0;
  int b = 0;
  int oe = 0;
  int os = 0;
  while ((l = (e - s) >> p) != 0) {
    int ms = std::max(s + os - l, 0x00);
    int ps = std::min(s + os + l, 0xFF);
    int me = std::max(e - oe - l, 0x00);
    int pe = std::min(e - oe + l, 0xFF);

    int error0 = 0;
    int error1 = 0;
    int error2 = 0;
    int error3 = 0;
    int error4 = 0;
    for (int v = 0; v < 16; v++) {
      // find the closest code
      int dist0 = 0xFFFF;
      int dist1 = 0xFFFF;
      int dist2 = 0xFFFF;
      int dist3 = 0xFFFF;
      int dist4 = 0xFFFF;

      // for all possible codebook-entries
      for (int f = 0; f < 8; f++) {
        int cb0 = (((7 - f) *  s + f *  e) / 7);
        int cb1 = (((7 - f) *  s + f * me) / 7);
        int cb2 = (((7 - f) *  s + f * pe) / 7);
        int cb3 = (((7 - f) * ms + f *  e) / 7);
        int cb4 = (((7 - f) * ps + f *  e) / 7);

        int d0 = std::abs(vals[v] - cb0); d0 *= d0;
        int d1 = std::abs(vals[v] - cb1); d1 *= d1;
        int d2 = std::abs(vals[v] - cb2); d2 *= d2;
        int d3 = std::abs(vals[v] - cb3); d3 *= d3;
        int d4 = std::abs(vals[v] - cb4); d4 *= d4;

        if (dist0 > d0) dist0 = d0;
        if (dist1 > d1) dist1 = d1;
        if (dist2 > d2) dist2 = d2;
        if (dist3 > d3) dist3 = d3;
        if (dist4 > d4) dist4 = d4;
      }

      // accumulate the error
      error0 += dist0;
      error1 += dist1;
      error2 += dist2;
      error3 += dist3;
      error4 += dist4;
    }

    int merr = 0xFFFFFF;
    if (merr > error0) merr = error0;
    if (merr > error1) merr = error1;
    if (merr > error4) merr = error4;
    if (merr > error2) merr = error2;
    if (merr > error3) merr = error3;

    if (merr == error0) b = 0, p++;	// half range
    if (merr == error1) b = 1, oe += l;	// range up
    if (merr == error4) b = 2, os += l;	// range dn
    if (merr == error2) b = 3, oe -= l;	// range up
    if (merr == error3) b = 4, os -= l;	// range dn
  }

  // final match
  s += os;
  e -= oe;
#endif

#if 0
template<class DataType = unsigned short, const int overlap = 1, const int payload0 = 4>
struct distributed1x1D {

public:
  static const int fragments    = 1;
  static const int size		= sizeof(DataType) * 8;
  static const int precision    = -(overlap * (fragments - 1)) + payload0;

  /* ...XXXXXXXXXXXXXX	  original number
   * ............00000?	  bits of fragment 0, with overlap bit
   * ............00000	  combined number
   */

};

template<class DataType = unsigned short, const int overlap = 1, const int payload0 = 8, const int payload1 = 8>
struct distributed2x1D {

public:
  static const int fragments    = 2;
  static const int size		= sizeof(DataType) * 8;
  static const int precision    = -(overlap * (fragments - 1)) + payload0 + payload1;

  /* ...XXXXXXXXXXXXXX	  original number
   * ........1111?....	  bits of fragment 1, with overlap bit
   * ............00000?	  bits of fragment 0, with overlap bit
   * ........111100000	  combined number
   */

};

template<class DataType = unsigned short, const int overlap = 1, const int payload0 = 5, const int payload1 = 6, const int payload2 = 5>
struct distributed3x1D {

public:
  static const int fragments    = 3;
  static const int size		= sizeof(DataType) * 8;
  static const int precision    = -(overlap * (fragments - 1)) + payload0 + payload1 + payload2;

  /* ...XXXXXXXXXXXXXX	  original number
   * ...22222?........	  bits of fragment 2, with overlap bit
   * ........1111?....	  bits of fragment 1, with overlap bit
   * ............00000?	  bits of fragment 0, with overlap bit
   * ...22222111100000	  combined number
   */

};

template<class DataType = unsigned short, const int overlap = 1, const int payload0 = 4, const int payload1 = 4, const int payload2 = 4, const int payload3 = 4>
struct distributed4x1D {

public:
  static const int fragments    = 4;
  static const int size		= sizeof(DataType) * 8;
  static const int precision    = -(overlap * (fragments - 1)) + payload0 + payload1 + payload2 + payload3;

  /* ....XXXXXXXXXXXXX	  original number
   * ....333?.........	  bits of fragment 3, with overlap bit
   * .......222?......	  bits of fragment 2, with overlap bit
   * ..........111?...	  bits of fragment 1, with overlap bit
   * .............0000?	  bits of fragment 0, with overlap bit
   * ....3332221110000	  combined number
   */

private:
  /* ---------------------------------------------------------------------------
   */
  inline void import(const DataType tht) {
  }
};
#endif

#if 0
  unsigned int ming = 0xFFFFFFFF, maxg = 0x00000000;
  unsigned int minb = 0xFFFFFFFF, maxb = 0x00000000;
  unsigned int minr = 0xFFFFFFFF, maxr = 0x00000000;

#define NB  14
  for (int i = 0; i < (1 << NB); i++) {
    // can be 0 or 1, this makes you see
    int testofs = 0;

    // ...ggggg?........ 14 - 6
    // ........bbbb?.... 14 - 5 - 5
    // ............rrrrr 14 - 0
    // ...gggggbbbbrrrrr
    int e = i;
    unsigned int ee, ge, be, re;

    // top 6/5 bit
    assert(e < (1 << (NB - 0)));
    int gp = (i >> (NB - 6)) & (~0); gp += testofs;
    bool gt = false;

    // reconstruction
//  int gr = (gp * ((1 << 12) + (1 << 6) + (1 << 0))) >> ((3 * 6) - (NB - 0));
    int gr = ((gp << 2) + (gp >> 4)) << (NB - 8);

    int tt = qerror::reconstruct<NB,6,8>(gp);
    assert(gr == tt);
    ee = e; ge = qerror::code<NB,6,8>(ee);

    // ensure truncation to next lower reconstruction value
    if (gt = gr > e) gp -= 1,
//    gr = (gp * ((1 << 12) + (1 << 6) + (1 << 0))) >> ((3 * 6) - (NB - 0));
      gr = ((gp << 2) + (gp >> 4)) << (NB - 8);

    // error
    e = e - gr;
    assert(ee == e);
    assert(ge == gp);
    ming = std::min(ming, ee);
    maxg = std::max(maxg, ee);
//  fprintf(stderr, "0x%04x - (0x%04x                  ) => 0x%04x\n", i, gr, e);

    // middle 5/4 bit
    assert(e < (1 << (NB - 5 - 0)));
    int bp = (e >> (NB - 5 - 5)) & (~0); bp += testofs;
    bool bt = false;

    // reconstruction
//  int br = (bp * ((1 << 5) + (1 << 0))) >> ((2 * 5) - (NB - 5));
    int br = ((bp << 3) + (bp >> 2)) << (NB - 5 - 8);

    int gg = qerror::reconstruct<NB-5,5,8>(bp);
    assert(br == gg);
    ee = e; be = qerror::code<NB-5,5,8>(ee);

    // ensure truncation to next lower reconstruction value
    if (bt = br > e) bp -= 1,
//    br = (bp * ((1 << 5) + (1 << 0))) >> ((2 * 5) - (NB - 5));
      br = ((bp << 3) + (bp >> 2)) << (NB - 5 - 8);

    // error
    e = e - br;
    assert(ee == e);
    assert(be == bp);
    minb = std::min(minb, ee);
    maxb = std::max(maxb, ee);
//  fprintf(stderr, "0x%04x - (0x%04x + 0x%04x          [---]) => 0x%04x\n", i, gr, br, e);

    // bottom 4/5 bit
    assert(e < (1 << (NB - 5 - 4 - 0)));
    int rp = (e >> (NB - 5 - 4 - 5)) & (~0); rp += testofs;
    bool rt = false;

    // reconstruction
//  int rr = (rp * ((1 << 0))) >> ((1 * 5) - (NB - 4 - 4));
    int rr = ((rp << 3) + (rp >> 2)) >> 3;

    int uu = qerror::reconstruct<NB-5-4,5,8>(rp);
    assert(rr == uu);
    ee = e; re = qerror::code<NB-5-4,5,8>(ee);

    // ensure truncation to next lower reconstruction value
    if (rt = rr > e) rp -= 1,
//    rr = (rp * ((1 << 0))) >> ((1 * 5) - (NB - 5 - 4));
      rr = ((rp << 3) + (rp >> 2)) >> 3;

    // error
    e = e - rr;
    assert(ee == e);
    assert(re == rp);
    minr = std::min(minr, ee);
    maxr = std::max(maxr, ee);

    ge = ((ge << 2) + (ge >> 4));
    re = ((re << 3) + (re >> 2));
    be = ((be << 3) + (be >> 2));

#define GB	6
#define BB	5
#define RB	5
    fprintf(stderr, "0x%04x - (0x%04x + 0x%04x + 0x%04x [%c%c%c]) => 0x%04x + 0x%04x + 0x%04x = 0x%04x [%8.4f%] => 0x%04x\n",
      i, gr, br, rr, gt ? 't' : ' ', bt ? 't' : ' ', rt ? 't' : ' ',
      ge, be, re,
      ((
        (ge << (0 + RB + BB - 1)) +
        (be << (0 + RB      - 1)) +
        (re << (0              ))
       ) >> (8 - RB)),
      ge * (1.0f * (1 << (0 + RB + BB - 1)) / (1 << (NB + (8 - RB)))) +
      be * (1.0f * (1 << (0 + RB      - 1)) / (1 << (NB + (8 - RB)))) +
      re * (1.0f * (1 << (0              )) / (1 << (NB + (8 - RB)))),
      e);
  }

  fprintf(stderr, "g [0x%04x,0x%04x] remaining corrective capacity is 0x%04x\n", ming, maxg, 0x1FF);
  fprintf(stderr, "b [0x%04x,0x%04x] remaining corrective capacity is 0x%04x\n", minb, maxb, 0x01F);
  fprintf(stderr, "r [0x%04x,0x%04x] remaining corrective capacity is 0x%04x\n", minr, maxr, 0x000);
#endif

#if 0
  unsigned int ming = 0xFFFFFFFF, maxg = 0x00000000;
  unsigned int minb = 0xFFFFFFFF, maxb = 0x00000000;
  unsigned int minr = 0xFFFFFFFF, maxr = 0x00000000;

#define NB  14
  for (int i = 0; i < (1 << NB); i++) {
    unsigned int ee = i, ge, be, re;

    ge = qerror::code<NB,6,8>(ee);

    ming = std::min(ming, ee);
    maxg = std::max(maxg, ee);

    be = qerror::code<NB-5,5,8>(ee);

    minb = std::min(minb, ee);
    maxb = std::max(maxb, ee);

    re = qerror::code<NB-5-4,5,8>(ee);

    minr = std::min(minr, ee);
    maxr = std::max(maxr, ee);

    ge = ((ge << 2) + (ge >> 4));
    be = ((be << 3) + (be >> 2));
    re = ((re << 3) + (re >> 2));

#define GB	6
#define BB	5
#define RB	5
    fprintf(stderr, "0x%04x => 0x%04x + 0x%04x + 0x%04x = 0x%04x/%8.4f% [%8.4f%] => 0x%04x\n",
      i,
      ge, be, re,
      ((
        (ge << (0 + RB + BB - 1)) +
        (be << (0 + RB      - 1)) +
        (re << (0              ))
       ) >> (8 - RB)),
      ((
        (ge << (0 + RB + BB - 1)) +
        (be << (0 + RB      - 1)) +
        (re << (0              ))
       ) >> (8 - RB)) / (1.0f * ((1 << NB) - 1)),

       /*
      (ge / 255.0f) * ((float)0xFF / (1 << (8 +               0))) +
      (be / 255.0f) * ((float)0xFF / (1 << (8 + GB          - 1))) +
      (re / 255.0f) * ((float)0xFF / (1 << (8 + GB - 1 + BB - 1))),
        */

      (ge / 255.0f) * ((float)((0xFF << (RB + BB - 1)) >> (8 - RB)) / ((1 << NB) - 1)) +
      (be / 255.0f) * ((float)((0xFF << (RB      - 1)) >> (8 - RB)) / ((1 << NB) - 1)) +
      (re / 255.0f) * ((float)((0xFF << (0          )) >> (8 - RB)) / ((1 << NB) - 1)),

       /*
        ((ge / 255.0f) * ((255.0f * 255.0f) / (256.0f * ((1 << ((8 - RB) + NB - (0 + RB + BB - 1))) - 1)))) +
        ((be / 255.0f) * ((255.0f * 255.0f) / (256.0f * ((1 << ((8 - RB) + NB - (0 + RB      - 1))) - 1)))) +
        ((re / 255.0f) * ((255.0f * 255.0f) / (256.0f * ((1 << ((8 - RB) + NB - (0              ))) - 1)))),
	*/
      ee);
  }

  fprintf(stderr, "g [0x%04x,0x%04x] remaining corrective capacity is 0x%04x\n", ming, maxg, 0x1FF);
  fprintf(stderr, "b [0x%04x,0x%04x] remaining corrective capacity is 0x%04x\n", minb, maxb, 0x01F);
  fprintf(stderr, "r [0x%04x,0x%04x] remaining corrective capacity is 0x%04x\n", minr, maxr, 0x000);
#endif

#if 0
#define S1(x)  int(x)
  const int weights_int[5][16] = {
    {S1(0)                                                                                                                        },  // 0
    {S1(0),                                                                                                                 S1(64)},  // 1
    {S1(0),                                 S1(21),                                 S1(43),                                 S1(64)},  // 2
    {S1(0),          S1(9),         S1(18),         S1(27),                 S1(37),         S1(46),         S1(55),         S1(64)},  // 3
    {S1(0),  S1(4),  S1(9), S1(13), S1(17), S1(21), S1(26), S1(30), S1(34), S1(38), S1(43), S1(47), S1(51), S1(55), S1(60), S1(64)}   // 4
  };

  for (int type = 0; type <= 1; type++) {
    int mini = 3, maxi = 4;
    int minb = 5, maxb = 6;
    int mins = 0, maxs = 0;
    if (type == 1) {
      mini = 4, maxi = 16;
      minb = 4, maxb = 8;
      mins = 0, maxs = 2;
    }

  for (int ipol = mini; ipol <= maxi; (ipol == 3 ? ipol++ : ipol *= 2))
  for (int bits = minb; bits <= maxb; bits++)
  for (int shrd = mins; shrd <= maxs; shrd++) {
    bool hit[16][1 << 8] = {{0}};

    for (int st = 0; st < (1 << bits); st++)
    for (int en = 0; en < (1 << bits); en++) {
      int sst = st;
      int sen = en;

      // clear shared bit
      if (shrd == 1)
      	sst &= ~1, sen &= ~1;
      // set shared bit
      else if (shrd == 2)
      	sst |=  1, sen |=  1;

      int rst = (sst * (1 << bits) + sst) >> (2 * bits - 8);
      int ren = (sen * (1 << bits) + sen) >> (2 * bits - 8);

      for (int i = 0; i < ipol; i++) {
	int ist = rst * (           i);
	int ien = ren * (ipol - 1 - i);
	int imd = (ist + ien) / (ipol - 1);

        if (type == 1) {
          int             lut = 0;
          if (ipol >=  2) lut = 1;
          if (ipol >=  4) lut = 2;
          if (ipol >=  8) lut = 3;
          if (ipol >= 16) lut = 4;

	  ist = rst * weights_int[lut][           i];
	  ien = ren * weights_int[lut][ipol - 1 - i];
	  imd = (ist + ien) / 64;
        }

	hit[i][imd] = true;
      }
    }

    fprintf(stderr, "%s: %d bits (shared %s) with %d intermediates:\n", !type ? "BC1" : "BC7", bits, shrd == 0 ? "disabled" : (shrd == 1 ? "cleared" : "set"), ipol);
    for (int i = 0; i < ipol; i++) {
      float error = 0.0f;
      int num = 0;
      for (int p = 0; p < 256; p++) {
        num += hit[i][p];
        if (!hit[i][p]) {
          int s = 1;

          for (; s < 256; s++) {
            if ((p - s) >= 0) {
              if (hit[i][p - s])
              	break;
            }
            else if ((p + s) < 256) {
              if (hit[i][p + s])
              	break;
            }
          }

          error += s * s;
        }
      }

      fprintf(stderr, " %2d: %7.3f%%, %8.4f RMSE\n", i, num * 100.0f / 256.0f, sqrtf(error / 256.0f));
    }
  }
  }
#endif

#if 0
  for (int m = 0; m <= 7; m++)
  for (int i = 0; i <= 127; i++) {
    int ipol = i * m * 32;
    int codes0 =  ipol / 7;
    int codes1 = -ipol / 7;
//  int codes1 = ((((ipol * 0x2493 + (ipol << 16)) >> 16)) >> 3);
//  int codes2 = ((((ipol * 0x2493 + (ipol << 16) + (0xFFFF << 3)) >> 16)) >> 3);
//  int codes2 = ((((ipol * 0x2493) >> 16) + ipol + 0x7) >> 3);
    int codes2 = (((( ipol * 0x2493) >> 16) +  ipol) >> 3);
    int codes3 = ((((-ipol * 0x2493) >> 16) + -ipol) >> 3) + (ipol != 0);

//  fprintf(stderr, "%d %4d 0x%08x 0x%08x%c 0x%08x%c\n", m, i, codes0, codes1, codes1 == codes0 ? '+' : '-', codes2, codes2 == codes0 ? '+' : '-');
    fprintf(stderr, "%d %4d 0x%08x 0x%08x %d | 0x%08x 0x%08x %d\n", m, i, codes0, codes2, codes0 - codes2, codes1, codes3, codes1 - codes3);
  }
#endif
  
#if 0
  for (int m = 0; m <= 5; m++)
  for (int i = -127; i <= 127; i++) {
    int ipol = i * m * 32;
    int codes0 = ipol / 5;
//  int codes2 = ((((ipol * 0x2493) >> 16) + ipol) >> 3) + (ipol < 0);
    int codes2 = ((((ipol * 0xCCCD) >> 16)       ) >> 2) + (ipol < 0);

//  fprintf(stderr, "%d %4d 0x%08x 0x%08x%c 0x%08x%c\n", m, i, codes0, codes1, codes1 == codes0 ? '+' : '-', codes2, codes2 == codes0 ? '+' : '-');
    fprintf(stderr, "%d %4d 0x%08x 0x%08x %d\n", m, i, codes0, codes2, codes0 - codes2);
  }
#endif
  
#if 0
  for (unsigned int S = 0; S <= 3; S++) {
  for (unsigned int A = 0; A <= 1; A++) {
  for (unsigned int O = 0; O <= 1; O++) {
  for (unsigned int N = 0; N <= 0xFFFF; N++) {

  bool good = true;
  for (int m = 0; m <= 5; m++)
  for (int i = 0; i <= 255; i++) {
    int xpol = (i * m * 32);
    int ipol = (i * m * 32);

    int codes0 = xpol / 5;
    int codes2 = (int((int(ipol * N) >> 16) + (A * ipol)) >> S) + (O * (ipol < 0));
    
    if (codes2 != codes0) {
    //fprintf(stderr, "((((ipol * 0x%04x) >> 16) + (%d * ipol >> 1) + %d) >> %d): %c\n", N, A, O, S, good ? '#' : '-');
      good = false; break; }
  }
  
  // ((((ipol * 0x3334) >> 16) + (0 * ipol)) >> 0) + (1 * (ipol < 0)): # divide by 5 from -127 to 127
  // ((((ipol * 0x4925) >> 16) + (0 * ipol)) >> 1) + (1 * (ipol < 0)): # divide by 7 from -127 to 127
  //
  // ((((ipol * 0xcccd) >> 16) + (0 * ipol)) >> 2) + (1 * (ipol < 0)): # divide by 5 from -127 to 127
  // ((((ipol * 0x2493) >> 16) + (1 * ipol)) >> 3) + (1 * (ipol < 0)): # divide by 7 from -127 to 127
  if (good) {
    fprintf(stderr, "((((ipol * 0x%04x) >> 16) + (%d * ipol)) >> %d) + (%d * (ipol < 0)): %c\n", N, A, S, O, good ? '#' : '-');
    fflush(stderr);
  }

  }
  }
  }
  }
#endif

  for (float X = -1.0f; X <= 1.0f; X += 0.125f)
  for (float Y = -1.0f; Y <= 1.0f; Y += 0.125f)
  for (float Z =  0.0f; Z <= 1.0f; Z += 0.125f)
  {
    float x = X, y = Y, z = Z;
    float len = sqrtf(x * x + y * y + z * z);
    x /= len; y /= len; z /= len;
    float u, v, w;

    float rx = fabsf(x);
    float ry = fabsf(y);
    float nx, ny, nz;
      nz = 1.0f;
    if (rx > ry)
      nx = 1.0f, ny = 0.0f;
    else
      ny = 1.0f, nx = 0.0f;
    if (x < 0.0f)
      nx = -nx;
    if (y < 0.0f)
      ny = -ny;
    
    // 0.70710678118654752440084436210485
//  float dist = sqrtf(0.5f * 0.5f + 0.5f * 0.5f);
//  float dist = (0.5f * 0.5f + 0.5f * 0.5f);
    float dist = 0.70710678118654752440084436210485;
//  float dot = (x * nx / dist + y * ny / dist + z * nz / dist);
//  float dot = (x * nx + y * ny + z * nz) / (dist);
    float dot = (x * nx + y * ny + z * nz);
    
    dot = (rx > ry ? rx : ry) + z;

    u = x / dot;
    v = y / dot;
    w = z / dot;
    
    float ru = fabsf(u);
    float rv = fabsf(v);
    float p = (ru > rv ? ru : rv) + w;

    fprintf(stdout, "%f,%f,%f | %f,%f,%f: %f | %f,%f,%f = %f\n", x, y, z, nx, ny, nz, dot, u, v, w, p);
  }
}
