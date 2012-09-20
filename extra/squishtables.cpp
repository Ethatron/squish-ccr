#include <stdio.h>
#include <math.h>

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
}
