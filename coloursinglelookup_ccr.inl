#ifndef SQUISH_CONFIG_H
#include "config.h"
#endif

/* different kind of LUTs
 *
 *  flat: won't fit into the immediate constant buffer
 *  packed: will put one value in each ICB-vector, set the rest to 0 -> {LUT-val, 0, 0, 0}
 *  vector: will put 4 values in each ICB-vector, smallest ICB of these
 *
 * the files will take care to make the correct code load, you just need to
 * change the include, rest is automatic
 */

//nclude "coloursinglelookup_ccr_flat.inl"
//nclude "coloursinglelookup_ccr_packed.inl"
#include "coloursinglelookup_ccr_vector.inl"
