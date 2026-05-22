/*
 * hamming.c - Hamming(7,4) single-error-correcting code.
 *
 * Bit positions are 1-indexed. Positions that are powers of two
 * (1, 2, 4) hold parity bits; the rest (3, 5, 6, 7) hold data.
 * Each parity bit covers the positions whose index has a particular
 * bit set, so the syndrome of a corrupted word is exactly the
 * 1-indexed position of the flipped bit.
 */
#include "coding/correct.h"

/* Extract bit at 1-indexed position p from a 7-bit codeword. */
static int getbit(uint8_t word, int p) {
    return (word >> (p - 1)) & 1;
}

/* Set bit at 1-indexed position p to value v. */
static uint8_t setbit(uint8_t word, int p, int v) {
    uint8_t mask = (uint8_t)(1u << (p - 1));
    return v ? (word | mask) : (word & (uint8_t)~mask);
}

uint8_t dsp_hamming74_encode(uint8_t nibble) {
    int d1 = (nibble >> 0) & 1;
    int d2 = (nibble >> 1) & 1;
    int d3 = (nibble >> 2) & 1;
    int d4 = (nibble >> 3) & 1;

    uint8_t c = 0;
    /* Data bits at positions 3, 5, 6, 7. */
    c = setbit(c, 3, d1);
    c = setbit(c, 5, d2);
    c = setbit(c, 6, d3);
    c = setbit(c, 7, d4);

    /* Parity bits cover positions sharing an index bit:
     *   p1 (pos 1) covers 3,5,7
     *   p2 (pos 2) covers 3,6,7
     *   p4 (pos 4) covers 5,6,7  */
    c = setbit(c, 1, d1 ^ d2 ^ d4);
    c = setbit(c, 2, d1 ^ d3 ^ d4);
    c = setbit(c, 4, d2 ^ d3 ^ d4);
    return c;
}

int dsp_hamming74_syndrome(uint8_t codeword) {
    /* Recompute each parity check; a failing check contributes its
     * position weight, and the sum is the flipped bit's position. */
    int s1 = getbit(codeword, 1) ^ getbit(codeword, 3)
           ^ getbit(codeword, 5) ^ getbit(codeword, 7);
    int s2 = getbit(codeword, 2) ^ getbit(codeword, 3)
           ^ getbit(codeword, 6) ^ getbit(codeword, 7);
    int s4 = getbit(codeword, 4) ^ getbit(codeword, 5)
           ^ getbit(codeword, 6) ^ getbit(codeword, 7);
    return s1 * 1 + s2 * 2 + s4 * 4;
}

uint8_t dsp_hamming74_decode(uint8_t codeword, uint8_t *corrected_out) {
    int syn = dsp_hamming74_syndrome(codeword);

    /* A non-zero syndrome points straight at the corrupt bit. */
    if (syn != 0)
        codeword ^= (uint8_t)(1u << (syn - 1));

    if (corrected_out)
        *corrected_out = codeword & 0x7F;

    /* Reassemble the 4 data bits from positions 3, 5, 6, 7. */
    int d1 = getbit(codeword, 3);
    int d2 = getbit(codeword, 5);
    int d3 = getbit(codeword, 6);
    int d4 = getbit(codeword, 7);
    return (uint8_t)(d1 | (d2 << 1) | (d3 << 2) | (d4 << 3));
}
