/*
 * correct.h - Forward Error CORRECTION (FEC)
 *
 * PROBLEM SOLVED
 *   FEC adds enough structured redundancy that the receiver can
 *   REPAIR corrupted bits with no retransmission. This is essential
 *   where a round-trip is impossible or too slow: deep-space links,
 *   broadcast, stored media, real-time streaming.
 *
 * SCHEMES PROVIDED
 *   Hamming(7,4)  - 4 data bits + 3 parity bits. Corrects any single-
 *                   bit error and detects (with the extra overall
 *                   parity bit, SECDED) double-bit errors. Used in
 *                   ECC memory and satellite links.
 *   Reed-Solomon  - operates on SYMBOLS (groups of bits) over a finite
 *                   field, so a damaged symbol counts as one error
 *                   however many of its bits flipped. This makes RS
 *                   excellent for BURST errors - CDs/DVDs, QR codes,
 *                   DVB. RS(n,k) corrects up to (n-k)/2 symbol errors.
 *   Convolutional - encodes a continuous bit STREAM; each output
 *                   depends on the current and several past input
 *                   bits. Decoded with the Viterbi algorithm, which
 *                   finds the most likely transmitted sequence.
 *
 * NOT IMPLEMENTED HERE (documented for completeness)
 *   LDPC codes  - sparse parity-check matrices decoded by iterative
 *                 belief propagation on a bipartite graph. Near
 *                 Shannon-limit performance; used in 5G and Wi-Fi 6.
 *   Turbo codes - two convolutional codes joined by an interleaver,
 *                 decoded iteratively. Used in 4G LTE and deep space.
 *   Both need large matrices / iterative soft decoders beyond the
 *   scope of this guide; see docs/ALGORITHMS.md for the comparison.
 *
 * SOFT-DECISION DECODING
 *   The Viterbi decoder here also accepts soft inputs - real-valued
 *   channel confidences rather than hard 0/1 bits. Using the analog
 *   information the demodulator already has typically buys ~2 dB of
 *   coding gain over hard-decision decoding.
 */
#ifndef DSP_CORRECT_H
#define DSP_CORRECT_H

#include <stddef.h>
#include <stdint.h>

/* ===================================================================
 * Hamming(7,4) - single-error-correcting block code
 * =================================================================== */

/*
 * Encode the low 4 bits of `nibble` into a 7-bit Hamming codeword.
 * Bit layout (1-indexed): positions 1,2,4 are parity; 3,5,6,7 data.
 */
uint8_t dsp_hamming74_encode(uint8_t nibble);

/*
 * Decode a 7-bit Hamming codeword.
 *   codeword     : received 7-bit value (possibly corrupted)
 *   corrected_out: if non-NULL, receives the corrected 7-bit codeword
 * Returns the 4-bit data nibble. Any single-bit error is corrected
 * automatically; the syndrome identifies the flipped position.
 */
uint8_t dsp_hamming74_decode(uint8_t codeword, uint8_t *corrected_out);

/*
 * Syndrome of a received codeword: 0 means no detected error,
 * otherwise the 1-indexed bit position that is in error.
 */
int dsp_hamming74_syndrome(uint8_t codeword);

/* ===================================================================
 * Reed-Solomon over GF(2^8) - burst-error-correcting block code
 * =================================================================== */

/* Maximum codeword length for GF(2^8): 2^8 - 1. */
#define DSP_RS_NN 255

/*
 * Reed-Solomon codec context. RS(n,k) with n = k + 2t parity symbols
 * corrects up to t symbol errors. Call dsp_rs_init once before use.
 */
typedef struct {
    int nroots;          /* number of parity symbols (n - k) = 2t */
    uint8_t gen[DSP_RS_NN + 1];   /* generator polynomial coefficients */
} dsp_rs;

/*
 * Initialise an RS codec with `nroots` parity symbols (must be even
 * and < 255). Returns 0 on success, -1 on bad parameters.
 */
int dsp_rs_init(dsp_rs *rs, int nroots);

/*
 * Encode: append rs->nroots parity symbols to `data` (length k).
 *   data   : k message symbols
 *   k      : message length in symbols
 *   parity : output buffer, length rs->nroots
 */
void dsp_rs_encode(const dsp_rs *rs, const uint8_t *data, size_t k,
                   uint8_t *parity);

/*
 * Decode a full codeword in place (data followed by parity).
 *   code : n = k + nroots symbols, corrected in place
 *   n    : total codeword length
 * Returns the number of symbol errors corrected, or -1 if the error
 * count exceeded the code's capability (uncorrectable).
 */
int dsp_rs_decode(const dsp_rs *rs, uint8_t *code, size_t n);

/* ===================================================================
 * Convolutional code (rate 1/2, K=3) with Viterbi decoding
 * =================================================================== */

/*
 * Encode a bit stream with the standard rate-1/2, constraint-length-3
 * convolutional code (generator polynomials 0o7 and 0o5). Each input
 * bit produces two output bits.
 *   bits  : input bits (one bit per byte, value 0/1), length nbits
 *   out   : output bits, length 2*nbits + tail; caller-allocated with
 *           room for 2*(nbits + 2) bits (2 flush bits return to state 0)
 * Returns the number of output bits written.
 */
size_t dsp_conv_encode(const uint8_t *bits, size_t nbits, uint8_t *out);

/* Output length produced by dsp_conv_encode for nbits input bits. */
static inline size_t dsp_conv_encoded_len(size_t nbits) {
    return 2 * (nbits + 2);   /* +2 tail bits to flush the register */
}

/*
 * Hard-decision Viterbi decode.
 *   recv   : received coded bits (0/1), length nsym
 *   nsym   : number of received bits (must be even)
 *   out    : decoded message bits, caller-allocated; length
 *            nsym/2 - 2 (the tail bits are dropped)
 * Returns the number of decoded message bits.
 */
size_t dsp_viterbi_decode(const uint8_t *recv, size_t nsym, uint8_t *out);

/*
 * Soft-decision Viterbi decode. Identical to the hard version but the
 * input is real-valued per-bit confidence in [-1, +1], where +1 means
 * "very likely a 1" and -1 means "very likely a 0". Using the analog
 * confidence rather than a thresholded bit improves correction.
 *   recv : soft values, length nsym
 */
size_t dsp_viterbi_decode_soft(const double *recv, size_t nsym,
                               uint8_t *out);

#endif /* DSP_CORRECT_H */
