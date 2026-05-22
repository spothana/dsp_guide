/*
 * qam.h - Quadrature Amplitude Modulation
 *
 * PROBLEM SOLVED
 *   A radio sends complex-valued SYMBOLS, not bits. A QAM mapper packs
 *   a group of bits into one point of a 2-D constellation - an (I, Q)
 *   pair, the real and imaginary parts of a complex symbol. The
 *   receiver does the reverse: it slices a noisy received point back
 *   to the nearest constellation point and reads off the bits.
 *
 *   bits/symbol = log2(M):
 *     QPSK    M = 4   2 bits/symbol
 *     16-QAM  M = 16  4 bits/symbol
 *     64-QAM  M = 64  6 bits/symbol
 *   Higher M carries more data per symbol but packs points closer
 *   together, so it needs a cleaner channel - the rate/robustness
 *   trade-off at the heart of adaptive modulation in Wi-Fi and 5G.
 *
 * GRAY CODING
 *   Constellation points are labelled so ADJACENT points differ in
 *   exactly ONE bit. Noise most often pushes a symbol to a neighbour,
 *   and Gray coding makes that common error cost a single bit flip
 *   rather than several.
 *
 * NORMALISATION
 *   Each constellation is scaled to unit average symbol energy, so
 *   different QAM orders can be compared at the same signal power.
 */
#ifndef DSP_QAM_H
#define DSP_QAM_H

#include "../common.h"
#include <stdint.h>

/* Supported constellations. */
typedef enum {
    DSP_QAM_QPSK   = 4,    /* 2 bits/symbol */
    DSP_QAM_16     = 16,   /* 4 bits/symbol */
    DSP_QAM_64     = 64    /* 6 bits/symbol */
} dsp_qam_order;

/* Bits carried by one symbol of the given order: log2(M). */
size_t dsp_qam_bits_per_symbol(dsp_qam_order order);

/*
 * Map a bit stream to QAM symbols.
 *   bits   : input bits (one bit per byte, value 0/1)
 *   nbits  : number of bits - MUST be a multiple of bits_per_symbol
 *   syms   : output complex symbols, length nbits / bits_per_symbol
 * Returns the number of symbols written, or 0 on a length mismatch.
 */
size_t dsp_qam_modulate(dsp_qam_order order, const uint8_t *bits,
                        size_t nbits, cplx *syms);

/* Number of symbols dsp_qam_modulate produces for nbits input bits. */
size_t dsp_qam_symbol_count(dsp_qam_order order, size_t nbits);

/*
 * Hard-decision demap: slice each received symbol to the nearest
 * constellation point and emit its bits.
 *   syms   : received complex symbols, length nsym
 *   bits   : output bits, length nsym * bits_per_symbol
 * Returns the number of bits written.
 */
size_t dsp_qam_demodulate(dsp_qam_order order, const cplx *syms,
                          size_t nsym, uint8_t *bits);

/*
 * Soft-decision demap: instead of hard bits, emit a log-likelihood
 * ratio per bit (positive favours 0, negative favours 1), suitable
 * for feeding an LDPC or Viterbi soft decoder.
 *   syms      : received symbols, length nsym
 *   noise_var : channel noise variance, used to scale the LLRs
 *   llr       : output LLRs, length nsym * bits_per_symbol
 * Returns the number of LLRs written.
 */
size_t dsp_qam_demodulate_soft(dsp_qam_order order, const cplx *syms,
                               size_t nsym, double noise_var,
                               double *llr);

#endif /* DSP_QAM_H */
