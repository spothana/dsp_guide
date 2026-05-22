/*
 * ldpc.h - Low-Density Parity-Check codes
 *
 * PROBLEM SOLVED
 *   LDPC codes get within a fraction of a decibel of the Shannon limit
 *   - the theoretical best any code can do. They are the FEC behind
 *   5G, Wi-Fi 6, and modern satellite and storage standards.
 *
 * THE IDEA
 *   An LDPC code is defined by a parity-check matrix H that is SPARSE
 *   (mostly zeros - the "low-density" in the name). A vector c is a
 *   valid codeword exactly when H * c = 0 over GF(2).
 *
 *   H is best pictured as a bipartite TANNER GRAPH:
 *     - one VARIABLE node per code bit       (columns of H)
 *     - one CHECK node per parity equation   (rows of H)
 *     - an edge wherever H has a 1.
 *   Sparsity means few edges, which is what makes iterative decoding
 *   both fast and effective.
 *
 * DECODING - why this module has two decoders
 *   Bit-flipping (hard decision): each check votes on whether its bits
 *     look consistent; the variable bits in the most unhappy checks
 *     get flipped. Simple, integer-only, good for intuition. This is
 *     the LDPC analogue of hard-decision Viterbi.
 *   Sum-product / belief propagation (soft decision): variable and
 *     check nodes exchange PROBABILITIES (as log-likelihood ratios)
 *     along the graph edges, refining their belief each iteration.
 *     This is the real-world decoder and the source of LDPC's
 *     near-Shannon performance - the same soft-information principle
 *     as soft-decision Viterbi, run iteratively on a graph.
 *
 * vs OTHER CODES IN THIS GUIDE
 *   Reed-Solomon corrects a guaranteed t errors with exact algebra.
 *   LDPC has no such hard guarantee - it is probabilistic - but on a
 *   noisy channel it corrects far more, far closer to the limit. The
 *   cost is iterative decoding: more computation and some latency, the
 *   classic trade-off the guide keeps returning to.
 *
 * SCOPE
 *   This is a compact, teaching-oriented implementation: small codes,
 *   a regular-code generator, both decoders. Production LDPC uses
 *   large structured matrices (e.g. quasi-cyclic) and fixed-point
 *   min-sum decoding; the principles are exactly the ones here.
 */
#ifndef DSP_LDPC_H
#define DSP_LDPC_H

#include <stddef.h>
#include <stdint.h>

/*
 * An LDPC code: a sparse m x n parity-check matrix H.
 *   n = code length (variable nodes), m = parity checks (check nodes).
 *   The code carries k = n - m information bits (when H is full rank).
 *
 * H is stored twice, as adjacency lists, for fast iteration in both
 * directions - this is the natural sparse representation of a Tanner
 * graph:
 *   row_var[i]  : variable nodes touching check i
 *   col_chk[j]  : check nodes touching variable j
 */
typedef struct {
    size_t n;            /* number of code bits / variable nodes */
    size_t m;            /* number of parity checks / check nodes */

    size_t **row_var;    /* row_var[i] = vars in check i */
    size_t  *row_deg;    /* row_deg[i] = number of those vars */

    size_t **col_chk;    /* col_chk[j] = checks on variable j */
    size_t  *col_deg;    /* col_deg[j] = number of those checks */
} dsp_ldpc;

/*
 * Build an LDPC code from a dense 0/1 parity-check matrix.
 *   H   : row-major m*n array of 0/1 bytes
 *   m,n : dimensions
 * Returns 0 on success, -1 on allocation failure. Pair with
 * dsp_ldpc_free.
 */
int dsp_ldpc_from_matrix(dsp_ldpc *code, const uint8_t *H,
                         size_t m, size_t n);

/*
 * Generate a regular (wc, wr) LDPC code: every variable node has
 * degree wc, every check node has degree wr. Requires n*wc == m*wr.
 * This is the classic Gallager-style random construction.
 *   m,n  : dimensions
 *   wc   : variable-node (column) degree
 *   wr   : check-node (row) degree
 *   seed : RNG seed for reproducibility
 * Returns 0 on success, -1 on bad parameters or allocation failure.
 */
int dsp_ldpc_make_regular(dsp_ldpc *code, size_t m, size_t n,
                          size_t wc, size_t wr, unsigned seed);

/* Release all memory owned by an LDPC code. */
void dsp_ldpc_free(dsp_ldpc *code);

/*
 * Syndrome check: returns 1 if `bits` (length n) satisfies every
 * parity check (H * bits == 0), 0 otherwise. A valid codeword and a
 * correctly decoded word both return 1.
 */
int dsp_ldpc_check(const dsp_ldpc *code, const uint8_t *bits);

/*
 * Number of unsatisfied parity checks for `bits` - the syndrome
 * weight. Zero means a valid codeword; larger means more corruption.
 */
size_t dsp_ldpc_syndrome_weight(const dsp_ldpc *code,
                                const uint8_t *bits);

/*
 * Bit-flipping (hard-decision) decoder.
 *   recv     : received hard bits, length n - corrected in place
 *   max_iter : iteration cap
 * Returns the iteration count on success (syndrome reached zero), or
 * -1 if still unsatisfied after max_iter passes.
 */
int dsp_ldpc_decode_bitflip(const dsp_ldpc *code, uint8_t *recv,
                            int max_iter);

/*
 * Sum-product (belief-propagation) soft-decision decoder.
 *   llr      : channel log-likelihood ratios, length n. Convention:
 *              llr[j] > 0 favours bit 0, llr[j] < 0 favours bit 1.
 *   out      : decoded hard bits, length n
 *   max_iter : iteration cap
 * Returns the iteration count when the syndrome reaches zero, or -1
 * if max_iter passes without a valid codeword (out still holds the
 * best estimate).
 */
int dsp_ldpc_decode_sumproduct(const dsp_ldpc *code, const double *llr,
                               uint8_t *out, int max_iter);

/*
 * Min-sum soft-decision decoder.
 *
 * Same message-passing schedule as sum-product, but the check-node
 * update replaces the exact tanh/atanh rule with its standard
 * approximation: the outgoing magnitude is the MINIMUM of the incoming
 * magnitudes, and the sign is their product. This drops all
 * transcendental functions, which is why min-sum (with a scaling
 * correction) is what real LDPC hardware - 5G, Wi-Fi 6 - actually
 * runs. It is faster and fixed-point friendly, at a small cost in
 * correction performance versus exact sum-product.
 *
 *   llr      : channel log-likelihood ratios, length n
 *   out      : decoded hard bits, length n
 *   scale    : attenuation applied to each check message, in (0, 1].
 *              Min-sum overestimates magnitudes; ~0.75 (normalised
 *              min-sum) recovers most of the gap to sum-product. Pass
 *              1.0 for plain, unscaled min-sum.
 *   max_iter : iteration cap
 * Returns the iteration count at convergence, or -1 on failure.
 */
int dsp_ldpc_decode_minsum(const dsp_ldpc *code, const double *llr,
                           uint8_t *out, double scale, int max_iter);

/*
 * Helper: convert a received antipodal sample (nominally +1 for bit 0,
 * -1 for bit 1) plus a noise variance into a channel LLR for an AWGN
 * channel. llr = 2 * sample / variance.
 */
double dsp_ldpc_awgn_llr(double sample, double noise_var);

/* Decoder selector for the BER sweep below. */
typedef enum {
    DSP_LDPC_BITFLIP,      /* hard-decision bit-flipping */
    DSP_LDPC_SUMPRODUCT,   /* exact belief propagation */
    DSP_LDPC_MINSUM        /* min-sum approximation */
} dsp_ldpc_decoder;

/*
 * Measure decoder performance: transmit the all-zero codeword over a
 * simulated AWGN channel `trials` times and return the resulting
 * bit-error rate (fraction of decoded bits still wrong).
 *
 *   decoder   : which decoder to exercise
 *   noise_std : channel noise standard deviation (higher = worse SNR)
 *   trials    : number of independent codewords to simulate
 *   max_iter  : per-codeword iteration cap
 *   seed      : RNG seed for reproducibility
 *
 * This is how LDPC codes are characterised in practice - sweeping
 * noise_std and plotting BER reveals the steep "waterfall" curve.
 * Returns the bit-error rate in [0, 1], or -1.0 on allocation failure.
 */
double dsp_ldpc_ber_sweep(const dsp_ldpc *code, dsp_ldpc_decoder decoder,
                          double noise_std, int trials, int max_iter,
                          unsigned seed);

#endif /* DSP_LDPC_H */
