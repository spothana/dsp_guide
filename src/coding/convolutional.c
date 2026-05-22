/*
 * convolutional.c - Rate-1/2, K=3 convolutional code with Viterbi.
 *
 * The encoder runs the input bit stream through a 2-bit shift register
 * and produces two output bits per input bit using the generator
 * polynomials G1 = 0o7 (111) and G2 = 0o5 (101). After the data, two
 * zero "tail" bits flush the register back to state 0.
 *
 * The Viterbi decoder finds the maximum-likelihood input sequence by
 * tracking, for every encoder state at every time step, the lowest-
 * cost path that reaches it (a trellis), then tracing the best path
 * back. Hard-decision uses Hamming distance; soft-decision uses a
 * correlation metric against real-valued channel confidences.
 */
#include "coding/correct.h"
#include <stdlib.h>
#include <string.h>

#define K        3            /* constraint length */
#define NSTATES  4            /* 2^(K-1) encoder states */
#define G1       7            /* generator polynomial 1 (octal 7) */
#define G2       5            /* generator polynomial 2 (octal 5) */

/* Parity (mod-2 sum of set bits) of a small integer. */
static int parity_of(unsigned v) {
    int p = 0;
    while (v) { p ^= (v & 1); v >>= 1; }
    return p;
}

/*
 * Encoder output for a given (state, input bit).
 * The 3-bit register content is (input << 2) | state.
 * Returns the two output bits packed as (o1 << 1) | o2.
 */
static int encode_step(int state, int bit) {
    unsigned reg = ((unsigned)bit << 2) | (unsigned)state;
    int o1 = parity_of(reg & G1);
    int o2 = parity_of(reg & G2);
    return (o1 << 1) | o2;
}

/* Next state after shifting `bit` into the register. */
static int next_state(int state, int bit) {
    return ((state >> 1) | (bit << 1)) & (NSTATES - 1);
}

size_t dsp_conv_encode(const uint8_t *bits, size_t nbits, uint8_t *out) {
    int state = 0;
    size_t o = 0;

    /* Data bits followed by 2 zero flush bits. */
    for (size_t i = 0; i < nbits + 2; ++i) {
        int bit = (i < nbits) ? (bits[i] & 1) : 0;
        int sym = encode_step(state, bit);
        out[o++] = (uint8_t)((sym >> 1) & 1);
        out[o++] = (uint8_t)(sym & 1);
        state = next_state(state, bit);
    }
    return o;
}

/*
 * Shared Viterbi engine. The branch metric is supplied as a callback
 * so the same trellis code serves both hard and soft decoding.
 *   recv     : received data (uint8_t* for hard, double* for soft)
 *   nsym     : number of received bits (must be even)
 *   out      : decoded message bits
 *   branch   : cost of expecting symbol bits (e1,e2) at trellis step s
 * Returns the number of decoded message bits.
 */
static size_t viterbi_run(size_t nsym, uint8_t *out,
                          double (*branch)(size_t step, int e1, int e2,
                                           const void *ctx),
                          const void *ctx) {
    if (nsym < 4 || (nsym & 1))
        return 0;

    size_t steps = nsym / 2;            /* trellis depth */
    const double BIG = 1e18;

    /* path_metric[state], plus a back-pointer trellis. */
    double *pm   = malloc(NSTATES * sizeof(double));
    double *npm  = malloc(NSTATES * sizeof(double));
    /* prev[step][state] = predecessor state; bit[step][state] = input. */
    int *prev = malloc(steps * NSTATES * sizeof(int));
    int *ibit = malloc(steps * NSTATES * sizeof(int));
    if (!pm || !npm || !prev || !ibit) {
        free(pm); free(npm); free(prev); free(ibit);
        return 0;
    }

    /* The encoder starts in state 0. */
    for (int s = 0; s < NSTATES; ++s)
        pm[s] = (s == 0) ? 0.0 : BIG;

    for (size_t step = 0; step < steps; ++step) {
        for (int s = 0; s < NSTATES; ++s)
            npm[s] = BIG;

        for (int s = 0; s < NSTATES; ++s) {
            if (pm[s] >= BIG) continue;
            for (int bit = 0; bit < 2; ++bit) {
                int sym = encode_step(s, bit);
                int e1  = (sym >> 1) & 1;
                int e2  = sym & 1;
                int ns  = next_state(s, bit);

                double cost = pm[s] + branch(step, e1, e2, ctx);
                if (cost < npm[ns]) {
                    npm[ns] = cost;
                    prev[step * NSTATES + ns] = s;
                    ibit[step * NSTATES + ns] = bit;
                }
            }
        }
        double *t = pm; pm = npm; npm = t;
    }

    /* The encoder was flushed to state 0 - trace back from there. */
    int state = 0;
    size_t ndata = steps - 2;            /* drop the 2 tail bits */
    for (size_t step = steps; step-- > 0; ) {
        int bit = ibit[step * NSTATES + state];
        int ps  = prev[step * NSTATES + state];
        if (step < ndata)
            out[step] = (uint8_t)bit;
        state = ps;
    }

    free(pm); free(npm); free(prev); free(ibit);
    return ndata;
}

/* ---- Hard-decision branch metric: Hamming distance ------------------ */

static double branch_hard(size_t step, int e1, int e2, const void *ctx) {
    const uint8_t *recv = (const uint8_t *)ctx;
    int r1 = recv[2 * step]     & 1;
    int r2 = recv[2 * step + 1] & 1;
    return (double)((r1 ^ e1) + (r2 ^ e2));
}

size_t dsp_viterbi_decode(const uint8_t *recv, size_t nsym, uint8_t *out) {
    return viterbi_run(nsym, out, branch_hard, recv);
}

/* ---- Soft-decision branch metric: correlation ----------------------- */

static double branch_soft(size_t step, int e1, int e2, const void *ctx) {
    const double *recv = (const double *)ctx;
    /* Map expected bit {0,1} to the antipodal symbol {+1,-1}. A larger
     * agreement with the soft value lowers the cost. */
    double exp1 = e1 ? -1.0 : 1.0;
    double exp2 = e2 ? -1.0 : 1.0;
    return -(recv[2 * step] * exp1 + recv[2 * step + 1] * exp2);
}

size_t dsp_viterbi_decode_soft(const double *recv, size_t nsym,
                               uint8_t *out) {
    return viterbi_run(nsym, out, branch_soft, recv);
}
