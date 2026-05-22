/*
 * ldpc_decode.c - LDPC decoders: bit-flipping and sum-product.
 *
 * Both decoders pass messages along the edges of the Tanner graph.
 * The difference is what the messages carry: bit-flipping passes hard
 * votes, sum-product passes probabilities (log-likelihood ratios).
 */
#include "coding/ldpc.h"
#include <stdlib.h>
#include <math.h>

/* Parity of check i under the current bit estimate. */
static int check_parity(const dsp_ldpc *code, size_t i,
                        const uint8_t *bits) {
    int p = 0;
    for (size_t d = 0; d < code->row_deg[i]; ++d)
        p ^= (bits[code->row_var[i][d]] & 1);
    return p;
}

/* ===================================================================
 * Bit-flipping decoder (hard decision)
 *
 * Each pass: compute every check's parity. A failing check is
 * "unhappy" with all the bits on it. Each bit tallies how many of its
 * checks are unhappy; the single bit with the most unhappy checks is
 * flipped. Repeat until every check is satisfied.
 *
 * Flipping exactly one bit per iteration (Gallager's Algorithm A in
 * spirit) is the stable choice - flipping every maximal-score bit at
 * once tends to overshoot and oscillate on short codes.
 * =================================================================== */

int dsp_ldpc_decode_bitflip(const dsp_ldpc *code, uint8_t *recv,
                            int max_iter) {
    size_t *votes = malloc(code->n * sizeof(size_t));
    uint8_t *fail = malloc(code->m * sizeof(uint8_t));
    if (!votes || !fail) {
        free(votes); free(fail);
        return -1;
    }

    int result = -1;

    for (int iter = 1; iter <= max_iter; ++iter) {
        /* Evaluate every parity check. */
        size_t nfail = 0;
        for (size_t i = 0; i < code->m; ++i) {
            fail[i] = (uint8_t)check_parity(code, i, recv);
            if (fail[i]) ++nfail;
        }
        if (nfail == 0) {           /* all checks satisfied -> done */
            result = iter - 1;
            break;
        }

        /* Each bit counts how many of its checks currently fail. */
        size_t worst = 0, worst_bit = 0;
        for (size_t j = 0; j < code->n; ++j) {
            size_t v = 0;
            for (size_t d = 0; d < code->col_deg[j]; ++d)
                if (fail[code->col_chk[j][d]])
                    ++v;
            votes[j] = v;
            if (v > worst) { worst = v; worst_bit = j; }
        }

        if (worst == 0)             /* no bit can improve -> stuck */
            break;

        /* Flip the single most-suspect bit. */
        recv[worst_bit] ^= 1;
    }

    free(votes);
    free(fail);
    return result;
}

/* ===================================================================
 * Soft-decision decoders: sum-product and min-sum
 *
 * Both share the same message-passing schedule on the Tanner graph.
 * Messages are log-likelihood ratios (LLRs): a real number whose sign
 * is the believed bit (positive -> 0, negative -> 1) and whose
 * magnitude is the confidence. Two message types flow each iteration:
 *
 *   var->check  M[j][i] : variable j tells check i its belief,
 *                         excluding what check i last said (so a
 *                         check never just hears its own echo).
 *   check->var  E[i][j] : check i tells variable j what the parity
 *                         constraint implies, from the OTHER bits.
 *
 * The ONLY difference between the two decoders is the check-node rule:
 *   sum-product : exact - combine LLRs via tanh(x/2), multiply, then
 *                 invert with 2*atanh. Best performance.
 *   min-sum     : approximate - outgoing magnitude is the minimum of
 *                 the incoming magnitudes, sign is their product. No
 *                 transcendental functions; this is what real LDPC
 *                 hardware uses. A scale factor < 1 corrects min-sum's
 *                 tendency to overestimate (normalised min-sum).
 * =================================================================== */

/* Shared decoder core. `use_minsum` selects the check-node rule;
 * `scale` attenuates the check messages (min-sum only; 1.0 otherwise). */
static int ldpc_soft_decode(const dsp_ldpc *code, const double *llr,
                            uint8_t *out, int max_iter,
                            int use_minsum, double scale) {
    size_t m = code->m, n = code->n;

    /* Edge-indexed message stores. Each variable j owns col_deg[j]
     * edges; we lay them out contiguously and index by an offset. */
    size_t *voff = malloc((n + 1) * sizeof(size_t));
    size_t *coff = malloc((m + 1) * sizeof(size_t));
    if (!voff || !coff) { free(voff); free(coff); return -1; }

    voff[0] = 0;
    for (size_t j = 0; j < n; ++j)
        voff[j + 1] = voff[j] + code->col_deg[j];
    coff[0] = 0;
    for (size_t i = 0; i < m; ++i)
        coff[i + 1] = coff[i] + code->row_deg[i];

    size_t nedge = voff[n];

    /* M : var->check messages, indexed by the variable's edge slot.
     * E : check->var messages, indexed by the check's edge slot.
     * To move values between the two layouts we precompute both
     * directions of the edge correspondence. */
    double *M = malloc(nedge * sizeof(double));
    double *E = malloc(nedge * sizeof(double));
    size_t *v2c = malloc(nedge * sizeof(size_t));  /* var slot -> chk slot */
    size_t *c2v = malloc(nedge * sizeof(size_t));  /* chk slot -> var slot */
    if (!M || !E || !v2c || !c2v) {
        free(voff); free(coff); free(M); free(E); free(v2c); free(c2v);
        return -1;
    }

    /* Build the var-edge <-> check-edge correspondence. For variable
     * j's d-th check i, find which slot of check i points back at j,
     * then record the mapping in both directions. */
    for (size_t j = 0; j < n; ++j) {
        for (size_t d = 0; d < code->col_deg[j]; ++d) {
            size_t i = code->col_chk[j][d];
            size_t slot = 0;
            for (size_t e = 0; e < code->row_deg[i]; ++e)
                if (code->row_var[i][e] == j) { slot = e; break; }
            v2c[voff[j] + d] = coff[i] + slot;
            c2v[coff[i] + slot] = voff[j] + d;
        }
    }

    /* Initialise: every var->check message is just the channel LLR. */
    for (size_t j = 0; j < n; ++j)
        for (size_t d = 0; d < code->col_deg[j]; ++d)
            M[voff[j] + d] = llr[j];

    int result = -1;

    for (int iter = 1; iter <= max_iter; ++iter) {
        /* --- Check-node update ---
         * c2v maps each check edge directly to its var-side M slot. */
        for (size_t i = 0; i < m; ++i) {
            size_t deg = code->row_deg[i];
            for (size_t d = 0; d < deg; ++d) {
                if (use_minsum) {
                    /* Min-sum: sign = product of signs, magnitude =
                     * minimum magnitude, over all other edges. */
                    double min_mag = 1e300;
                    int    sign    = 1;
                    for (size_t e = 0; e < deg; ++e) {
                        if (e == d) continue;
                        double v = M[c2v[coff[i] + e]];
                        double a = (v < 0.0) ? -v : v;
                        if (a < min_mag) min_mag = a;
                        if (v < 0.0) sign = -sign;
                    }
                    E[coff[i] + d] = scale * sign * min_mag;
                } else {
                    /* Sum-product: exact tanh rule.
                     * E = 2 * atanh( product of tanh(M/2) ). */
                    double prod = 1.0;
                    for (size_t e = 0; e < deg; ++e) {
                        if (e == d) continue;
                        prod *= tanh(M[c2v[coff[i] + e]] / 2.0);
                    }
                    /* Clamp to keep atanh finite. */
                    if (prod >  0.999999999) prod =  0.999999999;
                    if (prod < -0.999999999) prod = -0.999999999;
                    E[coff[i] + d] = 2.0 * atanh(prod);
                }
            }
        }

        /* --- Variable-node update ---
         * M[j][i] = channel LLR + sum of incoming E from OTHER checks.
         * Total belief = channel LLR + sum of ALL incoming E. */
        for (size_t j = 0; j < n; ++j) {
            size_t deg = code->col_deg[j];

            double total = llr[j];
            for (size_t d = 0; d < deg; ++d)
                total += E[v2c[voff[j] + d]];

            /* Hard decision from the total belief. */
            out[j] = (total < 0.0) ? 1 : 0;

            /* Extrinsic message excludes the recipient check's own E. */
            for (size_t d = 0; d < deg; ++d)
                M[voff[j] + d] = total - E[v2c[voff[j] + d]];
        }

        /* Converged once the hard decision is a valid codeword. */
        if (dsp_ldpc_check(code, out)) {
            result = iter;
            break;
        }
    }

    free(voff); free(coff);
    free(M); free(E); free(v2c); free(c2v);
    return result;
}

/* Public wrappers: exact sum-product and approximate min-sum. */

int dsp_ldpc_decode_sumproduct(const dsp_ldpc *code, const double *llr,
                               uint8_t *out, int max_iter) {
    return ldpc_soft_decode(code, llr, out, max_iter,
                            /*use_minsum=*/0, /*scale=*/1.0);
}

int dsp_ldpc_decode_minsum(const dsp_ldpc *code, const double *llr,
                           uint8_t *out, double scale, int max_iter) {
    if (scale <= 0.0 || scale > 1.0)
        scale = 1.0;
    return ldpc_soft_decode(code, llr, out, max_iter,
                            /*use_minsum=*/1, scale);
}

double dsp_ldpc_awgn_llr(double sample, double noise_var) {
    /* For antipodal signalling (bit 0 -> +1, bit 1 -> -1) on an AWGN
     * channel, the channel LLR is 2*y/sigma^2. */
    if (noise_var <= 0.0)
        noise_var = 1e-9;
    return 2.0 * sample / noise_var;
}

/* ===================================================================
 * BER sweep - characterise a decoder over a simulated AWGN channel
 * =================================================================== */

/* xorshift RNG, local to this file. */
static unsigned ber_rng(unsigned *s) {
    unsigned x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}

/* Uniform double in [0, 1). */
static double ber_uniform(unsigned *s) {
    return (ber_rng(s) >> 8) / 16777216.0;   /* 24 random bits */
}

/* One standard-normal sample via the Box-Muller transform. */
static double ber_gauss(unsigned *s) {
    static const double TWO_PI = 6.28318530717958647692;
    double u1 = ber_uniform(s);
    double u2 = ber_uniform(s);
    if (u1 < 1e-12) u1 = 1e-12;              /* avoid log(0) */
    return sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
}

double dsp_ldpc_ber_sweep(const dsp_ldpc *code, dsp_ldpc_decoder decoder,
                          double noise_std, int trials, int max_iter,
                          unsigned seed) {
    size_t n = code->n;
    uint8_t *recv = malloc(n * sizeof(uint8_t));
    uint8_t *dec  = malloc(n * sizeof(uint8_t));
    double  *llr  = malloc(n * sizeof(double));
    if (!recv || !dec || !llr) {
        free(recv); free(dec); free(llr);
        return -1.0;
    }

    unsigned state = seed ? seed : 0xC0FFEEu;
    double   var   = noise_std * noise_std;
    size_t   total_bits = 0, bit_errors = 0;

    for (int t = 0; t < trials; ++t) {
        /* Transmit the all-zero codeword as antipodal +1 symbols,
         * add Gaussian noise, and form per-bit channel LLRs. */
        for (size_t j = 0; j < n; ++j) {
            double sample = 1.0 + noise_std * ber_gauss(&state);
            llr[j]  = dsp_ldpc_awgn_llr(sample, var);
            recv[j] = (sample < 0.0) ? 1 : 0;   /* hard slice */
        }

        switch (decoder) {
            case DSP_LDPC_BITFLIP:
                dsp_ldpc_decode_bitflip(code, recv, max_iter);
                for (size_t j = 0; j < n; ++j) dec[j] = recv[j];
                break;
            case DSP_LDPC_SUMPRODUCT:
                dsp_ldpc_decode_sumproduct(code, llr, dec, max_iter);
                break;
            case DSP_LDPC_MINSUM:
                dsp_ldpc_decode_minsum(code, llr, dec, 0.75, max_iter);
                break;
        }

        /* The transmitted codeword was all-zero, so any 1 is an error. */
        for (size_t j = 0; j < n; ++j)
            if (dec[j]) ++bit_errors;
        total_bits += n;
    }

    free(recv); free(dec); free(llr);
    return total_bits ? (double)bit_errors / (double)total_bits : 0.0;
}
