/*
 * qam.c - QPSK, 16-QAM, 64-QAM mapping and demapping.
 *
 * All three are SQUARE constellations: the bits of a symbol split
 * into two halves, one choosing the I (real) coordinate and one the
 * Q (imaginary). Each half is a Gray-coded PAM level, so the 2-D
 * labelling is automatically Gray on both axes.
 *
 * For M-QAM there are L = sqrt(M) levels per axis:
 *   QPSK   -> L = 2,  levels {-1, +1}
 *   16-QAM -> L = 4,  levels {-3,-1,+1,+3}
 *   64-QAM -> L = 8,  levels {-7,-5,-3,-1,+1,+3,+5,+7}
 * before normalisation to unit average energy.
 */
#include "modulation/qam.h"
#include <stdlib.h>

size_t dsp_qam_bits_per_symbol(dsp_qam_order order) {
    switch (order) {
        case DSP_QAM_QPSK: return 2;
        case DSP_QAM_16:   return 4;
        case DSP_QAM_64:   return 6;
        default:           return 0;
    }
}

/* Levels per axis: sqrt(M). */
static size_t axis_levels(dsp_qam_order order) {
    switch (order) {
        case DSP_QAM_QPSK: return 2;
        case DSP_QAM_16:   return 4;
        case DSP_QAM_64:   return 8;
        default:           return 0;
    }
}

/*
 * Normalisation factor: 1/sqrt(mean symbol energy). For a square
 * M-QAM with L levels per axis the mean energy is 2*(M-1)/3.
 */
static double norm_factor(dsp_qam_order order) {
    double M = (double)order;
    return 1.0 / sqrt(2.0 * (M - 1.0) / 3.0);
}

/* Reflected binary (Gray) code of an integer, and its inverse. */
static unsigned to_gray(unsigned x) {
    return x ^ (x >> 1);
}
static unsigned from_gray(unsigned g) {
    unsigned x = 0;
    for (; g; g >>= 1)
        x ^= g;
    return x;
}

/*
 * Map a group of `b` bits (per axis) to a PAM coordinate.
 * The bits select a Gray-code index; the index maps to an odd-integer
 * level centred on zero: {-(L-1), ..., -1, +1, ..., +(L-1)}.
 */
static double bits_to_level(const uint8_t *bits, size_t b) {
    unsigned v = 0;
    for (size_t i = 0; i < b; ++i)
        v = (v << 1) | (bits[i] & 1);
    unsigned idx = from_gray(v);                 /* 0 .. L-1 */
    return 2.0 * (double)idx - ((double)(1u << b) - 1.0);
}

/* Inverse: slice a coordinate to the nearest level, emit its b bits. */
static void level_to_bits(double coord, size_t b, uint8_t *bits) {
    size_t L = (size_t)1u << b;
    /* Nearest odd-integer level -> index 0 .. L-1. */
    long idx = (long)floor((coord + (double)(L - 1)) / 2.0 + 0.5);
    if (idx < 0) idx = 0;
    if (idx > (long)L - 1) idx = (long)L - 1;
    unsigned g = to_gray((unsigned)idx);
    for (size_t i = 0; i < b; ++i)
        bits[i] = (uint8_t)((g >> (b - 1 - i)) & 1);
}

size_t dsp_qam_symbol_count(dsp_qam_order order, size_t nbits) {
    size_t bps = dsp_qam_bits_per_symbol(order);
    return bps ? nbits / bps : 0;
}

size_t dsp_qam_modulate(dsp_qam_order order, const uint8_t *bits,
                        size_t nbits, cplx *syms) {
    size_t bps = dsp_qam_bits_per_symbol(order);
    if (bps == 0 || nbits % bps != 0)
        return 0;

    size_t half = bps / 2;                       /* bits per axis */
    double scale = norm_factor(order);
    size_t nsym = nbits / bps;

    for (size_t s = 0; s < nsym; ++s) {
        const uint8_t *sym_bits = bits + s * bps;
        /* First half of the bits -> I axis, second half -> Q axis. */
        double i = bits_to_level(sym_bits,        half);
        double q = bits_to_level(sym_bits + half, half);
        syms[s] = dsp_cplx(i * scale, q * scale);
    }
    return nsym;
}

size_t dsp_qam_demodulate(dsp_qam_order order, const cplx *syms,
                          size_t nsym, uint8_t *bits) {
    size_t bps = dsp_qam_bits_per_symbol(order);
    if (bps == 0)
        return 0;

    size_t half = bps / 2;
    double scale = norm_factor(order);

    for (size_t s = 0; s < nsym; ++s) {
        /* Undo the energy normalisation, then slice each axis. */
        double i = creal(syms[s]) / scale;
        double q = cimag(syms[s]) / scale;
        level_to_bits(i, half, bits + s * bps);
        level_to_bits(q, half, bits + s * bps + half);
    }
    return nsym * bps;
}

size_t dsp_qam_demodulate_soft(dsp_qam_order order, const cplx *syms,
                               size_t nsym, double noise_var,
                               double *llr) {
    size_t bps = dsp_qam_bits_per_symbol(order);
    if (bps == 0)
        return 0;
    if (noise_var <= 0.0)
        noise_var = 1e-9;

    size_t half  = bps / 2;
    size_t L     = axis_levels(order);
    double scale = norm_factor(order);

    /* For each axis and each bit position, the LLR is the log-ratio
     * of the summed likelihoods of the levels with that bit = 0
     * versus = 1 (the max-log approximation keeps only the nearest
     * level in each set). */
    for (size_t s = 0; s < nsym; ++s) {
        double axis[2];
        axis[0] = creal(syms[s]) / scale;        /* I */
        axis[1] = cimag(syms[s]) / scale;        /* Q */

        for (int ax = 0; ax < 2; ++ax) {
            double y = axis[ax];
            for (size_t bit = 0; bit < half; ++bit) {
                double best0 = 1e300, best1 = 1e300;
                /* Walk every level on this axis. */
                for (size_t idx = 0; idx < L; ++idx) {
                    double level = 2.0 * (double)idx - (double)(L - 1);
                    double d = y - level;
                    double dist = d * d;          /* squared distance */
                    unsigned g = to_gray((unsigned)idx);
                    int b = (g >> (half - 1 - bit)) & 1;
                    if (b == 0) {
                        if (dist < best0) best0 = dist;
                    } else {
                        if (dist < best1) best1 = dist;
                    }
                }
                /* Max-log LLR: (d1^2 - d0^2) / (2*noise_var). Positive
                 * favours bit 0, matching the LDPC/Viterbi convention. */
                double L_val = (best1 - best0) / (2.0 * noise_var);
                size_t pos = s * bps
                           + (size_t)ax * half + bit;
                llr[pos] = L_val;
            }
        }
    }
    return nsym * bps;
}
