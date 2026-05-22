/*
 * fft.c - Iterative radix-2 Cooley-Tukey FFT.
 *
 * Two phases:
 *   1. Bit-reversal permutation. Repeatedly splitting samples by
 *      even/odd index is equivalent to sorting them by the bit-reversed
 *      value of their index. Doing this permutation up front lets every
 *      butterfly read and write the same two slots -> the FFT runs
 *      fully in place, needing no scratch array.
 *   2. log2(N) stages of butterflies. Stage s combines sub-DFTs of
 *      size `len`, doubling len each stage (2, 4, 8, ... N).
 */
#include "transforms/fft.h"

/* Reverse the low `bits` bits of x. Used to build the input permutation. */
static size_t reverse_bits(size_t x, int bits) {
    size_t r = 0;
    for (int i = 0; i < bits; ++i) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

/* Core engine: sign = -1 for forward transform, +1 for inverse. */
static int fft_core(cplx *data, size_t n, int sign) {
    if (!dsp_is_pow2(n)) return -1;
    if (n < 2) return 0;

    /* How many bits index n elements. */
    int bits = 0;
    for (size_t t = n; t > 1; t >>= 1) ++bits;

    /* --- Phase 1: bit-reversal permutation (in place) --- */
    for (size_t i = 0; i < n; ++i) {
        size_t j = reverse_bits(i, bits);
        if (j > i) {
            cplx tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
        }
    }

    /* --- Phase 2: log2(N) butterfly stages --- */
    for (size_t len = 2; len <= n; len <<= 1) {
        /* Twiddle step for this stage: exp(sign * -j*2*pi/len). */
        double theta = (double)sign * -2.0 * M_PI / (double)len;
        cplx wstep = cexp(theta * I);

        for (size_t start = 0; start < n; start += len) {
            cplx w = 1.0;
            for (size_t k = 0; k < len / 2; ++k) {
                /* One butterfly on the pair (a, b). */
                cplx a = data[start + k];
                cplx b = data[start + k + len / 2] * w;
                data[start + k]           = a + b;   /* + output */
                data[start + k + len / 2] = a - b;   /* - output */
                w *= wstep;
            }
        }
    }
    return 0;
}

int dsp_fft(cplx *data, size_t n) {
    return fft_core(data, n, -1);
}

int dsp_ifft(cplx *data, size_t n) {
    int rc = fft_core(data, n, +1);
    if (rc != 0) return rc;
    for (size_t i = 0; i < n; ++i)
        data[i] /= (double)n;   /* 1/N normalisation */
    return 0;
}
