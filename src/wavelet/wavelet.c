/*
 * wavelet.c - Haar discrete wavelet transform (Mallat pyramid).
 *
 * The Haar wavelet pairs adjacent samples:
 *   approximation = (a + b) / sqrt(2)   (a running average)
 *   detail        = (a - b) / sqrt(2)   (a local difference)
 * The 1/sqrt(2) factor makes the transform orthonormal, so it
 * preserves energy and the inverse is exact.
 */
#include "wavelet/wavelet.h"
#include <stdlib.h>

static const double INV_SQRT2 = 0.70710678118654752440;

void dsp_dwt_haar_step(const double *in, size_t n,
                       double *approx, double *detail) {
    size_t half = n / 2;
    for (size_t i = 0; i < half; ++i) {
        double a = in[2 * i];
        double b = in[2 * i + 1];
        approx[i] = (a + b) * INV_SQRT2;
        detail[i] = (a - b) * INV_SQRT2;
    }
}

void dsp_idwt_haar_step(const double *approx, const double *detail,
                        size_t half, double *out) {
    for (size_t i = 0; i < half; ++i) {
        double s = approx[i];
        double d = detail[i];
        /* Invert the 2x2 orthonormal Haar mixing matrix. */
        out[2 * i]     = (s + d) * INV_SQRT2;
        out[2 * i + 1] = (s - d) * INV_SQRT2;
    }
}

int dsp_dwt_haar(double *data, size_t n) {
    if (!dsp_is_pow2(n)) return -1;
    if (n < 2) return 0;

    double *tmp = malloc(n * sizeof(double));
    if (!tmp) return -1;

    int levels = 0;
    /* Each pass works on the leading `len` samples (the running
     * approximation) and shrinks len by half. */
    for (size_t len = n; len >= 2; len >>= 1) {
        size_t half = len / 2;
        dsp_dwt_haar_step(data, len, tmp, tmp + half);
        for (size_t i = 0; i < len; ++i)
            data[i] = tmp[i];
        ++levels;
    }

    free(tmp);
    return levels;
}

int dsp_idwt_haar(double *data, size_t n) {
    if (!dsp_is_pow2(n)) return -1;
    if (n < 2) return 0;

    double *tmp = malloc(n * sizeof(double));
    if (!tmp) return -1;

    /* Reverse the forward pyramid: start from the coarsest pair and
     * grow len back up to n. */
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len / 2;
        dsp_idwt_haar_step(data, data + half, half, tmp);
        for (size_t i = 0; i < len; ++i)
            data[i] = tmp[i];
    }

    free(tmp);
    return 0;
}
