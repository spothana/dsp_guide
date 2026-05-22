/*
 * dct.c - Orthonormal DCT-II (forward) and DCT-III (inverse).
 *
 * Orthonormal scaling is used so the transform is energy-preserving
 * and the inverse is simply the transpose. The alpha(k) factor below
 * is 1/sqrt(N) for the DC coefficient and sqrt(2/N) otherwise.
 */
#include "transforms/dct.h"

/* Orthonormal scale factor for coefficient index k. */
static double alpha(size_t k, size_t n) {
    return (k == 0) ? sqrt(1.0 / (double)n)
                    : sqrt(2.0 / (double)n);
}

void dsp_dct(const double *in, double *out, size_t n) {
    for (size_t k = 0; k < n; ++k) {
        double acc = 0.0;
        for (size_t t = 0; t < n; ++t) {
            acc += in[t] * cos(M_PI * (2.0 * (double)t + 1.0)
                                     * (double)k / (2.0 * (double)n));
        }
        out[k] = alpha(k, n) * acc;
    }
}

void dsp_idct(const double *in, double *out, size_t n) {
    for (size_t t = 0; t < n; ++t) {
        double acc = 0.0;
        for (size_t k = 0; k < n; ++k) {
            acc += alpha(k, n) * in[k]
                 * cos(M_PI * (2.0 * (double)t + 1.0)
                            * (double)k / (2.0 * (double)n));
        }
        out[t] = acc;
    }
}
