/*
 * convolution.c - Direct and FFT-based convolution.
 */
#include "operations/convolution.h"
#include "transforms/fft.h"
#include <stdlib.h>

void dsp_convolve(const double *x, size_t n,
                  const double *h, size_t m,
                  double *y) {
    size_t out = dsp_conv_len(n, m);
    for (size_t i = 0; i < out; ++i)
        y[i] = 0.0;

    for (size_t i = 0; i < n; ++i)
        for (size_t k = 0; k < m; ++k)
            y[i + k] += x[i] * h[k];
}

int dsp_convolve_fft(const double *x, size_t n,
                     const double *h, size_t m,
                     double *y) {
    size_t out = dsp_conv_len(n, m);
    if (out == 0) return 0;

    /* Pad to a power of two large enough to hold the linear result.
     * (Padding past the linear length avoids circular wraparound.) */
    size_t fftn = dsp_next_pow2(out);

    cplx *xa = calloc(fftn, sizeof(cplx));
    cplx *ha = calloc(fftn, sizeof(cplx));
    if (!xa || !ha) { free(xa); free(ha); return -1; }

    for (size_t i = 0; i < n; ++i) xa[i] = x[i];
    for (size_t i = 0; i < m; ++i) ha[i] = h[i];

    dsp_fft(xa, fftn);
    dsp_fft(ha, fftn);

    /* Convolution in time == pointwise multiply in frequency. */
    for (size_t i = 0; i < fftn; ++i)
        xa[i] *= ha[i];

    dsp_ifft(xa, fftn);

    for (size_t i = 0; i < out; ++i)
        y[i] = creal(xa[i]);

    free(xa);
    free(ha);
    return 0;
}
