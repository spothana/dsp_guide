/*
 * dft.c - Direct O(N^2) Discrete Fourier Transform.
 *
 * This is the textbook definition computed literally. The nested loop
 * makes the O(N^2) cost explicit: the outer loop walks the N output
 * bins, the inner loop sums over the N input samples.
 */
#include "transforms/dft.h"
#include <stdlib.h>
#include <string.h>

void dsp_dft(const cplx *in, cplx *out, size_t n) {
    /* Alias-safe: if out == in, writing out[k] would corrupt later
     * reads of in[t]. Snapshot the input into a scratch buffer first. */
    cplx *src = malloc(n * sizeof(cplx));
    if (!src) return;
    memcpy(src, in, n * sizeof(cplx));

    for (size_t k = 0; k < n; ++k) {
        cplx acc = 0.0;
        for (size_t t = 0; t < n; ++t) {
            /* Reference (twiddle) wave for frequency bin k at sample t. */
            double angle = -2.0 * M_PI * (double)k * (double)t / (double)n;
            acc += src[t] * cexp(angle * I);
        }
        out[k] = acc;
    }
    free(src);
}

void dsp_idft(const cplx *in, cplx *out, size_t n) {
    cplx *src = malloc(n * sizeof(cplx));
    if (!src) return;
    memcpy(src, in, n * sizeof(cplx));

    for (size_t t = 0; t < n; ++t) {
        cplx acc = 0.0;
        for (size_t k = 0; k < n; ++k) {
            /* Inverse uses +j (counter-clockwise) and a 1/N scale. */
            double angle = 2.0 * M_PI * (double)k * (double)t / (double)n;
            acc += src[k] * cexp(angle * I);
        }
        out[t] = acc / (double)n;
    }
    free(src);
}
