/*
 * resample.c - Decimation, interpolation, and rational resampling.
 *
 * Both primitives reuse the FIR module: dsp_fir_design_lowpass builds
 * the anti-alias / anti-image filter, dsp_fir_apply runs it.
 */
#include "sampling/resample.h"
#include "filtering/fir.h"
#include <stdlib.h>

/* Number of taps for the internal low-pass filters. Odd -> linear
 * phase Type-I FIR. Long enough for a reasonably sharp transition. */
#define RESAMPLE_TAPS 65

size_t dsp_decimate(const double *x, size_t n, size_t m, double *y) {
    if (m == 0 || n == 0) return 0;
    if (m == 1) {                       /* nothing to do */
        for (size_t i = 0; i < n; ++i) y[i] = x[i];
        return n;
    }

    double taps[RESAMPLE_TAPS];
    double *filtered = malloc(n * sizeof(double));
    if (!filtered) return 0;

    /* Anti-aliasing low-pass: cutoff just below the new Nyquist. */
    dsp_fir_design_lowpass(taps, RESAMPLE_TAPS, 0.5 / (double)m);
    dsp_fir_apply(x, filtered, n, taps, RESAMPLE_TAPS);

    /* Keep every m-th sample. */
    size_t out = dsp_decimate_len(n, m);
    for (size_t i = 0; i < out; ++i)
        y[i] = filtered[i * m];

    free(filtered);
    return out;
}

size_t dsp_interpolate(const double *x, size_t n, size_t l, double *y) {
    if (l == 0 || n == 0) return 0;
    if (l == 1) {
        for (size_t i = 0; i < n; ++i) y[i] = x[i];
        return n;
    }

    size_t out = dsp_interpolate_len(n, l);

    /* Step 1: zero-stuff - insert l-1 zeros between samples. */
    double *stuffed = calloc(out, sizeof(double));
    if (!stuffed) return 0;
    for (size_t i = 0; i < n; ++i)
        stuffed[i * l] = x[i];

    /* Step 2: anti-imaging low-pass, cutoff 0.5/l, gain l. */
    double taps[RESAMPLE_TAPS];
    dsp_fir_design_lowpass(taps, RESAMPLE_TAPS, 0.5 / (double)l);
    for (size_t k = 0; k < RESAMPLE_TAPS; ++k)
        taps[k] *= (double)l;           /* restore amplitude */

    dsp_fir_apply(stuffed, y, out, taps, RESAMPLE_TAPS);

    free(stuffed);
    return out;
}

size_t dsp_resample(const double *x, size_t n,
                    size_t l, size_t m, double *y) {
    if (l == 0 || m == 0 || n == 0) return 0;

    /* Upsample by l into a temporary buffer ... */
    size_t up_len = dsp_interpolate_len(n, l);
    double *up = malloc(up_len * sizeof(double));
    if (!up) return 0;

    if (dsp_interpolate(x, n, l, up) == 0) { free(up); return 0; }

    /* ... then decimate that by m straight into y. */
    size_t out = dsp_decimate(up, up_len, m, y);

    free(up);
    return out;
}
