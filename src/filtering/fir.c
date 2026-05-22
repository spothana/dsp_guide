/*
 * fir.c - Direct-form FIR filtering and windowed-sinc low-pass design.
 */
#include "filtering/fir.h"

void dsp_fir_apply(const double *x, double *y, size_t n,
                   const double *taps, size_t ntaps) {
    for (size_t i = 0; i < n; ++i) {
        double acc = 0.0;
        for (size_t k = 0; k < ntaps; ++k) {
            /* x[i-k], with samples before the start treated as zero. */
            if (i >= k)
                acc += taps[k] * x[i - k];
        }
        y[i] = acc;
    }
}

void dsp_fir_design_lowpass(double *taps, size_t ntaps, double fc) {
    /* Centre index - the filter is symmetric about this point. */
    double mid = (double)(ntaps - 1) / 2.0;
    double sum = 0.0;

    for (size_t k = 0; k < ntaps; ++k) {
        double m = (double)k - mid;

        /* Ideal low-pass impulse response: a sinc function. */
        double sinc;
        if (m == 0.0)
            sinc = 2.0 * fc;
        else
            sinc = sin(2.0 * M_PI * fc * m) / (M_PI * m);

        /* Hamming window to suppress truncation ripple. */
        double w = 0.54 - 0.46 * cos(2.0 * M_PI * (double)k
                                          / (double)(ntaps - 1));

        taps[k] = sinc * w;
        sum += taps[k];
    }

    /* Normalise to unity DC gain (passband gain == 1). */
    for (size_t k = 0; k < ntaps; ++k)
        taps[k] /= sum;
}
