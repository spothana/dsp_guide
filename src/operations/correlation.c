/*
 * correlation.c - Cross-correlation, auto-correlation, delay estimation.
 */
#include "operations/correlation.h"

void dsp_cross_correlate(const double *x, size_t n,
                         const double *y, size_t m,
                         double *r) {
    size_t out = n + m - 1;
    for (size_t i = 0; i < out; ++i)
        r[i] = 0.0;

    /* lag runs from -(m-1) .. (n-1); index = lag + (m-1). */
    for (long lag = -(long)(m - 1); lag <= (long)(n - 1); ++lag) {
        double acc = 0.0;
        for (size_t i = 0; i < n; ++i) {
            long j = (long)i - lag;          /* index into y */
            if (j >= 0 && j < (long)m)
                acc += x[i] * y[j];
        }
        r[lag + (long)(m - 1)] = acc;
    }
}

void dsp_auto_correlate(const double *x, size_t n, double *r) {
    for (size_t lag = 0; lag < n; ++lag) {
        double acc = 0.0;
        for (size_t i = 0; i + lag < n; ++i)
            acc += x[i] * x[i + lag];
        r[lag] = acc;
    }
}

long dsp_estimate_delay(const double *x, size_t n,
                        const double *y, size_t m) {
    double best_val = -1.0;
    long   best_lag = 0;
    int    seen     = 0;

    /* R_xy[tau] = sum_n x[n] * y[n+tau]. The peak tau is the shift
     * that best aligns y onto x; a positive value means y is delayed
     * relative to x (y's features appear tau samples later). */
    for (long tau = -(long)(n - 1); tau <= (long)(m - 1); ++tau) {
        double acc = 0.0;
        for (size_t i = 0; i < n; ++i) {
            long j = (long)i + tau;      /* index into y */
            if (j >= 0 && j < (long)m)
                acc += x[i] * y[j];
        }
        if (!seen || acc > best_val) {
            best_val = acc;
            best_lag = tau;
            seen = 1;
        }
    }
    return best_lag;
}
