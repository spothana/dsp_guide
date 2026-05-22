/*
 * correlation.h - Correlation
 *
 * PROBLEM SOLVED
 *   Measures the similarity between two signals as a function of a
 *   time lag tau:
 *
 *   R_xy[tau] = sum_n x[n] * y[n + tau]
 *
 *   It is mathematically close to convolution but does NOT time-reverse
 *   either signal - convolution flips one signal, correlation does not.
 *
 * USES
 *   - CROSS-correlation finds the time delay between related signals.
 *     Radar, sonar, and GPS correlate a transmitted pattern with the
 *     received echo; the lag at the peak gives the round-trip delay.
 *   - AUTO-correlation (a signal with itself) reveals periodicity and
 *     underpins pitch detection and power-spectrum estimation.
 *
 * FREQUENCY-DOMAIN VIEW
 *   Cross-correlation corresponds to X . conj(Y), i.e. multiply by the
 *   complex conjugate (convolution is just X . H).
 */
#ifndef DSP_CORRELATION_H
#define DSP_CORRELATION_H

#include "../common.h"

/*
 * Full cross-correlation of x (length n) with y (length m).
 *
 * Output `r` has length n + m - 1. Lag tau runs from -(m-1) to (n-1);
 * r[0] is the most negative lag and r[n+m-2] the most positive. The
 * helper dsp_xcorr_zero_index() gives the index of lag tau = 0.
 */
void dsp_cross_correlate(const double *x, size_t n,
                         const double *y, size_t m,
                         double *r);

/* Index within the cross-correlation output that corresponds to lag 0. */
static inline size_t dsp_xcorr_zero_index(size_t m) {
    return m - 1;
}

/*
 * Auto-correlation of x with itself. Output `r` has length n; r[k] is
 * the correlation at non-negative lag k (r[0] is the signal energy).
 */
void dsp_auto_correlate(const double *x, size_t n, double *r);

/*
 * Find the lag (which may be negative) at which the cross-correlation
 * of x and y is maximised - i.e. the estimated delay of y relative
 * to x. This is the core operation behind radar/sonar ranging.
 */
long dsp_estimate_delay(const double *x, size_t n,
                        const double *y, size_t m);

#endif /* DSP_CORRELATION_H */
