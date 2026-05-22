/*
 * convolution.h - Convolution
 *
 * PROBLEM SOLVED
 *   Convolution is how filtering works mathematically. For an LTI
 *   system with impulse response h, the output is x convolved with h:
 *
 *   y[n] = sum_k x[k] * h[n-k]
 *
 * KEY IDENTITY
 *   Convolution in time == multiplication in frequency:
 *       Y = X . H
 *   So a signal can be filtered by transforming to the frequency
 *   domain, multiplying, and transforming back.
 *
 * METHODS PROVIDED
 *   - dsp_convolve      : direct O(N*M) sum. Clear, good for short h.
 *   - dsp_convolve_fft  : fast convolution via the FFT, O(N log N).
 *                         Pads both inputs to a common power of two,
 *                         multiplies the spectra, and inverse-transforms.
 *
 * Output length for inputs of length n and m is always n + m - 1.
 */
#ifndef DSP_CONVOLUTION_H
#define DSP_CONVOLUTION_H

#include "../common.h"

/* Output length of a convolution of inputs sized n and m. */
static inline size_t dsp_conv_len(size_t n, size_t m) {
    return (n == 0 || m == 0) ? 0 : n + m - 1;
}

/*
 * Direct (time-domain) convolution.
 *   x,n : first signal and its length
 *   h,m : second signal and its length
 *   y   : output, length dsp_conv_len(n, m), provided by caller
 */
void dsp_convolve(const double *x, size_t n,
                  const double *h, size_t m,
                  double *y);

/*
 * FFT-based fast convolution. Same arguments and output length as
 * dsp_convolve, but O((n+m) log(n+m)) instead of O(n*m).
 * Returns 0 on success, -1 on allocation failure.
 */
int dsp_convolve_fft(const double *x, size_t n,
                     const double *h, size_t m,
                     double *y);

#endif /* DSP_CONVOLUTION_H */
