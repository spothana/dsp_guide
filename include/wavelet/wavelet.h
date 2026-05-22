/*
 * wavelet.h - Discrete Wavelet Transform (multi-resolution analysis)
 *
 * PROBLEM SOLVED
 *   The Fourier transform tells you WHAT frequencies are present, but
 *   not WHEN. The wavelet transform decomposes a signal using scaled
 *   and shifted copies of a localised "mother wavelet", giving
 *   ADAPTIVE time-frequency resolution:
 *     - high frequencies: narrow wavelet -> fine time resolution
 *     - low  frequencies: wide   wavelet -> fine frequency resolution
 *
 * WHY IT MATTERS
 *   Ideal for NON-STATIONARY signals whose spectral content changes
 *   over time - speech, music transients, ECG, seismic data. Signal
 *   energy concentrates in a few large coefficients while noise
 *   spreads thin, which makes wavelets excellent for compression and
 *   denoising (JPEG 2000 uses wavelets instead of the DCT).
 *
 * THIS IMPLEMENTATION
 *   The fast DWT via Mallat's pyramid algorithm with Haar filters -
 *   the simplest orthonormal wavelet. Each level splits the signal
 *   into an "approximation" (low-pass, averages) half and a "detail"
 *   (high-pass, differences) half, then recurses on the approximation.
 *
 * CONSTRAINT
 *   Signal length must be a power of two.
 */
#ifndef DSP_WAVELET_H
#define DSP_WAVELET_H

#include "../common.h"

/*
 * Single-level Haar DWT.
 *   in    : signal, length n (n even)
 *   approx: output low-pass coefficients,  length n/2
 *   detail: output high-pass coefficients, length n/2
 */
void dsp_dwt_haar_step(const double *in, size_t n,
                       double *approx, double *detail);

/*
 * Single-level inverse Haar DWT - reconstructs `out` (length n) from
 * approximation and detail halves (each length n/2).
 */
void dsp_idwt_haar_step(const double *approx, const double *detail,
                        size_t half, double *out);

/*
 * Full multi-level Haar DWT (Mallat pyramid), computed in place.
 *   data : signal of length n (MUST be a power of two); on return it
 *          holds, front to back, the coarsest approximation followed
 *          by detail coefficients from coarsest to finest level.
 * Returns the number of decomposition levels, or -1 if n is not a
 * power of two.
 */
int dsp_dwt_haar(double *data, size_t n);

/*
 * Full inverse multi-level Haar DWT - inverts dsp_dwt_haar in place.
 * Returns 0 on success, -1 if n is not a power of two.
 */
int dsp_idwt_haar(double *data, size_t n);

#endif /* DSP_WAVELET_H */
