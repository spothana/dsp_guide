/*
 * dft.h - Discrete Fourier Transform
 *
 * PROBLEM SOLVED
 *   Reveals the frequency content hidden in N time-domain samples.
 *   The DFT is the mathematical foundation of all spectral analysis
 *   on discrete signals.
 *
 *   X[k] = sum_{n=0}^{N-1} x[n] * exp(-j*2*pi*k*n/N)
 *
 * COMPLEXITY
 *   O(N^2) - every one of the N outputs sums over all N inputs.
 *   This is why the FFT exists (see fft.h). Use the DFT here only as
 *   a correctness reference; never in production for large N.
 *
 * NOTES
 *   - Assumes the N samples are one period of an infinitely periodic
 *     signal. If the signal does not fit a whole number of cycles in
 *     the window, energy smears across bins (spectral leakage).
 *   - dsp_idft recovers the time-domain signal exactly (round-trip).
 */
#ifndef DSP_DFT_H
#define DSP_DFT_H

#include "../common.h"

/* Forward DFT: time-domain `in` (length n) -> frequency-domain `out`. */
void dsp_dft(const cplx *in, cplx *out, size_t n);

/* Inverse DFT: frequency-domain `in` -> time-domain `out`. */
void dsp_idft(const cplx *in, cplx *out, size_t n);

#endif /* DSP_DFT_H */
