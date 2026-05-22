/*
 * fft.h - Fast Fourier Transform (radix-2 Cooley-Tukey)
 *
 * PROBLEM SOLVED
 *   Computes the EXACT same result as the DFT, far more efficiently.
 *   This is what makes real-time spectral analysis possible.
 *
 * COMPLEXITY
 *   O(N log N) instead of the DFT's O(N^2). For N = 1024 that is
 *   roughly a 100x speedup.
 *
 * HOW
 *   Divide and conquer: an N-point DFT splits into two N/2-point DFTs
 *   (even-indexed and odd-indexed samples), recombined with twiddle
 *   factors. The recursion bottoms out at size-1 DFTs.
 *
 *   The butterfly is the atomic step: because W_N^(k+N/2) = -W_N^k,
 *   two outputs are produced from two inputs with a single twiddle
 *   multiply (one uses +, one uses -).
 *
 * CONSTRAINT
 *   This radix-2 implementation requires n to be a power of two.
 *   Callers with other lengths should zero-pad to dsp_next_pow2(n).
 */
#ifndef DSP_FFT_H
#define DSP_FFT_H

#include "../common.h"

/*
 * In-place iterative radix-2 FFT.
 *   data : array of n complex samples, overwritten with the spectrum.
 *   n    : MUST be a power of two.
 * Returns 0 on success, -1 if n is not a power of two.
 */
int dsp_fft(cplx *data, size_t n);

/*
 * In-place inverse FFT. Same constraint on n. Output is scaled by 1/n
 * so that ifft(fft(x)) == x.
 */
int dsp_ifft(cplx *data, size_t n);

#endif /* DSP_FFT_H */
