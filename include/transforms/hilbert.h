/*
 * hilbert.h - Hilbert transform and the analytic signal
 *
 * PROBLEM SOLVED
 *   The Hilbert transform shifts every frequency component of a real
 *   signal by -90 degrees. On its own that is a curiosity; its value
 *   is what it builds - the ANALYTIC SIGNAL:
 *
 *       z(t) = x(t) + j * H{x(t)}
 *
 *   a complex signal whose real part is the original and whose
 *   imaginary part is its Hilbert transform. The analytic signal has
 *   only positive frequencies, and that makes two quantities fall
 *   straight out of it:
 *
 *     ENVELOPE              |z(t)|   - the instantaneous amplitude,
 *                                      the smooth outline of the
 *                                      signal (AM demodulation, onset
 *                                      detection, vibration analysis)
 *     INSTANTANEOUS PHASE   arg z(t) - and its time derivative, the
 *                                      INSTANTANEOUS FREQUENCY (FM
 *                                      demodulation, chirp analysis)
 *
 * HOW IT IS COMPUTED
 *   Frequency-domain method: take the FFT, then keep DC and Nyquist,
 *   double the positive-frequency bins, and zero the negative ones.
 *   The inverse FFT is the analytic signal. The Hilbert transform
 *   itself is its imaginary part. This is exact (no filter design),
 *   which is why it is the standard approach.
 *
 * RELATION TO THE REST OF THE GUIDE
 *   The Wigner-Ville distribution in spectral/timefreq.h uses the
 *   analytic signal to halve its cross terms; it now calls the public
 *   functions declared here.
 *
 * CONSTRAINT
 *   The transform length n must be a power of two (it uses the
 *   radix-2 FFT).
 */
#ifndef DSP_HILBERT_H
#define DSP_HILBERT_H

#include "../common.h"

/*
 * Analytic signal of a real input.
 *   x : real input samples, length n
 *   z : output analytic signal, length n; real part is x, imaginary
 *       part is the Hilbert transform of x
 * n must be a power of two. Returns 0 on success, -1 otherwise.
 */
int dsp_analytic_signal(const double *x, size_t n, cplx *z);

/*
 * Hilbert transform of a real signal (every frequency shifted by
 * -90 degrees).
 *   x : real input samples, length n
 *   h : output Hilbert transform, length n
 * n must be a power of two. Equivalent to the imaginary part of the
 * analytic signal. Returns 0 on success, -1 otherwise.
 */
int dsp_hilbert(const double *x, size_t n, double *h);

/*
 * Envelope (instantaneous amplitude) of a real signal: the magnitude
 * of its analytic signal, |z(t)|.
 *   x   : real input samples, length n
 *   env : output envelope, length n
 * n must be a power of two. Returns 0 on success, -1 otherwise.
 */
int dsp_envelope(const double *x, size_t n, double *env);

/*
 * Instantaneous frequency of a real signal, in cycles per sample.
 *   x    : real input samples, length n
 *   freq : output instantaneous frequency, length n
 * Computed as the sample-to-sample difference of the unwrapped
 * instantaneous phase, divided by 2*pi. freq[0] repeats freq[1] so
 * the output has the same length as the input. n must be a power of
 * two. Returns 0 on success, -1 otherwise.
 */
int dsp_instantaneous_frequency(const double *x, size_t n,
                                double *freq);

#endif /* DSP_HILBERT_H */
