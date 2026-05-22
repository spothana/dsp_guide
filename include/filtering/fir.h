/*
 * fir.h - Finite Impulse Response filtering
 *
 * PROBLEM SOLVED
 *   Frequency-selective filtering with guaranteed stability and the
 *   option of exactly linear phase.
 *
 *   y[n] = sum_{k=0}^{M} b[k] * x[n-k]            (no feedback)
 *
 * KEY PROPERTIES
 *   - ALWAYS STABLE. With no feedback path there are no poles away
 *     from the origin, so the filter cannot blow up for any choice
 *     of coefficients.
 *   - LINEAR PHASE when the coefficients are symmetric
 *     (b[k] == b[M-k]). Linear phase delays every frequency by the
 *     same amount, preserving waveform shape - critical for audio,
 *     communications, and biomedical signals.
 *   - COST: a sharp cutoff can need hundreds of taps, far more
 *     arithmetic than an IIR filter of equivalent selectivity.
 *
 * This module provides the direct-form FIR filter plus a windowed-sinc
 * design routine for low-pass filters.
 */
#ifndef DSP_FIR_H
#define DSP_FIR_H

#include "../common.h"

/*
 * Apply an FIR filter (direct form).
 *   x       : input signal, length n
 *   y       : output signal, length n (may not alias x)
 *   n       : signal length
 *   taps    : filter coefficients b[0..ntaps-1]
 *   ntaps   : number of coefficients
 * Output uses zero-state initial conditions (x[<0] treated as 0).
 */
void dsp_fir_apply(const double *x, double *y, size_t n,
                   const double *taps, size_t ntaps);

/*
 * Design a low-pass FIR filter by the windowed-sinc method.
 *   taps   : output buffer, length ntaps
 *   ntaps  : number of taps (odd recommended for a true linear-phase
 *            Type-I filter)
 *   fc     : normalised cutoff in (0, 0.5), i.e. fraction of the
 *            sample rate
 * The ideal sinc impulse response is multiplied by a Hamming window
 * to tame the truncation ripple. Coefficients are symmetric, so the
 * resulting filter has exactly linear phase.
 */
void dsp_fir_design_lowpass(double *taps, size_t ntaps, double fc);

#endif /* DSP_FIR_H */
