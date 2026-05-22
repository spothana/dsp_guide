/*
 * iir.h - Infinite Impulse Response filtering
 *
 * PROBLEM SOLVED
 *   Frequency-selective filtering that is very cheap per output sample
 *   because it reuses past OUTPUTS via feedback.
 *
 *   y[n] = sum_{k=0}^{M} b[k]*x[n-k] - sum_{k=1}^{L} a[k]*y[n-k]
 *
 * KEY PROPERTIES
 *   - EFFICIENT. A 4th-order IIR can match a magnitude response that
 *     would take 50+ FIR taps.
 *   - NONLINEAR PHASE. Feedback distorts waveform shape. Offline code
 *     can cancel this with forward-backward (zero-phase) filtering.
 *   - STABILITY NOT GUARANTEED. Every pole of the transfer function
 *     must lie strictly inside the unit circle. Coefficient
 *     quantisation can push a pole out and make the filter diverge.
 *
 * This module provides a direct-form-I biquad (2nd-order section), the
 * standard building block for IIR filters, plus an RBJ-cookbook
 * low-pass biquad designer and a simple stability check.
 */
#ifndef DSP_IIR_H
#define DSP_IIR_H

#include "../common.h"

/*
 * A biquad: second-order section with feedforward b0,b1,b2 and
 * feedback a1,a2 (a0 normalised to 1).
 *   H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2)
 */
typedef struct {
    double b0, b1, b2;   /* feedforward coefficients */
    double a1, a2;       /* feedback coefficients (a0 == 1) */
} dsp_biquad;

/*
 * Design a low-pass biquad (RBJ audio-EQ cookbook formulas).
 *   bq : output coefficients
 *   fc : normalised cutoff in (0, 0.5)
 *   q  : quality factor (0.707 == Butterworth, maximally flat)
 */
void dsp_iir_design_lowpass(dsp_biquad *bq, double fc, double q);

/*
 * Apply a biquad to a whole signal (direct form I, zero initial state).
 *   x, y : input and output, length n (y may not alias x)
 */
void dsp_iir_apply(const dsp_biquad *bq,
                   const double *x, double *y, size_t n);

/*
 * Stability check: returns 1 if both poles lie strictly inside the
 * unit circle, 0 otherwise. A stable biquad needs |a2| < 1 and
 * |a1| < 1 + a2.
 */
int dsp_iir_is_stable(const dsp_biquad *bq);

#endif /* DSP_IIR_H */
