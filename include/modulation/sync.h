/*
 * sync.h - Carrier and timing recovery
 *
 * PROBLEM SOLVED
 *   A receiver's local oscillator and sample clock never exactly
 *   match the transmitter's. Two impairments result, and a real
 *   modem must track and cancel both before it can demap symbols:
 *
 *   CARRIER OFFSET - a residual frequency/phase error makes the whole
 *     constellation spin. Left uncorrected, every symbol decision is
 *     wrong within a few samples. CARRIER RECOVERY is a feedback loop
 *     (a phase-locked loop, PLL) that estimates the rotation and
 *     counter-rotates the signal to lock the constellation in place.
 *
 *   TIMING OFFSET - the receiver does not know exactly when each
 *     symbol is centred, so naive sampling catches symbols off-peak,
 *     where pulse tails from neighbours still bleed in. TIMING
 *     RECOVERY finds the best sampling instant within each symbol
 *     period.
 *
 * THESE ARE FEEDBACK LOOPS
 *   Both are control loops, conceptually the same family as the LMS
 *   equalizer: measure an error, nudge a parameter to shrink it,
 *   repeat. A small loop gain converges slowly but tracks smoothly; a
 *   large gain locks fast but is jittery - the same speed/stability
 *   trade-off seen throughout the guide.
 *
 * SCOPE
 *   Teaching-oriented building blocks: a decision-directed carrier
 *   PLL for a QAM/PSK constellation, and a Gardner timing-error
 *   detector with an interpolating sampler. These show the structure
 *   of synchronisation without the full complexity of a production
 *   receiver.
 */
#ifndef DSP_SYNC_H
#define DSP_SYNC_H

#include "../common.h"

/* ===================================================================
 * Carrier recovery - decision-directed phase-locked loop
 * =================================================================== */

/*
 * Carrier-recovery PLL state.
 *   phase     : current estimated phase offset (radians)
 *   freq      : current estimated frequency offset (radians/sample)
 *   alpha     : proportional loop gain (phase correction)
 *   beta      : integral loop gain (frequency tracking)
 * alpha and beta set the loop bandwidth; beta << alpha for stability.
 */
typedef struct {
    double phase;
    double freq;
    double alpha;
    double beta;
} dsp_carrier_pll;

/*
 * Initialise a carrier-recovery PLL.
 *   alpha : proportional gain (try ~0.05)
 *   beta  : integral gain (try ~alpha/20 for a stable loop)
 */
void dsp_carrier_pll_init(dsp_carrier_pll *pll,
                          double alpha, double beta);

/*
 * Process one received sample.
 *   in  : the rotated received symbol
 * Returns the de-rotated symbol. The loop's phase-error detector is
 * decision-directed: it slices the corrected sample to the nearest
 * QPSK quadrant and uses the angle between the two as the error,
 * then advances phase and frequency estimates toward cancelling it.
 */
cplx dsp_carrier_pll_process(dsp_carrier_pll *pll, cplx in);

/*
 * Convenience: run the PLL over a whole block.
 *   in, out : input and output symbol arrays, length n
 */
void dsp_carrier_recover(dsp_carrier_pll *pll,
                         const cplx *in, cplx *out, size_t n);

/* ===================================================================
 * Timing recovery - Gardner timing-error detector
 * =================================================================== */

/*
 * Estimate the fractional symbol-timing offset of an oversampled
 * signal using the Gardner detector.
 *
 *   samples : received samples, length n, oversampled at `sps`
 *             samples per symbol (sps must be even - Gardner needs
 *             the mid-symbol sample)
 *   n       : sample count
 *   sps     : samples per symbol
 *
 * The Gardner error for symbol k uses three samples - the previous
 * symbol peak, the mid-point, and the current peak:
 *   e = real( mid * conj(curr - prev) )
 * averaged over the block. The returned value is the estimated
 * fractional offset in samples, in roughly (-sps/2, +sps/2); feed it
 * to dsp_timing_resample to correct the timing.
 */
double dsp_timing_error_gardner(const cplx *samples, size_t n,
                                size_t sps);

/*
 * Resample a signal at a shifted timing phase using linear
 * interpolation, then decimate to one sample per symbol.
 *   samples : oversampled input, length n
 *   n       : input sample count
 *   sps     : samples per symbol
 *   offset  : fractional sample offset to apply (from the Gardner
 *             detector); the sampler picks instants k*sps + offset
 *   syms    : output symbol-rate samples, length n / sps
 * Returns the number of output symbols.
 */
size_t dsp_timing_resample(const cplx *samples, size_t n, size_t sps,
                           double offset, cplx *syms);

#endif /* DSP_SYNC_H */
