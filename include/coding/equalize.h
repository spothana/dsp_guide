/*
 * equalize.h - Channel equalization (the DSP side of error handling)
 *
 * PROBLEM SOLVED
 *   Before any decoder can correct bit errors, the DSP front end must
 *   undo the distortion the physical channel imposed on the waveform:
 *   multipath echoes, bandwidth limiting, and similar linear effects
 *   that smear one symbol into the next (inter-symbol interference).
 *
 *   A channel acts like an unknown filter. An EQUALIZER is an adaptive
 *   filter that approximates the channel's inverse, so that channel
 *   followed by equalizer is roughly flat. Cleaning the signal here
 *   is what makes the downstream FEC (correct.h) effective.
 *
 * THIS MODULE
 *   An LMS (Least Mean Squares) adaptive FIR equalizer. The LMS rule
 *   nudges each tap a little, every sample, in the direction that
 *   reduces the squared error between the equalizer output and a known
 *   reference (a training sequence). It is the workhorse adaptive
 *   algorithm: cheap, stable, and simple. RLS converges faster at
 *   higher cost; the trade-off mirrors FIR-vs-IIR elsewhere in this
 *   guide.
 *
 * RELATED RECEIVER FUNCTIONS (context, not implemented here)
 *   Carrier recovery - tracks and removes residual frequency/phase
 *     offset so the constellation does not rotate.
 *   Timing recovery  - locks the sampling instant to the symbol
 *     centre so each symbol is sampled at its strongest point.
 *   Both are feedback loops conceptually similar to the LMS update.
 */
#ifndef DSP_EQUALIZE_H
#define DSP_EQUALIZE_H

#include <stddef.h>

/*
 * LMS adaptive FIR equalizer state.
 *   weights : the adaptive filter taps
 *   line    : delay line of recent inputs
 *   ntaps   : filter length
 *   mu      : step size (convergence speed vs stability trade-off)
 */
typedef struct {
    double *weights;
    double *line;
    size_t  ntaps;
    double  mu;
} dsp_lms;

/*
 * Initialise an LMS equalizer with `ntaps` taps and step size `mu`.
 * Allocates internal buffers. Returns 0 on success, -1 on allocation
 * failure. Call dsp_lms_free when done.
 */
int dsp_lms_init(dsp_lms *eq, size_t ntaps, double mu);

/* Release the buffers owned by an LMS equalizer. */
void dsp_lms_free(dsp_lms *eq);

/*
 * Process one input sample.
 *   eq      : equalizer state
 *   input   : the distorted sample arriving from the channel
 *   desired : the known reference (training) sample
 *   err_out : if non-NULL, receives the instantaneous error
 * Returns the equalizer output for this sample. The taps are updated
 * by the LMS rule:  w += mu * error * line.
 */
double dsp_lms_update(dsp_lms *eq, double input, double desired,
                      double *err_out);

/*
 * Run the equalizer over a training block, returning the final
 * mean-squared error. A decreasing MSE indicates convergence.
 */
double dsp_lms_train(dsp_lms *eq, const double *input,
                     const double *desired, size_t n);

#endif /* DSP_EQUALIZE_H */
