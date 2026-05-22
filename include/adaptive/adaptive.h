/*
 * adaptive.h - Adaptive filters
 *
 * PROBLEM SOLVED
 *   A fixed FIR/IIR filter has constant coefficients. An ADAPTIVE
 *   filter adjusts its coefficients automatically, sample by sample,
 *   to minimise the error between its output and a desired reference.
 *   This lets one filter cope with an unknown or changing system.
 *
 *   Every adaptive filter has the same skeleton:
 *     input    x(n)  -> filter -> output y(n)
 *     desired  d(n)
 *     error    e(n) = d(n) - y(n)
 *   and a learning rule that nudges the weights to shrink e(n).
 *
 * APPLICATIONS
 *   - Channel equalization  - invert an unknown channel's distortion
 *   - System identification - learn an unknown system's impulse response
 *   - Noise cancellation    - subtract noise correlated with a reference
 *   - Echo cancellation, prediction, beamforming
 *
 * THREE ALGORITHMS - a convergence-vs-cost spectrum
 *   LMS  : stochastic-gradient update. O(L) per sample, simple,
 *          robust, but slow and sensitive to the step size mu and to
 *          input signal power.
 *   NLMS : LMS with the step size normalised by the input power.
 *          Still O(L), but mu is easy to choose and convergence no
 *          longer degrades when the input level swings - the usual
 *          choice for echo and noise cancellation.
 *   RLS  : recursively solves the least-squares problem over all past
 *          samples. O(L^2) per sample - far costlier - but converges
 *          much faster and tracks rapidly changing signals better.
 *
 *   The same trade-off the guide meets in FIR-vs-IIR and
 *   direct-vs-FFT: cheap and steady, or expensive and fast.
 *
 * L is the filter length (number of taps) throughout.
 */
#ifndef DSP_ADAPTIVE_H
#define DSP_ADAPTIVE_H

#include <stddef.h>

/* ===================================================================
 * LMS - Least Mean Squares
 * =================================================================== */

/*
 * LMS adaptive FIR filter state.
 *   weights : the adaptive filter taps
 *   line    : delay line of recent inputs
 *   ntaps   : filter length L
 *   mu      : step size (learning rate) - the critical parameter
 */
typedef struct {
    double *weights;
    double *line;
    size_t  ntaps;
    double  mu;
} dsp_lms;

/*
 * Initialise an LMS filter with `ntaps` taps and step size `mu`.
 * Allocates internal buffers. Returns 0 on success, -1 on allocation
 * failure. Pair with dsp_lms_free.
 */
int dsp_lms_init(dsp_lms *f, size_t ntaps, double mu);

/* Release the buffers owned by an LMS filter. */
void dsp_lms_free(dsp_lms *f);

/*
 * Process one sample.
 *   input   : the sample arriving at the filter input
 *   desired : the reference sample to match
 *   err_out : if non-NULL, receives the instantaneous error
 * Returns the filter output. The weights are updated by the LMS rule
 *   w += mu * error * x.
 */
double dsp_lms_update(dsp_lms *f, double input, double desired,
                      double *err_out);

/*
 * Run the filter over a block, returning the final mean-squared
 * error. A decreasing MSE indicates convergence.
 */
double dsp_lms_train(dsp_lms *f, const double *input,
                     const double *desired, size_t n);

/* ===================================================================
 * NLMS - Normalized Least Mean Squares
 * =================================================================== */

/*
 * NLMS adaptive FIR filter state. Same structure as LMS plus eps, the
 * small constant that keeps the power normalisation from dividing by
 * zero on a silent input.
 */
typedef struct {
    double *weights;
    double *line;
    size_t  ntaps;
    double  mu;          /* normalised step size, typically 0.1 .. 1.0 */
    double  eps;         /* regularisation constant */
} dsp_nlms;

/*
 * Initialise an NLMS filter.
 *   ntaps : filter length
 *   mu    : normalised step size in (0, 2); ~0.5 is a safe default
 *   eps   : small positive regularisation constant (e.g. 1e-6)
 * Returns 0 on success, -1 on allocation failure.
 */
int dsp_nlms_init(dsp_nlms *f, size_t ntaps, double mu, double eps);

/* Release the buffers owned by an NLMS filter. */
void dsp_nlms_free(dsp_nlms *f);

/*
 * Process one sample. Identical interface to dsp_lms_update, but the
 * weight update divides the step by the input power:
 *   w += (mu / (||x||^2 + eps)) * error * x.
 * This makes convergence independent of the input signal level.
 */
double dsp_nlms_update(dsp_nlms *f, double input, double desired,
                       double *err_out);

/* Run the NLMS filter over a block; returns the final MSE. */
double dsp_nlms_train(dsp_nlms *f, const double *input,
                      const double *desired, size_t n);

/* ===================================================================
 * RLS - Recursive Least Squares
 * =================================================================== */

/*
 * RLS adaptive FIR filter state.
 *   weights : the adaptive filter taps
 *   line    : delay line of recent inputs
 *   P       : the L x L inverse-correlation matrix, updated each step
 *   gain    : scratch Kalman-gain vector, length L
 *   tmp     : scratch vector, length L
 *   ntaps   : filter length L
 *   lambda  : forgetting factor in (0, 1] - discounts old samples;
 *             1.0 weights all history equally, < 1 tracks change
 *   delta   : initial diagonal loading of P (sets early-step gain)
 */
typedef struct {
    double *weights;
    double *line;
    double *P;
    double *gain;
    double *tmp;
    size_t  ntaps;
    double  lambda;
    double  delta;
} dsp_rls;

/*
 * Initialise an RLS filter.
 *   ntaps  : filter length
 *   lambda : forgetting factor, typically 0.95 .. 0.9999
 *   delta  : initial value loaded onto the diagonal of P; a larger
 *            delta gives a gentler, more stable start (e.g. 100.0)
 * Returns 0 on success, -1 on allocation failure.
 */
int dsp_rls_init(dsp_rls *f, size_t ntaps, double lambda, double delta);

/* Release the buffers owned by an RLS filter. */
void dsp_rls_free(dsp_rls *f);

/*
 * Process one sample. Same interface as the LMS/NLMS updates, but the
 * weights are updated by the recursive least-squares rule, which
 * maintains an inverse-correlation matrix for fast convergence.
 * Cost is O(L^2) per sample.
 */
double dsp_rls_update(dsp_rls *f, double input, double desired,
                      double *err_out);

/* Run the RLS filter over a block; returns the final MSE. */
double dsp_rls_train(dsp_rls *f, const double *input,
                     const double *desired, size_t n);

#endif /* DSP_ADAPTIVE_H */
