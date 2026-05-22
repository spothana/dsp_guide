/*
 * estimation.h - Advanced spectral estimation
 *
 * PROBLEM SOLVED
 *   The FFT periodogram (the window functions in window.h prepare a
 *   signal for it) estimates a spectrum directly from the data. Its
 *   frequency resolution is fixed by the record length: two tones
 *   closer than about 1/N in normalised frequency blur into one peak,
 *   and a short or noisy record gives a ragged spectrum.
 *
 *   ADVANCED spectral estimators do better by assuming structure:
 *
 *   PARAMETRIC methods fit a model and read the spectrum off the
 *   model parameters:
 *     AR   - models the signal as an all-pole filter driven by white
 *            noise. Sharp spectral peaks, good for short records.
 *     ARMA - adds zeros (a moving-average part) for more general
 *            spectra that have both peaks and notches.
 *
 *   SUBSPACE methods split the data's covariance matrix into a signal
 *   and a noise subspace via an eigen-decomposition, then locate
 *   sinusoids with "super-resolution" - far finer than 1/N:
 *     MUSIC  - searches a pseudospectrum whose peaks are the tone
 *              frequencies.
 *     ESPRIT - solves for the frequencies directly, with no search.
 *
 * These shine in radar/sonar, array processing (direction finding),
 * tone detection, and speech analysis - wherever the FFT's resolution
 * is not enough.
 */
#ifndef DSP_ESTIMATION_H
#define DSP_ESTIMATION_H

#include <stddef.h>

/* ===================================================================
 * Autocorrelation - the common starting point
 * =================================================================== */

/*
 * Biased autocorrelation estimate of a real signal.
 *   x       : input samples, length n
 *   maxlag  : highest lag to compute
 *   r       : output, length maxlag+1; r[k] is the lag-k correlation
 * The biased form (dividing by n, not n-k) is used because it yields
 * a positive-semidefinite sequence, which the model fits require.
 */
void dsp_autocorr(const double *x, size_t n, size_t maxlag, double *r);

/* ===================================================================
 * AR - autoregressive spectral estimation
 * =================================================================== */

/*
 * Fit AR(p) coefficients by the Yule-Walker method, solved with the
 * Levinson-Durbin recursion.
 *   r      : autocorrelation, length order+1 (from dsp_autocorr)
 *   order  : AR model order p
 *   a      : output AR coefficients, length order; the model is
 *            x(n) = -sum a[k] x(n-1-k) + w(n)
 *   sigma2 : if non-NULL, receives the driving white-noise variance
 * Returns 0 on success, -1 on a degenerate (singular) correlation.
 */
int dsp_ar_yule_walker(const double *r, size_t order,
                       double *a, double *sigma2);

/*
 * Fit AR(p) coefficients by Burg's method (forward-backward linear
 * prediction, the maximum-entropy estimator). Often more accurate
 * than Yule-Walker for short records.
 *   x      : input samples, length n
 *   n      : sample count
 *   order  : AR model order p
 *   a      : output AR coefficients, length order
 *   sigma2 : if non-NULL, receives the driving-noise variance
 * Returns 0 on success, -1 on bad parameters.
 */
int dsp_ar_burg(const double *x, size_t n, size_t order,
                double *a, double *sigma2);

/*
 * Evaluate the AR power spectral density at `nf` frequencies evenly
 * spaced over [0, 0.5) (normalised frequency).
 *   a      : AR coefficients, length order
 *   order  : model order
 *   sigma2 : driving-noise variance
 *   psd    : output PSD, length nf
 */
void dsp_ar_psd(const double *a, size_t order, double sigma2,
                double *psd, size_t nf);

/* ===================================================================
 * ARMA - autoregressive moving-average estimation
 * =================================================================== */

/*
 * Fit an ARMA(p, q) model by the modified Yule-Walker method: the AR
 * part is estimated from the high-lag autocorrelations, then the MA
 * spectrum is obtained from the AR residual.
 *   x      : input samples, length n
 *   n      : sample count
 *   p      : AR order
 *   q      : MA order
 *   a      : output AR coefficients, length p
 *   b      : output MA coefficients, length q+1 (b[0] is the gain)
 * Returns 0 on success, -1 on bad parameters or a singular fit.
 */
int dsp_arma_estimate(const double *x, size_t n, size_t p, size_t q,
                      double *a, double *b);

/*
 * Evaluate the ARMA power spectral density at `nf` frequencies over
 * [0, 0.5).
 *   a   : AR coefficients, length p
 *   b   : MA coefficients, length q+1
 *   psd : output PSD, length nf
 */
void dsp_arma_psd(const double *a, size_t p,
                  const double *b, size_t q,
                  double *psd, size_t nf);

/* ===================================================================
 * Subspace methods - MUSIC and ESPRIT
 * =================================================================== */

/*
 * MUSIC pseudospectrum.
 *   x        : input samples, length n
 *   n        : sample count
 *   nsources : number of sinusoids to resolve (the signal-subspace
 *              dimension)
 *   msize    : covariance-matrix size; must exceed 2*nsources, and a
 *              larger value sharpens the estimate
 *   pseudo   : output pseudospectrum, length nf, over [0, 0.5);
 *              its PEAKS mark the estimated tone frequencies
 *   nf       : number of frequency grid points
 * Returns 0 on success, -1 on bad parameters.
 */
int dsp_music(const double *x, size_t n, size_t nsources,
              size_t msize, double *pseudo, size_t nf);

/*
 * ESPRIT frequency estimation - returns the frequencies directly,
 * with no spectral search.
 *   x        : input samples, length n
 *   n        : sample count
 *   nsources : number of sinusoids to resolve
 *   msize    : covariance-matrix size (> 2*nsources)
 *   freqs    : output normalised frequencies in [0, 0.5), length
 *              nsources; not sorted in any particular order
 * Returns the number of frequencies found (<= nsources), or -1 on
 * bad parameters.
 */
int dsp_esprit(const double *x, size_t n, size_t nsources,
               size_t msize, double *freqs);

#endif /* DSP_ESTIMATION_H */
