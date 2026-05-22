/*
 * timefreq.h - Time-frequency analysis
 *
 * PROBLEM SOLVED
 *   The FFT, and the parametric/subspace estimators in estimation.h,
 *   each produce ONE spectrum for the whole signal. That is fine for a
 *   stationary signal, but speech, music, radar chirps, and vibration
 *   signals are NON-STATIONARY - their frequency content changes over
 *   time. A single spectrum says which frequencies are present but not
 *   WHEN. Time-frequency analysis recovers that missing time axis.
 *
 * THREE APPROACHES
 *   STFT  - slide a short window along the signal and FFT each piece.
 *           Simple, invertible, the basis of the spectrogram. Its one
 *           limitation is a FIXED resolution: a short window localises
 *           time well but frequency poorly, and vice versa (the
 *           time-frequency uncertainty principle).
 *
 *   Filter bank - split the signal into frequency subbands with a set
 *           of bandpass filters, process each band, then reconstruct.
 *           A two-channel quadrature-mirror-filter (QMF) bank with the
 *           perfect-reconstruction property is the teaching example
 *           here; it is the same idea as the wavelet transform's
 *           low/high split, and underlies audio codecs (MP3, AAC).
 *
 *   Wigner-Ville - a quadratic distribution: the Fourier transform of
 *           the signal's instantaneous autocorrelation. It gives the
 *           sharpest possible time-frequency resolution with no window
 *           at all, but multiple components create spurious CROSS
 *           TERMS - oscillating artefacts halfway between real
 *           components - which is its defining drawback.
 *
 * The discrete wavelet transform (wavelet/) is a fourth time-frequency
 * tool, already in the guide, and is itself a wavelet filter bank.
 */
#ifndef DSP_TIMEFREQ_H
#define DSP_TIMEFREQ_H

#include "../common.h"
#include "spectral/window.h"

/* ===================================================================
 * STFT - Short-Time Fourier Transform
 * =================================================================== */

/*
 * The result of an STFT: a time-frequency grid of complex values.
 *   frames   : number of time frames
 *   bins     : FFT size (frequency bins per frame)
 *   data     : frames * bins complex values, row-major; data[t*bins+f]
 *              is the content of frequency bin f in time frame t
 *   hop      : hop size used between frames, in samples
 *   win_len  : analysis window length used
 */
typedef struct {
    size_t  frames;
    size_t  bins;
    size_t  hop;
    size_t  win_len;
    cplx   *data;
} dsp_stft;

/*
 * Compute the STFT of a real signal.
 *   x        : input samples, length n
 *   n        : sample count
 *   win_type : analysis window applied to each frame
 *   win_len  : window / FFT length; MUST be a power of two
 *   hop      : samples advanced between consecutive frames (hop <
 *              win_len gives overlap; win_len/4 is a common choice)
 *   out      : receives the STFT; free it with dsp_stft_free
 * Returns 0 on success, -1 on bad parameters or allocation failure.
 */
int dsp_stft_forward(const double *x, size_t n,
                     dsp_window_type win_type, size_t win_len,
                     size_t hop, dsp_stft *out);

/* Release the buffer held by an STFT result. */
void dsp_stft_free(dsp_stft *s);

/*
 * Invert an STFT back to a time signal by the overlap-add method.
 *   s   : an STFT produced by dsp_stft_forward
 *   x   : output samples, length at least dsp_stft_signal_len(s)
 * Returns the number of samples written. Reconstruction is exact (up
 * to edge effects) when the window and hop give constant overlap-add.
 */
size_t dsp_stft_inverse(const dsp_stft *s, double *x);

/* Length of the signal an STFT inverts to. */
size_t dsp_stft_signal_len(const dsp_stft *s);

/*
 * Spectrogram: the squared magnitude of the STFT, |STFT|^2.
 *   s    : an STFT result
 *   spec : output power grid, frames * bins, row-major
 * This is the standard visual time-frequency representation.
 */
void dsp_spectrogram(const dsp_stft *s, double *spec);

/* ===================================================================
 * Two-channel QMF filter bank
 * =================================================================== */

/*
 * A two-channel quadrature-mirror-filter analysis/synthesis bank.
 * The analysis stage splits a signal into a low-pass and a high-pass
 * subband, each decimated by two (so the total sample count is
 * unchanged - critical sampling). The synthesis stage reconstructs.
 * With the QMF relationship between the filters the aliasing from the
 * two bands cancels, giving (near) perfect reconstruction.
 *
 *   taps   : low-pass prototype filter coefficients
 *   ntaps  : prototype length
 * The high-pass filter and the synthesis filters are derived from the
 * prototype internally by the QMF mirror relations.
 */
typedef struct {
    double *taps;
    size_t  ntaps;
} dsp_qmf_bank;

/*
 * Initialise a QMF bank with a built-in low-pass prototype.
 * Returns 0 on success, -1 on allocation failure. Pair with
 * dsp_qmf_free.
 */
int dsp_qmf_init(dsp_qmf_bank *bank);

/* Release a QMF bank's filter buffer. */
void dsp_qmf_free(dsp_qmf_bank *bank);

/*
 * Analysis: split a signal into low and high subbands.
 *   x    : input samples, length n
 *   n    : sample count (even)
 *   lo   : low-pass subband, length n/2
 *   hi   : high-pass subband, length n/2
 * Each subband is filtered then decimated by two.
 */
void dsp_qmf_analyze(const dsp_qmf_bank *bank, const double *x,
                     size_t n, double *lo, double *hi);

/*
 * Synthesis: reconstruct a signal from its two subbands.
 *   lo, hi : subbands, length half each
 *   half   : subband length
 *   out    : reconstructed signal, length 2*half
 * Each subband is upsampled, filtered by the synthesis filters, and
 * summed; QMF design makes the aliasing terms cancel.
 */
void dsp_qmf_synthesize(const dsp_qmf_bank *bank,
                        const double *lo, const double *hi,
                        size_t half, double *out);

/* ===================================================================
 * Wigner-Ville distribution
 * =================================================================== */

/*
 * Discrete Wigner-Ville distribution of a real signal.
 *
 *   x      : input samples, length n
 *   n      : sample count; MUST be a power of two
 *   wvd    : output time-frequency grid, n * n doubles, row-major;
 *            wvd[t*n + f] is the energy density at time t, bin f
 *
 * The WVD is the Fourier transform, over a lag variable, of the
 * instantaneous autocorrelation x(t+tau/2) x*(t-tau/2). It achieves
 * the sharpest time-frequency localisation of any method here, but
 * for a multi-component signal it shows CROSS TERMS: oscillating
 * artefacts midway between every pair of true components. The signal
 * is internally converted to its analytic form to halve the cross
 * terms and remove the frequency aliasing of a real input.
 *
 * Returns 0 on success, -1 if n is not a power of two.
 */
int dsp_wigner_ville(const double *x, size_t n, double *wvd);

/*
 * Pseudo Wigner-Ville distribution: the WVD with a smoothing window
 * applied along the lag axis. The window suppresses the oscillating
 * cross terms at the cost of some resolution - the usual practical
 * compromise.
 *   x        : input samples, length n (power of two)
 *   win_type : smoothing window applied across lag
 *   wvd      : output grid, n * n doubles, row-major
 * Returns 0 on success, -1 if n is not a power of two.
 */
int dsp_pseudo_wigner_ville(const double *x, size_t n,
                            dsp_window_type win_type, double *wvd);

#endif /* DSP_TIMEFREQ_H */
