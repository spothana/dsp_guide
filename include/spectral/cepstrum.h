/*
 * cepstrum.h - Cepstral analysis: the cepstrum and MFCCs
 *
 * THE IDEA
 *   Ordinary spectral analysis transforms a signal once, into the
 *   frequency domain. Cepstral analysis goes a step further: it
 *   treats the LOG SPECTRUM itself as a signal and transforms THAT.
 *   The playful name says it all - "cepstrum" is "spectrum" with the
 *   first syllable reversed, and its independent variable is called
 *   QUEFRENCY (an anagram of "frequency"), measured in samples.
 *
 * WHY IT IS USEFUL - SOURCE/FILTER SEPARATION
 *   Many signals are a source convolved with a filter: voiced speech
 *   is a buzzy glottal source shaped by the vocal-tract filter; an
 *   echo is a signal convolved with a delayed impulse. Convolution in
 *   time is MULTIPLICATION of spectra - and taking the LOG turns that
 *   product into a SUM. A final transform then sends the slowly
 *   varying part (the filter / spectral envelope) to LOW quefrency
 *   and the fast-varying part (the source / pitch periodicity) to a
 *   sharp peak at HIGH quefrency. Two things that were tangled
 *   together in the signal sit at different quefrencies in the
 *   cepstrum, where they can be read off or separated.
 *
 *   This is the basis of pitch detection (find the high-quefrency
 *   peak), echo detection, and homomorphic deconvolution.
 *
 * REAL vs COMPLEX CEPSTRUM
 *   The REAL cepstrum uses only the log MAGNITUDE spectrum; it
 *   discards phase, so it is not invertible, but it is all that pitch
 *   detection and MFCCs need. The COMPLEX cepstrum keeps the phase
 *   (via the unwrapped phase spectrum), so the original signal can be
 *   reconstructed - which is what makes homomorphic deconvolution
 *   possible.
 *
 * MFCC - MEL-FREQUENCY CEPSTRAL COEFFICIENTS
 *   The dominant feature for speech and audio recognition. It is a
 *   cepstrum with two changes that mimic human hearing:
 *     1. the power spectrum is grouped into bands by a MEL-SCALE
 *        filterbank - fine resolution at low frequency, coarse at
 *        high, as the ear does;
 *     2. the final transform is a DCT (not an inverse FFT), which
 *        also decorrelates the coefficients.
 *   The first dozen-or-so MFCCs capture the spectral envelope - the
 *   shape that distinguishes one sound from another - compactly.
 *
 * CONSTRAINT
 *   FFT-based routines here require a power-of-two length.
 */
#ifndef DSP_CEPSTRUM_H
#define DSP_CEPSTRUM_H

#include "../common.h"

/* ===================================================================
 * Cepstrum
 * =================================================================== */

/*
 * Real cepstrum of a signal: IFFT( log |FFT(x)| ).
 *   x   : input samples, length n
 *   n   : length; MUST be a power of two
 *   cep : output real cepstrum, length n, indexed by quefrency
 * Uses only the log-magnitude spectrum (phase discarded), so it is
 * not invertible. Returns 0 on success, -1 on bad parameters.
 */
int dsp_cepstrum_real(const double *x, size_t n, double *cep);

/*
 * Complex cepstrum of a signal: IFFT( log(FFT(x)) ), where the
 * complex log combines log-magnitude with the UNWRAPPED phase.
 *   x   : input samples, length n
 *   n   : length; MUST be a power of two
 *   cep : output complex cepstrum, length n
 * Keeps phase information, so the signal can be reconstructed with
 * dsp_icepstrum_complex - the basis of homomorphic deconvolution.
 * Returns 0 on success, -1 on bad parameters.
 */
int dsp_cepstrum_complex(const double *x, size_t n, cplx *cep);

/*
 * Inverse complex cepstrum: reconstruct a signal from its complex
 * cepstrum (the inverse of dsp_cepstrum_complex).
 *   cep : complex cepstrum, length n
 *   n   : length; MUST be a power of two
 *   x   : output reconstructed signal, length n
 * Returns 0 on success, -1 on bad parameters.
 */
int dsp_icepstrum_complex(const cplx *cep, size_t n, double *x);

/*
 * Estimate the pitch period of a (assumed voiced) signal from its
 * real cepstrum. Voiced speech and other periodic signals produce a
 * sharp cepstral peak at a quefrency equal to the pitch period; this
 * finds the highest peak within a plausible quefrency range.
 *   x        : input samples, length n (power of two)
 *   n        : length
 *   min_q    : lowest quefrency (in samples) to search
 *   max_q    : highest quefrency to search
 * Returns the quefrency of the strongest peak - the estimated pitch
 * period in samples - or 0 if no clear peak is found. The pitch in
 * Hz is sample_rate / period.
 */
size_t dsp_cepstrum_pitch(const double *x, size_t n,
                          size_t min_q, size_t max_q);

/* ===================================================================
 * MFCC - Mel-Frequency Cepstral Coefficients
 * =================================================================== */

/*
 * A mel-scale triangular filterbank: the front end of MFCC.
 *   nfilters   : number of mel bands
 *   fft_size   : FFT size the filters will be applied to
 *   sample_rate: sampling rate, in Hz (sets the frequency mapping)
 *   weights    : output filter weights, nfilters * (fft_size/2 + 1),
 *                row-major; row f is the triangular weighting of
 *                filter f across the FFT magnitude bins
 * The filters are spaced equally on the mel scale - which is roughly
 * linear below 1 kHz and logarithmic above - so they are narrow at
 * low frequency and wide at high, matching the ear's resolution.
 */
typedef struct {
    size_t  nfilters;
    size_t  fft_size;
    double  sample_rate;
    double *weights;          /* nfilters x (fft_size/2 + 1) */
} dsp_mel_filterbank;

/*
 * Build a mel filterbank. Returns 0 on success, -1 on bad parameters
 * or allocation failure. Pair with dsp_mel_filterbank_free.
 */
int dsp_mel_filterbank_init(dsp_mel_filterbank *fb, size_t nfilters,
                            size_t fft_size, double sample_rate);

/* Release a mel filterbank's weight buffer. */
void dsp_mel_filterbank_free(dsp_mel_filterbank *fb);

/*
 * Compute the MFCCs of one signal frame.
 *   frame      : input samples, length fft_size (a power of two);
 *                normally one short, windowed frame of audio
 *   fft_size   : frame / FFT length
 *   fb         : a mel filterbank built for this fft_size
 *   ncoeffs    : number of cepstral coefficients to return
 *                (typically 12-13; must not exceed fb->nfilters)
 *   mfcc       : output coefficients, length ncoeffs
 *
 * The pipeline: power spectrum -> mel filterbank -> log -> DCT, with
 * the first `ncoeffs` DCT outputs kept. For a whole signal, the
 * caller frames it (see the STFT in timefreq.h) and calls this once
 * per frame.
 *
 * Returns 0 on success, -1 on bad parameters.
 */
int dsp_mfcc(const double *frame, size_t fft_size,
             const dsp_mel_filterbank *fb,
             size_t ncoeffs, double *mfcc);

#endif /* DSP_CEPSTRUM_H */
