/*
 * cepstrum.c - Real and complex cepstrum, and MFCC.
 *
 * Structure:
 *   - the real cepstrum   (IFFT of the log-magnitude spectrum)
 *   - the complex cepstrum (IFFT of the complex log; needs phase
 *     unwrapping) and its inverse
 *   - cepstral pitch detection
 *   - the mel-scale filterbank and the MFCC pipeline
 */
#include "spectral/cepstrum.h"
#include "transforms/fft.h"
#include "transforms/dct.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* A floor on the magnitude before taking a log, so a spectral null
 * gives a large-but-finite value instead of -infinity. */
#define CEP_LOG_FLOOR 1e-12

/* ===================================================================
 * Real cepstrum
 * =================================================================== */

int dsp_cepstrum_real(const double *x, size_t n, double *cep) {
    if (n == 0 || !dsp_is_pow2(n))
        return -1;

    cplx *spec = malloc(n * sizeof(cplx));
    if (!spec)
        return -1;

    /* Forward FFT. */
    for (size_t i = 0; i < n; ++i)
        spec[i] = x[i];
    if (dsp_fft(spec, n) != 0) {
        free(spec);
        return -1;
    }

    /* Replace each bin with the log of its magnitude (a real value),
     * then inverse-transform: the result is the real cepstrum. */
    for (size_t k = 0; k < n; ++k) {
        double mag = cabs(spec[k]);
        if (mag < CEP_LOG_FLOOR)
            mag = CEP_LOG_FLOOR;
        spec[k] = log(mag);
    }
    if (dsp_ifft(spec, n) != 0) {
        free(spec);
        return -1;
    }
    for (size_t i = 0; i < n; ++i)
        cep[i] = creal(spec[i]);

    free(spec);
    return 0;
}

/* ===================================================================
 * Complex cepstrum
 *
 * The complex log of a spectrum is  log|X| + j*arg(X). The magnitude
 * part is easy; the phase must be UNWRAPPED - the bin-to-bin jumps of
 * +/-2*pi that arg() introduces have to be removed so the phase is a
 * smooth function of frequency. Without unwrapping the inverse
 * transform is meaningless.
 * =================================================================== */

/* Unwrap a phase array in place: add multiples of 2*pi so that no
 * step between adjacent samples exceeds pi in magnitude. */
static void unwrap_phase(double *phase, size_t n) {
    double offset = 0.0;
    for (size_t k = 1; k < n; ++k) {
        double diff = phase[k] + offset - phase[k - 1];
        /* Pull the step into (-pi, pi]. */
        while (diff > M_PI)  { offset -= 2.0 * M_PI; diff -= 2.0 * M_PI; }
        while (diff <= -M_PI){ offset += 2.0 * M_PI; diff += 2.0 * M_PI; }
        phase[k] += offset;
    }
}

int dsp_cepstrum_complex(const double *x, size_t n, cplx *cep) {
    if (n == 0 || !dsp_is_pow2(n))
        return -1;

    cplx   *spec  = malloc(n * sizeof(cplx));
    double *phase = malloc(n * sizeof(double));
    if (!spec || !phase) {
        free(spec); free(phase);
        return -1;
    }

    for (size_t i = 0; i < n; ++i)
        spec[i] = x[i];
    if (dsp_fft(spec, n) != 0) {
        free(spec); free(phase);
        return -1;
    }

    /* Pull out the raw (wrapped) phase, then unwrap it. */
    for (size_t k = 0; k < n; ++k)
        phase[k] = carg(spec[k]);
    unwrap_phase(phase, n);

    /* Complex log: real part log|X|, imaginary part the unwrapped
     * phase. Inverse-transform to get the complex cepstrum. */
    for (size_t k = 0; k < n; ++k) {
        double mag = cabs(spec[k]);
        if (mag < CEP_LOG_FLOOR)
            mag = CEP_LOG_FLOOR;
        spec[k] = dsp_cplx(log(mag), phase[k]);
    }
    if (dsp_ifft(spec, n) != 0) {
        free(spec); free(phase);
        return -1;
    }
    memcpy(cep, spec, n * sizeof(cplx));

    free(spec); free(phase);
    return 0;
}

int dsp_icepstrum_complex(const cplx *cep, size_t n, double *x) {
    if (n == 0 || !dsp_is_pow2(n))
        return -1;

    cplx *buf = malloc(n * sizeof(cplx));
    if (!buf)
        return -1;

    /* Reverse the forward pipeline: FFT the cepstrum to recover the
     * complex log spectrum, exponentiate to undo the log, then
     * inverse-FFT back to the time domain. */
    memcpy(buf, cep, n * sizeof(cplx));
    if (dsp_fft(buf, n) != 0) {
        free(buf);
        return -1;
    }
    for (size_t k = 0; k < n; ++k)
        buf[k] = cexp(buf[k]);          /* exp of the complex log */
    if (dsp_ifft(buf, n) != 0) {
        free(buf);
        return -1;
    }
    for (size_t i = 0; i < n; ++i)
        x[i] = creal(buf[i]);

    free(buf);
    return 0;
}

size_t dsp_cepstrum_pitch(const double *x, size_t n,
                          size_t min_q, size_t max_q) {
    if (max_q >= n)
        max_q = n - 1;
    if (min_q < 1)
        min_q = 1;
    if (min_q >= max_q)
        return 0;

    double *cep = malloc(n * sizeof(double));
    if (!cep)
        return 0;
    if (dsp_cepstrum_real(x, n, cep) != 0) {
        free(cep);
        return 0;
    }

    /* A periodic signal puts a sharp peak in the real cepstrum at a
     * quefrency equal to its period. Find the largest value in the
     * plausible quefrency window. */
    size_t best_q = 0;
    double best_v = -1e300;
    for (size_t q = min_q; q <= max_q; ++q) {
        if (cep[q] > best_v) {
            best_v = cep[q];
            best_q = q;
        }
    }

    free(cep);
    return best_q;
}

/* ===================================================================
 * MFCC - mel filterbank and the cepstral pipeline
 * =================================================================== */

/* Hz <-> mel conversions (the common O'Shaughnessy formula). */
static double hz_to_mel(double hz) {
    return 2595.0 * log10(1.0 + hz / 700.0);
}
static double mel_to_hz(double mel) {
    return 700.0 * (pow(10.0, mel / 2595.0) - 1.0);
}

int dsp_mel_filterbank_init(dsp_mel_filterbank *fb, size_t nfilters,
                            size_t fft_size, double sample_rate) {
    if (nfilters == 0 || fft_size < 2 || sample_rate <= 0.0)
        return -1;

    size_t nbins = fft_size / 2 + 1;
    fb->nfilters    = nfilters;
    fb->fft_size    = fft_size;
    fb->sample_rate = sample_rate;
    fb->weights     = calloc(nfilters * nbins, sizeof(double));
    if (!fb->weights)
        return -1;

    /* Place nfilters+2 points equally on the mel scale, spanning 0 to
     * the Nyquist frequency. Each filter is a triangle whose base
     * runs from point i to point i+2 and whose apex is at point i+1. */
    double mel_lo = hz_to_mel(0.0);
    double mel_hi = hz_to_mel(sample_rate / 2.0);

    size_t npts = nfilters + 2;
    double *bin = malloc(npts * sizeof(double));
    if (!bin) {
        free(fb->weights);
        fb->weights = NULL;
        return -1;
    }
    /* Convert each mel point to an FFT bin index (fractional). */
    for (size_t i = 0; i < npts; ++i) {
        double mel = mel_lo + (mel_hi - mel_lo) * (double)i
                                              / (double)(npts - 1);
        double hz  = mel_to_hz(mel);
        bin[i] = hz * (double)fft_size / sample_rate;
    }

    /* Fill each filter's triangular weights across the magnitude
     * bins: rising from bin[f] to bin[f+1], falling to bin[f+2]. */
    for (size_t f = 0; f < nfilters; ++f) {
        double left   = bin[f];
        double centre = bin[f + 1];
        double right  = bin[f + 2];
        for (size_t k = 0; k < nbins; ++k) {
            double w = 0.0;
            double kk = (double)k;
            if (kk >= left && kk <= centre && centre > left)
                w = (kk - left) / (centre - left);
            else if (kk > centre && kk <= right && right > centre)
                w = (right - kk) / (right - centre);
            fb->weights[f * nbins + k] = w;
        }
    }

    free(bin);
    return 0;
}

void dsp_mel_filterbank_free(dsp_mel_filterbank *fb) {
    free(fb->weights);
    fb->weights = NULL;
    fb->nfilters = fb->fft_size = 0;
}

int dsp_mfcc(const double *frame, size_t fft_size,
             const dsp_mel_filterbank *fb,
             size_t ncoeffs, double *mfcc) {
    if (fft_size == 0 || !dsp_is_pow2(fft_size) || ncoeffs == 0)
        return -1;
    if (fb->fft_size != fft_size || ncoeffs > fb->nfilters)
        return -1;

    size_t nbins = fft_size / 2 + 1;

    cplx   *spec   = malloc(fft_size * sizeof(cplx));
    double *power  = malloc(nbins * sizeof(double));
    double *energy = malloc(fb->nfilters * sizeof(double));
    double *dctout = malloc(fb->nfilters * sizeof(double));
    if (!spec || !power || !energy || !dctout) {
        free(spec); free(power); free(energy); free(dctout);
        return -1;
    }

    /* 1. Power spectrum of the frame. */
    for (size_t i = 0; i < fft_size; ++i)
        spec[i] = frame[i];
    if (dsp_fft(spec, fft_size) != 0) {
        free(spec); free(power); free(energy); free(dctout);
        return -1;
    }
    for (size_t k = 0; k < nbins; ++k) {
        double re = creal(spec[k]);
        double im = cimag(spec[k]);
        power[k] = re * re + im * im;
    }

    /* 2. Mel filterbank: each band's energy is the power spectrum
     *    weighted by that triangular filter and summed.
     * 3. Log of each band energy. */
    for (size_t f = 0; f < fb->nfilters; ++f) {
        double acc = 0.0;
        const double *w = fb->weights + f * nbins;
        for (size_t k = 0; k < nbins; ++k)
            acc += w[k] * power[k];
        if (acc < CEP_LOG_FLOOR)
            acc = CEP_LOG_FLOOR;
        energy[f] = log(acc);
    }

    /* 4. DCT of the log-mel energies; keep the first ncoeffs. The DCT
     *    both performs the cepstral transform and decorrelates the
     *    coefficients. */
    dsp_dct(energy, dctout, fb->nfilters);
    for (size_t c = 0; c < ncoeffs; ++c)
        mfcc[c] = dctout[c];

    free(spec); free(power); free(energy); free(dctout);
    return 0;
}
