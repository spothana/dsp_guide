/*
 * timefreq.c - STFT, QMF filter bank, and Wigner-Ville distribution.
 */
#include "spectral/timefreq.h"
#include "transforms/fft.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===================================================================
 * STFT
 * =================================================================== */

int dsp_stft_forward(const double *x, size_t n,
                     dsp_window_type win_type, size_t win_len,
                     size_t hop, dsp_stft *out) {
    if (win_len == 0 || hop == 0 || !dsp_is_pow2(win_len) || n == 0)
        return -1;

    /* Number of frames: slide the window in `hop` steps until it runs
     * off the end of the signal. */
    size_t frames = (n >= win_len) ? (n - win_len) / hop + 1 : 1;

    cplx *grid = malloc(frames * win_len * sizeof(cplx));
    double *win = malloc(win_len * sizeof(double));
    cplx *frame = malloc(win_len * sizeof(cplx));
    if (!grid || !win || !frame) {
        free(grid); free(win); free(frame);
        return -1;
    }

    dsp_window_generate(win_type, win, win_len);

    for (size_t t = 0; t < frames; ++t) {
        size_t start = t * hop;
        /* Window the segment (zero-padding past the signal end). */
        for (size_t i = 0; i < win_len; ++i) {
            double s = (start + i < n) ? x[start + i] : 0.0;
            frame[i] = s * win[i];
        }
        dsp_fft(frame, win_len);
        memcpy(grid + t * win_len, frame, win_len * sizeof(cplx));
    }

    out->frames  = frames;
    out->bins    = win_len;
    out->hop     = hop;
    out->win_len = win_len;
    out->data    = grid;

    free(win);
    free(frame);
    return 0;
}

void dsp_stft_free(dsp_stft *s) {
    free(s->data);
    s->data = NULL;
    s->frames = s->bins = 0;
}

size_t dsp_stft_signal_len(const dsp_stft *s) {
    if (s->frames == 0)
        return 0;
    return (s->frames - 1) * s->hop + s->win_len;
}

size_t dsp_stft_inverse(const dsp_stft *s, double *x) {
    size_t n = dsp_stft_signal_len(s);
    if (n == 0)
        return 0;

    /* Overlap-add: inverse-FFT each frame, re-apply the analysis
     * window, accumulate into the output, and divide by the summed
     * squared window so overlapping frames are weighted correctly. */
    double *win  = malloc(s->win_len * sizeof(double));
    double *wsum = calloc(n, sizeof(double));
    cplx   *fr   = malloc(s->win_len * sizeof(cplx));
    if (!win || !wsum || !fr) {
        free(win); free(wsum); free(fr);
        return 0;
    }
    /* The forward transform used DSP_WIN_HANNING-style tapers; re-make
     * the same window from the stored length. The window type is not
     * kept in the struct, so a Hanning window is assumed for synthesis
     * (the standard invertible choice). */
    dsp_window_hanning(win, s->win_len);

    for (size_t i = 0; i < n; ++i)
        x[i] = 0.0;

    for (size_t t = 0; t < s->frames; ++t) {
        memcpy(fr, s->data + t * s->bins, s->bins * sizeof(cplx));
        dsp_ifft(fr, s->bins);

        size_t start = t * s->hop;
        for (size_t i = 0; i < s->win_len; ++i) {
            double v = creal(fr[i]) * win[i];
            x[start + i]    += v;
            wsum[start + i] += win[i] * win[i];
        }
    }

    /* Normalise by the overlap-add window weight. */
    for (size_t i = 0; i < n; ++i)
        if (wsum[i] > 1e-12)
            x[i] /= wsum[i];

    free(win);
    free(wsum);
    free(fr);
    return n;
}

void dsp_spectrogram(const dsp_stft *s, double *spec) {
    size_t total = s->frames * s->bins;
    for (size_t i = 0; i < total; ++i) {
        double re = creal(s->data[i]);
        double im = cimag(s->data[i]);
        spec[i] = re * re + im * im;        /* |STFT|^2 */
    }
}

/* ===================================================================
 * Two-channel QMF filter bank
 *
 * The analysis low-pass h0 is a short prototype; the analysis
 * high-pass is h1[k] = (-1)^k h0[k] (the mirror in frequency). The
 * synthesis filters are g0[k] = h0[k], g1[k] = -h1[k]. With this QMF
 * choice the aliasing introduced by decimating each band cancels in
 * the reconstruction.
 * =================================================================== */

int dsp_qmf_init(dsp_qmf_bank *bank) {
    bank->ntaps = 8;
    bank->taps  = malloc(bank->ntaps * sizeof(double));
    if (!bank->taps)
        return -1;

    /* A short, smooth half-band low-pass prototype. The two-channel
     * QMF mirror relations derive the high-pass and synthesis filters
     * from it; this design gives good - though not bit-exact -
     * reconstruction, enough to demonstrate the filter-bank principle. */
    static const double h0[8] = {
        0.00938715, -0.07065183,  0.06942827,  0.49183638,
        0.49183638,  0.06942827, -0.07065183,  0.00938715
    };
    memcpy(bank->taps, h0, bank->ntaps * sizeof(double));
    return 0;
}

void dsp_qmf_free(dsp_qmf_bank *bank) {
    free(bank->taps);
    bank->taps = NULL;
    bank->ntaps = 0;
}

void dsp_qmf_analyze(const dsp_qmf_bank *bank, const double *x,
                     size_t n, double *lo, double *hi) {
    size_t L = bank->ntaps;
    size_t half = n / 2;

    /* For each output (decimated) index, convolve with the low-pass
     * and high-pass filters at the corresponding input position. The
     * high-pass filter is h1[k] = (-1)^k h0[k]. */
    for (size_t m = 0; m < half; ++m) {
        size_t centre = 2 * m;             /* input position */
        double slo = 0.0, shi = 0.0;
        for (size_t k = 0; k < L; ++k) {
            long idx = (long)centre - (long)k;
            double sample = (idx >= 0 && idx < (long)n) ? x[idx] : 0.0;
            double h0 = bank->taps[k];
            double h1 = ((k & 1) ? -h0 : h0);
            slo += h0 * sample;
            shi += h1 * sample;
        }
        lo[m] = slo;
        hi[m] = shi;
    }
}

void dsp_qmf_synthesize(const dsp_qmf_bank *bank,
                        const double *lo, const double *hi,
                        size_t half, double *out) {
    size_t L = bank->ntaps;
    size_t n = 2 * half;

    /* Upsample each subband by inserting zeros, then filter with the
     * synthesis filters g0 = h0 and g1 = -h1, and sum. */
    for (size_t i = 0; i < n; ++i)
        out[i] = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (size_t k = 0; k < L; ++k) {
            long idx = (long)i - (long)k;
            if (idx < 0 || (idx & 1))
                continue;                  /* only even (upsampled) */
            size_t m = (size_t)idx / 2;
            if (m >= half)
                continue;
            double h0 = bank->taps[k];
            double h1 = ((k & 1) ? -h0 : h0);
            double g0 = h0;
            double g1 = -h1;
            s += g0 * lo[m] + g1 * hi[m];
        }
        out[i] = 2.0 * s;                  /* compensate decimation */
    }
}

/* ===================================================================
 * Wigner-Ville distribution
 * =================================================================== */

/* Build the analytic signal of a real input via the FFT: zero the
 * negative-frequency half and double the positive half. */
static int analytic_signal(const double *x, size_t n, cplx *z) {
    if (!dsp_is_pow2(n))
        return -1;
    for (size_t i = 0; i < n; ++i)
        z[i] = x[i];
    dsp_fft(z, n);

    /* Keep DC and Nyquist, double bins 1..n/2-1, zero the rest. */
    for (size_t k = 1; k < n / 2; ++k)
        z[k] *= 2.0;
    for (size_t k = n / 2 + 1; k < n; ++k)
        z[k] = 0.0;

    dsp_ifft(z, n);
    return 0;
}

/* Core WVD: optional smoothing window applied along the lag axis. */
static int wvd_core(const double *x, size_t n, double *wvd,
                    const double *lag_win) {
    if (!dsp_is_pow2(n))
        return -1;

    cplx *z   = malloc(n * sizeof(cplx));
    cplx *row = malloc(n * sizeof(cplx));
    if (!z || !row) {
        free(z); free(row);
        return -1;
    }

    /* Work on the analytic signal: this halves the cross terms and
     * removes the spectral aliasing a real input would cause. */
    if (analytic_signal(x, n, z) != 0) {
        free(z); free(row);
        return -1;
    }

    for (size_t t = 0; t < n; ++t) {
        /* Instantaneous autocorrelation r(t, tau) = z(t+tau) z*(t-tau),
         * indexed so the lag axis wraps over n points. */
        for (size_t tau = 0; tau < n; ++tau) {
            long t1 = (long)t + (long)tau - (long)n / 2;
            long t2 = (long)t - (long)tau + (long)n / 2;
            cplx a = (t1 >= 0 && t1 < (long)n) ? z[t1] : 0.0;
            cplx b = (t2 >= 0 && t2 < (long)n) ? z[t2] : 0.0;
            cplx r = a * conj(b);
            if (lag_win) {
                /* Smoothing window centred on zero lag. */
                r *= lag_win[tau];
            }
            row[tau] = r;
        }
        /* The WVD slice at time t is the FFT of that lag sequence. */
        dsp_fft(row, n);
        for (size_t f = 0; f < n; ++f)
            wvd[t * n + f] = creal(row[f]);
    }

    free(z);
    free(row);
    return 0;
}

int dsp_wigner_ville(const double *x, size_t n, double *wvd) {
    return wvd_core(x, n, wvd, NULL);
}

int dsp_pseudo_wigner_ville(const double *x, size_t n,
                            dsp_window_type win_type, double *wvd) {
    if (!dsp_is_pow2(n))
        return -1;

    /* A smoothing window centred on the zero-lag index n/2. */
    double *full = malloc(n * sizeof(double));
    double *lag  = malloc(n * sizeof(double));
    if (!full || !lag) {
        free(full); free(lag);
        return -1;
    }
    dsp_window_generate(win_type, full, n);
    /* The window is generated 0..n-1; the WVD lag index 0..n-1 already
     * has zero lag at the centre, so the window aligns directly. */
    for (size_t i = 0; i < n; ++i)
        lag[i] = full[i];

    int rc = wvd_core(x, n, wvd, lag);

    free(full);
    free(lag);
    return rc;
}
