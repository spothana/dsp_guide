/*
 * window.c - Rectangular, Hamming, Hanning, and Blackman windows.
 *
 * Each window is defined over k = 0 .. n-1 with the standard
 * "symmetric" normalisation, denominator (n - 1).
 */
#include "spectral/window.h"

void dsp_window_rectangular(double *w, size_t n) {
    for (size_t k = 0; k < n; ++k)
        w[k] = 1.0;
}

void dsp_window_hamming(double *w, size_t n) {
    if (n == 1) { w[0] = 1.0; return; }
    for (size_t k = 0; k < n; ++k)
        w[k] = 0.54 - 0.46 * cos(2.0 * M_PI * (double)k
                                      / (double)(n - 1));
}

void dsp_window_hanning(double *w, size_t n) {
    if (n == 1) { w[0] = 1.0; return; }
    for (size_t k = 0; k < n; ++k)
        w[k] = 0.5 - 0.5 * cos(2.0 * M_PI * (double)k
                                    / (double)(n - 1));
}

void dsp_window_blackman(double *w, size_t n) {
    if (n == 1) { w[0] = 1.0; return; }
    for (size_t k = 0; k < n; ++k) {
        double t = 2.0 * M_PI * (double)k / (double)(n - 1);
        w[k] = 0.42 - 0.5 * cos(t) + 0.08 * cos(2.0 * t);
    }
}

void dsp_window_generate(dsp_window_type type, double *w, size_t n) {
    switch (type) {
        case DSP_WIN_RECTANGULAR: dsp_window_rectangular(w, n); break;
        case DSP_WIN_HAMMING:     dsp_window_hamming(w, n);     break;
        case DSP_WIN_HANNING:     dsp_window_hanning(w, n);     break;
        case DSP_WIN_BLACKMAN:    dsp_window_blackman(w, n);    break;
        default:                  dsp_window_rectangular(w, n); break;
    }
}

void dsp_window_apply(const double *x, const double *w,
                      double *out, size_t n) {
    for (size_t k = 0; k < n; ++k)
        out[k] = x[k] * w[k];
}
