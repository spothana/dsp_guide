/*
 * window.h - Window functions for spectral analysis
 *
 * PROBLEM SOLVED
 *   Taking the DFT of a finite signal implicitly multiplies it by a
 *   rectangular window. The abrupt edges create artificial
 *   high-frequency content - SPECTRAL LEAKAGE. Window functions taper
 *   the signal smoothly to zero at its edges to reduce that leakage.
 *
 * THE CORE TRADE-OFF
 *   main-lobe width  vs  side-lobe level.
 *   A narrow main lobe resolves closely spaced frequencies; low side
 *   lobes keep a weak tone from being buried next to a strong one.
 *   You cannot maximise both.
 *
 * WINDOWS PROVIDED (peak side-lobe level, roughly)
 *   Rectangular : -13 dB  narrowest main lobe, best raw resolution
 *   Hamming     : -43 dB  good general-purpose choice
 *   Hanning     : -31 dB  similar, better side-lobe roll-off
 *   Blackman    : -58 dB  lowest side lobes, widest main lobe
 *
 * Each generator fills w[0..n-1] with the window coefficients; apply
 * a window by multiplying it pointwise with the signal before the FFT.
 */
#ifndef DSP_WINDOW_H
#define DSP_WINDOW_H

#include "../common.h"

/* Window identifiers. */
typedef enum {
    DSP_WIN_RECTANGULAR,
    DSP_WIN_HAMMING,
    DSP_WIN_HANNING,
    DSP_WIN_BLACKMAN
} dsp_window_type;

/* Rectangular window: all ones (i.e. no tapering). */
void dsp_window_rectangular(double *w, size_t n);

/* Hamming window: 0.54 - 0.46 cos(...). */
void dsp_window_hamming(double *w, size_t n);

/* Hanning (Hann) window: 0.5 - 0.5 cos(...). */
void dsp_window_hanning(double *w, size_t n);

/* Blackman window: 0.42 - 0.5 cos(...) + 0.08 cos(2...). */
void dsp_window_blackman(double *w, size_t n);

/* Dispatch by type - fills w[0..n-1] with the chosen window. */
void dsp_window_generate(dsp_window_type type, double *w, size_t n);

/* Multiply signal x by window w pointwise, writing into out. */
void dsp_window_apply(const double *x, const double *w,
                      double *out, size_t n);

#endif /* DSP_WINDOW_H */
