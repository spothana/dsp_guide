/*
 * resample.h - Sample rate conversion
 *
 * PROBLEM SOLVED
 *   Changing a signal's sample rate without introducing aliasing
 *   (when slowing down) or imaging (when speeding up) artifacts.
 *
 * DECIMATION (downsample by M)
 *   1. Low-pass anti-aliasing filter, cutoff pi/M.
 *   2. Keep every M-th sample, discard the rest.
 *   The filter MUST precede the discard - aliasing cannot be undone.
 *
 * INTERPOLATION (upsample by L)
 *   1. Insert L-1 zeros between every sample (creates spectral images).
 *   2. Low-pass anti-imaging filter, cutoff pi/L, gain L.
 *   The filter goes AFTER the zero-stuffing.
 *
 * MULTI-RATE / RATIONAL RESAMPLING (by L/M)
 *   Upsample by L, then downsample by M. Example: 44.1 kHz -> 48 kHz
 *   is L/M = 160/147. The two anti-alias/anti-image filters can be
 *   merged into one running at the high intermediate rate.
 *
 * MNEMONIC
 *   The low-pass filter always lives on the HIGH-rate side of the
 *   operation.
 */
#ifndef DSP_RESAMPLE_H
#define DSP_RESAMPLE_H

#include "../common.h"

/*
 * Decimate x (length n) by integer factor m.
 *   Output length is (n + m - 1) / m, written to y (caller-allocated).
 * An internal windowed-sinc low-pass filter (cutoff ~0.5/m) is applied
 * before discarding samples. Returns the output length, or 0 on error.
 */
size_t dsp_decimate(const double *x, size_t n, size_t m, double *y);

/* Output length produced by dsp_decimate for inputs n, m. */
static inline size_t dsp_decimate_len(size_t n, size_t m) {
    return (m == 0) ? 0 : (n + m - 1) / m;
}

/*
 * Interpolate x (length n) by integer factor l.
 *   Output length is n * l, written to y (caller-allocated).
 * Zero-stuffs then applies an internal low-pass anti-imaging filter
 * (cutoff ~0.5/l, gain l). Returns the output length, or 0 on error.
 */
size_t dsp_interpolate(const double *x, size_t n, size_t l, double *y);

/* Output length produced by dsp_interpolate for inputs n, l. */
static inline size_t dsp_interpolate_len(size_t n, size_t l) {
    return n * l;
}

/*
 * Rational resample x (length n) by the factor l/m: upsample by l,
 * then decimate by m. Output length is dsp_resample_len(n, l, m),
 * written to y. Returns the output length, or 0 on error.
 */
size_t dsp_resample(const double *x, size_t n,
                    size_t l, size_t m, double *y);

/* Output length produced by dsp_resample for inputs n, l, m. */
static inline size_t dsp_resample_len(size_t n, size_t l, size_t m) {
    return dsp_decimate_len(dsp_interpolate_len(n, l), m);
}

#endif /* DSP_RESAMPLE_H */
