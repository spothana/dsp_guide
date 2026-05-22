/*
 * common.h - Shared types and helpers for the DSP Guide
 *
 * A minimal complex-number type and a few utilities. We use C99's
 * native <complex.h> double _Complex for arithmetic, but wrap it in
 * a `cplx` typedef so the intent is obvious at call sites.
 */
#ifndef DSP_COMMON_H
#define DSP_COMMON_H

#include <complex.h>
#include <math.h>
#include <stddef.h>

/* Native C99 complex double. All transforms operate on this type. */
typedef double _Complex cplx;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Construct a complex number from real and imaginary parts. */
static inline cplx dsp_cplx(double re, double im) {
    return re + im * I;
}

/* Magnitude (absolute value) of a complex number. */
static inline double dsp_mag(cplx z) {
    return cabs(z);
}

/* Phase angle in radians, in (-pi, pi]. */
static inline double dsp_phase(cplx z) {
    return carg(z);
}

/* True when n is a power of two (n > 0). Used by the radix-2 FFT. */
static inline int dsp_is_pow2(size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

/* Smallest power of two >= n. Used for zero-padding before an FFT. */
static inline size_t dsp_next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

#endif /* DSP_COMMON_H */
