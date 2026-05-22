/*
 * hilbert.c - Hilbert transform, analytic signal, and the quantities
 *             that fall out of it (envelope, instantaneous frequency).
 */
#include "transforms/hilbert.h"
#include "transforms/fft.h"
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int dsp_analytic_signal(const double *x, size_t n, cplx *z) {
    if (n == 0 || !dsp_is_pow2(n))
        return -1;

    /* Start with the real signal in the complex buffer and FFT it. */
    for (size_t i = 0; i < n; ++i)
        z[i] = x[i];
    if (dsp_fft(z, n) != 0)
        return -1;

    /* Build the one-sided spectrum: keep DC and Nyquist, double the
     * positive-frequency bins, zero the negative-frequency bins. The
     * inverse FFT of that is the analytic signal.
     *
     * Note the sign convention of this FFT: a positive-frequency
     * exponential e^{+j...} lands in the UPPER half of the spectrum
     * (bins n/2+1 .. n-1), so it is those bins that are doubled and
     * the lower half that is zeroed. */
    for (size_t k = 1; k < n / 2; ++k)
        z[k] = 0.0;
    for (size_t k = n / 2 + 1; k < n; ++k)
        z[k] *= 2.0;
    /* z[0] (DC) and z[n/2] (Nyquist) are left unchanged. */

    if (dsp_ifft(z, n) != 0)
        return -1;
    return 0;
}

int dsp_hilbert(const double *x, size_t n, double *h) {
    /* The Hilbert transform is the imaginary part of the analytic
     * signal. */
    cplx *z = malloc(n * sizeof(cplx));
    if (!z)
        return -1;
    if (dsp_analytic_signal(x, n, z) != 0) {
        free(z);
        return -1;
    }
    for (size_t i = 0; i < n; ++i)
        h[i] = cimag(z[i]);
    free(z);
    return 0;
}

int dsp_envelope(const double *x, size_t n, double *env) {
    /* The envelope is the magnitude of the analytic signal. */
    cplx *z = malloc(n * sizeof(cplx));
    if (!z)
        return -1;
    if (dsp_analytic_signal(x, n, z) != 0) {
        free(z);
        return -1;
    }
    for (size_t i = 0; i < n; ++i)
        env[i] = cabs(z[i]);
    free(z);
    return 0;
}

int dsp_instantaneous_frequency(const double *x, size_t n,
                                double *freq) {
    if (n < 2)
        return -1;

    cplx *z = malloc(n * sizeof(cplx));
    if (!z)
        return -1;
    if (dsp_analytic_signal(x, n, z) != 0) {
        free(z);
        return -1;
    }

    /* Instantaneous frequency is the rate of change of the
     * instantaneous phase arg(z). Rather than unwrap the phase
     * explicitly, take the phase DIFFERENCE between adjacent samples
     * as the argument of z[i] * conj(z[i-1]) - this lands in
     * (-pi, pi] automatically, which is the unwrapping. Divide by
     * 2*pi to get cycles per sample. */
    for (size_t i = 1; i < n; ++i) {
        cplx ratio = z[i] * conj(z[i - 1]);
        double dphase = carg(ratio);          /* already in (-pi, pi] */
        freq[i] = dphase / (2.0 * M_PI);
    }
    /* freq[0] has no preceding sample; repeat the first valid value. */
    freq[0] = freq[1];

    free(z);
    return 0;
}
