/*
 * pulse.c - Root-raised-cosine pulse shaping.
 */
#include "modulation/pulse.h"
#include <stdlib.h>

void dsp_rrc_design(double *taps, size_t ntaps, size_t sps, double beta) {
    /* Centre index - the filter is symmetric about this point. */
    double mid = (double)(ntaps - 1) / 2.0;

    for (size_t k = 0; k < ntaps; ++k) {
        /* Time in symbol units: t = (k - mid) / sps. */
        double t = ((double)k - mid) / (double)sps;
        double h;

        if (fabs(t) < 1e-9) {
            /* Closed-form limit at t = 0. */
            h = 1.0 - beta + 4.0 * beta / M_PI;
        } else if (beta > 0.0 &&
                   fabs(fabs(t) - 1.0 / (4.0 * beta)) < 1e-9) {
            /* Closed-form limit at the removable singularity
             * t = +/- 1/(4*beta). */
            double a = (1.0 + 2.0 / M_PI)
                     * sin(M_PI / (4.0 * beta));
            double b = (1.0 - 2.0 / M_PI)
                     * cos(M_PI / (4.0 * beta));
            h = (beta / sqrt(2.0)) * (a + b);
        } else {
            /* General RRC impulse response. */
            double num = sin(M_PI * t * (1.0 - beta))
                       + 4.0 * beta * t
                         * cos(M_PI * t * (1.0 + beta));
            double den = M_PI * t
                       * (1.0 - (4.0 * beta * t) * (4.0 * beta * t));
            h = num / den;
        }
        taps[k] = h;
    }

    /* Normalise to unit energy so a TX/RX RRC pair has unity gain. */
    double energy = 0.0;
    for (size_t k = 0; k < ntaps; ++k)
        energy += taps[k] * taps[k];
    double scale = (energy > 0.0) ? 1.0 / sqrt(energy) : 1.0;
    for (size_t k = 0; k < ntaps; ++k)
        taps[k] *= scale;
}

size_t dsp_pulse_shape(const cplx *syms, size_t nsym, size_t sps,
                       const double *taps, size_t ntaps, cplx *out) {
    size_t n = nsym * sps;

    /* Step 1: upsample into a scratch buffer - place each symbol at
     * index k*sps with zeros between (the zero-stuffing the RRC then
     * interpolates over). */
    cplx *up = malloc(n * sizeof(cplx));
    if (!up)
        return 0;
    for (size_t i = 0; i < n; ++i)
        up[i] = 0.0;
    for (size_t k = 0; k < nsym; ++k)
        up[k * sps] = syms[k];

    /* Step 2: filter the upsampled stream with the RRC taps
     * (direct-form FIR, zero initial state). */
    for (size_t i = 0; i < n; ++i) {
        cplx acc = 0.0;
        for (size_t j = 0; j < ntaps; ++j) {
            if (i >= j)
                acc += taps[j] * up[i - j];
        }
        out[i] = acc;
    }

    free(up);
    return n;
}

void dsp_matched_filter(const cplx *in, size_t n,
                        const double *taps, size_t ntaps, cplx *out) {
    /* Convolve the received samples with the RRC taps. The RRC is
     * symmetric, so it is its own matched filter. */
    for (size_t i = n; i-- > 0; ) {
        cplx acc = 0.0;
        for (size_t j = 0; j < ntaps; ++j) {
            if (i >= j)
                acc += taps[j] * in[i - j];
        }
        out[i] = acc;
    }
}
