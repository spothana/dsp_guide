/*
 * sync.c - Carrier-recovery PLL and Gardner timing recovery.
 */
#include "modulation/sync.h"

/* ===================================================================
 * Carrier recovery - decision-directed PLL
 * =================================================================== */

void dsp_carrier_pll_init(dsp_carrier_pll *pll,
                          double alpha, double beta) {
    pll->phase = 0.0;
    pll->freq  = 0.0;
    pll->alpha = alpha;
    pll->beta  = beta;
}

cplx dsp_carrier_pll_process(dsp_carrier_pll *pll, cplx in) {
    /* De-rotate by the current phase estimate. */
    cplx corrected = in * cexp(-pll->phase * I);

    /* Decision-directed phase error: slice the corrected sample to
     * the nearest QPSK quadrant (sign of I and Q), then measure the
     * angle between the sample and that ideal point. */
    double i = creal(corrected);
    double q = cimag(corrected);
    cplx decision = dsp_cplx((i >= 0.0) ? 1.0 : -1.0,
                             (q >= 0.0) ? 1.0 : -1.0);

    /* The phase error is the argument of corrected * conj(decision):
     * how far the sample sits from its assumed-correct quadrant. */
    cplx err_vec = corrected * conj(decision);
    double error = carg(err_vec);

    /* Second-order loop filter: the proportional term corrects phase
     * directly, the integral term accumulates into a frequency
     * estimate so a constant frequency offset is tracked out. */
    pll->freq  += pll->beta * error;
    pll->phase += pll->freq + pll->alpha * error;

    /* Keep the phase bounded for numerical tidiness. */
    while (pll->phase >  M_PI) pll->phase -= 2.0 * M_PI;
    while (pll->phase < -M_PI) pll->phase += 2.0 * M_PI;

    return corrected;
}

void dsp_carrier_recover(dsp_carrier_pll *pll,
                         const cplx *in, cplx *out, size_t n) {
    for (size_t k = 0; k < n; ++k)
        out[k] = dsp_carrier_pll_process(pll, in[k]);
}

/* ===================================================================
 * Timing recovery - Gardner timing-error detector
 * =================================================================== */

double dsp_timing_error_gardner(const cplx *samples, size_t n,
                                size_t sps) {
    if (sps < 2 || (sps & 1))
        return 0.0;                 /* Gardner needs an even sps */

    size_t half = sps / 2;
    double acc = 0.0;
    size_t count = 0;

    /* For each pair of adjacent symbol peaks, the Gardner error uses
     * the mid-point sample between them:
     *   e = real( mid * conj(curr_peak - prev_peak) )
     * It is zero when the peaks are correctly timed. */
    for (size_t k = sps; k + sps <= n; k += sps) {
        cplx prev = samples[k - sps];
        cplx curr = samples[k];
        cplx mid  = samples[k - half];

        cplx diff = curr - prev;
        acc += creal(mid) * creal(diff) + cimag(mid) * cimag(diff);
        ++count;
    }

    if (count == 0)
        return 0.0;

    /* The averaged error, scaled into a fractional-sample offset.
     * The sign points toward the correction direction. */
    return (acc / (double)count) * (double)half;
}

size_t dsp_timing_resample(const cplx *samples, size_t n, size_t sps,
                           double offset, cplx *syms) {
    size_t nsym = n / sps;
    size_t out  = 0;

    for (size_t k = 0; k < nsym; ++k) {
        /* Ideal sampling instant for symbol k, shifted by `offset`. */
        double pos = (double)(k * sps) + offset;

        /* Linear interpolation between the two straddling samples. */
        long   base = (long)floor(pos);
        double frac = pos - (double)base;

        if (base < 0 || (size_t)(base + 1) >= n)
            continue;               /* instant falls outside the data */

        cplx a = samples[base];
        cplx b = samples[base + 1];
        syms[out++] = a + (b - a) * frac;
    }
    return out;
}
