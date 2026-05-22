/*
 * channel.c - Multipath + AWGN channel model.
 */
#include "modulation/channel.h"
#include "transforms/fft.h"
#include <stdlib.h>
#include <string.h>

/* xorshift RNG for reproducible noise. */
static unsigned ch_rng(unsigned *s) {
    unsigned x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}

/* Uniform double in [0, 1). */
static double ch_uniform(unsigned *s) {
    return (ch_rng(s) >> 8) / 16777216.0;
}

/* One standard-normal sample (Box-Muller). */
static double ch_gauss(unsigned *s) {
    static const double TWO_PI = 6.28318530717958647692;
    double u1 = ch_uniform(s);
    double u2 = ch_uniform(s);
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
}

int dsp_channel_init(dsp_channel *ch, const cplx *taps, size_t ntaps,
                     double noise_std, unsigned seed) {
    if (ntaps == 0)
        return -1;
    ch->taps = malloc(ntaps * sizeof(cplx));
    if (!ch->taps)
        return -1;
    memcpy(ch->taps, taps, ntaps * sizeof(cplx));
    ch->ntaps     = ntaps;
    ch->noise_std = (noise_std < 0.0) ? 0.0 : noise_std;
    ch->rng_state = seed ? seed : 0xA5A5A5A5u;
    return 0;
}

int dsp_channel_init_awgn(dsp_channel *ch, double noise_std,
                          unsigned seed) {
    cplx unit = dsp_cplx(1.0, 0.0);
    return dsp_channel_init(ch, &unit, 1, noise_std, seed);
}

void dsp_channel_free(dsp_channel *ch) {
    free(ch->taps);
    ch->taps = NULL;
    ch->ntaps = 0;
}

void dsp_channel_apply(dsp_channel *ch, const cplx *in, cplx *out,
                       size_t n) {
    /* Step 1: multipath - convolve the input with the channel taps.
     * out[i] = sum_k taps[k] * in[i-k], zero state before the start. */
    for (size_t i = 0; i < n; ++i) {
        cplx acc = 0.0;
        for (size_t k = 0; k < ch->ntaps && k <= i; ++k)
            acc += ch->taps[k] * in[i - k];
        out[i] = acc;
    }

    /* Step 2: add complex AWGN. The noise power is split equally
     * between the real and imaginary parts, so each gets
     * noise_std / sqrt(2). */
    if (ch->noise_std > 0.0) {
        double sigma = ch->noise_std / sqrt(2.0);
        for (size_t i = 0; i < n; ++i) {
            double nr = sigma * ch_gauss(&ch->rng_state);
            double ni = sigma * ch_gauss(&ch->rng_state);
            out[i] += dsp_cplx(nr, ni);
        }
    }
}

void dsp_channel_frequency_response(const dsp_channel *ch,
                                    cplx *resp, size_t nfft) {
    /* Zero-pad the tap vector to nfft and take its DFT. The result is
     * the channel's gain at each of the nfft subcarrier frequencies. */
    for (size_t i = 0; i < nfft; ++i)
        resp[i] = (i < ch->ntaps) ? ch->taps[i] : 0.0;

    if (dsp_is_pow2(nfft)) {
        dsp_fft(resp, nfft);
    } else {
        /* Fall back to a direct DFT for non-power-of-two sizes. */
        cplx *tmp = malloc(nfft * sizeof(cplx));
        if (!tmp) return;
        for (size_t k = 0; k < nfft; ++k) {
            cplx s = 0.0;
            for (size_t t = 0; t < nfft; ++t) {
                double ang = -2.0 * 3.14159265358979323846
                           * (double)k * (double)t / (double)nfft;
                s += resp[t] * cexp(ang * I);
            }
            tmp[k] = s;
        }
        memcpy(resp, tmp, nfft * sizeof(cplx));
        free(tmp);
    }
}
