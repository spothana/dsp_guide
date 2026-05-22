/*
 * ofdm.c - OFDM modulator, demodulator, and per-subcarrier equalizer.
 */
#include "modulation/ofdm.h"
#include "transforms/fft.h"
#include <stdlib.h>
#include <string.h>

int dsp_ofdm_init(dsp_ofdm *o, size_t nfft, size_t cp_len) {
    if (!dsp_is_pow2(nfft))
        return -1;
    if (cp_len >= nfft)
        return -1;
    o->nfft   = nfft;
    o->cp_len = cp_len;
    return 0;
}

size_t dsp_ofdm_modulate(const dsp_ofdm *o, const cplx *freq,
                         cplx *time) {
    size_t nfft = o->nfft, cp = o->cp_len;

    /* Step 1: IFFT - turn the per-subcarrier symbols into a
     * time-domain block. The IFFT is the OFDM modulator. */
    cplx *block = malloc(nfft * sizeof(cplx));
    if (!block)
        return 0;
    memcpy(block, freq, nfft * sizeof(cplx));
    if (dsp_ifft(block, nfft) != 0) {
        free(block);
        return 0;
    }

    /* The inverse FFT here carries a full 1/N normalisation, which
     * would leave each time sample with only 1/N of a subcarrier's
     * power. Multiply by sqrt(N) so the transmitted signal has the
     * same average power as the frequency-domain constellation - the
     * standard unitary OFDM convention. The demodulator divides it
     * back out. Without this, a time-domain noise level cannot be
     * related to the per-symbol SNR. */
    double scale = sqrt((double)nfft);
    for (size_t i = 0; i < nfft; ++i)
        block[i] *= scale;

    /* Step 2: cyclic prefix - copy the last cp samples to the front.
     * Output layout: [ block[nfft-cp .. nfft-1] | block[0 .. nfft-1] ]. */
    for (size_t i = 0; i < cp; ++i)
        time[i] = block[nfft - cp + i];
    for (size_t i = 0; i < nfft; ++i)
        time[cp + i] = block[i];

    free(block);
    return nfft + cp;
}

size_t dsp_ofdm_demodulate(const dsp_ofdm *o, const cplx *time,
                           cplx *freq) {
    size_t nfft = o->nfft, cp = o->cp_len;

    /* Step 1: discard the cyclic prefix, keep the nfft payload
     * samples. Because the prefix absorbed the multipath transient,
     * what remains is one clean circular-convolution period. */
    for (size_t i = 0; i < nfft; ++i)
        freq[i] = time[cp + i];

    /* Step 2: FFT - recover every subcarrier symbol at once. */
    if (dsp_fft(freq, nfft) != 0)
        return 0;

    /* Undo the sqrt(N) power scaling applied by the modulator, so a
     * clean channel returns the original constellation unchanged. */
    double scale = 1.0 / sqrt((double)nfft);
    for (size_t k = 0; k < nfft; ++k)
        freq[k] *= scale;

    return nfft;
}

void dsp_ofdm_equalize(const dsp_ofdm *o, cplx *freq,
                       const cplx *chan_fr) {
    /* Thanks to the cyclic prefix, the channel acts as a pointwise
     * multiply in the FFT domain: Y[k] = H[k] * X[k]. Zero-forcing
     * equalization simply inverts that - one complex divide per
     * subcarrier. (A guard against a near-zero gain keeps it finite.) */
    for (size_t k = 0; k < o->nfft; ++k) {
        cplx h = chan_fr[k];
        double mag2 = creal(h) * creal(h) + cimag(h) * cimag(h);
        if (mag2 < 1e-12)
            continue;                 /* deep fade - leave as-is */
        freq[k] = freq[k] / h;
    }
}
