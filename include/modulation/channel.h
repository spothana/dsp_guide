/*
 * channel.h - Communication channel model
 *
 * PROBLEM SOLVED
 *   To test a modem honestly you need a realistic channel between
 *   transmitter and receiver. A real radio channel does two things:
 *
 *   1. MULTIPATH - the signal reaches the receiver by several paths
 *      of different length, so delayed echoes add to the direct
 *      signal. Mathematically the channel is an FIR filter; its taps
 *      are the path gains. Multipath smears symbols into each other
 *      (inter-symbol interference) - exactly what an OFDM cyclic
 *      prefix and an equalizer are designed to undo.
 *
 *   2. NOISE - thermal noise adds a random complex Gaussian value to
 *      every sample (additive white Gaussian noise, AWGN).
 *
 *   channel output = (input convolved with multipath taps) + noise
 *
 * This module provides a reusable complex channel: set the multipath
 * taps and a noise level, then push samples through. It is the test
 * harness for the OFDM modem and pairs with the BER measurements
 * already used for LDPC.
 */
#ifndef DSP_CHANNEL_H
#define DSP_CHANNEL_H

#include "../common.h"

/*
 * A complex multipath + AWGN channel.
 *   taps      : complex FIR taps (path gains); taps[0] is the direct
 *               path. A single tap {1+0j} is a pure-AWGN channel.
 *   ntaps     : number of taps
 *   noise_std : standard deviation of the complex Gaussian noise
 *               (split evenly across the real and imaginary parts)
 *   rng_state : internal RNG state - do not modify directly
 */
typedef struct {
    cplx    *taps;
    size_t   ntaps;
    double   noise_std;
    unsigned rng_state;
} dsp_channel;

/*
 * Initialise a channel.
 *   taps      : array of ntaps complex path gains (copied internally)
 *   ntaps     : number of taps (>= 1)
 *   noise_std : complex AWGN standard deviation (0 for a noiseless
 *               channel)
 *   seed      : RNG seed for reproducible noise
 * Returns 0 on success, -1 on bad parameters or allocation failure.
 * Pair with dsp_channel_free.
 */
int dsp_channel_init(dsp_channel *ch, const cplx *taps, size_t ntaps,
                     double noise_std, unsigned seed);

/*
 * Convenience: initialise a pure-AWGN channel (single unit tap, no
 * multipath). Equivalent to dsp_channel_init with taps = {1+0j}.
 */
int dsp_channel_init_awgn(dsp_channel *ch, double noise_std,
                          unsigned seed);

/* Release the channel's internal buffers. */
void dsp_channel_free(dsp_channel *ch);

/*
 * Pass a block of complex samples through the channel.
 *   in  : transmitted samples, length n
 *   out : received samples, length n (may not alias in)
 * The multipath filter uses zero initial state, then AWGN is added.
 * To model continuous streaming, process the whole stream in one call.
 */
void dsp_channel_apply(dsp_channel *ch, const cplx *in, cplx *out,
                       size_t n);

/*
 * Frequency response of the multipath filter at `nfft` equally spaced
 * points (the DFT of the zero-padded tap vector). This is exactly the
 * per-subcarrier gain an OFDM receiver must equalize.
 *   resp : output array, length nfft
 */
void dsp_channel_frequency_response(const dsp_channel *ch,
                                    cplx *resp, size_t nfft);

#endif /* DSP_CHANNEL_H */
