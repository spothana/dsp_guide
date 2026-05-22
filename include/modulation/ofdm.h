/*
 * ofdm.h - Orthogonal Frequency-Division Multiplexing
 *
 * PROBLEM SOLVED
 *   A single high-rate carrier sent through a multipath channel
 *   smears symbols together (inter-symbol interference), and undoing
 *   that needs a long, expensive equalizer. OFDM sidesteps the
 *   problem: it splits the data across many SLOW, narrow subcarriers
 *   sent in parallel. Each subcarrier is narrow enough that the
 *   channel looks flat across it - so equalization becomes one
 *   complex multiply per subcarrier.
 *
 * THE ELEGANT PART - it is just an FFT
 *   Modulator : load one QAM symbol onto each subcarrier, then take
 *               an IFFT. The IFFT *is* the OFDM modulator - it
 *               produces the time-domain transmit signal in a single
 *               transform.
 *   Demodulator: take an FFT of the received block and read each
 *               subcarrier's symbol straight back out.
 *
 * THE CYCLIC PREFIX
 *   Before transmitting, copy the LAST cp_len samples of each IFFT
 *   block to its FRONT. This guard interval absorbs the multipath
 *   echo, and - the key trick - it makes the channel's LINEAR
 *   convolution look like a CIRCULAR convolution over the block. A
 *   circular convolution is a pointwise product in the FFT domain, so
 *   each subcarrier just picks up a single complex channel gain.
 *   Equalization is then one divide per subcarrier (see _equalize).
 *
 *   transmit:  QAM symbols -> IFFT -> add cyclic prefix -> channel
 *   receive:   channel -> remove cyclic prefix -> FFT -> equalize
 *
 * USED IN
 *   Wi-Fi, 4G LTE, 5G, DVB-T, ADSL - OFDM is the dominant modulation
 *   of modern broadband wireless and wireline systems.
 *
 * CONSTRAINT
 *   The number of subcarriers nfft must be a power of two (radix-2
 *   FFT). The cyclic prefix must be at least as long as the channel's
 *   multipath span for the circular-convolution trick to hold exactly.
 */
#ifndef DSP_OFDM_H
#define DSP_OFDM_H

#include "../common.h"

/*
 * OFDM system parameters.
 *   nfft   : number of subcarriers / FFT size (power of two)
 *   cp_len : cyclic-prefix length in samples (< nfft)
 */
typedef struct {
    size_t nfft;
    size_t cp_len;
} dsp_ofdm;

/* Initialise OFDM parameters. Returns 0 on success, -1 if nfft is not
 * a power of two or cp_len >= nfft. */
int dsp_ofdm_init(dsp_ofdm *o, size_t nfft, size_t cp_len);

/* Length of one transmitted OFDM symbol: nfft + cp_len samples. */
static inline size_t dsp_ofdm_symbol_len(const dsp_ofdm *o) {
    return o->nfft + o->cp_len;
}

/*
 * Modulate one OFDM symbol.
 *   freq : nfft QAM symbols, one per subcarrier
 *   time : output time-domain samples, length nfft + cp_len
 * Performs an IFFT, then prepends the cyclic prefix.
 * Returns the number of output samples, or 0 on error.
 */
size_t dsp_ofdm_modulate(const dsp_ofdm *o, const cplx *freq,
                         cplx *time);

/*
 * Demodulate one received OFDM symbol.
 *   time : received samples, length nfft + cp_len
 *   freq : output subcarrier symbols, length nfft
 * Removes the cyclic prefix, then performs an FFT.
 * Returns the number of subcarriers, or 0 on error.
 */
size_t dsp_ofdm_demodulate(const dsp_ofdm *o, const cplx *time,
                           cplx *freq);

/*
 * Per-subcarrier zero-forcing equalization.
 *   freq    : demodulated subcarrier symbols, length nfft - corrected
 *             in place
 *   chan_fr : channel frequency response, length nfft (the FFT of the
 *             channel taps - see dsp_channel_frequency_response)
 * Divides each subcarrier by its channel gain, undoing the multipath
 * distortion. This single complex divide per subcarrier replaces the
 * long time-domain equalizer a single-carrier system would need.
 */
void dsp_ofdm_equalize(const dsp_ofdm *o, cplx *freq,
                       const cplx *chan_fr);

#endif /* DSP_OFDM_H */
