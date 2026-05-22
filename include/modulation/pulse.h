/*
 * pulse.h - Pulse shaping (root-raised-cosine filtering)
 *
 * PROBLEM SOLVED
 *   A symbol stream cannot be transmitted as bare impulses: sharp
 *   edges need infinite bandwidth, and rectangular pulses smear into
 *   neighbouring symbol slots (inter-symbol interference, ISI).
 *   PULSE SHAPING replaces each symbol with a smooth pulse that is
 *   bandlimited yet still ISI-free at the sampling instants.
 *
 *   The raised-cosine (RC) pulse has the Nyquist property: it is zero
 *   at every symbol instant except its own, so neighbouring symbols
 *   do not interfere when sampled correctly. In practice the RC
 *   response is split into a ROOT-raised-cosine (RRC) filter at the
 *   transmitter and an identical one at the receiver - cascaded they
 *   form the full RC response, and the receive RRC also acts as the
 *   matched filter that maximises SNR.
 *
 * THE ROLL-OFF FACTOR (beta, 0..1)
 *   beta = 0   : narrowest bandwidth, slowest-decaying pulse tails
 *   beta = 1   : widest bandwidth (2x), compact tails, easy timing
 *   Typical systems use 0.2-0.35 - the bandwidth/robustness trade-off.
 *
 * USAGE
 *   transmit: upsample symbols by sps, then filter with the RRC taps
 *   receive : filter the received samples with the same RRC taps,
 *             then sample at the symbol instants
 */
#ifndef DSP_PULSE_H
#define DSP_PULSE_H

#include "../common.h"

/*
 * Design a root-raised-cosine FIR filter.
 *   taps   : output buffer, length ntaps (odd recommended -> linear
 *            phase, exact centre tap)
 *   ntaps  : number of taps; span/sps symbols of pulse history
 *   sps    : samples per symbol (the oversampling factor)
 *   beta   : roll-off factor in [0, 1]
 * The taps are normalised to unit energy so a transmit/receive RRC
 * pair has unity gain.
 */
void dsp_rrc_design(double *taps, size_t ntaps, size_t sps, double beta);

/*
 * Pulse-shape a complex symbol stream for transmission.
 *   syms   : input symbols, length nsym
 *   sps    : samples per symbol
 *   taps   : RRC filter taps
 *   ntaps  : number of taps
 *   out    : output samples, length nsym * sps
 * Upsamples by inserting sps-1 zeros between symbols, then applies the
 * RRC filter. Returns the number of output samples.
 */
size_t dsp_pulse_shape(const cplx *syms, size_t nsym, size_t sps,
                       const double *taps, size_t ntaps, cplx *out);

/* Output length produced by dsp_pulse_shape. */
static inline size_t dsp_pulse_shaped_len(size_t nsym, size_t sps) {
    return nsym * sps;
}

/*
 * Matched-filter the received samples with the same RRC taps.
 *   in     : received samples, length n
 *   n      : sample count
 *   taps   : RRC filter taps (the same ones used to transmit)
 *   ntaps  : number of taps
 *   out    : filtered output, length n
 */
void dsp_matched_filter(const cplx *in, size_t n,
                        const double *taps, size_t ntaps, cplx *out);

#endif /* DSP_PULSE_H */
