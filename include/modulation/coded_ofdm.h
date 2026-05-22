/*
 * coded_ofdm.h - Coded OFDM: the full transmit/receive chain
 *
 * PROBLEM SOLVED
 *   OFDM alone (ofdm.h) handles multipath, but it does not correct
 *   bit errors - a subcarrier in a deep fade is simply lost. Real
 *   systems therefore wrap OFDM around a forward-error-correction
 *   code: this is CODED OFDM, the scheme used by Wi-Fi, LTE, and 5G.
 *
 *   The two halves of this guide finally meet here:
 *     - the coding/ modules (LDPC, Reed-Solomon) supply error
 *       correction,
 *     - the modulation/ modules (QAM, OFDM, channel) supply the
 *       physical link.
 *
 * THE CHAIN
 *   transmit:  message bits
 *                -> FEC encode        (add redundancy)
 *                -> QAM map           (bits -> subcarrier symbols)
 *                -> OFDM modulate     (IFFT + cyclic prefix)
 *                -> channel           (multipath + noise)
 *   receive:   received samples
 *                -> OFDM demodulate   (remove prefix + FFT)
 *                -> equalize          (one divide per subcarrier)
 *                -> QAM soft-demap    (symbols -> per-bit LLRs)
 *                -> FEC decode        (correct the残 errors)
 *
 *   The FEC and the OFDM physical layer each handle what the other
 *   cannot: OFDM tames multipath, the code cleans up the residual
 *   noise errors - including the bits carried on faded subcarriers.
 *
 * WHY SOFT DECODING MATTERS HERE
 *   The QAM demapper produces a log-likelihood ratio per bit, and the
 *   equalizer's channel knowledge scales those LLRs - a bit from a
 *   weak subcarrier arrives with low confidence. Passing that soft
 *   information into the LDPC decoder is what makes coded OFDM work
 *   close to the channel capacity.
 *
 * SCOPE
 *   This module wires the existing pieces into one transceiver and a
 *   BER measurement. It uses the all-zero codeword convention (as the
 *   LDPC BER sweep does) so no separate systematic encoder is needed:
 *   the focus is the end-to-end signal chain, not key generation.
 */
#ifndef DSP_CODED_OFDM_H
#define DSP_CODED_OFDM_H

#include "../common.h"
#include "modulation/qam.h"
#include "modulation/ofdm.h"
#include "modulation/channel.h"
#include "coding/ldpc.h"

/*
 * Coded-OFDM transceiver configuration.
 *   ofdm  : OFDM parameters (FFT size, cyclic prefix)
 *   order : QAM constellation for every subcarrier
 *   code  : LDPC code supplying forward error correction
 */
typedef struct {
    dsp_ofdm        ofdm;
    dsp_qam_order   order;
    const dsp_ldpc *code;
} dsp_coded_ofdm;

/*
 * Run one full coded-OFDM frame through the chain and report errors.
 *
 * The transmitted codeword is the all-zero LDPC codeword (always
 * valid). It is QAM-mapped, OFDM-modulated, passed through `ch`,
 * demodulated, equalized against the known channel response, soft-
 * demapped to LLRs, and LDPC-decoded. The function compares the
 * decoded bits with the transmitted all-zero word.
 *
 *   cfg       : transceiver configuration
 *   ch        : channel to transmit through
 *   max_iter  : LDPC decoder iteration cap
 *   bit_errors: if non-NULL, receives the number of wrong bits
 *   raw_errors: if non-NULL, receives the pre-decoding bit-error count
 *               (errors at the QAM demapper, before FEC) - useful for
 *               showing how much work the code did
 * Returns 0 on success, -1 on a configuration or allocation error.
 *
 * Requires code->n == ofdm.nfft * bits_per_symbol(order), i.e. one
 * full codeword fills exactly one OFDM symbol.
 */
int dsp_coded_ofdm_run_frame(const dsp_coded_ofdm *cfg,
                             dsp_channel *ch, int max_iter,
                             size_t *bit_errors, size_t *raw_errors);

/*
 * Measure the coded-OFDM bit-error rate over many frames.
 *   cfg      : transceiver configuration
 *   ch       : channel (its RNG advances across frames)
 *   frames   : number of OFDM frames to simulate
 *   max_iter : LDPC iteration cap per frame
 * Returns the post-decoding bit-error rate in [0, 1], or -1.0 on error.
 */
double dsp_coded_ofdm_ber(const dsp_coded_ofdm *cfg, dsp_channel *ch,
                          int frames, int max_iter);

#endif /* DSP_CODED_OFDM_H */
