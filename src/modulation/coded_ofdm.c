/*
 * coded_ofdm.c - End-to-end coded-OFDM transceiver.
 *
 * This is where the coding/ and modulation/ halves of the guide join.
 * One frame: encode -> QAM map -> OFDM -> channel -> OFDM^-1 ->
 * equalize -> QAM soft-demap -> LDPC decode.
 */
#include "modulation/coded_ofdm.h"
#include "transforms/fft.h"
#include <stdlib.h>

int dsp_coded_ofdm_run_frame(const dsp_coded_ofdm *cfg,
                             dsp_channel *ch, int max_iter,
                             size_t *bit_errors, size_t *raw_errors) {
    size_t nfft = cfg->ofdm.nfft;
    size_t bps  = dsp_qam_bits_per_symbol(cfg->order);
    if (bps == 0 || cfg->code == NULL)
        return -1;

    /* One LDPC codeword must fill exactly one OFDM symbol. */
    size_t nbits = nfft * bps;
    if (cfg->code->n != nbits)
        return -1;

    size_t sym_len = dsp_ofdm_symbol_len(&cfg->ofdm);

    uint8_t *tx_bits  = malloc(nbits * sizeof(uint8_t));
    cplx    *tx_freq  = malloc(nfft  * sizeof(cplx));
    cplx    *tx_time  = malloc(sym_len * sizeof(cplx));
    cplx    *rx_time  = malloc(sym_len * sizeof(cplx));
    cplx    *rx_freq  = malloc(nfft  * sizeof(cplx));
    cplx    *chan_fr  = malloc(nfft  * sizeof(cplx));
    double  *llr      = malloc(nbits * sizeof(double));
    uint8_t *rx_bits  = malloc(nbits * sizeof(uint8_t));
    uint8_t *raw_bits = malloc(nbits * sizeof(uint8_t));
    if (!tx_bits || !tx_freq || !tx_time || !rx_time || !rx_freq ||
        !chan_fr || !llr || !rx_bits || !raw_bits) {
        free(tx_bits); free(tx_freq); free(tx_time); free(rx_time);
        free(rx_freq); free(chan_fr); free(llr); free(rx_bits);
        free(raw_bits);
        return -1;
    }

    /* --- Transmit ---------------------------------------------------
     * The all-zero codeword is always valid, so we use it directly:
     * no separate systematic encoder is needed to exercise the chain. */
    for (size_t i = 0; i < nbits; ++i)
        tx_bits[i] = 0;

    /* Bits -> QAM symbols, one per subcarrier. */
    dsp_qam_modulate(cfg->order, tx_bits, nbits, tx_freq);

    /* QAM symbols -> OFDM time-domain samples (IFFT + cyclic prefix). */
    dsp_ofdm_modulate(&cfg->ofdm, tx_freq, tx_time);

    /* --- Channel ---------------------------------------------------- */
    dsp_channel_apply(ch, tx_time, rx_time, sym_len);

    /* --- Receive ----------------------------------------------------
     * OFDM demodulate: remove the cyclic prefix and FFT. */
    dsp_ofdm_demodulate(&cfg->ofdm, rx_time, rx_freq);

    /* Per-subcarrier equalization against the known channel response.
     * (A real receiver estimates this from pilots; here we use the
     * channel's exact frequency response.) */
    dsp_channel_frequency_response(ch, chan_fr, nfft);
    dsp_ofdm_equalize(&cfg->ofdm, rx_freq, chan_fr);

    /* Pre-decoding ("raw") errors: hard-demap before the FEC runs,
     * so the caller can see how many errors the code had to fix. */
    if (raw_errors) {
        dsp_qam_demodulate(cfg->order, rx_freq, nfft, raw_bits);
        size_t re = 0;
        for (size_t i = 0; i < nbits; ++i)
            if (raw_bits[i] != tx_bits[i]) ++re;
        *raw_errors = re;
    }

    /* Soft-demap to per-bit LLRs, scaled by the noise variance, then
     * hand the soft information to the LDPC decoder.
     *
     * The OFDM demodulator divides by sqrt(nfft), so the channel noise
     * seen at the subcarriers is scaled down by the same factor: its
     * effective variance is noise_var / nfft. The soft demapper must
     * use this post-demodulation variance for the LLRs to be correctly
     * scaled. */
    double noise_var = ch->noise_std * ch->noise_std;
    double eff_var   = noise_var / (double)nfft;
    if (eff_var < 1e-12) eff_var = 1e-12;
    dsp_qam_demodulate_soft(cfg->order, rx_freq, nfft, eff_var, llr);

    dsp_ldpc_decode_sumproduct(cfg->code, llr, rx_bits, max_iter);

    /* --- Tally post-decoding errors --------------------------------- */
    if (bit_errors) {
        size_t be = 0;
        for (size_t i = 0; i < nbits; ++i)
            if (rx_bits[i] != tx_bits[i]) ++be;
        *bit_errors = be;
    }

    free(tx_bits); free(tx_freq); free(tx_time); free(rx_time);
    free(rx_freq); free(chan_fr); free(llr); free(rx_bits);
    free(raw_bits);
    return 0;
}

double dsp_coded_ofdm_ber(const dsp_coded_ofdm *cfg, dsp_channel *ch,
                          int frames, int max_iter) {
    size_t bps = dsp_qam_bits_per_symbol(cfg->order);
    if (bps == 0 || cfg->code == NULL || frames <= 0)
        return -1.0;

    size_t nbits = cfg->ofdm.nfft * bps;
    size_t total = 0, errors = 0;

    for (int f = 0; f < frames; ++f) {
        size_t be = 0;
        if (dsp_coded_ofdm_run_frame(cfg, ch, max_iter, &be, NULL) != 0)
            return -1.0;
        errors += be;
        total  += nbits;
    }
    return total ? (double)errors / (double)total : 0.0;
}
