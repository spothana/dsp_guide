/*
 * test_all.c - Unit tests for the DSP Guide.
 *
 * A tiny assert-and-count harness. Each test checks a numerical
 * property that must hold for a correct implementation. Exit status
 * is 0 when every test passes, 1 otherwise (so ctest reports it).
 */
#include "dsp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int tests_run = 0;
static int tests_failed = 0;

#define CHECK(cond, msg)                                        \
    do {                                                        \
        ++tests_run;                                            \
        if (!(cond)) {                                          \
            ++tests_failed;                                     \
            printf("  FAIL: %s\n", (msg));                      \
        } else {                                                \
            printf("  ok  : %s\n", (msg));                      \
        }                                                       \
    } while (0)

/* Approximate equality for floating-point comparisons. */
static int close(double a, double b, double tol) {
    double d = a - b;
    return (d < 0 ? -d : d) <= tol;
}

/* ---- transforms ----------------------------------------------------- */

static void test_dft_fft_agree(void) {
    printf("[transforms] DFT and FFT produce the same spectrum\n");
    enum { N = 16 };
    cplx sig[N], a[N], b[N];
    for (int i = 0; i < N; ++i)
        sig[i] = dsp_cplx(sin(i) + 0.3 * i, 0.0);

    memcpy(a, sig, sizeof sig);
    memcpy(b, sig, sizeof sig);
    dsp_dft(a, a, N);          /* note: dft tolerates in==out here */
    dsp_fft(b, N);

    int match = 1;
    for (int k = 0; k < N; ++k)
        if (!close(dsp_mag(a[k]), dsp_mag(b[k]), 1e-9))
            match = 0;
    CHECK(match, "FFT magnitude == DFT magnitude for all bins");
}

static void test_fft_roundtrip(void) {
    printf("[transforms] inverse FFT recovers the original signal\n");
    enum { N = 32 };
    cplx orig[N], buf[N];
    for (int i = 0; i < N; ++i)
        orig[i] = dsp_cplx((double)(i % 7) - 3.0, 0.0);
    memcpy(buf, orig, sizeof orig);

    dsp_fft(buf, N);
    dsp_ifft(buf, N);

    int ok = 1;
    for (int i = 0; i < N; ++i)
        if (!close(creal(buf[i]), creal(orig[i]), 1e-9))
            ok = 0;
    CHECK(ok, "ifft(fft(x)) == x");
}

static void test_fft_rejects_non_pow2(void) {
    printf("[transforms] radix-2 FFT rejects non-power-of-two lengths\n");
    cplx buf[6] = {0};
    CHECK(dsp_fft(buf, 6) == -1, "fft on length 6 returns error");
    CHECK(dsp_fft(buf, 4) == 0,  "fft on length 4 succeeds");
}

static void test_dct_roundtrip(void) {
    printf("[transforms] inverse DCT recovers the original signal\n");
    enum { N = 16 };
    double orig[N], buf[N], back[N];
    for (int i = 0; i < N; ++i) orig[i] = 1.0 + 0.5 * i - 0.1 * i * i;

    dsp_dct(orig, buf, N);
    dsp_idct(buf, back, N);

    int ok = 1;
    for (int i = 0; i < N; ++i)
        if (!close(orig[i], back[i], 1e-9)) ok = 0;
    CHECK(ok, "idct(dct(x)) == x");
}

static void test_dct_energy_compaction(void) {
    printf("[transforms] DCT compacts a smooth signal's energy\n");
    enum { N = 32 };
    double ramp[N], c[N];
    for (int i = 0; i < N; ++i) ramp[i] = (double)i;
    dsp_dct(ramp, c, N);

    /* Energy in the first 4 coefficients vs the total. */
    double head = 0.0, total = 0.0;
    for (int k = 0; k < N; ++k) {
        double e = c[k] * c[k];
        total += e;
        if (k < 4) head += e;
    }
    CHECK(head / total > 0.99,
          "first 4 of 32 DCT coefficients hold >99% of the energy");
}

/* ---- filtering ------------------------------------------------------ */

static void test_fir_linear_phase(void) {
    printf("[filtering] windowed-sinc FIR has symmetric coefficients\n");
    enum { TAPS = 31 };
    double taps[TAPS];
    dsp_fir_design_lowpass(taps, TAPS, 0.2);

    int sym = 1;
    for (int k = 0; k < TAPS; ++k)
        if (!close(taps[k], taps[TAPS - 1 - k], 1e-12))
            sym = 0;
    CHECK(sym, "taps[k] == taps[M-k] -> exactly linear phase");
}

static void test_fir_dc_gain(void) {
    printf("[filtering] FIR low-pass has unity DC gain\n");
    enum { TAPS = 31 };
    double taps[TAPS];
    dsp_fir_design_lowpass(taps, TAPS, 0.2);
    double sum = 0.0;
    for (int k = 0; k < TAPS; ++k) sum += taps[k];
    CHECK(close(sum, 1.0, 1e-9), "sum of taps == 1");
}

static void test_fir_attenuates_highfreq(void) {
    printf("[filtering] FIR low-pass attenuates a high-frequency tone\n");
    enum { N = 128, TAPS = 41 };
    double taps[TAPS], x[N], y[N];
    dsp_fir_design_lowpass(taps, TAPS, 0.10);
    for (int n = 0; n < N; ++n)
        x[n] = sin(2.0 * M_PI * 0.40 * n);   /* well above cutoff */
    dsp_fir_apply(x, y, N, taps, TAPS);

    /* Compare steady-state amplitude (skip the filter's transient). */
    double in_pk = 0.0, out_pk = 0.0;
    for (int n = TAPS; n < N; ++n) {
        if (fabs(x[n]) > in_pk)  in_pk  = fabs(x[n]);
        if (fabs(y[n]) > out_pk) out_pk = fabs(y[n]);
    }
    CHECK(out_pk < 0.25 * in_pk,
          "0.40 tone is attenuated to <25% amplitude");
}

static void test_iir_stability(void) {
    printf("[filtering] IIR stability test flags pole locations\n");
    dsp_biquad good;
    dsp_iir_design_lowpass(&good, 0.15, 0.707);
    CHECK(dsp_iir_is_stable(&good), "designed Butterworth biquad is stable");

    /* Hand-crafted unstable section: a pole outside the unit circle. */
    dsp_biquad bad = { 1.0, 0.0, 0.0, /*a1*/ 0.0, /*a2*/ -1.5 };
    CHECK(!dsp_iir_is_stable(&bad), "biquad with |a2|>1 flagged unstable");
}

static void test_iir_attenuates_highfreq(void) {
    printf("[filtering] IIR low-pass attenuates a high-frequency tone\n");
    enum { N = 256 };
    double x[N], y[N];
    dsp_biquad bq;
    dsp_iir_design_lowpass(&bq, 0.05, 0.707);
    for (int n = 0; n < N; ++n)
        x[n] = sin(2.0 * M_PI * 0.35 * n);
    dsp_iir_apply(&bq, x, y, N);

    double in_pk = 0.0, out_pk = 0.0;
    for (int n = N / 2; n < N; ++n) {       /* steady state */
        if (fabs(x[n]) > in_pk)  in_pk  = fabs(x[n]);
        if (fabs(y[n]) > out_pk) out_pk = fabs(y[n]);
    }
    CHECK(out_pk < 0.3 * in_pk,
          "0.35 tone attenuated to <30% amplitude");
}

/* ---- operations ----------------------------------------------------- */

static void test_convolution_methods_agree(void) {
    printf("[operations] direct and FFT convolution agree\n");
    double x[10], h[5], yd[14], yf[14];
    for (int i = 0; i < 10; ++i) x[i] = i - 4.0;
    for (int i = 0; i < 5;  ++i) h[i] = (i == 2) ? 1.0 : 0.25;

    dsp_convolve(x, 10, h, 5, yd);
    dsp_convolve_fft(x, 10, h, 5, yf);

    int ok = 1;
    for (int i = 0; i < 14; ++i)
        if (!close(yd[i], yf[i], 1e-9)) ok = 0;
    CHECK(ok, "direct convolution == FFT convolution");
}

static void test_convolution_identity(void) {
    printf("[operations] convolving with a unit impulse is identity\n");
    double x[5] = {2, -1, 3, 0, 4};
    double imp[1] = {1.0};
    double y[5];
    dsp_convolve(x, 5, imp, 1, y);
    int ok = 1;
    for (int i = 0; i < 5; ++i)
        if (!close(x[i], y[i], 1e-12)) ok = 0;
    CHECK(ok, "x convolved with delta == x");
}

static void test_correlation_delay(void) {
    printf("[operations] cross-correlation finds a known delay\n");
    enum { N = 32 };
    double a[N], b[N];
    for (int i = 0; i < N; ++i) {
        a[i] = (i == 6)  ? 1.0 : 0.0;
        b[i] = (i == 13) ? 1.0 : 0.0;        /* delayed by 7 */
    }
    CHECK(dsp_estimate_delay(a, N, b, N) == 7,
          "estimated lag == 7");
}

static void test_autocorrelation_peak(void) {
    printf("[operations] auto-correlation peaks at zero lag\n");
    enum { N = 20 };
    double x[N], r[N];
    for (int i = 0; i < N; ++i) x[i] = sin(0.5 * i);
    dsp_auto_correlate(x, N, r);
    int ok = 1;
    for (int k = 1; k < N; ++k)
        if (r[k] > r[0]) ok = 0;
    CHECK(ok, "R[0] is the maximum (signal energy)");
}

/* ---- spectral ------------------------------------------------------- */

static void test_window_endpoints(void) {
    printf("[spectral] tapered windows go to (near) zero at the edges\n");
    enum { N = 32 };
    double w[N];

    dsp_window_hanning(w, N);
    CHECK(close(w[0], 0.0, 1e-9) && close(w[N-1], 0.0, 1e-9),
          "Hanning window endpoints are 0");

    dsp_window_rectangular(w, N);
    CHECK(close(w[0], 1.0, 1e-12) && close(w[N-1], 1.0, 1e-12),
          "Rectangular window endpoints are 1 (no taper)");
}

static void test_window_symmetry(void) {
    printf("[spectral] window functions are symmetric\n");
    enum { N = 31 };
    double w[N];
    dsp_window_blackman(w, N);
    int sym = 1;
    for (int k = 0; k < N; ++k)
        if (!close(w[k], w[N-1-k], 1e-12)) sym = 0;
    CHECK(sym, "Blackman window is symmetric about its centre");
}

static void test_window_leakage_order(void) {
    printf("[spectral] tapering reduces total edge discontinuity\n");
    enum { N = 64 };
    double rect[N], hann[N], black[N];
    dsp_window_rectangular(rect, N);
    dsp_window_hanning(hann, N);
    dsp_window_blackman(black, N);
    /* Edge value is a proxy for leakage: lower edge -> less leakage. */
    CHECK(hann[0] < rect[0] && black[0] < rect[0],
          "Hanning and Blackman start lower than Rectangular");
}

/* ---- sampling ------------------------------------------------------- */

static void test_decimate_length(void) {
    printf("[sampling] decimation produces the expected length\n");
    enum { N = 100 };
    double x[N], y[N];
    for (int i = 0; i < N; ++i) x[i] = sin(0.1 * i);
    size_t out = dsp_decimate(x, N, 4, y);
    CHECK(out == dsp_decimate_len(N, 4), "decimate-by-4 length matches");
}

static void test_interpolate_length_and_zeros(void) {
    printf("[sampling] interpolation expands length by the factor\n");
    enum { N = 16, L = 3 };
    double x[N];
    double *y = malloc(N * L * sizeof(double));
    for (int i = 0; i < N; ++i) x[i] = 1.0;
    size_t out = dsp_interpolate(x, N, L, y);
    CHECK(out == (size_t)N * L, "interpolate-by-3 length == 3N");
    free(y);
}

static void test_resample_ratio(void) {
    printf("[sampling] rational resampling applies the L/M ratio\n");
    enum { N = 48 };
    double x[N];
    double *y = malloc(dsp_resample_len(N, 3, 2) * sizeof(double));
    for (int i = 0; i < N; ++i) x[i] = sin(0.05 * i);
    size_t out = dsp_resample(x, N, 3, 2, y);
    /* 48 -> *3 = 144 -> /2 ~= 72 */
    CHECK(out == dsp_resample_len(N, 3, 2) && out > N,
          "resample by 3/2 lengthens the signal");
    free(y);
}

/* ---- wavelet -------------------------------------------------------- */

static void test_dwt_roundtrip(void) {
    printf("[wavelet] inverse DWT reconstructs the original signal\n");
    enum { N = 16 };
    double orig[N], buf[N];
    for (int i = 0; i < N; ++i) orig[i] = (i * 7 % 11) - 5.0;
    memcpy(buf, orig, sizeof orig);

    int levels = dsp_dwt_haar(buf, N);
    dsp_idwt_haar(buf, N);

    int ok = (levels == 4);                  /* log2(16) */
    for (int i = 0; i < N; ++i)
        if (!close(buf[i], orig[i], 1e-9)) ok = 0;
    CHECK(ok, "idwt(dwt(x)) == x and level count == log2(N)");
}

static void test_dwt_rejects_non_pow2(void) {
    printf("[wavelet] DWT rejects non-power-of-two lengths\n");
    double buf[12] = {0};
    CHECK(dsp_dwt_haar(buf, 12) == -1, "dwt on length 12 returns error");
}

static void test_dwt_constant_signal(void) {
    printf("[wavelet] DWT of a constant signal has zero detail\n");
    enum { N = 8 };
    double buf[N];
    for (int i = 0; i < N; ++i) buf[i] = 5.0;
    dsp_dwt_haar(buf, N);
    /* All detail coefficients (indices 1..N-1) should vanish. */
    int ok = 1;
    for (int i = 1; i < N; ++i)
        if (!close(buf[i], 0.0, 1e-9)) ok = 0;
    CHECK(ok, "constant signal -> all detail coefficients are 0");
}

/* ---- coding: error detection ---------------------------------------- */

static void test_parity_detects_odd_flips(void) {
    printf("[coding] parity detects an odd number of bit errors\n");
    uint8_t data[4] = { 0xA5, 0x3C, 0xFF, 0x00 };
    int p = dsp_parity_compute(data, 4, DSP_PARITY_EVEN);
    CHECK(dsp_parity_check(data, 4, DSP_PARITY_EVEN, p),
          "correct parity bit passes the check");

    uint8_t flipped[4] = { 0xA5, 0x3C, 0xFF, 0x01 };  /* one bit flipped */
    CHECK(!dsp_parity_check(flipped, 4, DSP_PARITY_EVEN, p),
          "single-bit error fails the parity check");
}

static void test_checksum_roundtrip(void) {
    printf("[coding] Internet checksum verifies intact data\n");
    uint8_t data[10] = { 1,2,3,4,5,6,7,8,9,10 };
    uint16_t c = dsp_checksum16(data, 10);
    CHECK(dsp_checksum16_verify(data, 10, c),
          "intact buffer verifies against its checksum");

    data[4] ^= 0x80;
    CHECK(!dsp_checksum16_verify(data, 10, c),
          "corrupted buffer fails verification");
}

static void test_crc32_known_vector(void) {
    printf("[coding] CRC-32 matches the standard test vector\n");
    /* CRC-32 of the ASCII string "123456789" is 0xCBF43926. */
    const uint8_t v[9] = { '1','2','3','4','5','6','7','8','9' };
    CHECK(dsp_crc32(v, 9) == 0xCBF43926u,
          "CRC-32(\"123456789\") == 0xCBF43926");
}

static void test_crc32_detects_burst(void) {
    printf("[coding] CRC-32 detects a burst error\n");
    uint8_t data[16];
    for (int i = 0; i < 16; ++i) data[i] = (uint8_t)(i * 17 + 3);
    uint32_t good = dsp_crc32(data, 16);
    data[7] ^= 0xFF;            /* an 8-bit burst */
    CHECK(dsp_crc32(data, 16) != good,
          "burst error changes the CRC");
}

/* ---- coding: forward error correction ------------------------------- */

static void test_hamming_corrects_single_bit(void) {
    printf("[coding] Hamming(7,4) corrects any single-bit error\n");
    int all_ok = 1;
    for (uint8_t nib = 0; nib < 16; ++nib) {
        uint8_t code = dsp_hamming74_encode(nib);
        /* Flip each of the 7 bit positions in turn. */
        for (int b = 0; b < 7; ++b) {
            uint8_t bad = code ^ (uint8_t)(1u << b);
            uint8_t got = dsp_hamming74_decode(bad, NULL);
            if (got != nib) all_ok = 0;
        }
    }
    CHECK(all_ok, "every nibble recovers from every single-bit flip");
}

static void test_hamming_clean_syndrome(void) {
    printf("[coding] Hamming syndrome is zero for a clean codeword\n");
    int all_ok = 1;
    for (uint8_t nib = 0; nib < 16; ++nib)
        if (dsp_hamming74_syndrome(dsp_hamming74_encode(nib)) != 0)
            all_ok = 0;
    CHECK(all_ok, "uncorrupted codewords have syndrome 0");
}

static void test_rs_corrects_burst(void) {
    printf("[coding] Reed-Solomon corrects a burst of symbol errors\n");
    dsp_rs rs;
    int rc = dsp_rs_init(&rs, 8);          /* 8 parity -> t = 4 */
    CHECK(rc == 0, "RS codec initialises with 8 parity symbols");

    uint8_t data[12];
    for (int i = 0; i < 12; ++i) data[i] = (uint8_t)(i * 11 + 7);
    uint8_t parity[8];
    dsp_rs_encode(&rs, data, 12, parity);

    uint8_t cw[20];
    for (int i = 0; i < 12; ++i) cw[i] = data[i];
    for (int i = 0; i < 8;  ++i) cw[12 + i] = parity[i];

    /* Corrupt 4 symbols - exactly the correction limit. */
    cw[3] ^= 0x9D; cw[4] ^= 0x12; cw[5] ^= 0xC4; cw[6] ^= 0x77;

    int nfix = dsp_rs_decode(&rs, cw, 20);
    int ok = (nfix == 4);
    for (int i = 0; i < 12; ++i) if (cw[i] != data[i]) ok = 0;
    CHECK(ok, "4 corrupted symbols are located and repaired");
}

static void test_rs_clean_codeword(void) {
    printf("[coding] Reed-Solomon reports zero errors on clean data\n");
    dsp_rs rs;
    dsp_rs_init(&rs, 4);
    uint8_t data[8] = { 5,5,5,5,5,5,5,5 };
    uint8_t parity[4];
    dsp_rs_encode(&rs, data, 8, parity);
    uint8_t cw[12];
    for (int i = 0; i < 8; ++i) cw[i] = data[i];
    for (int i = 0; i < 4; ++i) cw[8 + i] = parity[i];
    CHECK(dsp_rs_decode(&rs, cw, 12) == 0,
          "an untouched codeword decodes with 0 corrections");
}

static void test_viterbi_hard_decoding(void) {
    printf("[coding] Viterbi (hard) corrects errors in a coded stream\n");
    uint8_t bits[12] = { 1,1,0,1,0,0,1,0,1,1,0,0 };
    uint8_t enc[2 * (12 + 2)];
    size_t elen = dsp_conv_encode(bits, 12, enc);

    uint8_t rx[2 * (12 + 2)];
    memcpy(rx, enc, elen);
    rx[3] ^= 1; rx[10] ^= 1;               /* two well-separated errors */

    uint8_t dec[12];
    size_t dn = dsp_viterbi_decode(rx, elen, dec);
    int ok = (dn == 12);
    for (size_t i = 0; i < dn; ++i) if (dec[i] != bits[i]) ok = 0;
    CHECK(ok, "hard-decision Viterbi recovers the original 12 bits");
}

static void test_viterbi_soft_decoding(void) {
    printf("[coding] Viterbi (soft) decodes from analog confidences\n");
    uint8_t bits[10] = { 0,1,1,0,1,0,0,1,1,0 };
    uint8_t enc[2 * (10 + 2)];
    size_t elen = dsp_conv_encode(bits, 10, enc);

    double soft[2 * (10 + 2)];
    for (size_t i = 0; i < elen; ++i)
        soft[i] = enc[i] ? -1.0 : 1.0;
    /* Two near-zero (ambiguous) samples instead of hard flips. */
    soft[4] = 0.05; soft[11] = -0.05;

    uint8_t dec[10];
    size_t dn = dsp_viterbi_decode_soft(soft, elen, dec);
    int ok = (dn == 10);
    for (size_t i = 0; i < dn; ++i) if (dec[i] != bits[i]) ok = 0;
    CHECK(ok, "soft-decision Viterbi recovers the original 10 bits");
}

static void test_conv_encoded_length(void) {
    printf("[coding] convolutional encoder output length is 2(n+2)\n");
    uint8_t bits[5] = { 1,0,1,0,1 };
    uint8_t enc[2 * (5 + 2)];
    size_t elen = dsp_conv_encode(bits, 5, enc);
    CHECK(elen == dsp_conv_encoded_len(5),
          "rate-1/2 code with 2 tail bits -> length 2(n+2)");
}

/* ---- coding: channel equalization ----------------------------------- */

static void test_lms_converges(void) {
    printf("[coding] LMS equalizer reduces error as it converges\n");
    enum { N = 600 };
    double tx[N], rx[N];
    unsigned seed = 999;
    for (int i = 0; i < N; ++i) {
        seed = seed * 1103515245u + 12345u;
        tx[i] = ((seed >> 16) & 1) ? 1.0 : -1.0;
    }
    /* 2-tap multipath channel. */
    rx[0] = tx[0];
    for (int i = 1; i < N; ++i)
        rx[i] = tx[i] + 0.5 * tx[i - 1];

    dsp_lms eq;
    int rc = dsp_lms_init(&eq, 11, 0.02);
    CHECK(rc == 0, "LMS equalizer initialises");

    double mse_early = dsp_lms_train(&eq, rx, tx, 100);
    dsp_lms_train(&eq, rx + 100, tx + 100, 400);
    double mse_late = dsp_lms_train(&eq, rx + 500, tx + 500, 100);

    CHECK(mse_late < mse_early,
          "mean-squared error decreases as the taps adapt");
    dsp_lms_free(&eq);
}

int main(void) {
    printf("DSP STUDY GUIDE - test suite\n");
    printf("============================\n\n");

    test_dft_fft_agree();
    test_fft_roundtrip();
    test_fft_rejects_non_pow2();
    test_dct_roundtrip();
    test_dct_energy_compaction();

    test_fir_linear_phase();
    test_fir_dc_gain();
    test_fir_attenuates_highfreq();
    test_iir_stability();
    test_iir_attenuates_highfreq();

    test_convolution_methods_agree();
    test_convolution_identity();
    test_correlation_delay();
    test_autocorrelation_peak();

    test_window_endpoints();
    test_window_symmetry();
    test_window_leakage_order();

    test_decimate_length();
    test_interpolate_length_and_zeros();
    test_resample_ratio();

    test_dwt_roundtrip();
    test_dwt_rejects_non_pow2();
    test_dwt_constant_signal();

    test_parity_detects_odd_flips();
    test_checksum_roundtrip();
    test_crc32_known_vector();
    test_crc32_detects_burst();

    test_hamming_corrects_single_bit();
    test_hamming_clean_syndrome();
    test_rs_corrects_burst();
    test_rs_clean_codeword();
    test_viterbi_hard_decoding();
    test_viterbi_soft_decoding();
    test_conv_encoded_length();

    test_lms_converges();

    printf("\n============================\n");
    printf("Tests run: %d   Passed: %d   Failed: %d\n",
           tests_run, tests_run - tests_failed, tests_failed);
    return tests_failed == 0 ? 0 : 1;
}
