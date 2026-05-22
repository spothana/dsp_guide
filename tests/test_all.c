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

/* ---- coding: interleaving -------------------------------------------- */

static void test_block_interleave_roundtrip(void) {
    printf("[coding] block interleave then deinterleave is identity\n");
    enum { R = 5, C = 7, LEN = R * C };
    uint8_t in[LEN], mid[LEN], out[LEN];
    for (int i = 0; i < LEN; ++i) in[i] = (uint8_t)(i * 3 + 1);

    int rc1 = dsp_block_interleave(in, mid, LEN, R, C);
    int rc2 = dsp_block_deinterleave(mid, out, LEN, R, C);
    int ok = (rc1 == 0 && rc2 == 0);
    for (int i = 0; i < LEN; ++i) if (out[i] != in[i]) ok = 0;
    CHECK(ok, "deinterleave(interleave(x)) == x");
}

static void test_block_interleave_size_check(void) {
    printf("[coding] block interleaver rejects a bad matrix size\n");
    uint8_t in[10], out[10];
    /* 10 != 3 * 4, so this must fail. */
    CHECK(dsp_block_interleave(in, out, 10, 3, 4) == -1,
          "len != rows*cols is rejected");
}

static void test_block_interleave_spreads_burst(void) {
    printf("[coding] block interleaving spreads a burst across rows\n");
    enum { R = 4, C = 8, LEN = R * C };
    uint8_t in[LEN], mid[LEN], out[LEN];
    for (int i = 0; i < LEN; ++i) in[i] = (uint8_t)i;

    dsp_block_interleave(in, mid, LEN, R, C);
    /* A burst of R consecutive errors on the interleaved stream. */
    for (int i = 10; i < 10 + R; ++i) mid[i] ^= 0xFF;
    dsp_block_deinterleave(mid, out, LEN, R, C);

    /* Count errors per row; each row models one codeword. */
    int max_per_row = 0;
    for (int r = 0; r < R; ++r) {
        int e = 0;
        for (int c = 0; c < C; ++c)
            if (out[r * C + c] != in[r * C + c]) ++e;
        if (e > max_per_row) max_per_row = e;
    }
    CHECK(max_per_row <= 1,
          "an R-symbol burst leaves at most 1 error per row");
}

static void test_block_interleave_rescues_rs(void) {
    printf("[coding] interleaving lets RS survive an over-limit burst\n");
    dsp_rs rs;
    dsp_rs_init(&rs, 4);                        /* t = 2 */
    enum { NCW = 4, N = 15, KK = 11, TOTAL = NCW * N };

    uint8_t cw[TOTAL];
    for (int w = 0; w < NCW; ++w) {
        uint8_t msg[KK], par[4];
        for (int i = 0; i < KK; ++i) msg[i] = (uint8_t)(w * 13 + i * 2);
        dsp_rs_encode(&rs, msg, KK, par);
        for (int i = 0; i < KK; ++i) cw[w * N + i] = msg[i];
        for (int i = 0; i < 4;  ++i) cw[w * N + KK + i] = par[i];
    }

    /* Without interleaving: a 6-symbol burst inside one codeword. */
    uint8_t plain[TOTAL];
    for (int i = 0; i < TOTAL; ++i) plain[i] = cw[i];
    for (int i = 3; i < 9; ++i) plain[i] ^= 0x7C;
    CHECK(dsp_rs_decode(&rs, plain, N) < 0,
          "6-symbol burst alone exceeds RS t=2 (uncorrectable)");

    /* With interleaving: the same burst spread across all codewords.
     * The matrix has NCW rows (one per codeword) and N columns, so
     * reading out by column makes consecutive symbols come from
     * consecutive codewords - a burst is dealt round-robin. */
    uint8_t tx[TOTAL], rx[TOTAL];
    dsp_block_interleave(cw, tx, TOTAL, NCW, N);
    for (int i = 3; i < 9; ++i) tx[i] ^= 0x7C;
    dsp_block_deinterleave(tx, rx, TOTAL, NCW, N);

    int ok = 1;
    for (int w = 0; w < NCW; ++w) {
        if (dsp_rs_decode(&rs, rx + w * N, N) < 0) ok = 0;
        for (int i = 0; i < KK; ++i)
            if (rx[w * N + i] != cw[w * N + i]) ok = 0;
    }
    CHECK(ok, "interleaved: every RS codeword recovers from the burst");
}

static void test_conv_interleaver_roundtrip(void) {
    printf("[coding] convolutional interleaver pair preserves the stream\n");
    dsp_conv_interleaver itl, dtl;
    int rc1 = dsp_conv_interleaver_init(&itl, 4, 3);
    int rc2 = dsp_conv_deinterleaver_init(&dtl, 4, 3);
    CHECK(rc1 == 0 && rc2 == 0, "interleaver and deinterleaver init");

    enum { N = 64 };
    uint8_t in[N];
    for (int i = 0; i < N; ++i) in[i] = (uint8_t)(i * 5 + 9);

    /* Push through interleaver then deinterleaver; the pair imposes a
     * fixed total latency, so output i matches input i once primed. */
    size_t lat = dsp_conv_interleaver_latency(4, 3);
    uint8_t out[N];
    for (int i = 0; i < N; ++i) {
        uint8_t a = dsp_conv_interleave_step(&itl, in[i]);
        out[i] = dsp_conv_interleave_step(&dtl, a);
    }
    int ok = 1;
    for (size_t i = lat; i < N; ++i)
        if (out[i] != in[i - lat]) ok = 0;
    CHECK(ok, "deinterleaved stream equals the input, delayed by latency");

    dsp_conv_interleaver_free(&itl);
    dsp_conv_interleaver_free(&dtl);
}

/* ---- coding: LDPC --------------------------------------------------- */

static void test_ldpc_from_matrix(void) {
    printf("[coding] LDPC builds from a dense parity-check matrix\n");
    /* A small 3x6 parity-check matrix. */
    uint8_t H[3 * 6] = {
        1,1,0,1,0,0,
        0,1,1,0,1,0,
        1,0,1,0,0,1
    };
    dsp_ldpc code;
    int rc = dsp_ldpc_from_matrix(&code, H, 3, 6);
    CHECK(rc == 0, "code constructed from matrix");

    /* The all-zero word satisfies every parity check. */
    uint8_t zero[6] = {0};
    CHECK(dsp_ldpc_check(&code, zero) == 1,
          "all-zero word is a valid codeword");

    /* Degrees read back from H: variable 1 is in checks 0 and 1. */
    int deg_ok = (code.col_deg[1] == 2 && code.row_deg[0] == 3);
    CHECK(deg_ok, "node degrees match the matrix");

    dsp_ldpc_free(&code);
}

static void test_ldpc_syndrome(void) {
    printf("[coding] LDPC syndrome weight flags corruption\n");
    uint8_t H[3 * 6] = {
        1,1,0,1,0,0,
        0,1,1,0,1,0,
        1,0,1,0,0,1
    };
    dsp_ldpc code;
    dsp_ldpc_from_matrix(&code, H, 3, 6);

    uint8_t word[6] = {0};
    CHECK(dsp_ldpc_syndrome_weight(&code, word) == 0,
          "clean codeword has syndrome weight 0");

    word[0] ^= 1;     /* bit 0 is in checks 0 and 2 */
    CHECK(dsp_ldpc_syndrome_weight(&code, word) == 2,
          "flipping bit 0 fails its 2 checks");

    dsp_ldpc_free(&code);
}

static void test_ldpc_regular_construction(void) {
    printf("[coding] regular LDPC generator yields exact node degrees\n");
    dsp_ldpc code;
    int rc = dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1);
    CHECK(rc == 0, "regular (wc=3,wr=4) code constructed");

    int cols_ok = 1, rows_ok = 1;
    for (size_t j = 0; j < code.n; ++j)
        if (code.col_deg[j] != 3) cols_ok = 0;
    for (size_t i = 0; i < code.m; ++i)
        if (code.row_deg[i] != 4) rows_ok = 0;
    CHECK(cols_ok, "every variable node has degree 3");
    CHECK(rows_ok, "every check node has degree 4");

    dsp_ldpc_free(&code);
}

static void test_ldpc_regular_rejects_bad_params(void) {
    printf("[coding] regular LDPC generator rejects inconsistent sizes\n");
    dsp_ldpc code;
    /* n*wc = 12*3 = 36, m*wr = 7*4 = 28 -> mismatch, must fail. */
    CHECK(dsp_ldpc_make_regular(&code, 7, 12, 3, 4, 1) == -1,
          "n*wc != m*wr is rejected");
}

static void test_ldpc_bitflip_single_error(void) {
    printf("[coding] LDPC bit-flipping corrects a single-bit error\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1);
    uint8_t cw[12] = {0};

    int all_ok = 1;
    for (int b = 0; b < 12; ++b) {
        uint8_t r[12];
        for (int i = 0; i < 12; ++i) r[i] = cw[i];
        r[b] ^= 1;
        int it = dsp_ldpc_decode_bitflip(&code, r, 40);
        if (it < 0) { all_ok = 0; continue; }
        for (int i = 0; i < 12; ++i)
            if (r[i] != cw[i]) all_ok = 0;
    }
    CHECK(all_ok, "every single-bit error is corrected");

    dsp_ldpc_free(&code);
}

static void test_ldpc_sumproduct_corrects(void) {
    printf("[coding] LDPC sum-product corrects multiple soft errors\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1);
    uint8_t cw[12] = {0};

    /* Try every pair of wrong-sign bits; sum-product should fix all. */
    int all_ok = 1;
    for (int a = 0; a < 12 && all_ok; ++a) {
        for (int b = a + 1; b < 12; ++b) {
            double llr[12];
            for (int i = 0; i < 12; ++i) llr[i] = 3.0;
            llr[a] = -1.0;
            llr[b] = -1.0;
            uint8_t out[12];
            int it = dsp_ldpc_decode_sumproduct(&code, llr, out, 50);
            if (it < 0) { all_ok = 0; break; }
            for (int i = 0; i < 12; ++i)
                if (out[i] != cw[i]) { all_ok = 0; break; }
        }
    }
    CHECK(all_ok, "all 2-bit wrong-sign patterns are corrected");

    dsp_ldpc_free(&code);
}

static void test_ldpc_sumproduct_clean(void) {
    printf("[coding] LDPC sum-product accepts a clean codeword fast\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1);

    double llr[12];
    for (int i = 0; i < 12; ++i) llr[i] = 4.0;   /* all strongly bit 0 */
    uint8_t out[12];
    int it = dsp_ldpc_decode_sumproduct(&code, llr, out, 50);

    int zeros = 1;
    for (int i = 0; i < 12; ++i) if (out[i] != 0) zeros = 0;
    CHECK(it >= 0 && it <= 1 && zeros,
          "clean input decodes to the all-zero codeword in <=1 pass");

    dsp_ldpc_free(&code);
}

static void test_ldpc_minsum_corrects(void) {
    printf("[coding] LDPC min-sum corrects multiple soft errors\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1);
    uint8_t cw[12] = {0};

    /* Same 2-bit wrong-sign sweep used for sum-product. */
    int all_ok = 1;
    for (int a = 0; a < 12 && all_ok; ++a) {
        for (int b = a + 1; b < 12; ++b) {
            double llr[12];
            for (int i = 0; i < 12; ++i) llr[i] = 3.0;
            llr[a] = -1.0;
            llr[b] = -1.0;
            uint8_t out[12];
            int it = dsp_ldpc_decode_minsum(&code, llr, out, 0.75, 50);
            if (it < 0) { all_ok = 0; break; }
            for (int i = 0; i < 12; ++i)
                if (out[i] != cw[i]) { all_ok = 0; break; }
        }
    }
    CHECK(all_ok, "min-sum corrects all 2-bit wrong-sign patterns");

    dsp_ldpc_free(&code);
}

static void test_ldpc_ber_monotonic(void) {
    printf("[coding] LDPC bit-error rate rises with channel noise\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1);

    /* The waterfall: more noise must not give a lower BER. */
    double low  = dsp_ldpc_ber_sweep(&code, DSP_LDPC_SUMPRODUCT,
                                     0.4, 1000, 50, 7);
    double high = dsp_ldpc_ber_sweep(&code, DSP_LDPC_SUMPRODUCT,
                                     1.0, 1000, 50, 7);
    CHECK(low >= 0.0 && high >= 0.0 && high > low,
          "BER at high noise exceeds BER at low noise");

    dsp_ldpc_free(&code);
}

static void test_ldpc_minsum_beats_bitflip(void) {
    printf("[coding] soft decoders outperform hard bit-flipping\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1);

    /* At a moderate noise level the soft decoders should have a
     * clearly lower bit-error rate than hard-decision bit-flipping. */
    double bf = dsp_ldpc_ber_sweep(&code, DSP_LDPC_BITFLIP,
                                   0.7, 2000, 50, 99);
    double sp = dsp_ldpc_ber_sweep(&code, DSP_LDPC_SUMPRODUCT,
                                   0.7, 2000, 50, 99);
    double mn = dsp_ldpc_ber_sweep(&code, DSP_LDPC_MINSUM,
                                   0.7, 2000, 50, 99);
    CHECK(sp < bf && mn < bf,
          "sum-product and min-sum both beat bit-flipping");

    dsp_ldpc_free(&code);
}

static void test_ldpc_awgn_llr_sign(void) {
    printf("[coding] AWGN LLR helper has the correct sign convention\n");
    /* A +1 sample favours bit 0 (positive LLR); -1 favours bit 1. */
    CHECK(dsp_ldpc_awgn_llr(1.0, 0.5) > 0.0,
          "+1 sample -> positive LLR (bit 0)");
    CHECK(dsp_ldpc_awgn_llr(-1.0, 0.5) < 0.0,
          "-1 sample -> negative LLR (bit 1)");
}

/* ---- modulation: QAM ------------------------------------------------ */

static void test_qam_roundtrip(void) {
    printf("[modulation] QAM map/demap round-trips for all orders\n");
    dsp_qam_order orders[3] = { DSP_QAM_QPSK, DSP_QAM_16, DSP_QAM_64 };
    int all_ok = 1;
    for (int o = 0; o < 3; ++o) {
        size_t bps = dsp_qam_bits_per_symbol(orders[o]);
        size_t nb  = bps * 100;
        uint8_t bits[6 * 100], back[6 * 100];
        cplx sym[100];
        unsigned s = 12345;
        for (size_t i = 0; i < nb; ++i) {
            s = s * 1103515245u + 12345u;
            bits[i] = (uint8_t)((s >> 16) & 1);
        }
        size_t ns = dsp_qam_modulate(orders[o], bits, nb, sym);
        dsp_qam_demodulate(orders[o], sym, ns, back);
        for (size_t i = 0; i < nb; ++i)
            if (bits[i] != back[i]) all_ok = 0;
    }
    CHECK(all_ok, "QPSK, 16-QAM, 64-QAM all round-trip noiselessly");
}

static void test_qam_unit_energy(void) {
    printf("[modulation] QAM constellations have unit average energy\n");
    dsp_qam_order orders[3] = { DSP_QAM_QPSK, DSP_QAM_16, DSP_QAM_64 };
    int all_ok = 1;
    for (int o = 0; o < 3; ++o) {
        size_t bps = dsp_qam_bits_per_symbol(orders[o]);
        size_t M   = (size_t)orders[o];
        double energy = 0.0;
        for (size_t v = 0; v < M; ++v) {
            uint8_t bits[6];
            for (size_t b = 0; b < bps; ++b)
                bits[b] = (uint8_t)((v >> (bps - 1 - b)) & 1);
            cplx sym;
            dsp_qam_modulate(orders[o], bits, bps, &sym);
            energy += creal(sym) * creal(sym) + cimag(sym) * cimag(sym);
        }
        if (!close(energy / (double)M, 1.0, 1e-9)) all_ok = 0;
    }
    CHECK(all_ok, "mean symbol energy is 1.0 for every order");
}

static void test_qam_soft_llr_sign(void) {
    printf("[modulation] QAM soft demap has the right LLR sign\n");
    /* Map all-zero bits, demap softly: every LLR must favour bit 0. */
    uint8_t bits[4] = {0, 0, 0, 0};
    cplx sym;
    dsp_qam_modulate(DSP_QAM_16, bits, 4, &sym);
    double llr[4];
    dsp_qam_demodulate_soft(DSP_QAM_16, &sym, 1, 0.1, llr);
    int ok = 1;
    for (int i = 0; i < 4; ++i) if (llr[i] <= 0.0) ok = 0;
    CHECK(ok, "all-zero symbol yields positive (bit-0) LLRs");
}

/* ---- modulation: channel -------------------------------------------- */

static void test_channel_awgn_noiseless(void) {
    printf("[modulation] noiseless AWGN channel passes samples through\n");
    dsp_channel ch;
    dsp_channel_init_awgn(&ch, 0.0, 1);
    cplx in[8], out[8];
    for (int i = 0; i < 8; ++i) in[i] = dsp_cplx(i - 4.0, 0.5 * i);
    dsp_channel_apply(&ch, in, out, 8);
    int ok = 1;
    for (int i = 0; i < 8; ++i)
        if (!close(creal(in[i]), creal(out[i]), 1e-12) ||
            !close(cimag(in[i]), cimag(out[i]), 1e-12)) ok = 0;
    CHECK(ok, "a unit-tap noiseless channel is an identity");
    dsp_channel_free(&ch);
}

static void test_channel_multipath(void) {
    printf("[modulation] multipath channel convolves with its taps\n");
    cplx taps[2] = { dsp_cplx(1.0, 0.0), dsp_cplx(0.5, 0.0) };
    dsp_channel ch;
    dsp_channel_init(&ch, taps, 2, 0.0, 1);
    cplx in[4]  = { dsp_cplx(1,0), dsp_cplx(0,0),
                    dsp_cplx(0,0), dsp_cplx(0,0) };
    cplx out[4];
    dsp_channel_apply(&ch, in, out, 4);
    /* An impulse in -> the tap vector out. */
    CHECK(close(creal(out[0]), 1.0, 1e-12) &&
          close(creal(out[1]), 0.5, 1e-12),
          "impulse response equals the channel taps");
    dsp_channel_free(&ch);
}

/* ---- modulation: OFDM ----------------------------------------------- */

static void test_ofdm_roundtrip(void) {
    printf("[modulation] OFDM modulate/demodulate round-trips\n");
    dsp_ofdm o;
    int rc = dsp_ofdm_init(&o, 64, 16);
    CHECK(rc == 0, "OFDM initialises with 64 subcarriers, 16-CP");

    cplx freq[64], time[80], back[64];
    unsigned s = 555;
    for (int i = 0; i < 64; ++i) {
        s = s * 1103515245u + 12345u;
        freq[i] = dsp_cplx(((s >> 16) & 1) ? 1.0 : -1.0,
                           ((s >> 17) & 1) ? 1.0 : -1.0);
    }
    dsp_ofdm_modulate(&o, freq, time);
    dsp_ofdm_demodulate(&o, time, back);

    int ok = 1;
    for (int i = 0; i < 64; ++i)
        if (cabs(back[i] - freq[i]) > 1e-9) ok = 0;
    CHECK(ok, "subcarriers survive the IFFT/CP/FFT round-trip");
}

static void test_ofdm_rejects_non_pow2(void) {
    printf("[modulation] OFDM rejects a non-power-of-two FFT size\n");
    dsp_ofdm o;
    CHECK(dsp_ofdm_init(&o, 48, 8) == -1, "nfft=48 is rejected");
    CHECK(dsp_ofdm_init(&o, 64, 80) == -1, "cp >= nfft is rejected");
}

static void test_ofdm_equalizes_multipath(void) {
    printf("[modulation] OFDM equalizer undoes multipath distortion\n");
    dsp_ofdm o;
    dsp_ofdm_init(&o, 64, 16);

    cplx freq[64], time[80], rx[80], rxf[64], cfr[64];
    unsigned s = 909;
    for (int i = 0; i < 64; ++i) {
        s = s * 1103515245u + 12345u;
        freq[i] = dsp_cplx(((s >> 16) & 1) ? 0.707 : -0.707,
                           ((s >> 17) & 1) ? 0.707 : -0.707);
    }
    dsp_ofdm_modulate(&o, freq, time);

    /* 3-tap multipath, noiseless: equalization should be near-exact. */
    cplx taps[3] = { dsp_cplx(1.0, 0.0),
                     dsp_cplx(0.3, 0.1),
                     dsp_cplx(0.1, -0.05) };
    dsp_channel ch;
    dsp_channel_init(&ch, taps, 3, 0.0, 1);
    dsp_channel_apply(&ch, time, rx, 80);
    dsp_ofdm_demodulate(&o, rx, rxf);
    dsp_channel_frequency_response(&ch, cfr, 64);
    dsp_ofdm_equalize(&o, rxf, cfr);

    double max_err = 0.0;
    for (int i = 0; i < 64; ++i) {
        double d = cabs(rxf[i] - freq[i]);
        if (d > max_err) max_err = d;
    }
    CHECK(max_err < 1e-6,
          "per-subcarrier equalization recovers the symbols");
    dsp_channel_free(&ch);
}

/* ---- modulation: coded OFDM ----------------------------------------- */

static void test_coded_ofdm_noiseless(void) {
    printf("[modulation] coded OFDM is error-free on a clean channel\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 64, 128, 3, 6, 1);
    dsp_ofdm o;
    dsp_ofdm_init(&o, 64, 16);
    dsp_coded_ofdm cfg = { .ofdm = o,
                           .order = DSP_QAM_QPSK,
                           .code = &code };

    /* Noiseless 3-tap multipath: OFDM + equalizer should be exact, so
     * the decoded frame has zero bit errors. */
    cplx taps[3] = { dsp_cplx(1.0, 0.0),
                     dsp_cplx(0.3, 0.1),
                     dsp_cplx(0.1, -0.05) };
    dsp_channel ch;
    dsp_channel_init(&ch, taps, 3, 0.0, 1);

    size_t be = 999, re = 999;
    int rc = dsp_coded_ofdm_run_frame(&cfg, &ch, 50, &be, &re);
    CHECK(rc == 0 && be == 0 && re == 0,
          "noiseless multipath frame decodes with 0 errors");

    dsp_channel_free(&ch);
    dsp_ldpc_free(&code);
}

static void test_coded_ofdm_dimension_check(void) {
    printf("[modulation] coded OFDM rejects a codeword/symbol mismatch\n");
    /* QPSK on 64 subcarriers needs a 128-bit codeword; give it 120. */
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 60, 120, 3, 6, 1);
    dsp_ofdm o;
    dsp_ofdm_init(&o, 64, 16);
    dsp_coded_ofdm cfg = { .ofdm = o,
                           .order = DSP_QAM_QPSK,
                           .code = &code };
    dsp_channel ch;
    dsp_channel_init_awgn(&ch, 0.1, 1);
    CHECK(dsp_coded_ofdm_run_frame(&cfg, &ch, 50, NULL, NULL) == -1,
          "code length != nfft*bits_per_symbol is rejected");
    dsp_channel_free(&ch);
    dsp_ldpc_free(&code);
}

static void test_coded_ofdm_corrects_errors(void) {
    printf("[modulation] coded OFDM: FEC lowers BER below the raw rate\n");
    dsp_ldpc code;
    dsp_ldpc_make_regular(&code, 64, 128, 3, 6, 1);
    dsp_ofdm o;
    dsp_ofdm_init(&o, 64, 16);
    dsp_coded_ofdm cfg = { .ofdm = o,
                           .order = DSP_QAM_QPSK,
                           .code = &code };

    cplx taps[3] = { dsp_cplx(1.0, 0.0),
                     dsp_cplx(0.25, 0.1),
                     dsp_cplx(0.1, -0.05) };
    dsp_channel ch;
    dsp_channel_init(&ch, taps, 3, 0.6, 4242);

    /* Accumulate raw (pre-FEC) and coded (post-FEC) errors. */
    size_t raw_tot = 0, cod_tot = 0;
    for (int f = 0; f < 150; ++f) {
        size_t be, re;
        dsp_coded_ofdm_run_frame(&cfg, &ch, 50, &be, &re);
        raw_tot += re;
        cod_tot += be;
    }
    CHECK(raw_tot > 0 && cod_tot < raw_tot,
          "LDPC reduces the bit-error count versus raw demapping");

    dsp_channel_free(&ch);
    dsp_ldpc_free(&code);
}

/* ---- modulation: pulse shaping & sync -------------------------------- */

static void test_rrc_energy_and_symmetry(void) {
    printf("[modulation] RRC filter is unit-energy and symmetric\n");
    enum { NT = 65 };
    double taps[NT];
    dsp_rrc_design(taps, NT, 8, 0.25);

    double energy = 0.0;
    for (int i = 0; i < NT; ++i) energy += taps[i] * taps[i];
    CHECK(close(energy, 1.0, 1e-9), "RRC tap energy normalised to 1");

    int sym = 1;
    for (int i = 0; i < NT; ++i)
        if (!close(taps[i], taps[NT - 1 - i], 1e-9)) sym = 0;
    CHECK(sym, "RRC taps are symmetric (linear phase)");
}

static void test_pulse_shape_length(void) {
    printf("[modulation] pulse shaping expands by samples-per-symbol\n");
    enum { NSYM = 20, SPS = 4, NT = 33 };
    cplx syms[NSYM], out[NSYM * SPS];
    double taps[NT];
    dsp_rrc_design(taps, NT, SPS, 0.3);
    for (int i = 0; i < NSYM; ++i) syms[i] = dsp_cplx(1.0, -1.0);
    size_t n = dsp_pulse_shape(syms, NSYM, SPS, taps, NT, out);
    CHECK(n == dsp_pulse_shaped_len(NSYM, SPS),
          "shaped length == nsym * sps");
}

static void test_carrier_pll_locks(void) {
    printf("[modulation] carrier PLL locks out a phase offset\n");
    cplx in[300], out[300];
    for (int i = 0; i < 300; ++i) {
        cplx sym = dsp_cplx((i & 1) ? 0.707 : -0.707,
                            (i & 2) ? 0.707 : -0.707);
        in[i] = sym * cexp(0.35 * I);        /* fixed 0.35 rad offset */
    }
    dsp_carrier_pll pll;
    dsp_carrier_pll_init(&pll, 0.05, 0.0025);
    dsp_carrier_recover(&pll, in, out, 300);

    /* After locking, the last symbols should match the un-rotated
     * constellation closely. */
    double err = 0.0;
    for (int i = 280; i < 300; ++i) {
        cplx sym = dsp_cplx((i & 1) ? 0.707 : -0.707,
                            (i & 2) ? 0.707 : -0.707);
        err += cabs(out[i] - sym);
    }
    CHECK(err / 20.0 < 0.05,
          "constellation error is near zero after the loop locks");
}

static void test_timing_resample_length(void) {
    printf("[modulation] timing resampler decimates to symbol rate\n");
    enum { N = 128, SPS = 4 };
    cplx samples[N], syms[N / SPS];
    for (int i = 0; i < N; ++i) samples[i] = dsp_cplx((double)i, 0.0);
    size_t ns = dsp_timing_resample(samples, N, SPS, 0.0, syms);
    CHECK(ns == N / SPS, "one output symbol per sps input samples");
}

/* ---- image: container ----------------------------------------------- */

static void test_image_alloc_and_access(void) {
    printf("[image] image allocation and clamped pixel access\n");
    dsp_image img;
    int rc = dsp_image_alloc(&img, 8, 6);
    CHECK(rc == 0 && img.width == 8 && img.height == 6,
          "image allocates with the requested dimensions");

    img.data[2 * 8 + 3] = 42.0;
    CHECK(close(dsp_image_at(&img, 3, 2), 42.0, 1e-12),
          "pixel read-back is correct");
    /* Out-of-range coordinates clamp to the edge, not crash. */
    CHECK(close(dsp_image_at(&img, -5, -5), img.data[0], 1e-12),
          "negative coordinates clamp to the corner pixel");

    dsp_image_free(&img);
}

/* ---- image: 2-D transforms ------------------------------------------ */

static void test_fft2d_roundtrip(void) {
    printf("[image] 2-D FFT round-trips back to the original image\n");
    dsp_image img;
    dsp_image_alloc(&img, 16, 16);
    for (size_t i = 0; i < 256; ++i)
        img.data[i] = (double)((i * 13 + 7) % 200);

    double re[256], im[256];
    dsp_image back;
    dsp_image_alloc(&back, 16, 16);
    int rc1 = dsp_fft2d(&img, re, im);
    int rc2 = dsp_ifft2d(re, im, &back);

    int ok = (rc1 == 0 && rc2 == 0);
    for (size_t i = 0; i < 256; ++i)
        if (!close(img.data[i], back.data[i], 1e-6)) ok = 0;
    CHECK(ok, "ifft2d(fft2d(image)) == image");

    dsp_image_free(&img);
    dsp_image_free(&back);
}

static void test_fft2d_rejects_non_pow2(void) {
    printf("[image] 2-D FFT rejects non-power-of-two dimensions\n");
    dsp_image img;
    dsp_image_alloc(&img, 20, 16);
    double re[320], im[320];
    CHECK(dsp_fft2d(&img, re, im) == -1,
          "a width of 20 is rejected");
    dsp_image_free(&img);
}

static void test_dct2d_roundtrip(void) {
    printf("[image] 2-D DCT round-trips back to the original image\n");
    dsp_image img;
    dsp_image_alloc(&img, 24, 24);
    for (size_t i = 0; i < 576; ++i)
        img.data[i] = (double)((i * 5 + 11) % 256);

    double coeffs[576];
    dsp_image back;
    dsp_image_alloc(&back, 24, 24);
    dsp_dct2d(&img, coeffs);
    dsp_idct2d(coeffs, &back);

    int ok = 1;
    for (size_t i = 0; i < 576; ++i)
        if (!close(img.data[i], back.data[i], 1e-6)) ok = 0;
    CHECK(ok, "idct2d(dct2d(image)) == image");

    dsp_image_free(&img);
    dsp_image_free(&back);
}

static void test_dct8x8_compaction(void) {
    printf("[image] 8x8 DCT compacts a smooth block's energy\n");
    /* A smooth gradient block - energy should pile into the corner. */
    double block[64], coeffs[64];
    for (int y = 0; y < 8; ++y)
        for (int x = 0; x < 8; ++x)
            block[y * 8 + x] = (double)(x + y) * 10.0;
    dsp_dct8x8(block, coeffs);

    double total = 0.0, low = 0.0;
    for (int i = 0; i < 64; ++i) {
        double e = coeffs[i] * coeffs[i];
        total += e;
        if ((i / 8) < 3 && (i % 8) < 3) low += e;
    }
    CHECK(low / total > 0.98,
          "the low-frequency 3x3 corner holds >98% of the energy");
}

/* ---- image: spatial filters ----------------------------------------- */

static void test_gaussian_preserves_flat(void) {
    printf("[image] Gaussian blur preserves a flat region\n");
    dsp_image img, out;
    dsp_image_alloc(&img, 32, 32);
    dsp_image_alloc(&out, 32, 32);
    for (size_t i = 0; i < 1024; ++i) img.data[i] = 90.0;
    dsp_gaussian_blur(&img, &out, 1.5);

    int ok = 1;
    for (size_t i = 0; i < 1024; ++i)
        if (!close(out.data[i], 90.0, 1e-6)) ok = 0;
    CHECK(ok, "blurring a constant image leaves it unchanged");

    dsp_image_free(&img);
    dsp_image_free(&out);
}

static void test_sobel_finds_edges(void) {
    printf("[image] Sobel responds to an edge, not to flat regions\n");
    dsp_image img, out;
    dsp_image_alloc(&img, 32, 32);
    dsp_image_alloc(&out, 32, 32);
    /* A vertical step edge down the middle. */
    for (size_t y = 0; y < 32; ++y)
        for (size_t x = 0; x < 32; ++x)
            img.data[y * 32 + x] = (x < 16) ? 0.0 : 255.0;
    dsp_sobel(&img, &out);

    double on_edge  = out.data[16 * 32 + 15];
    double off_edge = out.data[16 * 32 + 4];
    CHECK(on_edge > 100.0 && off_edge < 1.0,
          "strong response at the edge, near zero in flat areas");

    dsp_image_free(&img);
    dsp_image_free(&out);
}

static void test_median_removes_spike(void) {
    printf("[image] median filter removes an isolated spike\n");
    dsp_image img, out;
    dsp_image_alloc(&img, 16, 16);
    dsp_image_alloc(&out, 16, 16);
    for (size_t i = 0; i < 256; ++i) img.data[i] = 60.0;
    img.data[8 * 16 + 8] = 255.0;            /* salt noise */
    dsp_median_filter(&img, &out, 3);
    CHECK(close(out.data[8 * 16 + 8], 60.0, 1e-9),
          "the spike is replaced by the neighbourhood median");

    dsp_image_free(&img);
    dsp_image_free(&out);
}

/* ---- image: point operators & 2-D wavelet --------------------------- */

static void test_histogram_total(void) {
    printf("[image] histogram bins sum to the pixel count\n");
    dsp_image img;
    dsp_image_alloc(&img, 20, 15);
    for (size_t i = 0; i < 300; ++i) img.data[i] = (double)(i % 256);
    size_t hist[256];
    dsp_histogram(&img, hist);
    size_t total = 0;
    for (int i = 0; i < 256; ++i) total += hist[i];
    CHECK(total == 300, "histogram counts add up to width*height");
    dsp_image_free(&img);
}

static void test_otsu_threshold_bimodal(void) {
    printf("[image] Otsu threshold separates a bimodal image\n");
    dsp_image img;
    dsp_image_alloc(&img, 32, 32);
    for (size_t i = 0; i < 1024; ++i)
        img.data[i] = (i < 512) ? 40.0 : 210.0;
    double t = dsp_threshold_otsu(&img);
    CHECK(t >= 40.0 && t < 210.0,
          "threshold falls between the two intensity modes");
    dsp_image_free(&img);
}

static void test_threshold_binarizes(void) {
    printf("[image] thresholding produces a pure black/white image\n");
    dsp_image img, out;
    dsp_image_alloc(&img, 16, 16);
    dsp_image_alloc(&out, 16, 16);
    for (size_t i = 0; i < 256; ++i) img.data[i] = (double)i;
    dsp_threshold(&img, &out, 128.0);
    int ok = 1;
    for (size_t i = 0; i < 256; ++i)
        if (out.data[i] != 0.0 && out.data[i] != 255.0) ok = 0;
    CHECK(ok, "every output pixel is exactly 0 or 255");
    dsp_image_free(&img);
    dsp_image_free(&out);
}

static void test_dwt2d_roundtrip(void) {
    printf("[image] 2-D Haar wavelet round-trips back to the image\n");
    dsp_image img, work;
    dsp_image_alloc(&img, 16, 16);
    for (size_t i = 0; i < 256; ++i)
        img.data[i] = (double)((i * 9 + 4) % 220);
    dsp_image_copy(&work, &img);

    int rc1 = dsp_dwt2d_haar(&work);
    int rc2 = dsp_idwt2d_haar(&work);
    int ok = (rc1 == 0 && rc2 == 0);
    for (size_t i = 0; i < 256; ++i)
        if (!close(img.data[i], work.data[i], 1e-9)) ok = 0;
    CHECK(ok, "idwt2d(dwt2d(image)) == image");

    dsp_image_free(&img);
    dsp_image_free(&work);
}

static void test_dwt2d_rejects_odd(void) {
    printf("[image] 2-D wavelet rejects odd dimensions\n");
    dsp_image img;
    dsp_image_alloc(&img, 15, 16);
    CHECK(dsp_dwt2d_haar(&img) == -1, "an odd width is rejected");
    dsp_image_free(&img);
}

int main(void) {
    printf("DSP GUIDE - test suite\n");
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

    test_block_interleave_roundtrip();
    test_block_interleave_size_check();
    test_block_interleave_spreads_burst();
    test_block_interleave_rescues_rs();
    test_conv_interleaver_roundtrip();

    test_ldpc_from_matrix();
    test_ldpc_syndrome();
    test_ldpc_regular_construction();
    test_ldpc_regular_rejects_bad_params();
    test_ldpc_bitflip_single_error();
    test_ldpc_sumproduct_corrects();
    test_ldpc_sumproduct_clean();
    test_ldpc_minsum_corrects();
    test_ldpc_ber_monotonic();
    test_ldpc_minsum_beats_bitflip();
    test_ldpc_awgn_llr_sign();

    test_qam_roundtrip();
    test_qam_unit_energy();
    test_qam_soft_llr_sign();
    test_channel_awgn_noiseless();
    test_channel_multipath();
    test_ofdm_roundtrip();
    test_ofdm_rejects_non_pow2();
    test_ofdm_equalizes_multipath();
    test_coded_ofdm_noiseless();
    test_coded_ofdm_dimension_check();
    test_coded_ofdm_corrects_errors();

    test_rrc_energy_and_symmetry();
    test_pulse_shape_length();
    test_carrier_pll_locks();
    test_timing_resample_length();

    test_image_alloc_and_access();
    test_fft2d_roundtrip();
    test_fft2d_rejects_non_pow2();
    test_dct2d_roundtrip();
    test_dct8x8_compaction();
    test_gaussian_preserves_flat();
    test_sobel_finds_edges();
    test_median_removes_spike();
    test_histogram_total();
    test_otsu_threshold_bimodal();
    test_threshold_binarizes();
    test_dwt2d_roundtrip();
    test_dwt2d_rejects_odd();

    test_lms_converges();

    printf("\n============================\n");
    printf("Tests run: %d   Passed: %d   Failed: %d\n",
           tests_run, tests_run - tests_failed, tests_failed);
    return tests_failed == 0 ? 0 : 1;
}
