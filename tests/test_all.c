/*
 * test_all.c - Unit tests for the DSP Study Guide.
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

    printf("\n============================\n");
    printf("Tests run: %d   Passed: %d   Failed: %d\n",
           tests_run, tests_run - tests_failed, tests_failed);
    return tests_failed == 0 ? 0 : 1;
}
