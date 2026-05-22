/*
 * main.c - DSP Study Guide demo runner.
 *
 * Exercises every module with a small, readable example and prints
 * annotated output. Build target: dsp_demo.
 */
#include "dsp.h"
#include <stdio.h>
#include <stdlib.h>

static void section(const char *title) {
    printf("\n========================================\n");
    printf("  %s\n", title);
    printf("========================================\n");
}

/* ---- Frequency-domain transforms ------------------------------------ */

static void demo_transforms(void) {
    section("FREQUENCY-DOMAIN TRANSFORMS");

    /* A length-8 signal: a pure cosine at exactly 1 cycle per window.
     * Its energy should land entirely in bin 1 (and its mirror, 7). */
    enum { N = 8 };
    cplx sig[N];
    for (int t = 0; t < N; ++t)
        sig[t] = dsp_cplx(cos(2.0 * M_PI * 1.0 * t / N), 0.0);

    cplx dft_out[N], fft_buf[N];
    dsp_dft(sig, dft_out, N);
    for (int i = 0; i < N; ++i) fft_buf[i] = sig[i];
    dsp_fft(fft_buf, N);

    printf("Pure cosine, 1 cycle in an 8-sample window.\n");
    printf("DFT and FFT must agree (FFT is the fast route to the same answer):\n");
    printf("  bin |   DFT |mag   FFT |mag\n");
    for (int k = 0; k < N; ++k)
        printf("  %3d | %8.3f      %8.3f\n",
               k, dsp_mag(dft_out[k]), dsp_mag(fft_buf[k]));
    printf("-> energy concentrates in bins 1 and 7, as expected.\n");

    /* DCT energy compaction on a smooth ramp. */
    double ramp[N], coeff[N];
    for (int i = 0; i < N; ++i) ramp[i] = (double)i;
    dsp_dct(ramp, coeff, N);
    printf("\nDCT of a smooth ramp - note how energy compacts into the\n");
    printf("first few coefficients (the basis of JPEG compression):\n  ");
    for (int k = 0; k < N; ++k) printf("%7.3f ", coeff[k]);
    printf("\n");
}

/* ---- Digital filtering ---------------------------------------------- */

static void demo_filtering(void) {
    section("DIGITAL FILTERING (FIR vs IIR)");

    enum { N = 64 };
    double x[N], y_fir[N], y_iir[N];

    /* Input: low-frequency tone plus a fast tone we want to remove. */
    for (int n = 0; n < N; ++n)
        x[n] = sin(2.0 * M_PI * 0.05 * n)      /* keep this  */
             + 0.5 * sin(2.0 * M_PI * 0.40 * n); /* reject this */

    /* FIR low-pass, 31 symmetric taps -> exactly linear phase. */
    enum { TAPS = 31 };
    double taps[TAPS];
    dsp_fir_design_lowpass(taps, TAPS, 0.10);
    dsp_fir_apply(x, y_fir, N, taps, TAPS);

    /* IIR low-pass biquad -> same job, far fewer coefficients. */
    dsp_biquad bq;
    dsp_iir_design_lowpass(&bq, 0.10, 0.707);   /* Butterworth Q */
    dsp_iir_apply(&bq, x, y_iir, N);

    printf("FIR : 31 taps, symmetric  -> linear phase, always stable\n");
    printf("IIR : 5 coefficients      -> efficient, stable? %s\n",
           dsp_iir_is_stable(&bq) ? "yes" : "NO");
    printf("\n   n |   input |  FIR out |  IIR out\n");
    for (int n = 20; n < 30; ++n)
        printf("  %3d | %8.3f | %8.3f | %8.3f\n",
               n, x[n], y_fir[n], y_iir[n]);
    printf("-> both attenuate the fast 0.40 tone; IIR does it with\n");
    printf("   5 coefficients vs the FIR's 31 taps.\n");
}

/* ---- Signal operations ---------------------------------------------- */

static void demo_operations(void) {
    section("SIGNAL OPERATIONS (convolution & correlation)");

    /* Convolution: a 3-point moving-average filter applied to a step. */
    double x[6] = {1, 1, 1, 1, 1, 1};
    double h[3] = {1.0/3, 1.0/3, 1.0/3};
    size_t clen = dsp_conv_len(6, 3);
    double *yc  = malloc(clen * sizeof(double));
    double *yf  = malloc(clen * sizeof(double));
    dsp_convolve(x, 6, h, 3, yc);
    dsp_convolve_fft(x, 6, h, 3, yf);

    printf("Convolution of a step with a 3-point averager.\n");
    printf("Direct and FFT-based methods agree:\n  direct: ");
    for (size_t i = 0; i < clen; ++i) printf("%.3f ", yc[i]);
    printf("\n  fft   : ");
    for (size_t i = 0; i < clen; ++i) printf("%.3f ", yf[i]);
    printf("\n");
    free(yc); free(yf);

    /* Correlation: find a known delay between two signals. */
    enum { M = 16 };
    double a[M], b[M];
    for (int i = 0; i < M; ++i) {
        a[i] = (i == 4) ? 1.0 : 0.0;          /* pulse at index 4  */
        b[i] = (i == 9) ? 1.0 : 0.0;          /* same pulse, +5    */
    }
    long delay = dsp_estimate_delay(a, M, b, M);
    printf("\nCross-correlation delay estimation (radar/sonar principle):\n");
    printf("  pulse in A at 4, in B at 9  -> estimated lag = %ld (expect 5)\n",
           delay);
}

/* ---- Spectral analysis ---------------------------------------------- */

static void demo_spectral(void) {
    section("SPECTRAL ANALYSIS (windowing)");

    enum { N = 16 };
    double w[N];
    dsp_window_type types[4] = {
        DSP_WIN_RECTANGULAR, DSP_WIN_HAMMING,
        DSP_WIN_HANNING,     DSP_WIN_BLACKMAN
    };
    const char *names[4] = { "Rectangular", "Hamming",
                             "Hanning", "Blackman" };

    printf("Window taper shapes (16 points). Smoother edges -> less\n");
    printf("spectral leakage, at the cost of a wider main lobe.\n");
    for (int t = 0; t < 4; ++t) {
        dsp_window_generate(types[t], w, N);
        printf("  %-12s: ", names[t]);
        for (int k = 0; k < N; k += 2) printf("%.2f ", w[k]);
        printf("\n");
    }
}

/* ---- Sample rate conversion ----------------------------------------- */

static void demo_sampling(void) {
    section("SAMPLE RATE CONVERSION");

    enum { N = 32 };
    double x[N];
    for (int n = 0; n < N; ++n)
        x[n] = sin(2.0 * M_PI * 0.05 * n);

    /* Decimate by 2. */
    size_t dlen = dsp_decimate_len(N, 2);
    double *dec = malloc(dlen * sizeof(double));
    dsp_decimate(x, N, 2, dec);

    /* Interpolate by 2. */
    size_t ilen = dsp_interpolate_len(N, 2);
    double *itp = malloc(ilen * sizeof(double));
    dsp_interpolate(x, N, 2, itp);

    /* Rational resample by 3/2. */
    size_t rlen = dsp_resample_len(N, 3, 2);
    double *res = malloc(rlen * sizeof(double));
    dsp_resample(x, N, 3, 2, res);

    printf("Input length %d (sine wave).\n", N);
    printf("  decimate  by 2   -> length %zu  (lower rate)\n", dlen);
    printf("  interpolate by 2 -> length %zu  (higher rate)\n", ilen);
    printf("  resample by 3/2  -> length %zu  (rational change)\n", rlen);
    printf("-> the anti-alias / anti-image low-pass filter always sits\n");
    printf("   on the high-rate side of the operation.\n");

    free(dec); free(itp); free(res);
}

/* ---- Multi-resolution analysis -------------------------------------- */

static void demo_wavelet(void) {
    section("MULTI-RESOLUTION ANALYSIS (wavelet transform)");

    enum { N = 8 };
    double data[N]  = {4, 6, 10, 12, 8, 6, 5, 5};
    double orig[N];
    for (int i = 0; i < N; ++i) orig[i] = data[i];

    printf("Original signal : ");
    for (int i = 0; i < N; ++i) printf("%6.2f ", orig[i]);
    printf("\n");

    int levels = dsp_dwt_haar(data, N);
    printf("Haar DWT (%d levels): ", levels);
    for (int i = 0; i < N; ++i) printf("%6.2f ", data[i]);
    printf("\n-> leading value is the coarse average; the rest are\n");
    printf("   detail coefficients from coarse to fine.\n");

    dsp_idwt_haar(data, N);
    printf("Reconstructed   : ");
    for (int i = 0; i < N; ++i) printf("%6.2f ", data[i]);
    printf("\n-> inverse DWT recovers the original exactly.\n");
}

int main(void) {
    printf("DSP STUDY GUIDE - annotated demo\n");
    printf("C implementation of common digital signal processing algorithms.\n");

    demo_transforms();
    demo_filtering();
    demo_operations();
    demo_spectral();
    demo_sampling();
    demo_wavelet();

    printf("\nAll demos complete.\n");
    return 0;
}
