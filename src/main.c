/*
 * main.c - DSP Guide demo runner.
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

/* ---- Error detection ------------------------------------------------ */

static void demo_detection(void) {
    section("ERROR DETECTION (CRC, parity, checksum)");

    uint8_t msg[8] = { 'D','S','P','-','g','u','i','d' };

    /* Parity. */
    int par = dsp_parity_compute(msg, 8, DSP_PARITY_EVEN);
    printf("Even parity bit for the message: %d\n", par);
    printf("  check with correct bit  : %s\n",
           dsp_parity_check(msg, 8, DSP_PARITY_EVEN, par) ? "ok" : "FAIL");
    printf("  check with flipped bit  : %s (error detected)\n",
           dsp_parity_check(msg, 8, DSP_PARITY_EVEN, par ^ 1)
               ? "ok" : "mismatch");

    /* Internet checksum. */
    uint16_t sum = dsp_checksum16(msg, 8);
    printf("\n16-bit Internet checksum: 0x%04X\n", sum);
    printf("  verify intact message : %s\n",
           dsp_checksum16_verify(msg, 8, sum) ? "ok" : "FAIL");
    msg[2] ^= 0x20;                       /* corrupt one byte */
    printf("  verify after corruption: %s (corruption caught)\n",
           dsp_checksum16_verify(msg, 8, sum) ? "ok" : "mismatch");
    msg[2] ^= 0x20;                       /* restore */

    /* CRC-32. */
    uint32_t crc = dsp_crc32(msg, 8);
    printf("\nCRC-32 of the message: 0x%08X\n", crc);
    uint8_t burst[8];
    for (int i = 0; i < 8; ++i) burst[i] = msg[i];
    burst[3] ^= 0xFF;                     /* an 8-bit burst error */
    printf("  CRC after a burst error: 0x%08X %s\n",
           dsp_crc32(burst, 8),
           dsp_crc32(burst, 8) == crc ? "(missed!)" : "(burst detected)");
}

/* ---- Forward error correction --------------------------------------- */

static void demo_correction(void) {
    section("FORWARD ERROR CORRECTION (Hamming, RS, Viterbi)");

    /* Hamming(7,4): corrupt one bit, watch it get repaired. */
    uint8_t nibble = 0x0B;                /* 1011 */
    uint8_t code   = dsp_hamming74_encode(nibble);
    uint8_t bad    = code ^ 0x04;         /* flip bit 3 */
    uint8_t fixed;
    uint8_t recovered = dsp_hamming74_decode(bad, &fixed);
    printf("Hamming(7,4): sent nibble 0x%X, codeword 0x%02X\n",
           nibble, code);
    printf("  received 0x%02X (bit 3 flipped) -> syndrome %d\n",
           bad, dsp_hamming74_syndrome(bad));
    printf("  decoded nibble 0x%X  %s\n", recovered,
           recovered == nibble ? "(single-bit error corrected)" : "(FAIL)");

    /* Reed-Solomon: corrupt a burst of symbols, correct them. */
    dsp_rs rs;
    dsp_rs_init(&rs, 6);                  /* 6 parity symbols -> t = 3 */
    uint8_t data[10] = { 10,20,30,40,50,60,70,80,90,100 };
    uint8_t parity[6];
    dsp_rs_encode(&rs, data, 10, parity);

    uint8_t cw[16];
    for (int i = 0; i < 10; ++i) cw[i] = data[i];
    for (int i = 0; i < 6;  ++i) cw[10 + i] = parity[i];
    cw[2] ^= 0x55;  cw[3] ^= 0xAA;  cw[4] ^= 0x0F;   /* 3-symbol burst */

    int nfix = dsp_rs_decode(&rs, cw, 16);
    int rs_ok = 1;
    for (int i = 0; i < 10; ++i) if (cw[i] != data[i]) rs_ok = 0;
    printf("\nReed-Solomon RS(16,10): corrupted 3 symbols\n");
    printf("  decoder corrected %d symbol errors  %s\n", nfix,
           rs_ok ? "(message fully recovered)" : "(FAIL)");

    /* Convolutional + Viterbi: hard vs soft decision. */
    uint8_t bits[8] = { 1,0,1,1,0,0,1,0 };
    uint8_t enc[2 * (8 + 2)];
    size_t enc_len = dsp_conv_encode(bits, 8, enc);

    /* Inject two bit errors into the coded stream. */
    uint8_t rx[2 * (8 + 2)];
    for (size_t i = 0; i < enc_len; ++i) rx[i] = enc[i];
    rx[5] ^= 1;  rx[12] ^= 1;

    uint8_t dec[8];
    size_t dn = dsp_viterbi_decode(rx, enc_len, dec);
    int v_ok = (dn == 8);
    for (size_t i = 0; i < dn; ++i) if (dec[i] != bits[i]) v_ok = 0;
    printf("\nConvolutional code + Viterbi (hard decision):\n");
    printf("  8 bits encoded to %zu, 2 errors injected\n", enc_len);
    printf("  decoded %zu bits  %s\n", dn,
           v_ok ? "(both errors corrected)" : "(FAIL)");

    /* Soft-decision: same errors, but as weak confidences not flips. */
    double soft[2 * (8 + 2)];
    for (size_t i = 0; i < enc_len; ++i)
        soft[i] = enc[i] ? -1.0 : 1.0;    /* clean antipodal mapping */
    soft[5]  = -0.1;  soft[12] = 0.1;     /* two low-confidence samples */
    uint8_t decs[8];
    size_t dns = dsp_viterbi_decode_soft(soft, enc_len, decs);
    int s_ok = (dns == 8);
    for (size_t i = 0; i < dns; ++i) if (decs[i] != bits[i]) s_ok = 0;
    printf("  soft-decision decode %s\n",
           s_ok ? "(uses analog confidence for better accuracy)"
                : "(FAIL)");
}

/* ---- Channel equalization ------------------------------------------- */

static void demo_equalization(void) {
    section("CHANNEL EQUALIZATION (LMS adaptive filter)");

    /* A known training sequence (+1/-1 symbols). */
    enum { N = 400 };
    double tx[N], rx[N];
    unsigned seed = 12345;
    for (int i = 0; i < N; ++i) {
        seed = seed * 1103515245u + 12345u;
        tx[i] = ((seed >> 16) & 1) ? 1.0 : -1.0;
    }

    /* A simple 2-tap multipath channel: a delayed echo smears symbols. */
    rx[0] = tx[0];
    for (int i = 1; i < N; ++i)
        rx[i] = tx[i] + 0.6 * tx[i - 1];

    dsp_lms eq;
    dsp_lms_init(&eq, 9, 0.02);

    /* MSE over the first vs last quarter shows the filter converging. */
    double mse_early = dsp_lms_train(&eq, rx, tx, N / 4);
    dsp_lms_train(&eq, rx + N / 4, tx + N / 4, N / 2);
    double mse_late  = dsp_lms_train(&eq, rx + 3 * N / 4,
                                     tx + 3 * N / 4, N / 4);

    printf("2-tap multipath channel (echo at 0.6x, one symbol late).\n");
    printf("LMS equalizer, 9 taps, learns the channel inverse:\n");
    printf("  mean-squared error, first quarter : %.5f\n", mse_early);
    printf("  mean-squared error, last quarter  : %.5f\n", mse_late);
    printf("-> error drops as the taps converge; the cleaned signal\n");
    printf("   is what the FEC decoder then operates on.\n");

    dsp_lms_free(&eq);
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
    demo_detection();
    demo_correction();
    demo_equalization();

    printf("\nAll demos complete.\n");
    return 0;
}
