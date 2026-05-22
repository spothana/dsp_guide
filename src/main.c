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

/* ---- Advanced spectral estimation ----------------------------------- */

/* Locate the two highest interior peaks of a spectrum/pseudospectrum. */
static void find_two_peaks(const double *p, size_t nf,
                           double *f1, double *f2) {
    double best1 = -1.0, best2 = -1.0;
    size_t i1 = 0, i2 = 0;
    for (size_t i = 1; i + 1 < nf; ++i) {
        if (p[i] > p[i - 1] && p[i] > p[i + 1]) {
            if (p[i] > best1) {
                best2 = best1; i2 = i1;
                best1 = p[i];  i1 = i;
            } else if (p[i] > best2) {
                best2 = p[i];  i2 = i;
            }
        }
    }
    *f1 = 0.5 * (double)i1 / (double)nf;
    *f2 = 0.5 * (double)i2 / (double)nf;
}

static void demo_estimation(void) {
    section("ADVANCED SPECTRAL ESTIMATION (AR, ARMA, MUSIC, ESPRIT)");

    /* Two closely spaced tones in a SHORT, noisy record. The FFT bin
     * width here is 1/N = 1/64 ~ 0.016, and the tones sit only 0.03
     * apart - right at the edge of what the FFT can resolve. */
    enum { N = 64 };
    double x[N];
    unsigned seed = 20240;
    for (int n = 0; n < N; ++n) {
        seed = seed * 1103515245u + 12345u;
        double noise = ((double)((seed >> 16) & 0xFFFF) / 65535.0) - 0.5;
        x[n] = cos(2.0 * M_PI * 0.20 * n)
             + cos(2.0 * M_PI * 0.23 * n)
             + 0.10 * noise;
    }
    printf("Test signal: two tones at f = 0.20 and 0.23, in a short\n");
    printf("64-sample record. The FFT bin width (~0.016) barely\n");
    printf("separates them; these methods resolve them cleanly.\n");

    enum { ORDER = 14, NF = 1000 };

    /* AR via Yule-Walker (Levinson-Durbin). */
    double r[ORDER + 1];
    dsp_autocorr(x, N, ORDER, r);
    double a_yw[ORDER], sig_yw;
    double psd[NF], f1, f2;
    dsp_ar_yule_walker(r, ORDER, a_yw, &sig_yw);
    dsp_ar_psd(a_yw, ORDER, sig_yw, psd, NF);
    find_two_peaks(psd, NF, &f1, &f2);
    printf("\nAR - Yule-Walker (all-pole model, Levinson-Durbin):\n");
    printf("  spectral peaks at f = %.3f and %.3f\n", f1, f2);

    /* AR via Burg's method. */
    double a_bg[ORDER], sig_bg;
    dsp_ar_burg(x, N, ORDER, a_bg, &sig_bg);
    dsp_ar_psd(a_bg, ORDER, sig_bg, psd, NF);
    find_two_peaks(psd, NF, &f1, &f2);
    printf("\nAR - Burg (forward-backward, maximum entropy):\n");
    printf("  spectral peaks at f = %.3f and %.3f\n", f1, f2);

    /* ARMA via the modified Yule-Walker method. */
    double arma_a[6], arma_b[5];
    if (dsp_arma_estimate(x, N, 6, 4, arma_a, arma_b) == 0) {
        dsp_arma_psd(arma_a, 6, arma_b, 4, psd, NF);
        find_two_peaks(psd, NF, &f1, &f2);
        printf("\nARMA - poles and zeros (modified Yule-Walker):\n");
        printf("  spectral peaks at f = %.3f and %.3f\n", f1, f2);
    }

    /* MUSIC subspace pseudospectrum. */
    double pmusic[NF];
    if (dsp_music(x, N, 2, 16, pmusic, NF) == 0) {
        find_two_peaks(pmusic, NF, &f1, &f2);
        printf("\nMUSIC - subspace pseudospectrum (super-resolution):\n");
        printf("  pseudospectrum peaks at f = %.3f and %.3f\n", f1, f2);
    }

    /* ESPRIT - direct frequency estimates, no spectral search. */
    double freqs[2];
    int ne = dsp_esprit(x, N, 2, 16, freqs);
    if (ne > 0) {
        printf("\nESPRIT - direct estimate, no grid search:\n");
        printf("  frequencies:");
        for (int i = 0; i < ne; ++i) printf(" %.3f", freqs[i]);
        printf("\n");
    }

    printf("\n-> every method recovers the 0.20 / 0.23 pair that the\n");
    printf("   classical FFT periodogram would blur into one lobe.\n");
}

/* ---- Time-frequency analysis ---------------------------------------- */

static void demo_timefreq(void) {
    section("TIME-FREQUENCY ANALYSIS (STFT, filter bank, Wigner-Ville)");

    /* A linear chirp: a single tone whose frequency sweeps upward.
     * The FFT alone would just show a smear of all those frequencies;
     * time-frequency methods reveal the frequency RISING with time. */
    enum { N = 512 };
    double x[N];
    for (int n = 0; n < N; ++n) {
        /* Instantaneous frequency sweeps 0.05 -> 0.40. */
        double phase = 2.0 * M_PI
                     * (0.05 * n + 0.35 * 0.5 * n * n / N);
        x[n] = cos(phase);
    }
    printf("Test signal: a linear chirp sweeping f = 0.05 -> 0.40.\n");
    printf("A single FFT shows only that all those frequencies exist;\n");
    printf("these methods show WHEN each one occurs.\n");

    /* --- STFT / spectrogram --- */
    dsp_stft s;
    dsp_stft_forward(x, N, DSP_WIN_HANNING, 64, 16, &s);
    double *spec = malloc(s.frames * s.bins * sizeof(double));
    dsp_spectrogram(&s, spec);

    printf("\nSTFT: %zu time frames x %zu frequency bins "
           "(64-pt Hann window, hop 16)\n", s.frames, s.bins);
    printf("  dominant frequency per frame (rising = chirp tracked):\n   ");
    for (size_t t = 0; t < s.frames; t += 4) {
        size_t pk = 1;
        double mx = 0.0;
        for (size_t f = 1; f < s.bins / 2; ++f)
            if (spec[t * s.bins + f] > mx) {
                mx = spec[t * s.bins + f];
                pk = f;
            }
        printf(" %.2f", (double)pk / (double)s.bins);
    }
    printf("\n");

    /* STFT is invertible: overlap-add reconstructs the signal. */
    size_t slen = dsp_stft_signal_len(&s);
    double *recon = malloc(slen * sizeof(double));
    dsp_stft_inverse(&s, recon);
    double err = 0.0;
    int cnt = 0;
    for (size_t i = 64; i < slen - 64 && i < (size_t)N; ++i) {
        err += fabs(recon[i] - x[i]);
        ++cnt;
    }
    printf("  overlap-add inverse: mean reconstruction error %.2e\n",
           cnt ? err / cnt : 0.0);

    /* --- QMF filter bank --- */
    dsp_qmf_bank bank;
    dsp_qmf_init(&bank);
    double *lo = malloc((N / 2) * sizeof(double));
    double *hi = malloc((N / 2) * sizeof(double));
    double *qrec = malloc(N * sizeof(double));
    dsp_qmf_analyze(&bank, x, N, lo, hi);
    dsp_qmf_synthesize(&bank, lo, hi, N / 2, qrec);

    /* Subband energies show how the chirp's energy splits low/high. */
    double e_lo = 0.0, e_hi = 0.0;
    for (int i = 0; i < N / 2; ++i) {
        e_lo += lo[i] * lo[i];
        e_hi += hi[i] * hi[i];
    }
    /* Reconstruction matches the input delayed by ntaps-1 samples. */
    size_t delay = bank.ntaps - 1;
    double qerr = 0.0;
    int qcnt = 0;
    for (size_t i = delay + 16; i + 16 < (size_t)N; ++i) {
        qerr += fabs(qrec[i] - x[i - delay]);
        ++qcnt;
    }
    printf("\nQMF filter bank: splits the signal into 2 subbands, "
           "each\n  decimated by 2 (critical sampling - no extra "
           "samples).\n");
    printf("  subband energy: low = %.0f, high = %.0f\n", e_lo, e_hi);
    printf("  analysis + synthesis reconstruction error: %.4f\n",
           qcnt ? qerr / qcnt : 0.0);
    printf("  -> the QMF mirror design cancels the aliasing from\n");
    printf("     decimation, the basis of subband audio coding.\n");

    /* --- Wigner-Ville --- */
    enum { NW = 128 };
    double xw[NW];
    for (int n = 0; n < NW; ++n) {
        double phase = 2.0 * M_PI
                     * (0.10 * n + 0.30 * 0.5 * n * n / NW);
        xw[n] = cos(phase);
    }
    double *wvd = malloc(NW * NW * sizeof(double));
    dsp_wigner_ville(xw, NW, wvd);

    printf("\nWigner-Ville distribution: the sharpest time-frequency\n");
    printf("  resolution of any method here, with no window needed.\n");
    printf("  chirp ridge - instantaneous frequency per time slice:\n   ");
    for (int t = 16; t < NW; t += 16) {
        int pk = 0;
        double mx = -1e30;
        for (int f = 0; f < NW; ++f)
            if (wvd[t * NW + f] > mx) {
                mx = wvd[t * NW + f];
                pk = f;
            }
        /* The WVD frequency axis spans [0, 0.5) over NW bins. */
        printf(" %.2f", (double)pk / (2.0 * NW));
    }
    printf("\n");
    printf("  -> the ridge tracks the sweep tightly; for a signal with\n");
    printf("     several components the WVD also shows spurious cross\n");
    printf("     terms midway between them - its defining drawback,\n");
    printf("     which the pseudo-WVD smooths away.\n");

    dsp_stft_free(&s);
    dsp_qmf_free(&bank);
    free(spec); free(recon);
    free(lo); free(hi); free(qrec); free(wvd);
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

/* ---- Adaptive filters ------------------------------------------------ */

/* Fill an array with pseudo-random +/-1 (or uniform) values. */
static unsigned adapt_seed = 24680;
static double adapt_rand(void) {
    adapt_seed = adapt_seed * 1103515245u + 12345u;
    return ((double)((adapt_seed >> 16) & 0xFFFF) / 65535.0) * 2.0 - 1.0;
}

static void demo_adaptive(void) {
    section("ADAPTIVE FILTERS (LMS, NLMS, RLS)");

    printf("Three algorithms on a convergence-vs-cost spectrum:\n");
    printf("  LMS  - O(L), simple, slow, step-size sensitive\n");
    printf("  NLMS - O(L), power-normalised step, easy to tune\n");
    printf("  RLS  - O(L^2), fast convergence, tracks change well\n");

    enum { L = 8, N = 1500 };

    /* --- Application 1: system identification ----------------------
     * An unknown FIR system; each filter must learn its taps from the
     * input/output pair alone. */
    double sys[L] = { 0.8, -0.5, 0.3, 0.25, -0.15, 0.1, -0.05, 0.02 };
    double x[N], d[N];
    for (int n = 0; n < N; ++n) x[n] = adapt_rand();
    for (int n = 0; n < N; ++n) {
        double a = 0.0;
        for (int k = 0; k < L; ++k)
            if (n >= k) a += sys[k] * x[n - k];
        d[n] = a;
    }

    dsp_lms  lms;  dsp_lms_init(&lms,  L, 0.05);
    dsp_nlms nlms; dsp_nlms_init(&nlms, L, 0.5, 1e-6);
    dsp_rls  rls;  dsp_rls_init(&rls,  L, 0.99, 100.0);

    /* Early MSE (first 50 samples) measures convergence speed. */
    double l_early = dsp_lms_train(&lms,  x, d, 50);
    double n_early = dsp_nlms_train(&nlms, x, d, 50);
    double r_early = dsp_rls_train(&rls,  x, d, 50);
    /* Train on the rest to let every filter settle. */
    dsp_lms_train(&lms,  x + 50, d + 50, N - 50);
    dsp_nlms_train(&nlms, x + 50, d + 50, N - 50);
    dsp_rls_train(&rls,  x + 50, d + 50, N - 50);

    printf("\n1. System identification - learn an unknown 8-tap filter\n");
    printf("   early MSE (first 50 samples, lower = faster convergence):\n");
    printf("     LMS  %.5f   NLMS %.5f   RLS %.5f\n",
           l_early, n_early, r_early);
    double cerr = 0.0;
    for (int k = 0; k < L; ++k) cerr += fabs(rls.weights[k] - sys[k]);
    printf("   RLS recovered the system's taps, total error %.5f\n", cerr);
    printf("   -> RLS converges fastest, NLMS next, LMS slowest.\n");

    dsp_lms_free(&lms); dsp_nlms_free(&nlms); dsp_rls_free(&rls);

    /* --- Application 2: channel equalization -----------------------
     * A multipath channel smears symbols; the adaptive filter learns
     * its inverse so the equalized output matches the sent symbols. */
    double tx[N], rx[N];
    for (int n = 0; n < N; ++n)
        tx[n] = (adapt_rand() >= 0.0) ? 1.0 : -1.0;
    rx[0] = tx[0];
    for (int n = 1; n < N; ++n)
        rx[n] = tx[n] + 0.6 * tx[n - 1];      /* 2-tap echo */

    dsp_nlms eq;
    dsp_nlms_init(&eq, 9, 0.5, 1e-6);
    double eq_early = dsp_nlms_train(&eq, rx, tx, N / 6);
    dsp_nlms_train(&eq, rx + N / 6, tx + N / 6, N - N / 3);
    double eq_late  = dsp_nlms_train(&eq, rx + 5 * N / 6,
                                     tx + 5 * N / 6, N / 6);
    printf("\n2. Channel equalization - undo a 2-tap multipath echo\n");
    printf("   NLMS equalizer MSE: %.5f (early) -> %.5f (converged)\n",
           eq_early, eq_late);
    printf("   -> the cleaned signal is what a demapper/FEC then sees.\n");
    dsp_nlms_free(&eq);

    /* --- Application 3: adaptive noise cancellation ----------------
     * A clean signal is buried under noise. A reference input carries
     * noise correlated with (but not equal to) the contaminating
     * noise. The adaptive filter shapes the reference to match the
     * noise, and the error signal is the recovered clean signal. */
    enum { NC_N = 4000 };
    double clean[NC_N], noise_ref[NC_N], primary[NC_N], recovered[NC_N];
    for (int n = 0; n < NC_N; ++n) {
        /* The signal we want back: a weak slow sine, buried in noise. */
        clean[n] = 0.3 * sin(2.0 * M_PI * 0.01 * n);
        noise_ref[n] = adapt_rand();
    }
    /* The primary input = clean signal + a filtered version of the
     * reference noise (the path from noise source to the mic). */
    for (int n = 0; n < NC_N; ++n) {
        double leaked = 0.7 * noise_ref[n]
                      + 0.3 * (n > 0 ? noise_ref[n - 1] : 0.0);
        primary[n] = clean[n] + leaked;
    }

    /* The filter adapts noise_ref toward the leaked noise; its error
     * (primary - output) is the recovered clean signal. */
    dsp_nlms nc;
    dsp_nlms_init(&nc, 4, 0.2, 1e-6);
    for (int n = 0; n < NC_N; ++n) {
        double err;
        dsp_nlms_update(&nc, noise_ref[n], primary[n], &err);
        recovered[n] = err;          /* error == recovered clean signal */
    }

    /* Compare noise power before and after, over the settled tail. */
    double in_noise = 0.0, out_noise = 0.0;
    for (int n = 3 * NC_N / 4; n < NC_N; ++n) {
        double before = primary[n]   - clean[n];
        double after  = recovered[n] - clean[n];
        in_noise  += before * before;
        out_noise += after  * after;
    }
    printf("\n3. Adaptive noise cancellation - recover a buried signal\n");
    printf("   residual noise power: %.4f (before) -> %.4f (after)\n",
           in_noise / (NC_N / 4), out_noise / (NC_N / 4));
    printf("   -> the filter models the noise path; the error signal\n");
    printf("      is the cleaned output (headphones, sensors).\n");
    dsp_nlms_free(&nc);
}

/* ---- Kalman filtering & sensor fusion ------------------------------- */

static unsigned kal_seed = 13579;
static double kal_gauss(void) {
    kal_seed = kal_seed * 1103515245u + 12345u;
    double u1 = ((kal_seed >> 16) & 0xFFFF) / 65535.0;
    kal_seed = kal_seed * 1103515245u + 12345u;
    double u2 = ((kal_seed >> 16) & 0xFFFF) / 65535.0;
    if (u1 < 1e-9) u1 = 1e-9;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static void demo_kalman(void) {
    section("KALMAN FILTERING & SENSOR FUSION");

    printf("The Kalman filter estimates the hidden state of a known\n");
    printf("dynamic system from noisy measurements, via a predict /\n");
    printf("update cycle. It is the optimal linear estimator and the\n");
    printf("standard tool for tracking and sensor fusion.\n");

    /* --- 1. Constant-velocity tracking --- */
    enum { N = 80 };
    double dt = 1.0, true_v = 2.5;
    dsp_kalman kf;
    dsp_kalman_tracker_init(&kf, 1, dt, 0.1, 4.0);
    kf.x[0] = 0.0;  kf.x[1] = 0.0;          /* initial pos, vel */
    kf.P[0] = 100.0; kf.P[3] = 100.0;       /* large initial doubt */

    double err_raw = 0.0, err_kf = 0.0;
    int cnt = 0;
    for (int n = 0; n < N; ++n) {
        double true_pos = true_v * dt * n;
        double meas = true_pos + 4.0 * kal_gauss();
        dsp_kalman_predict(&kf);
        dsp_kalman_update(&kf, &meas);
        if (n > 20) {                       /* after convergence */
            err_raw += (meas - true_pos) * (meas - true_pos);
            err_kf  += (kf.x[0] - true_pos) * (kf.x[0] - true_pos);
            ++cnt;
        }
    }
    printf("\n1. Constant-velocity tracking (position sensor, noise 4.0)\n");
    printf("   RMS error - raw measurement: %.2f, Kalman estimate: %.2f\n",
           sqrt(err_raw / cnt), sqrt(err_kf / cnt));
    printf("   recovered velocity %.2f (true %.1f) - never measured\n",
           kf.x[1], true_v);
    printf("   -> the filter smooths position AND infers the hidden\n");
    printf("      velocity from the position stream alone.\n");
    dsp_kalman_free(&kf);

    /* --- 2. Static sensor fusion --- */
    double meas[3]  = { 10.4, 9.6, 10.1 };
    double var[3]   = { 4.00, 1.00, 0.25 };   /* sensor 3 the best */
    double fused_var;
    double fused = dsp_kalman_fuse(meas, var, 3, &fused_var);
    printf("\n2. Sensor fusion - 3 sensors of differing accuracy\n");
    printf("   readings 10.4 / 9.6 / 10.1, variances 4.0 / 1.0 / 0.25\n");
    printf("   fused estimate %.3f, fused variance %.3f\n",
           fused, fused_var);
    printf("   -> inverse-variance weighting: the fused variance "
           "(%.3f)\n", fused_var);
    printf("      is below even the best single sensor (0.25).\n");

    /* --- 3. Kalman fusion of two sensors over time --- */
    dsp_kalman kf2;
    dsp_kalman_init(&kf2, 2, 2);             /* state [pos,vel], 2 meas */
    kf2.F[0] = 1.0; kf2.F[1] = 1.0; kf2.F[3] = 1.0;
    kf2.Q[0] = 0.01; kf2.Q[3] = 0.01;
    kf2.H[0] = 1.0; kf2.H[2] = 1.0;          /* both sense position */
    kf2.R[0] = 4.0; kf2.R[3] = 0.5;          /* sensor B 8x sharper */
    kf2.x[0] = 0.0; kf2.x[1] = 0.0;
    kf2.P[0] = 100.0; kf2.P[3] = 100.0;

    double err_a = 0.0, err_b = 0.0, err_f = 0.0;
    cnt = 0;
    for (int n = 0; n < N; ++n) {
        double tp = 2.0 * n;
        double za = tp + 2.00 * kal_gauss();   /* noisy sensor A */
        double zb = tp + 0.71 * kal_gauss();   /* sharp sensor B */
        double z[2] = { za, zb };
        dsp_kalman_predict(&kf2);
        dsp_kalman_update(&kf2, z);
        if (n > 20) {
            err_a += (za - tp) * (za - tp);
            err_b += (zb - tp) * (zb - tp);
            err_f += (kf2.x[0] - tp) * (kf2.x[0] - tp);
            ++cnt;
        }
    }
    printf("\n3. Kalman fusion of two position sensors over time\n");
    printf("   RMS error - sensor A: %.2f, sensor B: %.2f, fused: %.2f\n",
           sqrt(err_a / cnt), sqrt(err_b / cnt), sqrt(err_f / cnt));
    printf("   -> the fused track beats even the better sensor: the\n");
    printf("      update step weights each reading by its precision.\n");
    dsp_kalman_free(&kf2);

    printf("\nThe Extended Kalman Filter (dsp_ekf) handles nonlinear\n");
    printf("systems - radar range/bearing, GPS - by linearising the\n");
    printf("model at each step through caller-supplied Jacobians.\n");
}

/* ---- Interleaving ---------------------------------------------------- */

static void demo_interleaving(void) {
    section("INTERLEAVING (burst-error resilience)");

    /* Block interleaver: show how a burst gets spread out. */
    enum { R = 4, C = 6, LEN = R * C };
    uint8_t data[LEN], inter[LEN], deinter[LEN];
    for (int i = 0; i < LEN; ++i) data[i] = (uint8_t)i;

    dsp_block_interleave(data, inter, LEN, R, C);
    printf("Block interleaver, %dx%d matrix.\n", R, C);
    printf("  original order : ");
    for (int i = 0; i < LEN; ++i) printf("%2d ", data[i]);
    printf("\n  interleaved    : ");
    for (int i = 0; i < LEN; ++i) printf("%2d ", inter[i]);

    /* A 4-symbol burst hits the interleaved stream... */
    for (int i = 8; i < 12; ++i) inter[i] = 0xEE;
    dsp_block_deinterleave(inter, deinter, LEN, R, C);
    printf("\n  after a 4-symbol burst, de-interleaved error map:\n  ");
    for (int i = 0; i < LEN; ++i)
        printf("%s ", deinter[i] == data[i] ? " ." : " X");
    printf("\n-> the 4 errors are scattered, not clustered: each row\n");
    printf("   (codeword) now holds at most one.\n");

    /* End-to-end: interleaving rescues Reed-Solomon from a long burst. */
    printf("\nEnd-to-end test: Reed-Solomon RS(15,11), corrects t=2.\n");
    dsp_rs rs;
    dsp_rs_init(&rs, 4);                       /* 4 parity -> t = 2 */

    /* Build 4 codewords of 15 symbols each = 60 symbols. */
    enum { NCW = 4, N = 15, KK = 11, TOTAL = NCW * N };
    uint8_t cw[TOTAL];
    for (int w = 0; w < NCW; ++w) {
        uint8_t msg[KK], par[4];
        for (int i = 0; i < KK; ++i) msg[i] = (uint8_t)(w * 20 + i);
        dsp_rs_encode(&rs, msg, KK, par);
        for (int i = 0; i < KK; ++i) cw[w * N + i] = msg[i];
        for (int i = 0; i < 4;  ++i) cw[w * N + KK + i] = par[i];
    }

    /* Case A: a 6-symbol burst with NO interleaving. It all lands in
     * one codeword -> 6 errors > t=2 -> that codeword is lost. */
    uint8_t plain[TOTAL];
    for (int i = 0; i < TOTAL; ++i) plain[i] = cw[i];
    for (int i = 16; i < 22; ++i) plain[i] ^= 0x5A;   /* burst in cw 1 */
    int lost = 0;
    for (int w = 0; w < NCW; ++w)
        if (dsp_rs_decode(&rs, plain + w * N, N) < 0) ++lost;
    printf("  without interleaving: 6-symbol burst -> %d codeword(s) "
           "uncorrectable\n", lost);

    /* Case B: same burst, but interleave before the channel and
     * de-interleave after. The matrix has NCW rows (one per codeword)
     * and N columns, so the 6 errors spread one-per-codeword. */
    uint8_t tx[TOTAL], rx[TOTAL];
    dsp_block_interleave(cw, tx, TOTAL, NCW, N);   /* 4x15 matrix */
    for (int i = 16; i < 22; ++i) tx[i] ^= 0x5A;   /* same 6-burst */
    dsp_block_deinterleave(tx, rx, TOTAL, NCW, N);
    int recovered = 1;
    for (int w = 0; w < NCW; ++w)
        if (dsp_rs_decode(&rs, rx + w * N, N) < 0) recovered = 0;
    int data_ok = 1;
    for (int w = 0; w < NCW; ++w)
        for (int i = 0; i < KK; ++i)
            if (rx[w * N + i] != cw[w * N + i]) data_ok = 0;
    printf("  with interleaving   : same burst -> %s\n",
           (recovered && data_ok)
               ? "all 4 codewords recovered"
               : "decode failed");
    printf("-> interleaving spreads the burst so each RS codeword\n");
    printf("   sees <= 2 errors, within its correction limit.\n");

    /* Convolutional interleaver: continuous-stream alternative. */
    dsp_conv_interleaver itl, dtl;
    dsp_conv_interleaver_init(&itl, 4, 2);
    dsp_conv_deinterleaver_init(&dtl, 4, 2);
    printf("\nConvolutional interleaver (4 branches, depth 2):\n");
    printf("  end-to-end latency: %zu symbols\n",
           dsp_conv_interleaver_latency(4, 2));
    printf("-> streams continuously, no block boundary; about half the\n");
    printf("   memory and latency of an equivalent block interleaver.\n");
    dsp_conv_interleaver_free(&itl);
    dsp_conv_interleaver_free(&dtl);
}

/* ---- LDPC -------------------------------------------------------------- */

static void demo_ldpc(void) {
    section("LDPC CODES (near-Shannon-limit FEC)");

    /* A small regular (wc=3, wr=4) LDPC code: n=12 bits, m=9 checks. */
    dsp_ldpc code;
    if (dsp_ldpc_make_regular(&code, 9, 12, 3, 4, 1) != 0) {
        printf("  (code construction failed)\n");
        return;
    }
    printf("Regular LDPC code: 12 variable nodes, 9 check nodes.\n");
    printf("  every bit is in 3 checks; every check covers 4 bits.\n");
    printf("  the sparse parity-check matrix H is a bipartite Tanner "
           "graph.\n");

    /* The all-zero word is always a valid codeword - use it as the
     * transmitted codeword so the demo stays self-contained. */
    uint8_t cw[12] = {0};
    printf("\nTransmitted codeword satisfies every parity check: %s\n",
           dsp_ldpc_check(&code, cw) ? "yes" : "no");

    /* --- Bit-flipping decoder (hard decision) --- */
    uint8_t hard[12];
    for (int i = 0; i < 12; ++i) hard[i] = cw[i];
    hard[5] ^= 1;                         /* one received bit flipped */
    printf("\nBit-flipping decoder (hard decision):\n");
    printf("  received word has syndrome weight %zu (a check fails)\n",
           dsp_ldpc_syndrome_weight(&code, hard));
    int bf_iter = dsp_ldpc_decode_bitflip(&code, hard, 40);
    int bf_ok = 1;
    for (int i = 0; i < 12; ++i) if (hard[i] != cw[i]) bf_ok = 0;
    printf("  flipped the most-suspect bit each pass -> converged in "
           "%d iteration(s)  %s\n", bf_iter,
           bf_ok ? "(error corrected)" : "(FAIL)");

    /* --- Sum-product decoder (soft decision) --- */
    /* Channel LLRs: positive favours bit 0. Two bits arrive with the
     * wrong sign - the demodulator was fooled - but with low
     * confidence. Bit-flipping would see two hard errors; the
     * sum-product decoder weighs the analog confidence instead. */
    double llr[12];
    for (int i = 0; i < 12; ++i) llr[i] = 3.0;   /* clean, favour 0 */
    llr[2] = -1.0;                                /* wrong-sign, weak */
    llr[8] = -0.8;                                /* wrong-sign, weak */

    uint8_t soft[12];
    int sp_iter = dsp_ldpc_decode_sumproduct(&code, llr, soft, 50);
    int sp_ok = 1;
    for (int i = 0; i < 12; ++i) if (soft[i] != cw[i]) sp_ok = 0;
    printf("\nSum-product decoder (soft decision, belief propagation):\n");
    printf("  variable and check nodes exchange probabilities (LLRs)\n");
    printf("  along the Tanner-graph edges, refining belief each pass\n");
    printf("  2 wrong-sign inputs -> converged in %d iteration(s)  %s\n",
           sp_iter, sp_ok ? "(both errors corrected)" : "(FAIL)");

    /* --- Min-sum decoder (the hardware approximation) --- */
    uint8_t ms[12];
    int ms_iter = dsp_ldpc_decode_minsum(&code, llr, ms, 0.75, 50);
    int ms_ok = 1;
    for (int i = 0; i < 12; ++i) if (ms[i] != cw[i]) ms_ok = 0;
    printf("\nMin-sum decoder (normalised, scale 0.75):\n");
    printf("  replaces the tanh rule with a plain minimum-of-magnitudes\n");
    printf("  -> no transcendental functions, fixed-point friendly\n");
    printf("  same 2 wrong-sign inputs -> converged in %d iteration(s)  %s\n",
           ms_iter, ms_ok ? "(both errors corrected)" : "(FAIL)");

    /* --- BER sweep: the three decoders over a noisy channel --- */
    printf("\nBit-error rate over a simulated AWGN channel "
           "(2000 codewords each):\n");
    printf("  noise |  bit-flip | sum-product |   min-sum\n");
    double noise_levels[3] = { 0.5, 0.7, 0.9 };
    for (int k = 0; k < 3; ++k) {
        double ns = noise_levels[k];
        double bf = dsp_ldpc_ber_sweep(&code, DSP_LDPC_BITFLIP,
                                       ns, 2000, 50, 2024);
        double sp = dsp_ldpc_ber_sweep(&code, DSP_LDPC_SUMPRODUCT,
                                       ns, 2000, 50, 2024);
        double mn = dsp_ldpc_ber_sweep(&code, DSP_LDPC_MINSUM,
                                       ns, 2000, 50, 2024);
        printf("   %.1f  | %9.4f | %11.4f | %9.4f\n", ns, bf, sp, mn);
    }
    printf("-> sum-product and min-sum track each other closely, and\n");
    printf("   both crush the hard-decision bit-flipper. Min-sum gets\n");
    printf("   that performance without a single tanh - which is why\n");
    printf("   it is the decoder in real 5G and Wi-Fi 6 silicon.\n");

    dsp_ldpc_free(&code);
}

/* ---- OFDM and coded OFDM ---------------------------------------------- */

static void demo_ofdm(void) {
    section("OFDM (multicarrier modulation)");

    /* QAM: map bits onto a constellation. */
    printf("QAM constellations carry log2(M) bits per symbol:\n");
    dsp_qam_order orders[3] = { DSP_QAM_QPSK, DSP_QAM_16, DSP_QAM_64 };
    const char *qnames[3]   = { "QPSK", "16-QAM", "64-QAM" };
    for (int k = 0; k < 3; ++k)
        printf("  %-7s : %zu bits/symbol\n", qnames[k],
               dsp_qam_bits_per_symbol(orders[k]));

    /* OFDM round-trip through a multipath channel. */
    enum { NFFT = 64, CP = 16 };
    dsp_ofdm ofdm;
    dsp_ofdm_init(&ofdm, NFFT, CP);

    cplx freq[NFFT], tx[NFFT + CP], rx[NFFT + CP], rxf[NFFT], cfr[NFFT];
    unsigned seed = 4321;
    for (int i = 0; i < NFFT; ++i) {
        seed = seed * 1103515245u + 12345u;
        double re = ((seed >> 16) & 1) ? 0.707 : -0.707;
        double im = ((seed >> 17) & 1) ? 0.707 : -0.707;
        freq[i] = dsp_cplx(re, im);          /* one QPSK symbol/carrier */
    }

    printf("\nOFDM symbol: %d subcarriers, %d-sample cyclic prefix.\n",
           NFFT, CP);
    printf("  modulator = IFFT + cyclic prefix; demodulator = "
           "remove prefix + FFT\n");

    /* A 3-tap multipath channel, noiseless, to isolate the equalizer. */
    cplx taps[3] = { dsp_cplx(1.0, 0.0),
                     dsp_cplx(0.3, 0.1),
                     dsp_cplx(0.1, -0.05) };
    dsp_channel ch;
    dsp_channel_init(&ch, taps, 3, 0.0, 1);

    dsp_ofdm_modulate(&ofdm, freq, tx);
    dsp_channel_apply(&ch, tx, rx, NFFT + CP);
    dsp_ofdm_demodulate(&ofdm, rx, rxf);

    /* Before equalization the multipath has rotated every subcarrier. */
    double err_before = 0.0;
    for (int i = 0; i < NFFT; ++i) {
        double d = cabs(rxf[i] - freq[i]);
        if (d > err_before) err_before = d;
    }

    /* One complex divide per subcarrier undoes the channel. */
    dsp_channel_frequency_response(&ch, cfr, NFFT);
    dsp_ofdm_equalize(&ofdm, rxf, cfr);
    double err_after = 0.0;
    for (int i = 0; i < NFFT; ++i) {
        double d = cabs(rxf[i] - freq[i]);
        if (d > err_after) err_after = d;
    }

    printf("  max subcarrier error before equalization: %.3f\n",
           err_before);
    printf("  max subcarrier error after  equalization: %.2e\n",
           err_after);
    printf("-> the cyclic prefix turns multipath into a per-subcarrier\n");
    printf("   complex gain, so equalization is one divide per carrier.\n");

    dsp_channel_free(&ch);
}

static void demo_coded_ofdm(void) {
    section("CODED OFDM (FEC + OFDM end-to-end)");

    /* nfft=64 subcarriers, QPSK -> 128 coded bits per OFDM symbol, so
     * the LDPC codeword length must be 128. */
    enum { NFFT = 64, CP = 16 };
    dsp_ldpc code;
    if (dsp_ldpc_make_regular(&code, 64, 128, 3, 6, 1) != 0) {
        printf("  (LDPC construction failed)\n");
        return;
    }

    dsp_ofdm ofdm;
    dsp_ofdm_init(&ofdm, NFFT, CP);
    dsp_coded_ofdm cfg = { .ofdm = ofdm,
                           .order = DSP_QAM_QPSK,
                           .code = &code };

    printf("Full chain: LDPC encode -> QPSK map -> OFDM (IFFT+CP)\n");
    printf("            -> multipath+AWGN channel -> OFDM demod\n");
    printf("            -> equalize -> soft demap -> LDPC decode\n\n");

    cplx taps[3] = { dsp_cplx(1.0, 0.0),
                     dsp_cplx(0.25, 0.1),
                     dsp_cplx(0.1, -0.05) };

    printf("Bit-error rate over a 3-tap multipath channel "
           "(200 OFDM frames):\n");
    printf("  noise |  raw BER  | coded BER\n");
    double levels[4] = { 0.3, 0.5, 0.7, 0.9 };
    for (int k = 0; k < 4; ++k) {
        dsp_channel ch;
        dsp_channel_init(&ch, taps, 3, levels[k], 777);

        size_t raw_tot = 0, cod_tot = 0;
        size_t nbits = NFFT * 2;             /* QPSK: 2 bits/subcarrier */
        for (int f = 0; f < 200; ++f) {
            size_t be, re;
            dsp_coded_ofdm_run_frame(&cfg, &ch, 50, &be, &re);
            raw_tot += re;
            cod_tot += be;
        }
        printf("   %.1f  | %9.4f | %9.4f\n", levels[k],
               (double)raw_tot / (200.0 * nbits),
               (double)cod_tot / (200.0 * nbits));

        dsp_channel_free(&ch);
    }
    printf("-> 'raw' is the error rate straight off the QAM demapper;\n");
    printf("   'coded' is after LDPC. The code cleans up the residual\n");
    printf("   noise errors - this FEC + OFDM pairing is exactly what\n");
    printf("   carries data in Wi-Fi, LTE, and 5G.\n");

    dsp_ldpc_free(&code);
}

/* ---- Pulse shaping and synchronization -------------------------------- */

static void demo_sync(void) {
    section("PULSE SHAPING & SYNCHRONIZATION");

    /* Root-raised-cosine pulse shaping. */
    enum { NSYM = 32, SPS = 8, NTAPS = 65 };
    double rrc[NTAPS];
    dsp_rrc_design(rrc, NTAPS, SPS, 0.25);
    double energy = 0.0;
    for (int i = 0; i < NTAPS; ++i) energy += rrc[i] * rrc[i];
    printf("Root-raised-cosine filter: %d taps, %d samples/symbol, "
           "roll-off 0.25\n", NTAPS, SPS);
    printf("  tap energy = %.4f (normalised to 1)\n", energy);
    printf("  shapes each symbol into a bandlimited, ISI-free pulse.\n");

    /* Shape a QPSK symbol stream, then matched-filter it back. */
    cplx syms[NSYM], shaped[NSYM * SPS], matched[NSYM * SPS];
    unsigned s = 31;
    for (int k = 0; k < NSYM; ++k) {
        s = s * 1103515245u + 12345u;
        syms[k] = dsp_cplx(((s >> 16) & 1) ? 0.707 : -0.707,
                           ((s >> 17) & 1) ? 0.707 : -0.707);
    }
    dsp_pulse_shape(syms, NSYM, SPS, rrc, NTAPS, shaped);
    dsp_matched_filter(shaped, NSYM * SPS, rrc, NTAPS, matched);
    printf("  %d symbols -> %d shaped samples -> matched filtered\n",
           NSYM, NSYM * SPS);

    /* Carrier recovery: lock out a fixed phase offset. */
    cplx rotated[200], recovered[200];
    for (int i = 0; i < 200; ++i) {
        cplx sym = dsp_cplx((i & 1) ? 0.707 : -0.707,
                            (i & 2) ? 0.707 : -0.707);
        rotated[i] = sym * cexp(0.4 * I);    /* 0.4 rad carrier offset */
    }
    dsp_carrier_pll pll;
    dsp_carrier_pll_init(&pll, 0.05, 0.0025);
    dsp_carrier_recover(&pll, rotated, recovered, 200);

    double err_start = 0.0, err_end = 0.0;
    for (int i = 0; i < 20; ++i) {
        cplx sym = dsp_cplx((i & 1) ? 0.707 : -0.707,
                            (i & 2) ? 0.707 : -0.707);
        err_start += cabs(recovered[i] - sym);
    }
    for (int i = 180; i < 200; ++i) {
        cplx sym = dsp_cplx((i & 1) ? 0.707 : -0.707,
                            (i & 2) ? 0.707 : -0.707);
        err_end += cabs(recovered[i] - sym);
    }
    printf("\nCarrier-recovery PLL locking out a 0.4 rad offset:\n");
    printf("  constellation error, first 20 symbols : %.3f\n",
           err_start / 20.0);
    printf("  constellation error, last 20 symbols  : %.3f\n",
           err_end / 20.0);
    printf("-> the loop tracks the offset and pulls the error to ~0;\n");
    printf("   timing recovery (Gardner detector) locks the sampling\n");
    printf("   instant the same way, with a feedback loop.\n");
}

/* ---- Image processing ------------------------------------------------- */

static void demo_image(void) {
    section("IMAGE PROCESSING (2-D DSP)");

    /* Build a synthetic test image: a smooth gradient with a bright
     * square, so the filters and transforms have structure to act on. */
    enum { W = 64, H = 64 };
    dsp_image img;
    dsp_image_alloc(&img, W, H);
    for (size_t y = 0; y < H; ++y)
        for (size_t x = 0; x < W; ++x) {
            double v = 40.0 + 120.0 * ((double)x / W);
            if (x > 20 && x < 44 && y > 20 && y < 44)
                v = 220.0;                   /* a bright square */
            img.data[y * W + x] = v;
        }
    printf("Synthetic %dx%d test image: gradient + a bright square.\n",
           W, H);

    /* 2-D FFT round-trip. */
    double *re = malloc(W * H * sizeof(double));
    double *im = malloc(W * H * sizeof(double));
    dsp_image recon;
    dsp_image_alloc(&recon, W, H);
    dsp_fft2d(&img, re, im);
    dsp_ifft2d(re, im, &recon);
    double fft_err = 0.0;
    for (size_t i = 0; i < W * H; ++i) {
        double d = fabs(img.data[i] - recon.data[i]);
        if (d > fft_err) fft_err = d;
    }
    printf("\n2-D FFT (separable: 1-D FFT on rows, then columns):\n");
    printf("  round-trip max error: %.2e\n", fft_err);

    /* 2-D DCT energy compaction. */
    double *coeffs = malloc(W * H * sizeof(double));
    dsp_dct2d(&img, coeffs);
    double total = 0.0, corner = 0.0;
    for (size_t y = 0; y < H; ++y)
        for (size_t x = 0; x < W; ++x) {
            double e = coeffs[y * W + x] * coeffs[y * W + x];
            total += e;
            if (x < 8 && y < 8) corner += e;     /* low-frequency 8x8 */
        }
    printf("\n2-D DCT (the JPEG transform):\n");
    printf("  energy in the low-frequency 8x8 corner: %.2f%%\n",
           100.0 * corner / total);
    printf("  -> compaction into few coefficients is why JPEG works.\n");

    /* Spatial filters. */
    dsp_image blurred, edges;
    dsp_image_alloc(&blurred, W, H);
    dsp_image_alloc(&edges, W, H);
    dsp_gaussian_blur(&img, &blurred, 1.5);
    dsp_sobel(&img, &edges);

    /* The square's border should light up under Sobel. */
    double edge_on  = edges.data[32 * W + 20];   /* on the square edge */
    double edge_off = edges.data[10 * W + 10];   /* smooth gradient    */
    printf("\nSpatial filters (2-D convolution):\n");
    printf("  Gaussian blur: smooths detail (unit-gain kernel)\n");
    printf("  Sobel edge magnitude at the square's border: %.0f\n",
           edge_on);
    printf("  Sobel response in the smooth gradient      : %.1f\n",
           edge_off);

    /* Median filter removes salt noise a blur cannot. */
    dsp_image noisy, cleaned;
    dsp_image_copy(&noisy, &img);
    noisy.data[30 * W + 30] = 255.0;             /* an impulse */
    noisy.data[40 * W + 15] = 0.0;
    dsp_image_alloc(&cleaned, W, H);
    dsp_median_filter(&noisy, &cleaned, 3);
    printf("\nMedian filter (non-linear): erases salt-and-pepper\n");
    printf("  spike pixel was 255, after median: %.0f\n",
           cleaned.data[30 * W + 30]);

    /* Otsu automatic thresholding. */
    double t = dsp_threshold_otsu(&img);
    printf("\nOtsu automatic threshold: %.0f (no manual tuning)\n", t);

    /* 2-D wavelet: split into four resolution quadrants. */
    dsp_image wav;
    dsp_image_copy(&wav, &img);
    dsp_dwt2d_haar(&wav);
    dsp_idwt2d_haar(&wav);
    double wav_err = 0.0;
    for (size_t i = 0; i < W * H; ++i) {
        double d = fabs(img.data[i] - wav.data[i]);
        if (d > wav_err) wav_err = d;
    }
    printf("\n2-D Haar wavelet (the JPEG 2000 transform):\n");
    printf("  splits the image into approximation + 3 detail "
           "quadrants\n");
    printf("  round-trip max error: %.2e\n", wav_err);

    free(re); free(im); free(coeffs);
    dsp_image_free(&img);   dsp_image_free(&recon);
    dsp_image_free(&blurred); dsp_image_free(&edges);
    dsp_image_free(&noisy); dsp_image_free(&cleaned);
    dsp_image_free(&wav);
}

/* ---- Array signal processing ---------------------------------------- */

static void demo_array(void) {
    section("ARRAY PROCESSING (beamforming & DOA)");

    /* An 8-element uniform linear array, half-wavelength spacing.
     * Two plane-wave sources arrive from -20 and +30 degrees. */
    enum { M = 8, T = 300, NANG = 361 };
    double deg = M_PI / 180.0;
    double src[2] = { -20.0 * deg, 30.0 * deg };

    cplx *snap = malloc(T * M * sizeof(cplx));
    cplx *R    = malloc(M * M * sizeof(cplx));
    dsp_array_synthesize(M, 0.5, T, src, 2, 0.2, 2024, snap);
    dsp_array_covariance(snap, M, T, R);

    printf("8-sensor uniform linear array, half-wavelength spacing.\n");
    printf("Two sources arrive from -20 and +30 degrees.\n");
    printf("The array's job: find those directions from the sensor\n");
    printf("covariance alone.\n");

    /* Helper to report the strongest peaks of an angle spectrum. */
    double *pc = malloc(NANG * sizeof(double));
    double *pm = malloc(NANG * sizeof(double));
    double *pu = malloc(NANG * sizeof(double));

    /* --- Conventional (delay-and-sum) beamformer --- */
    dsp_beamform_conventional(R, M, 0.5, pc, NANG);
    size_t cpk = 1;
    for (size_t i = 1; i + 1 < NANG; ++i)
        if (pc[i] > pc[cpk]) cpk = i;
    printf("\nConventional (delay-and-sum) beamformer:\n");
    printf("  steers the array across all angles, measures power\n");
    printf("  strongest peak at %.0f degrees\n",
           -90.0 + 180.0 * (double)cpk / (NANG - 1));
    printf("  -> finds a source, but its wide beam blurs the two.\n");

    /* --- MVDR / Capon beamformer --- */
    dsp_beamform_mvdr(R, M, 0.5, pm, NANG);
    double mmax = 0.0;
    for (size_t i = 0; i < NANG; ++i)
        if (pm[i] > mmax) mmax = pm[i];
    printf("\nMVDR / Capon beamformer (minimum variance):\n");
    printf("  keeps unit gain on the look angle, nulls everything\n");
    printf("  else - far sharper. Peaks at:");
    for (size_t i = 1; i + 1 < NANG; ++i)
        if (pm[i] > pm[i - 1] && pm[i] > pm[i + 1] && pm[i] > 0.2 * mmax)
            printf(" %.0f", -90.0 + 180.0 * (double)i / (NANG - 1));
    printf(" degrees\n");

    /* --- Spatial MUSIC --- */
    dsp_doa_music(R, M, 0.5, 2, pu, NANG);
    double umax = 0.0;
    for (size_t i = 0; i < NANG; ++i)
        if (pu[i] > umax) umax = pu[i];
    printf("\nSpatial MUSIC (subspace super-resolution):\n");
    printf("  eigen-splits the covariance into signal + noise\n");
    printf("  subspaces. Pseudospectrum peaks at:");
    for (size_t i = 1; i + 1 < NANG; ++i)
        if (pu[i] > pu[i - 1] && pu[i] > pu[i + 1] && pu[i] > 0.1 * umax)
            printf(" %.0f", -90.0 + 180.0 * (double)i / (NANG - 1));
    printf(" degrees\n");

    /* --- Spatial ESPRIT --- */
    double doa[2];
    int ne = dsp_doa_esprit(R, M, 0.5, 2, doa);
    printf("\nSpatial ESPRIT (direct, no angle search):\n");
    printf("  exploits the array's shift-invariance. Angles:");
    for (int i = 0; i < ne; ++i)
        printf(" %.0f", doa[i] / deg);
    printf(" degrees\n");

    printf("\n-> MUSIC and ESPRIT are the same algorithms used for\n");
    printf("   temporal frequency estimation (see spectral/): a\n");
    printf("   direction and a frequency are both the phase slope of\n");
    printf("   a complex exponential.\n");

    free(snap); free(R);
    free(pc); free(pm); free(pu);
}

int main(void) {
    printf("DSP GUIDE - annotated demo\n");
    printf("C implementation of common digital signal processing algorithms.\n");

    demo_transforms();
    demo_filtering();
    demo_operations();
    demo_spectral();
    demo_estimation();
    demo_timefreq();
    demo_sampling();
    demo_wavelet();
    demo_detection();
    demo_correction();
    demo_ldpc();
    demo_interleaving();
    demo_adaptive();
    demo_kalman();
    demo_ofdm();
    demo_coded_ofdm();
    demo_sync();
    demo_image();
    demo_array();

    printf("\nAll demos complete.\n");
    return 0;
}
