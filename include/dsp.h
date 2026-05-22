/*
 * dsp.h - Master header for the DSP Guide
 *
 * Include this single header to access every transform, filter,
 * signal operation, window, sample-rate converter, and the wavelet
 * transform implemented in this project.
 *
 * Organisation mirrors the Fiveable "Common Digital Signal Processing
 * Algorithms" guide:
 *
 *   transforms/  - DFT, FFT, DCT          (time <-> frequency)
 *   filtering/   - FIR, IIR               (frequency-selective filters)
 *   operations/  - convolution, correlation
 *   spectral/    - window functions, advanced spectral estimation
 *                  (AR, ARMA, MUSIC, ESPRIT), and time-frequency
 *                  analysis (STFT, QMF filter bank, Wigner-Ville)
 *   sampling/    - decimation, interpolation, rational resampling
 *   wavelet/     - discrete wavelet transform (multi-resolution)
 *   coding/      - error detection (CRC, parity, checksum) and
 *                  forward error correction (Hamming, Reed-Solomon,
 *                  convolutional/Viterbi, LDPC), with interleaving
 *   adaptive/    - adaptive filters (LMS, NLMS, RLS) for equalization,
 *                  system identification, and noise cancellation
 *   modulation/  - QAM, OFDM, channel model, the coded-OFDM
 *                  transceiver, pulse shaping, and carrier/timing
 *                  recovery
 *   image/       - grayscale image type, 2-D FFT and DCT, spatial
 *                  filtering, point operators, and the 2-D wavelet
 *   array/       - array signal processing: beamforming and
 *                  direction-of-arrival estimation
 *
 * Every module is annotated with: what problem it solves, its
 * computational complexity, and the trade-offs versus its alternatives.
 */
#ifndef DSP_H
#define DSP_H

#include "common.h"

#include "transforms/dft.h"
#include "transforms/fft.h"
#include "transforms/dct.h"

#include "filtering/fir.h"
#include "filtering/iir.h"

#include "operations/convolution.h"
#include "operations/correlation.h"

#include "spectral/window.h"
#include "spectral/estimation.h"
#include "spectral/timefreq.h"

#include "sampling/resample.h"

#include "wavelet/wavelet.h"

#include "coding/detect.h"
#include "coding/correct.h"
#include "coding/ldpc.h"
#include "coding/interleave.h"

#include "adaptive/adaptive.h"

#include "modulation/qam.h"
#include "modulation/channel.h"
#include "modulation/ofdm.h"
#include "modulation/coded_ofdm.h"
#include "modulation/pulse.h"
#include "modulation/sync.h"

#include "image/image_all.h"

#include "array/array.h"

#endif /* DSP_H */
