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
 *   spectral/    - window functions       (leakage reduction)
 *   sampling/    - decimation, interpolation, rational resampling
 *   wavelet/     - discrete wavelet transform (multi-resolution)
 *   coding/      - error detection (CRC, parity, checksum), forward
 *                  error correction (Hamming, Reed-Solomon,
 *                  convolutional/Viterbi), and channel equalization
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

#include "sampling/resample.h"

#include "wavelet/wavelet.h"

#include "coding/detect.h"
#include "coding/correct.h"
#include "coding/equalize.h"

#endif /* DSP_H */
