/*
 * dct.h - Discrete Cosine Transform (DCT-II) and its inverse (DCT-III)
 *
 * PROBLEM SOLVED
 *   Decomposes a real signal using only cosine basis functions,
 *   producing purely real coefficients. Its defining property is
 *   strong ENERGY COMPACTION: most of a typical signal's energy lands
 *   in a few low-frequency coefficients.
 *
 * WHY IT MATTERS
 *   That compaction is exactly what compression wants. JPEG transforms
 *   each 8x8 pixel block with a 2-D DCT, then quantises away the many
 *   near-zero high-frequency coefficients. MPEG video works similarly.
 *
 * vs DFT
 *   Both move a signal to the frequency domain. The DCT uses cosines
 *   only and stays real for real input. It also handles block edges
 *   better - the cosine basis implies a symmetric extension, avoiding
 *   the discontinuity that causes DFT leakage.
 *
 * COMPLEXITY
 *   This is the direct O(N^2) definition (clear, not fast). Fast
 *   O(N log N) DCTs exist via an FFT of a reordered/extended sequence.
 */
#ifndef DSP_DCT_H
#define DSP_DCT_H

#include "../common.h"

/* Forward DCT-II: real `in` (length n) -> real coefficients `out`. */
void dsp_dct(const double *in, double *out, size_t n);

/* Inverse DCT (DCT-III): coefficients `in` -> reconstructed `out`. */
void dsp_idct(const double *in, double *out, size_t n);

#endif /* DSP_DCT_H */
