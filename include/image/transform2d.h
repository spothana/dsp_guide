/*
 * transform2d.h - 2-D frequency transforms for images
 *
 * PROBLEM SOLVED
 *   The transforms in transforms/ are 1-D. Images are 2-D, so image
 *   processing needs 2-D versions. The key fact that makes this easy:
 *   the 2-D FFT and 2-D DCT are SEPARABLE - a 2-D transform is just a
 *   1-D transform applied to every row, then to every column. No new
 *   mathematics, only bookkeeping.
 *
 * 2-D FFT
 *   Reveals the spatial-frequency content of an image: smooth regions
 *   concentrate energy near the centre (low frequency), edges and
 *   texture spread it outward (high frequency). The basis of
 *   frequency-domain filtering, and of fast 2-D convolution.
 *
 * 2-D DCT
 *   The transform behind JPEG. JPEG splits an image into 8x8 blocks
 *   and applies a 2-D DCT to each; the strong energy compaction packs
 *   most of a block's energy into a few low-frequency coefficients,
 *   which can then be quantised cheaply. This module provides both a
 *   whole-image 2-D DCT and the 8x8 block transform JPEG uses.
 */
#ifndef DSP_TRANSFORM2D_H
#define DSP_TRANSFORM2D_H

#include "../common.h"
#include "image/image.h"

/*
 * 2-D FFT of a real image.
 *   img  : input image (width and height must both be powers of two)
 *   re   : output real part, width*height, row-major
 *   im   : output imaginary part, width*height, row-major
 * Returns 0 on success, -1 if a dimension is not a power of two.
 */
int dsp_fft2d(const dsp_image *img, double *re, double *im);

/*
 * Inverse 2-D FFT back to a real image.
 *   re, im : input spectrum, width*height each
 *   img    : output image (must be pre-allocated to width x height)
 * Returns 0 on success, -1 on a dimension error.
 */
int dsp_ifft2d(const double *re, const double *im, dsp_image *img);

/*
 * 2-D DCT of a whole image (separable DCT-II).
 *   img    : input image
 *   coeffs : output coefficients, width*height, row-major
 * Any dimensions are accepted (this uses the direct O(N^2) 1-D DCT
 * per row and column).
 */
void dsp_dct2d(const dsp_image *img, double *coeffs);

/*
 * Inverse 2-D DCT.
 *   coeffs : input coefficients, width*height
 *   img    : output image (pre-allocated to width x height)
 */
void dsp_idct2d(const double *coeffs, dsp_image *img);

/*
 * Forward 8x8 block DCT - the exact transform JPEG applies to each
 * pixel block. `block` and `out` are 64-element row-major arrays.
 */
void dsp_dct8x8(const double *block, double *out);

/* Inverse 8x8 block DCT. */
void dsp_idct8x8(const double *coeffs, double *block);

#endif /* DSP_TRANSFORM2D_H */
