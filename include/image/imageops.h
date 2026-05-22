/*
 * imageops.h - Point operators and 2-D wavelet for images
 *
 * This module collects the image operations that are not spatial
 * convolutions:
 *
 *   POINT OPERATORS map each pixel through a function of its own
 *   value alone (no neighbourhood):
 *     - histogram equalization spreads the pixel values to use the
 *       full tonal range, improving contrast in a washed-out image;
 *     - thresholding converts a grayscale image to pure black/white.
 *
 *   The 2-D WAVELET TRANSFORM is the multi-resolution analysis of
 *   wavelet/ extended to images. It is separable - a 1-D wavelet step
 *   on every row, then on every column - and splits an image into
 *   four quadrants: a half-size approximation plus horizontal,
 *   vertical, and diagonal detail. JPEG 2000 uses the 2-D wavelet in
 *   place of the DCT, because its multi-resolution structure
 *   compresses and scales better.
 */
#ifndef DSP_IMAGEOPS_H
#define DSP_IMAGEOPS_H

#include "image/image.h"

/*
 * Histogram equalization.
 *   img : image to equalize, modified in place
 * Builds the 256-bin intensity histogram, forms its cumulative
 * distribution, and remaps each pixel through that CDF so the output
 * intensities are spread as evenly as possible across 0..255.
 */
void dsp_histogram_equalize(dsp_image *img);

/*
 * Compute the 256-bin intensity histogram of an image.
 *   img  : input image
 *   hist : output array of 256 counts
 */
void dsp_histogram(const dsp_image *img, size_t hist[256]);

/*
 * Fixed-threshold binarization.
 *   src    : input image
 *   dst    : output image (pre-allocated, same size)
 *   thresh : pixels >= thresh become 255 (white), the rest 0 (black)
 */
void dsp_threshold(const dsp_image *src, dsp_image *dst,
                   double thresh);

/*
 * Otsu's automatic threshold. Returns the threshold value (0..255)
 * that best separates the histogram into two classes by maximising
 * the between-class variance - no manual tuning needed.
 */
double dsp_threshold_otsu(const dsp_image *img);

/*
 * One-level 2-D Haar wavelet transform.
 *   img : image to transform, modified in place. Width and height
 *         must both be even.
 * After the call the image holds four quadrants:
 *   top-left     = approximation (half-resolution image)
 *   top-right    = horizontal detail
 *   bottom-left  = vertical detail
 *   bottom-right = diagonal detail
 * Returns 0 on success, -1 if a dimension is odd.
 */
int dsp_dwt2d_haar(dsp_image *img);

/*
 * Inverse of dsp_dwt2d_haar - reconstructs the image in place from
 * its four wavelet quadrants. Returns 0 on success, -1 on an odd
 * dimension.
 */
int dsp_idwt2d_haar(dsp_image *img);

#endif /* DSP_IMAGEOPS_H */
