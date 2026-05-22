/*
 * filter2d.h - Spatial filtering for images
 *
 * PROBLEM SOLVED
 *   Most image processing - blurring, sharpening, edge detection - is
 *   2-D CONVOLUTION: slide a small kernel over the image and replace
 *   each pixel with a weighted sum of its neighbourhood. This is the
 *   image-domain counterpart of the 1-D convolution in operations/.
 *
 * KERNELS PROVIDED
 *   Gaussian blur  - smooths noise and detail; weight falls off with
 *                    distance from the centre.
 *   Box blur       - the crude uniform-average blur.
 *   Sharpen        - boosts local contrast (centre-weighted kernel).
 *   Sobel          - gradient operator, the classic edge detector;
 *                    separate horizontal and vertical kernels whose
 *                    magnitude gives edge strength.
 *   Laplacian      - second-derivative operator; responds to edges
 *                    and fine detail regardless of orientation.
 *
 * NON-LINEAR FILTERING
 *   The median filter is included here too. It is NOT a convolution -
 *   it replaces each pixel with the median of its neighbourhood - and
 *   cannot be expressed as a kernel. Median filtering removes
 *   salt-and-pepper noise while preserving edges, which a linear blur
 *   cannot do.
 *
 * BORDER HANDLING
 *   All filters use clamp-to-edge borders (via dsp_image_at), so a
 *   window overhanging the image reuses the nearest edge pixels.
 */
#ifndef DSP_FILTER2D_H
#define DSP_FILTER2D_H

#include "image/image.h"

/*
 * General 2-D convolution.
 *   src    : input image
 *   dst    : output image (pre-allocated to the same size as src)
 *   kernel : kw*kh kernel, row-major
 *   kw, kh : kernel dimensions (odd sizes recommended)
 * Each output pixel is the kernel-weighted sum of the neighbourhood
 * centred on it. Clamp-to-edge borders are used.
 */
void dsp_convolve2d(const dsp_image *src, dsp_image *dst,
                    const double *kernel, size_t kw, size_t kh);

/*
 * Gaussian blur.
 *   src, dst : input and output images
 *   sigma    : standard deviation of the Gaussian; larger = blurrier
 * Builds a (truncated) Gaussian kernel internally and convolves.
 */
void dsp_gaussian_blur(const dsp_image *src, dsp_image *dst,
                       double sigma);

/* Uniform box blur with a `size` x `size` averaging kernel. */
void dsp_box_blur(const dsp_image *src, dsp_image *dst, size_t size);

/* Sharpen: a centre-weighted 3x3 kernel that boosts local contrast. */
void dsp_sharpen(const dsp_image *src, dsp_image *dst);

/*
 * Sobel edge detection.
 *   src : input image
 *   dst : output edge-magnitude image (pre-allocated, same size)
 * Computes the horizontal and vertical Sobel gradients and writes
 * their magnitude sqrt(gx^2 + gy^2) - bright where the image has a
 * strong edge, dark in flat regions.
 */
void dsp_sobel(const dsp_image *src, dsp_image *dst);

/* Laplacian: an orientation-independent second-derivative edge/detail
 * operator (3x3 kernel). */
void dsp_laplacian(const dsp_image *src, dsp_image *dst);

/*
 * Median filter (non-linear). Replaces each pixel with the median of
 * its `size` x `size` neighbourhood - removes salt-and-pepper noise
 * while keeping edges sharp. `size` should be odd.
 */
void dsp_median_filter(const dsp_image *src, dsp_image *dst,
                       size_t size);

#endif /* DSP_FILTER2D_H */
