/*
 * filter2d.c - 2-D convolution and the common image filters.
 */
#include "image/filter2d.h"
#include <stdlib.h>
#include <math.h>

void dsp_convolve2d(const dsp_image *src, dsp_image *dst,
                    const double *kernel, size_t kw, size_t kh) {
    long cx = (long)kw / 2;        /* kernel centre offsets */
    long cy = (long)kh / 2;

    for (size_t y = 0; y < src->height; ++y) {
        for (size_t x = 0; x < src->width; ++x) {
            double acc = 0.0;
            /* Weighted sum over the kernel neighbourhood. */
            for (size_t ky = 0; ky < kh; ++ky) {
                for (size_t kx = 0; kx < kw; ++kx) {
                    long sx = (long)x + (long)kx - cx;
                    long sy = (long)y + (long)ky - cy;
                    acc += kernel[ky * kw + kx]
                         * dsp_image_at(src, sx, sy);
                }
            }
            dst->data[y * src->width + x] = acc;
        }
    }
}

void dsp_gaussian_blur(const dsp_image *src, dsp_image *dst,
                       double sigma) {
    if (sigma <= 0.0) sigma = 0.5;

    /* Kernel radius ~3 sigma captures essentially all the weight. */
    long radius = (long)ceil(3.0 * sigma);
    size_t k = (size_t)(2 * radius + 1);

    double *kernel = malloc(k * k * sizeof(double));
    if (!kernel) return;

    /* Build a normalised 2-D Gaussian. */
    double sum = 0.0;
    for (long dy = -radius; dy <= radius; ++dy) {
        for (long dx = -radius; dx <= radius; ++dx) {
            double v = exp(-(double)(dx * dx + dy * dy)
                           / (2.0 * sigma * sigma));
            kernel[(dy + radius) * k + (dx + radius)] = v;
            sum += v;
        }
    }
    for (size_t i = 0; i < k * k; ++i)
        kernel[i] /= sum;          /* unit DC gain -> no brightness shift */

    dsp_convolve2d(src, dst, kernel, k, k);
    free(kernel);
}

void dsp_box_blur(const dsp_image *src, dsp_image *dst, size_t size) {
    if (size == 0) size = 1;
    double w = 1.0 / (double)(size * size);

    double *kernel = malloc(size * size * sizeof(double));
    if (!kernel) return;
    for (size_t i = 0; i < size * size; ++i)
        kernel[i] = w;

    dsp_convolve2d(src, dst, kernel, size, size);
    free(kernel);
}

void dsp_sharpen(const dsp_image *src, dsp_image *dst) {
    /* Centre-weighted kernel: original plus its difference from a
     * blurred version (an unsharp-mask style 3x3). */
    static const double kernel[9] = {
         0.0, -1.0,  0.0,
        -1.0,  5.0, -1.0,
         0.0, -1.0,  0.0
    };
    dsp_convolve2d(src, dst, kernel, 3, 3);
}

void dsp_sobel(const dsp_image *src, dsp_image *dst) {
    /* Horizontal and vertical Sobel gradient kernels. */
    static const double gx[9] = {
        -1.0, 0.0, 1.0,
        -2.0, 0.0, 2.0,
        -1.0, 0.0, 1.0
    };
    static const double gy[9] = {
        -1.0, -2.0, -1.0,
         0.0,  0.0,  0.0,
         1.0,  2.0,  1.0
    };

    /* Convolve with each kernel, then take the gradient magnitude. */
    dsp_image hx, hy;
    if (dsp_image_alloc(&hx, src->width, src->height) != 0) return;
    if (dsp_image_alloc(&hy, src->width, src->height) != 0) {
        dsp_image_free(&hx); return;
    }
    dsp_convolve2d(src, &hx, gx, 3, 3);
    dsp_convolve2d(src, &hy, gy, 3, 3);

    size_t n = dsp_image_size(src);
    for (size_t i = 0; i < n; ++i) {
        double a = hx.data[i], b = hy.data[i];
        dst->data[i] = sqrt(a * a + b * b);
    }

    dsp_image_free(&hx);
    dsp_image_free(&hy);
}

void dsp_laplacian(const dsp_image *src, dsp_image *dst) {
    /* 3x3 Laplacian: sum of second derivatives, orientation-free. */
    static const double kernel[9] = {
         0.0,  1.0,  0.0,
         1.0, -4.0,  1.0,
         0.0,  1.0,  0.0
    };
    dsp_convolve2d(src, dst, kernel, 3, 3);
}

/* Insertion sort - fine for the small windows a median filter uses. */
static void sort_doubles(double *a, size_t n) {
    for (size_t i = 1; i < n; ++i) {
        double key = a[i];
        size_t j = i;
        while (j > 0 && a[j - 1] > key) {
            a[j] = a[j - 1];
            --j;
        }
        a[j] = key;
    }
}

void dsp_median_filter(const dsp_image *src, dsp_image *dst,
                       size_t size) {
    if (size == 0) size = 1;
    long c = (long)size / 2;

    double *window = malloc(size * size * sizeof(double));
    if (!window) return;

    for (size_t y = 0; y < src->height; ++y) {
        for (size_t x = 0; x < src->width; ++x) {
            /* Gather the neighbourhood, sort it, take the middle. */
            size_t cnt = 0;
            for (long dy = -c; dy <= c; ++dy)
                for (long dx = -c; dx <= c; ++dx)
                    window[cnt++] = dsp_image_at(src,
                                                 (long)x + dx,
                                                 (long)y + dy);
            sort_doubles(window, cnt);
            dst->data[y * src->width + x] = window[cnt / 2];
        }
    }

    free(window);
}
