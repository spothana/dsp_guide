/*
 * transform2d.c - Separable 2-D FFT and 2-D DCT.
 *
 * Every 2-D transform here is separable: apply the 1-D transform to
 * each row, then to each column of the result. The 1-D FFT and 1-D
 * DCT from the transforms/ modules do the actual work.
 */
#include "image/transform2d.h"
#include "transforms/fft.h"
#include "transforms/dct.h"
#include <stdlib.h>
#include <string.h>

/* ---- 2-D FFT -------------------------------------------------------- */

int dsp_fft2d(const dsp_image *img, double *re, double *im) {
    size_t w = img->width, h = img->height;
    if (!dsp_is_pow2(w) || !dsp_is_pow2(h))
        return -1;

    /* Pack the real image into a complex working buffer. */
    cplx *buf = malloc(w * h * sizeof(cplx));
    cplx *line = malloc((w > h ? w : h) * sizeof(cplx));
    if (!buf || !line) { free(buf); free(line); return -1; }

    for (size_t i = 0; i < w * h; ++i)
        buf[i] = img->data[i];

    /* Transform each row. */
    for (size_t y = 0; y < h; ++y) {
        for (size_t x = 0; x < w; ++x) line[x] = buf[y * w + x];
        dsp_fft(line, w);
        for (size_t x = 0; x < w; ++x) buf[y * w + x] = line[x];
    }

    /* Transform each column. */
    for (size_t x = 0; x < w; ++x) {
        for (size_t y = 0; y < h; ++y) line[y] = buf[y * w + x];
        dsp_fft(line, h);
        for (size_t y = 0; y < h; ++y) buf[y * w + x] = line[y];
    }

    for (size_t i = 0; i < w * h; ++i) {
        re[i] = creal(buf[i]);
        im[i] = cimag(buf[i]);
    }

    free(buf);
    free(line);
    return 0;
}

int dsp_ifft2d(const double *re, const double *im, dsp_image *img) {
    size_t w = img->width, h = img->height;
    if (!dsp_is_pow2(w) || !dsp_is_pow2(h))
        return -1;

    cplx *buf = malloc(w * h * sizeof(cplx));
    cplx *line = malloc((w > h ? w : h) * sizeof(cplx));
    if (!buf || !line) { free(buf); free(line); return -1; }

    for (size_t i = 0; i < w * h; ++i)
        buf[i] = re[i] + im[i] * I;

    /* Inverse-transform columns, then rows (order does not matter). */
    for (size_t x = 0; x < w; ++x) {
        for (size_t y = 0; y < h; ++y) line[y] = buf[y * w + x];
        dsp_ifft(line, h);
        for (size_t y = 0; y < h; ++y) buf[y * w + x] = line[y];
    }
    for (size_t y = 0; y < h; ++y) {
        for (size_t x = 0; x < w; ++x) line[x] = buf[y * w + x];
        dsp_ifft(line, w);
        for (size_t x = 0; x < w; ++x) buf[y * w + x] = line[x];
    }

    for (size_t i = 0; i < w * h; ++i)
        img->data[i] = creal(buf[i]);

    free(buf);
    free(line);
    return 0;
}

/* ---- 2-D DCT -------------------------------------------------------- */

/* Apply the separable 2-D DCT-II (forward) or DCT-III (inverse). */
static void dct2d_core(const double *in, double *out,
                       size_t w, size_t h, int inverse) {
    double *tmp  = malloc(w * h * sizeof(double));
    double *line = malloc((w > h ? w : h) * sizeof(double));
    double *res  = malloc((w > h ? w : h) * sizeof(double));
    if (!tmp || !line || !res) { free(tmp); free(line); free(res); return; }

    /* Rows. */
    for (size_t y = 0; y < h; ++y) {
        for (size_t x = 0; x < w; ++x) line[x] = in[y * w + x];
        if (inverse) dsp_idct(line, res, w);
        else         dsp_dct(line, res, w);
        for (size_t x = 0; x < w; ++x) tmp[y * w + x] = res[x];
    }

    /* Columns. */
    for (size_t x = 0; x < w; ++x) {
        for (size_t y = 0; y < h; ++y) line[y] = tmp[y * w + x];
        if (inverse) dsp_idct(line, res, h);
        else         dsp_dct(line, res, h);
        for (size_t y = 0; y < h; ++y) out[y * w + x] = res[y];
    }

    free(tmp);
    free(line);
    free(res);
}

void dsp_dct2d(const dsp_image *img, double *coeffs) {
    dct2d_core(img->data, coeffs, img->width, img->height, 0);
}

void dsp_idct2d(const double *coeffs, dsp_image *img) {
    dct2d_core(coeffs, img->data, img->width, img->height, 1);
}

void dsp_dct8x8(const double *block, double *out) {
    dct2d_core(block, out, 8, 8, 0);
}

void dsp_idct8x8(const double *coeffs, double *block) {
    dct2d_core(coeffs, block, 8, 8, 1);
}
