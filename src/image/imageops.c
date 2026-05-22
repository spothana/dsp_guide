/*
 * imageops.c - Point operators and the 2-D Haar wavelet transform.
 */
#include "image/imageops.h"
#include <stdlib.h>
#include <math.h>

/* ---- Histogram ------------------------------------------------------ */

void dsp_histogram(const dsp_image *img, size_t hist[256]) {
    for (int i = 0; i < 256; ++i)
        hist[i] = 0;

    size_t n = dsp_image_size(img);
    for (size_t i = 0; i < n; ++i) {
        long v = (long)floor(img->data[i] + 0.5);
        if (v < 0)   v = 0;
        if (v > 255) v = 255;
        hist[v]++;
    }
}

void dsp_histogram_equalize(dsp_image *img) {
    size_t hist[256];
    dsp_histogram(img, hist);

    size_t n = dsp_image_size(img);
    if (n == 0) return;

    /* Cumulative distribution function of the histogram. */
    size_t cdf[256];
    size_t running = 0;
    for (int i = 0; i < 256; ++i) {
        running += hist[i];
        cdf[i] = running;
    }

    /* Find the first non-zero CDF entry - the equalization reference. */
    size_t cdf_min = 0;
    for (int i = 0; i < 256; ++i) {
        if (cdf[i] != 0) { cdf_min = cdf[i]; break; }
    }

    /* Build the remap table: stretch the CDF across the full 0..255
     * range so output intensities are as uniform as possible. */
    double denom = (double)(n - cdf_min);
    if (denom <= 0.0) denom = 1.0;

    double lut[256];
    for (int i = 0; i < 256; ++i) {
        double v = ((double)cdf[i] - (double)cdf_min) / denom * 255.0;
        lut[i] = (v < 0.0) ? 0.0 : (v > 255.0 ? 255.0 : v);
    }

    /* Apply the remap to every pixel. */
    for (size_t i = 0; i < n; ++i) {
        long v = (long)floor(img->data[i] + 0.5);
        if (v < 0)   v = 0;
        if (v > 255) v = 255;
        img->data[i] = lut[v];
    }
}

/* ---- Thresholding --------------------------------------------------- */

void dsp_threshold(const dsp_image *src, dsp_image *dst,
                   double thresh) {
    size_t n = dsp_image_size(src);
    for (size_t i = 0; i < n; ++i)
        dst->data[i] = (src->data[i] >= thresh) ? 255.0 : 0.0;
}

double dsp_threshold_otsu(const dsp_image *img) {
    size_t hist[256];
    dsp_histogram(img, hist);

    size_t n = dsp_image_size(img);
    if (n == 0) return 128.0;

    /* Total intensity-weighted sum, for the class means. */
    double total_sum = 0.0;
    for (int i = 0; i < 256; ++i)
        total_sum += (double)i * (double)hist[i];

    double sum_bg   = 0.0;       /* weighted sum of the background class */
    size_t count_bg = 0;         /* pixel count of the background class */
    double best_var = -1.0;
    int    best_t   = 128;

    /* Sweep every possible threshold; keep the one that maximises the
     * between-class variance. */
    for (int t = 0; t < 256; ++t) {
        count_bg += hist[t];
        if (count_bg == 0)
            continue;
        size_t count_fg = n - count_bg;
        if (count_fg == 0)
            break;

        sum_bg += (double)t * (double)hist[t];
        double mean_bg = sum_bg / (double)count_bg;
        double mean_fg = (total_sum - sum_bg) / (double)count_fg;

        double diff = mean_bg - mean_fg;
        double between = (double)count_bg * (double)count_fg
                       * diff * diff;
        if (between > best_var) {
            best_var = between;
            best_t   = t;
        }
    }
    return (double)best_t;
}

/* ---- 2-D Haar wavelet ----------------------------------------------- */

static const double INV_SQRT2 = 0.70710678118654752440;

/* One-level 1-D Haar transform of a length-`len` vector (len even).
 * Output: [ approximation (len/2) | detail (len/2) ]. */
static void haar_fwd(const double *in, double *out, size_t len) {
    size_t half = len / 2;
    for (size_t i = 0; i < half; ++i) {
        double a = in[2 * i], b = in[2 * i + 1];
        out[i]        = (a + b) * INV_SQRT2;   /* average  */
        out[half + i] = (a - b) * INV_SQRT2;   /* difference */
    }
}

/* Inverse of haar_fwd. */
static void haar_inv(const double *in, double *out, size_t len) {
    size_t half = len / 2;
    for (size_t i = 0; i < half; ++i) {
        double s = in[i], d = in[half + i];
        out[2 * i]     = (s + d) * INV_SQRT2;
        out[2 * i + 1] = (s - d) * INV_SQRT2;
    }
}

int dsp_dwt2d_haar(dsp_image *img) {
    size_t w = img->width, h = img->height;
    if ((w & 1) || (h & 1))
        return -1;

    double *line = malloc((w > h ? w : h) * sizeof(double));
    double *res  = malloc((w > h ? w : h) * sizeof(double));
    if (!line || !res) { free(line); free(res); return -1; }

    /* Transform each row, then each column. The row pass leaves the
     * left half holding averages and the right half details; the
     * column pass then splits those into the four wavelet quadrants. */
    for (size_t y = 0; y < h; ++y) {
        for (size_t x = 0; x < w; ++x) line[x] = img->data[y * w + x];
        haar_fwd(line, res, w);
        for (size_t x = 0; x < w; ++x) img->data[y * w + x] = res[x];
    }
    for (size_t x = 0; x < w; ++x) {
        for (size_t y = 0; y < h; ++y) line[y] = img->data[y * w + x];
        haar_fwd(line, res, h);
        for (size_t y = 0; y < h; ++y) img->data[y * w + x] = res[y];
    }

    free(line);
    free(res);
    return 0;
}

int dsp_idwt2d_haar(dsp_image *img) {
    size_t w = img->width, h = img->height;
    if ((w & 1) || (h & 1))
        return -1;

    double *line = malloc((w > h ? w : h) * sizeof(double));
    double *res  = malloc((w > h ? w : h) * sizeof(double));
    if (!line || !res) { free(line); free(res); return -1; }

    /* Reverse the forward order: columns first, then rows. */
    for (size_t x = 0; x < w; ++x) {
        for (size_t y = 0; y < h; ++y) line[y] = img->data[y * w + x];
        haar_inv(line, res, h);
        for (size_t y = 0; y < h; ++y) img->data[y * w + x] = res[y];
    }
    for (size_t y = 0; y < h; ++y) {
        for (size_t x = 0; x < w; ++x) line[x] = img->data[y * w + x];
        haar_inv(line, res, w);
        for (size_t x = 0; x < w; ++x) img->data[y * w + x] = res[x];
    }

    free(line);
    free(res);
    return 0;
}
