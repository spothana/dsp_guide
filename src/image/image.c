/*
 * image.c - Grayscale image container and PGM file I/O.
 */
#include "image/image.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

int dsp_image_alloc(dsp_image *img, size_t width, size_t height) {
    img->width  = width;
    img->height = height;
    img->data   = calloc(width * height, sizeof(double));
    return img->data ? 0 : -1;
}

void dsp_image_free(dsp_image *img) {
    free(img->data);
    img->data = NULL;
    img->width = img->height = 0;
}

int dsp_image_copy(dsp_image *dst, const dsp_image *src) {
    if (dsp_image_alloc(dst, src->width, src->height) != 0)
        return -1;
    memcpy(dst->data, src->data,
           dsp_image_size(src) * sizeof(double));
    return 0;
}

double dsp_image_at(const dsp_image *img, long x, long y) {
    /* Clamp-to-edge border handling: out-of-range coordinates take
     * the value of the nearest edge pixel. This avoids dark borders
     * when a filter window overhangs the image. */
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x >= (long)img->width)  x = (long)img->width  - 1;
    if (y >= (long)img->height) y = (long)img->height - 1;
    return img->data[(size_t)y * img->width + (size_t)x];
}

void dsp_image_clamp(dsp_image *img) {
    size_t n = dsp_image_size(img);
    for (size_t i = 0; i < n; ++i) {
        double v = img->data[i];
        img->data[i] = (v < 0.0) ? 0.0 : (v > 255.0 ? 255.0 : v);
    }
}

/* Skip whitespace and #-comments in a PGM header. */
static int pgm_skip(FILE *f) {
    int c;
    for (;;) {
        c = fgetc(f);
        if (c == EOF) return EOF;
        if (c == '#') {                 /* comment to end of line */
            while (c != '\n' && c != EOF) c = fgetc(f);
        } else if (!isspace(c)) {
            ungetc(c, f);
            return 0;
        }
    }
}

int dsp_image_load_pgm(dsp_image *img, const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;

    char magic[3] = {0};
    if (fscanf(f, "%2s", magic) != 1) { fclose(f); return -1; }

    int ascii;
    if (strcmp(magic, "P2") == 0)      ascii = 1;
    else if (strcmp(magic, "P5") == 0) ascii = 0;
    else { fclose(f); return -1; }

    long w, h, maxval;
    if (pgm_skip(f) == EOF || fscanf(f, "%ld", &w) != 1 ||
        pgm_skip(f) == EOF || fscanf(f, "%ld", &h) != 1 ||
        pgm_skip(f) == EOF || fscanf(f, "%ld", &maxval) != 1 ||
        w <= 0 || h <= 0 || maxval <= 0) {
        fclose(f);
        return -1;
    }
    fgetc(f);                           /* single whitespace after header */

    if (dsp_image_alloc(img, (size_t)w, (size_t)h) != 0) {
        fclose(f);
        return -1;
    }

    size_t n = (size_t)w * (size_t)h;
    if (ascii) {
        for (size_t i = 0; i < n; ++i) {
            long v;
            if (fscanf(f, "%ld", &v) != 1) {
                dsp_image_free(img); fclose(f); return -1;
            }
            img->data[i] = (double)v;
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            int v = fgetc(f);
            if (v == EOF) {
                dsp_image_free(img); fclose(f); return -1;
            }
            img->data[i] = (double)v;
        }
    }

    fclose(f);
    return 0;
}

int dsp_image_save_pgm(const dsp_image *img, const char *path) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    fprintf(f, "P5\n%zu %zu\n255\n", img->width, img->height);

    size_t n = dsp_image_size(img);
    for (size_t i = 0; i < n; ++i) {
        double v = img->data[i];
        long   q = (long)floor(v + 0.5);    /* round to nearest */
        if (q < 0)   q = 0;
        if (q > 255) q = 255;
        fputc((int)q, f);
    }

    fclose(f);
    return 0;
}
