/*
 * image.h - Grayscale image type and PGM I/O
 *
 * The image-processing modules in this guide are 2-D. They all share
 * one simple container: a grayscale image stored as a flat array of
 * double-precision pixels in row-major order.
 *
 *   pixel (x, y) is at data[y * width + x]
 *
 * Pixels are doubles, not bytes, so transforms and filters can work
 * in their natural numeric range without repeated quantisation.
 * Convention: 0.0 = black, 255.0 = white, though intermediate results
 * may fall outside that range until clamped on output.
 *
 * I/O uses PGM (Portable GrayMap), the simplest standard image format
 * - a tiny text header followed by pixel values. No external library
 * is needed.
 */
#ifndef DSP_IMAGE_H
#define DSP_IMAGE_H

#include <stddef.h>

/*
 * A grayscale image.
 *   width, height : dimensions in pixels
 *   data          : width*height pixels, row-major, 0..255 nominal
 */
typedef struct {
    size_t  width;
    size_t  height;
    double *data;
} dsp_image;

/*
 * Allocate an image of the given size. Pixels are zero-initialised.
 * Returns 0 on success, -1 on allocation failure. Pair with
 * dsp_image_free.
 */
int dsp_image_alloc(dsp_image *img, size_t width, size_t height);

/* Release an image's pixel buffer. */
void dsp_image_free(dsp_image *img);

/* Deep-copy `src` into `dst` (dst is allocated by this call). */
int dsp_image_copy(dsp_image *dst, const dsp_image *src);

/* Total pixel count. */
static inline size_t dsp_image_size(const dsp_image *img) {
    return img->width * img->height;
}

/*
 * Read/write a pixel with edge handling. dsp_image_at clamps
 * out-of-range coordinates to the nearest edge pixel (the standard
 * "clamp" border mode for filtering); plain data[] indexing is used
 * when coordinates are known to be in range.
 */
double dsp_image_at(const dsp_image *img, long x, long y);

/* Clamp every pixel into [0, 255] - call before saving. */
void dsp_image_clamp(dsp_image *img);

/*
 * Load a binary (P5) or ASCII (P2) PGM file into `img`.
 * Returns 0 on success, -1 on a file or format error.
 */
int dsp_image_load_pgm(dsp_image *img, const char *path);

/*
 * Save `img` as a binary (P5) PGM file. Pixels are rounded and
 * clamped to 0..255. Returns 0 on success, -1 on a file error.
 */
int dsp_image_save_pgm(const dsp_image *img, const char *path);

#endif /* DSP_IMAGE_H */
