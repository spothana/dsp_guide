/*
 * interleave.c - Block and convolutional interleavers.
 */
#include "coding/interleave.h"
#include <stdlib.h>

/* ===================================================================
 * Block interleaver
 *
 * The matrix is filled row by row and emptied column by column:
 *   write index  (r, c) = r * cols + c
 *   read  index  (r, c) = c * rows + r
 * Symbols that were adjacent on input end up `rows` apart on output,
 * so a burst of up to `rows` consecutive errors is spread one-per-
 * codeword once the receiver reverses the permutation.
 * =================================================================== */

int dsp_block_interleave(const uint8_t *in, uint8_t *out, size_t len,
                         size_t rows, size_t cols) {
    if (len != rows * cols)
        return -1;

    /* Read out column by column. */
    size_t k = 0;
    for (size_t c = 0; c < cols; ++c)
        for (size_t r = 0; r < rows; ++r)
            out[k++] = in[r * cols + c];
    return 0;
}

int dsp_block_deinterleave(const uint8_t *in, uint8_t *out, size_t len,
                           size_t rows, size_t cols) {
    if (len != rows * cols)
        return -1;

    /* Exact inverse: the input arrived in column-major order, so
     * place each element back at its original row-major slot. */
    size_t k = 0;
    for (size_t c = 0; c < cols; ++c)
        for (size_t r = 0; r < rows; ++r)
            out[r * cols + c] = in[k++];
    return 0;
}

/* ===================================================================
 * Convolutional interleaver
 *
 * `branches` FIFO delay lines, branch i holding i * depth symbols.
 * Symbols are committed to branches in round-robin order. Branch 0
 * has zero delay (a pass-through); each later branch delays more, so
 * the output is a continuously scrambled version of the input.
 * =================================================================== */

/* Build the delay profile shared by interleaver and deinterleaver. */
static int conv_alloc(dsp_conv_interleaver *ci,
                      size_t branches, size_t depth, int reversed) {
    if (branches == 0 || depth == 0)
        return -1;

    ci->branches = branches;
    ci->depth    = depth;
    ci->cur      = 0;

    ci->pos      = calloc(branches, sizeof(size_t));
    ci->line_len = calloc(branches, sizeof(size_t));
    ci->offset   = calloc(branches, sizeof(size_t));
    if (!ci->pos || !ci->line_len || !ci->offset) {
        free(ci->pos); free(ci->line_len); free(ci->offset);
        ci->pos = ci->line_len = ci->offset = NULL;
        return -1;
    }

    /* Interleaver: branch i delays by i * depth.
     * Deinterleaver: mirrored, branch i delays by (B-1-i) * depth,
     * so the two stages together delay every branch equally. */
    size_t total = 0;
    for (size_t i = 0; i < branches; ++i) {
        size_t d = reversed ? (branches - 1 - i) : i;
        ci->line_len[i] = d * depth;
        ci->offset[i]   = total;
        total += ci->line_len[i];
    }

    /* One flat buffer backs every delay line (branch 0 may be empty). */
    ci->store = (total > 0) ? calloc(total, sizeof(uint8_t)) : NULL;
    if (total > 0 && !ci->store) {
        free(ci->pos); free(ci->line_len); free(ci->offset);
        ci->pos = ci->line_len = ci->offset = NULL;
        return -1;
    }
    return 0;
}

int dsp_conv_interleaver_init(dsp_conv_interleaver *ci,
                              size_t branches, size_t depth) {
    return conv_alloc(ci, branches, depth, 0);
}

int dsp_conv_deinterleaver_init(dsp_conv_interleaver *ci,
                                size_t branches, size_t depth) {
    return conv_alloc(ci, branches, depth, 1);
}

void dsp_conv_interleaver_free(dsp_conv_interleaver *ci) {
    free(ci->store);
    free(ci->pos);
    free(ci->line_len);
    free(ci->offset);
    ci->store = NULL;
    ci->pos = ci->line_len = ci->offset = NULL;
    ci->branches = ci->depth = ci->cur = 0;
}

uint8_t dsp_conv_interleave_step(dsp_conv_interleaver *ci, uint8_t sym) {
    size_t b = ci->cur;
    ci->cur = (ci->cur + 1) % ci->branches;

    /* A zero-length branch is a pass-through. */
    if (ci->line_len[b] == 0)
        return sym;

    /* FIFO ring buffer for this branch: read the oldest symbol out,
     * then overwrite that slot with the new one. */
    uint8_t *line = ci->store + ci->offset[b];
    size_t   p    = ci->pos[b];

    uint8_t out = line[p];
    line[p] = sym;
    ci->pos[b] = (p + 1) % ci->line_len[b];
    return out;
}
