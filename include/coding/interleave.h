/*
 * interleave.h - Interleaving for burst-error resilience
 *
 * PROBLEM SOLVED
 *   Every FEC code in correct.h has a limit: Hamming fixes 1 bit per
 *   codeword, Reed-Solomon t symbols, a convolutional code only a few
 *   scattered errors. A real channel, though, tends to corrupt data in
 *   BURSTS - a fade, an impulse, a scratch on a disc - dumping many
 *   errors into one stretch. A burst longer than a code's limit lands
 *   entirely inside one codeword and is uncorrectable.
 *
 *   An INTERLEAVER fixes this without touching the code. The
 *   transmitter shuffles symbols into a scrambled order; the receiver
 *   unshuffles them. A burst that hits the scrambled stream is spread
 *   thin once de-interleaved - its errors are scattered across many
 *   codewords, so each codeword sees only a few, within its limit.
 *
 *   transmit:  encode -> interleave -> channel
 *   receive:   channel -> deinterleave -> decode
 *
 * TWO TYPES PROVIDED
 *   Block interleaver - fill an R x C matrix row by row, read it out
 *     column by column. A burst of up to R consecutive errors is
 *     spread so that consecutive codewords get one error each. Simple,
 *     but needs the whole block before it can start, adding latency.
 *     Used in DVB-T and many storage systems.
 *   Convolutional interleaver - a bank of N delay lines of staggered
 *     length. It streams continuously (no block boundary) and reaches
 *     the same burst-spreading with about half the latency and memory
 *     of an equivalent block interleaver. Used in DVB and ADSL.
 *
 * KEY POINT
 *   Interleaving adds NO redundancy and corrects NOTHING on its own.
 *   It only rearranges symbols so the FEC stage that follows can do
 *   its job. Interleaver depth trades burst protection against latency.
 */
#ifndef DSP_INTERLEAVE_H
#define DSP_INTERLEAVE_H

#include <stddef.h>
#include <stdint.h>

/* ===================================================================
 * Block interleaver - R x C matrix, write rows / read columns
 * =================================================================== */

/*
 * Block-interleave `len` symbols using an `rows` x `cols` matrix.
 *   in   : input symbols, length len
 *   out  : output symbols, length len
 *   len  : symbol count - MUST equal rows * cols
 * Writes the input across matrix rows, reads it out down columns.
 * Returns 0 on success, -1 if len != rows * cols.
 */
int dsp_block_interleave(const uint8_t *in, uint8_t *out, size_t len,
                         size_t rows, size_t cols);

/*
 * Block-deinterleave: the exact inverse of dsp_block_interleave.
 * Writes down columns, reads across rows, restoring the original order.
 */
int dsp_block_deinterleave(const uint8_t *in, uint8_t *out, size_t len,
                           size_t rows, size_t cols);

/* ===================================================================
 * Convolutional interleaver - bank of staggered delay lines
 * =================================================================== */

/*
 * Convolutional interleaver / deinterleaver state.
 *
 * Branch i (i = 0 .. branches-1) holds a FIFO delay line of length
 * i * depth. Symbols are fed to branches in round-robin order. The
 * matching deinterleaver uses the mirrored delay profile so that
 * every symbol experiences the same total delay.
 */
typedef struct {
    size_t   branches;     /* number of delay lines */
    size_t   depth;        /* delay increment per branch */
    size_t   cur;          /* next branch to use (round-robin) */
    uint8_t *store;        /* flat backing buffer for all delay lines */
    size_t  *pos;          /* per-branch ring-buffer cursor */
    size_t  *line_len;     /* per-branch delay length */
    size_t  *offset;       /* per-branch start offset into `store` */
} dsp_conv_interleaver;

/*
 * Initialise a convolutional interleaver with `branches` delay lines
 * and a per-branch delay increment of `depth`.
 * Returns 0 on success, -1 on bad parameters or allocation failure.
 * Pair with dsp_conv_interleaver_free.
 */
int dsp_conv_interleaver_init(dsp_conv_interleaver *ci,
                              size_t branches, size_t depth);

/* Release the buffers owned by a convolutional interleaver. */
void dsp_conv_interleaver_free(dsp_conv_interleaver *ci);

/*
 * Push one symbol through the interleaver, returning one symbol.
 * Because each branch delays by a different amount, the output is a
 * time-scrambled version of the input. The deinterleaver is the same
 * structure with the branch order reversed (see _deint_init below).
 */
uint8_t dsp_conv_interleave_step(dsp_conv_interleaver *ci, uint8_t sym);

/*
 * Initialise a convolutional DEinterleaver matching an interleaver of
 * the same (branches, depth). It has the reversed delay profile, so
 * interleaver + deinterleaver gives every symbol an equal total delay.
 */
int dsp_conv_deinterleaver_init(dsp_conv_interleaver *ci,
                                size_t branches, size_t depth);

/*
 * Total end-to-end latency (in symbols) of an interleaver +
 * deinterleaver pair: branches * (branches - 1) * depth.
 */
static inline size_t dsp_conv_interleaver_latency(size_t branches,
                                                  size_t depth) {
    return branches * (branches - 1) * depth;
}

#endif /* DSP_INTERLEAVE_H */
