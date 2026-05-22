/*
 * detect.h - Error DETECTION techniques
 *
 * PROBLEM SOLVED
 *   A receiver needs to know whether a block of bits arrived intact.
 *   Detection codes append a small amount of redundancy that fails a
 *   recomputed check when the data has been corrupted. Detection only
 *   flags errors - it does not fix them (see correct.h for that).
 *
 * TECHNIQUES PROVIDED
 *   Parity   - one bit making the total number of 1s even (or odd).
 *              Detects any ODD number of bit flips; misses even ones.
 *              Cheapest possible check.
 *   Checksum - sum of the data bytes, appended for verification.
 *              Catches general corruption; weak against multi-bit
 *              errors that cancel in the sum. Used in TCP/IP.
 *   CRC      - treats the message as a polynomial and divides by a
 *              fixed generator polynomial; the remainder is appended.
 *              Strong against BURST errors, which dominate wireless
 *              and stored-media channels. CRC-32 is used in Ethernet.
 *
 * DETECTION vs CORRECTION
 *   Detection is cheap but only enables a retransmission request (ARQ).
 *   Forward error correction (correct.h) adds enough redundancy to
 *   repair errors in place, trading bandwidth for zero round-trips.
 */
#ifndef DSP_DETECT_H
#define DSP_DETECT_H

#include <stddef.h>
#include <stdint.h>

/* ---- Parity --------------------------------------------------------- */

typedef enum {
    DSP_PARITY_EVEN,   /* parity bit chosen so total 1-count is even */
    DSP_PARITY_ODD     /* ... so total 1-count is odd                */
} dsp_parity_type;

/*
 * Compute the parity bit for a byte buffer.
 * Returns 0 or 1 - the bit that, appended to the data, gives the
 * requested overall parity.
 */
int dsp_parity_compute(const uint8_t *data, size_t len,
                        dsp_parity_type type);

/*
 * Check data against a received parity bit.
 * Returns 1 if parity is consistent (no detected error), 0 otherwise.
 */
int dsp_parity_check(const uint8_t *data, size_t len,
                     dsp_parity_type type, int parity_bit);

/* ---- Checksum ------------------------------------------------------- */

/*
 * One's-complement 16-bit checksum (the Internet checksum used by
 * IP/TCP/UDP). Processes the buffer as big-endian 16-bit words.
 */
uint16_t dsp_checksum16(const uint8_t *data, size_t len);

/*
 * Verify a buffer plus its stored checksum. A correct buffer makes the
 * recomputed checksum (over data + checksum) equal zero.
 * Returns 1 if valid, 0 if corruption is detected.
 */
int dsp_checksum16_verify(const uint8_t *data, size_t len,
                          uint16_t checksum);

/* ---- CRC ------------------------------------------------------------ */

/*
 * CRC-32 (IEEE 802.3 / Ethernet polynomial 0xEDB8820, reflected).
 * Standard reflected algorithm with 0xFFFFFFFF init and final XOR.
 */
uint32_t dsp_crc32(const uint8_t *data, size_t len);

/*
 * Verify data against an expected CRC-32. Returns 1 if they match.
 */
int dsp_crc32_verify(const uint8_t *data, size_t len, uint32_t expected);

#endif /* DSP_DETECT_H */
