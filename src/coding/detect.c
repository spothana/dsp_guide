/*
 * detect.c - Parity, Internet checksum, and CRC-32.
 */
#include "coding/detect.h"

/* ---- Parity --------------------------------------------------------- */

/* Count of set bits across the whole buffer, taken modulo 2. */
static int ones_parity(const uint8_t *data, size_t len) {
    int p = 0;
    for (size_t i = 0; i < len; ++i) {
        uint8_t b = data[i];
        b ^= b >> 4;
        b ^= b >> 2;
        b ^= b >> 1;
        p ^= (b & 1);
    }
    return p;
}

int dsp_parity_compute(const uint8_t *data, size_t len,
                       dsp_parity_type type) {
    int p = ones_parity(data, len);
    /* Even parity: bit makes the count even -> bit == current parity.
     * Odd parity: invert that. */
    return (type == DSP_PARITY_EVEN) ? p : (p ^ 1);
}

int dsp_parity_check(const uint8_t *data, size_t len,
                     dsp_parity_type type, int parity_bit) {
    return dsp_parity_compute(data, len, type) == (parity_bit & 1);
}

/* ---- Checksum (16-bit one's complement, the Internet checksum) ------ */

/* Add with end-around carry, as one's-complement arithmetic requires. */
static uint32_t add_carry(uint32_t sum, uint16_t word) {
    sum += word;
    sum = (sum & 0xFFFF) + (sum >> 16);   /* fold the carry back in */
    return sum;
}

uint16_t dsp_checksum16(const uint8_t *data, size_t len) {
    uint32_t sum = 0;
    size_t i = 0;

    /* Sum big-endian 16-bit words. */
    for (; i + 1 < len; i += 2) {
        uint16_t word = (uint16_t)((data[i] << 8) | data[i + 1]);
        sum = add_carry(sum, word);
    }
    /* Odd trailing byte is padded with a zero low byte. */
    if (i < len)
        sum = add_carry(sum, (uint16_t)(data[i] << 8));

    /* Fold any residual carry, then one's-complement. */
    sum = (sum & 0xFFFF) + (sum >> 16);
    return (uint16_t)(~sum & 0xFFFF);
}

int dsp_checksum16_verify(const uint8_t *data, size_t len,
                          uint16_t checksum) {
    /* Recompute over the data and add the stored checksum: a clean
     * buffer yields all-ones, whose complement is zero. */
    uint32_t sum = 0;
    size_t i = 0;
    for (; i + 1 < len; i += 2) {
        uint16_t word = (uint16_t)((data[i] << 8) | data[i + 1]);
        sum = add_carry(sum, word);
    }
    if (i < len)
        sum = add_carry(sum, (uint16_t)(data[i] << 8));
    sum = add_carry(sum, checksum);
    sum = (sum & 0xFFFF) + (sum >> 16);
    return ((uint16_t)(~sum & 0xFFFF)) == 0;
}

/* ---- CRC-32 (reflected, Ethernet polynomial) ------------------------ */

static uint32_t crc32_table[256];
static int      crc32_table_ready = 0;

/* Build the 256-entry lookup table on first use (reflected form). */
static void crc32_init(void) {
    for (uint32_t i = 0; i < 256; ++i) {
        uint32_t c = i;
        for (int k = 0; k < 8; ++k)
            c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
        crc32_table[i] = c;
    }
    crc32_table_ready = 1;
}

uint32_t dsp_crc32(const uint8_t *data, size_t len) {
    if (!crc32_table_ready)
        crc32_init();

    uint32_t crc = 0xFFFFFFFFu;            /* standard initial value */
    for (size_t i = 0; i < len; ++i) {
        uint8_t idx = (uint8_t)(crc ^ data[i]);
        crc = crc32_table[idx] ^ (crc >> 8);
    }
    return crc ^ 0xFFFFFFFFu;              /* standard final XOR */
}

int dsp_crc32_verify(const uint8_t *data, size_t len, uint32_t expected) {
    return dsp_crc32(data, len) == expected;
}
