/*
 * ldpc_build.c - LDPC code construction and syndrome checks.
 *
 * The Tanner graph is stored as two sets of adjacency lists so that
 * both message-passing directions (check->var and var->check) are
 * cheap to walk.
 */
#include "coding/ldpc.h"
#include <stdlib.h>
#include <string.h>

/* Free any partially built code (safe on NULL fields). */
static void ldpc_free_partial(dsp_ldpc *code) {
    if (code->row_var) {
        for (size_t i = 0; i < code->m; ++i)
            free(code->row_var[i]);
        free(code->row_var);
    }
    if (code->col_chk) {
        for (size_t j = 0; j < code->n; ++j)
            free(code->col_chk[j]);
        free(code->col_chk);
    }
    free(code->row_deg);
    free(code->col_deg);
    code->row_var = code->col_chk = NULL;
    code->row_deg = code->col_deg = NULL;
}

void dsp_ldpc_free(dsp_ldpc *code) {
    ldpc_free_partial(code);
    code->m = code->n = 0;
}

/*
 * Build the adjacency lists from a list of (check, var) edges.
 * `edges` holds nedge pairs flattened as [chk0,var0, chk1,var1, ...].
 */
static int build_from_edges(dsp_ldpc *code, size_t m, size_t n,
                            const size_t *edges, size_t nedge) {
    code->m = m;
    code->n = n;
    code->row_var = calloc(m, sizeof(size_t *));
    code->col_chk = calloc(n, sizeof(size_t *));
    code->row_deg = calloc(m, sizeof(size_t));
    code->col_deg = calloc(n, sizeof(size_t));
    if (!code->row_var || !code->col_chk ||
        !code->row_deg || !code->col_deg) {
        ldpc_free_partial(code);
        return -1;
    }

    /* Pass 1: count the degree of every node. */
    for (size_t e = 0; e < nedge; ++e) {
        code->row_deg[edges[2 * e]]     += 1;
        code->col_deg[edges[2 * e + 1]] += 1;
    }

    /* Allocate each adjacency list to its exact degree. */
    for (size_t i = 0; i < m; ++i) {
        code->row_var[i] = code->row_deg[i]
            ? malloc(code->row_deg[i] * sizeof(size_t)) : NULL;
        if (code->row_deg[i] && !code->row_var[i]) {
            ldpc_free_partial(code); return -1;
        }
    }
    for (size_t j = 0; j < n; ++j) {
        code->col_chk[j] = code->col_deg[j]
            ? malloc(code->col_deg[j] * sizeof(size_t)) : NULL;
        if (code->col_deg[j] && !code->col_chk[j]) {
            ldpc_free_partial(code); return -1;
        }
    }

    /* Pass 2: fill the lists. Reuse the degree arrays as write
     * cursors, then restore them afterwards. */
    size_t *rfill = calloc(m, sizeof(size_t));
    size_t *cfill = calloc(n, sizeof(size_t));
    if (!rfill || !cfill) {
        free(rfill); free(cfill);
        ldpc_free_partial(code); return -1;
    }
    for (size_t e = 0; e < nedge; ++e) {
        size_t chk = edges[2 * e];
        size_t var = edges[2 * e + 1];
        code->row_var[chk][rfill[chk]++] = var;
        code->col_chk[var][cfill[var]++] = chk;
    }
    free(rfill);
    free(cfill);
    return 0;
}

int dsp_ldpc_from_matrix(dsp_ldpc *code, const uint8_t *H,
                         size_t m, size_t n) {
    memset(code, 0, sizeof *code);

    /* Collect the coordinates of every 1 in H. */
    size_t cap = 16, nedge = 0;
    size_t *edges = malloc(cap * 2 * sizeof(size_t));
    if (!edges) return -1;

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (H[i * n + j]) {
                if (nedge == cap) {
                    cap *= 2;
                    size_t *t = realloc(edges, cap * 2 * sizeof(size_t));
                    if (!t) { free(edges); return -1; }
                    edges = t;
                }
                edges[2 * nedge]     = i;
                edges[2 * nedge + 1] = j;
                ++nedge;
            }
        }
    }

    int rc = build_from_edges(code, m, n, edges, nedge);
    free(edges);
    return rc;
}

/* Small xorshift RNG - reproducible across platforms. */
static unsigned rng_next(unsigned *s) {
    unsigned x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    return x;
}

int dsp_ldpc_make_regular(dsp_ldpc *code, size_t m, size_t n,
                          size_t wc, size_t wr, unsigned seed) {
    memset(code, 0, sizeof *code);

    /* A regular code needs the total edge count to agree both ways. */
    if (m == 0 || n == 0 || wc == 0 || wr == 0)
        return -1;
    if (n * wc != m * wr)
        return -1;
    if (wr > n || wc > m)          /* a check/var cannot exceed the count */
        return -1;

    size_t nedge = n * wc;

    /*
     * Progressive edge construction. For each variable node in turn we
     * pick wc distinct check nodes, preferring checks that still have
     * spare capacity and are not already connected to this variable.
     * This avoids the all-or-nothing failure of a global socket
     * shuffle, which almost never produces a duplicate-free pairing
     * for dense small codes.
     */
    size_t *chk_rem = calloc(m, sizeof(size_t));   /* spare slots / check */
    uint8_t *dense  = calloc(m * n, sizeof(uint8_t));
    size_t  *edges  = malloc(nedge * 2 * sizeof(size_t));
    size_t  *cand   = malloc(m * sizeof(size_t));
    if (!chk_rem || !dense || !edges || !cand) {
        free(chk_rem); free(dense); free(edges); free(cand);
        return -1;
    }

    unsigned state = seed ? seed : 0x1234567u;
    int ok = 0;

    for (int attempt = 0; attempt < 200 && !ok; ++attempt) {
        for (size_t i = 0; i < m; ++i) chk_rem[i] = wr;
        memset(dense, 0, m * n);
        int feasible = 1;

        for (size_t j = 0; j < n && feasible; ++j) {
            for (size_t d = 0; d < wc; ++d) {
                /* Candidate checks: have capacity, not yet linked to j. */
                size_t ncand = 0;
                for (size_t i = 0; i < m; ++i)
                    if (chk_rem[i] > 0 && !dense[i * n + j])
                        cand[ncand++] = i;

                if (ncand == 0) { feasible = 0; break; }

                /* Pick the candidate(s) with the most remaining slots,
                 * breaking ties at random - keeps degrees balanced. */
                size_t best_rem = 0;
                for (size_t c = 0; c < ncand; ++c)
                    if (chk_rem[cand[c]] > best_rem)
                        best_rem = chk_rem[cand[c]];
                size_t top = 0;
                for (size_t c = 0; c < ncand; ++c)
                    if (chk_rem[cand[c]] == best_rem)
                        cand[top++] = cand[c];

                size_t pick = cand[rng_next(&state) % top];
                dense[pick * n + j] = 1;
                chk_rem[pick] -= 1;
            }
        }

        /* A valid attempt has used every check's slots exactly. */
        if (feasible) {
            int all_used = 1;
            for (size_t i = 0; i < m; ++i)
                if (chk_rem[i] != 0) { all_used = 0; break; }
            ok = all_used;
        }
    }

    int rc = -1;
    if (ok) {
        /* Flatten the dense matrix into an edge list. */
        size_t e = 0;
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                if (dense[i * n + j]) {
                    edges[2 * e]     = i;
                    edges[2 * e + 1] = j;
                    ++e;
                }
        rc = build_from_edges(code, m, n, edges, e);
    }

    free(chk_rem);
    free(dense);
    free(edges);
    free(cand);
    return rc;
}

/* Parity of one check: XOR of the bits on its variable nodes. */
static int check_parity(const dsp_ldpc *code, size_t i,
                         const uint8_t *bits) {
    int p = 0;
    for (size_t d = 0; d < code->row_deg[i]; ++d)
        p ^= (bits[code->row_var[i][d]] & 1);
    return p;
}

size_t dsp_ldpc_syndrome_weight(const dsp_ldpc *code,
                                const uint8_t *bits) {
    size_t w = 0;
    for (size_t i = 0; i < code->m; ++i)
        if (check_parity(code, i, bits))
            ++w;
    return w;
}

int dsp_ldpc_check(const dsp_ldpc *code, const uint8_t *bits) {
    return dsp_ldpc_syndrome_weight(code, bits) == 0;
}
