/*
 * array.c - Beamforming and DOA estimation for a uniform linear array.
 *
 * Structure:
 *   - steering vectors, snapshot synthesis, covariance
 *   - complex Hermitian linear solve (for MVDR)
 *   - complex Hermitian Jacobi eigensolver (for MUSIC / ESPRIT)
 *   - delay-and-sum and MVDR beamformers
 *   - spatial MUSIC and spatial ESPRIT
 */
#include "array/array.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===================================================================
 * Steering vectors, synthesis, covariance
 * =================================================================== */

void dsp_array_steering(size_t nsensors, double spacing,
                        double theta, cplx *a) {
    /* Sensor k sees a phase advance of 2*pi*d*k*sin(theta). */
    double s = sin(theta);
    for (size_t k = 0; k < nsensors; ++k) {
        double phase = 2.0 * M_PI * spacing * (double)k * s;
        a[k] = cexp(phase * I);
    }
}

/* xorshift RNG for reproducible noise. */
static unsigned arr_rng(unsigned *s) {
    unsigned x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}
static double arr_uniform(unsigned *s) {
    return (arr_rng(s) >> 8) / 16777216.0;
}
/* Standard-normal sample via Box-Muller. */
static double arr_gauss(unsigned *s) {
    static const double TWO_PI = 6.28318530717958647692;
    double u1 = arr_uniform(s);
    double u2 = arr_uniform(s);
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
}

void dsp_array_synthesize(size_t nsensors, double spacing,
                          size_t nsnapshots,
                          const double *angles, size_t nsources,
                          double noise_std, unsigned seed,
                          cplx *snapshots) {
    unsigned state = seed ? seed : 0x1234567u;

    /* Pre-compute each source's steering vector. */
    cplx *steer = malloc(nsources * nsensors * sizeof(cplx));
    if (!steer)
        return;
    for (size_t s = 0; s < nsources; ++s)
        dsp_array_steering(nsensors, spacing, angles[s],
                           steer + s * nsensors);

    double nsig = noise_std / sqrt(2.0);   /* per-component noise */

    for (size_t t = 0; t < nsnapshots; ++t) {
        cplx *row = snapshots + t * nsensors;
        for (size_t k = 0; k < nsensors; ++k)
            row[k] = 0.0;

        /* Each source has a random complex amplitude this snapshot
         * (uncorrelated sources), shared across all sensors. */
        for (size_t s = 0; s < nsources; ++s) {
            double ar = arr_gauss(&state);
            double ai = arr_gauss(&state);
            cplx amp = dsp_cplx(ar, ai);
            const cplx *a = steer + s * nsensors;
            for (size_t k = 0; k < nsensors; ++k)
                row[k] += amp * a[k];
        }
        /* Independent complex sensor noise. */
        if (noise_std > 0.0) {
            for (size_t k = 0; k < nsensors; ++k) {
                double nr = nsig * arr_gauss(&state);
                double ni = nsig * arr_gauss(&state);
                row[k] += dsp_cplx(nr, ni);
            }
        }
    }
    free(steer);
}

void dsp_array_covariance(const cplx *snapshots, size_t nsensors,
                          size_t nsnapshots, cplx *R) {
    for (size_t i = 0; i < nsensors * nsensors; ++i)
        R[i] = 0.0;

    /* R = (1/T) sum_t x(t) x(t)^H. */
    for (size_t t = 0; t < nsnapshots; ++t) {
        const cplx *x = snapshots + t * nsensors;
        for (size_t i = 0; i < nsensors; ++i)
            for (size_t j = 0; j < nsensors; ++j)
                R[i * nsensors + j] += x[i] * conj(x[j]);
    }
    double inv = (nsnapshots > 0) ? 1.0 / (double)nsnapshots : 1.0;
    for (size_t i = 0; i < nsensors * nsensors; ++i)
        R[i] *= inv;
}

/* ===================================================================
 * Complex Hermitian linear solve - for the MVDR beamformer
 *
 * Solves R w = b for a Hermitian positive-definite R by Gaussian
 * elimination with partial pivoting on the complex matrix.
 * =================================================================== */
static int hermitian_solve(const cplx *R, size_t m,
                           const cplx *b, cplx *w) {
    cplx *A = malloc(m * m * sizeof(cplx));
    cplx *y = malloc(m * sizeof(cplx));
    if (!A || !y) { free(A); free(y); return -1; }
    memcpy(A, R, m * m * sizeof(cplx));
    memcpy(y, b, m * sizeof(cplx));

    for (size_t col = 0; col < m; ++col) {
        /* Partial pivot on the largest-magnitude entry. */
        size_t piv = col;
        double best = cabs(A[col * m + col]);
        for (size_t row = col + 1; row < m; ++row) {
            double v = cabs(A[row * m + col]);
            if (v > best) { best = v; piv = row; }
        }
        if (best < 1e-12) { free(A); free(y); return -1; }
        if (piv != col) {
            for (size_t k = 0; k < m; ++k) {
                cplx t = A[col * m + k];
                A[col * m + k] = A[piv * m + k];
                A[piv * m + k] = t;
            }
            cplx t = y[col]; y[col] = y[piv]; y[piv] = t;
        }
        for (size_t row = col + 1; row < m; ++row) {
            cplx f = A[row * m + col] / A[col * m + col];
            for (size_t k = col; k < m; ++k)
                A[row * m + k] -= f * A[col * m + k];
            y[row] -= f * y[col];
        }
    }
    for (size_t i = m; i-- > 0; ) {
        cplx s = y[i];
        for (size_t k = i + 1; k < m; ++k)
            s -= A[i * m + k] * w[k];
        w[i] = s / A[i * m + i];
    }
    free(A); free(y);
    return 0;
}

/* ===================================================================
 * Complex Hermitian Jacobi eigensolver
 *
 * The Hermitian analogue of the real Jacobi method: repeated complex
 * plane rotations zero the largest off-diagonal element. A Hermitian
 * matrix has real eigenvalues; the eigenvectors are complex and are
 * accumulated in evec (columns). Suited to the small covariance
 * matrices array processing uses.
 * =================================================================== */
static void hermitian_eig(cplx *A, size_t m,
                          double *eval, cplx *evec) {
    /* evec starts as the identity. */
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < m; ++j)
            evec[i * m + j] = (i == j) ? 1.0 : 0.0;

    for (int sweep = 0; sweep < 200; ++sweep) {
        /* Largest off-diagonal magnitude. */
        double off = 0.0;
        size_t p = 0, q = 1;
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i + 1; j < m; ++j) {
                double v = cabs(A[i * m + j]);
                if (v > off) { off = v; p = i; q = j; }
            }
        if (off < 1e-14)
            break;

        /* Complex Jacobi rotation that zeros A[p][q].
         * App, Aqq are real (Hermitian diagonal); Apq is complex. */
        double app = creal(A[p * m + p]);
        double aqq = creal(A[q * m + q]);
        cplx   apq = A[p * m + q];
        double apq_abs = cabs(apq);
        if (apq_abs < 1e-300)
            continue;

        /* Phase of the off-diagonal element. */
        cplx phase = apq / apq_abs;

        /* Rotation angle (as in the real symmetric case, using the
         * magnitude of the off-diagonal term). */
        double theta = 0.5 * atan2(2.0 * apq_abs, aqq - app);
        double c = cos(theta);
        double s = sin(theta);

        /* The unitary rotation mixes rows/columns p and q:
         *   col p <-  c*col_p - conj(phase)*s*col_q
         *   col q <-  phase*s*col_p + c*col_q                       */
        cplx s_p = phase * s;          /* applied to p-side */
        cplx s_q = conj(phase) * s;    /* applied to q-side */

        for (size_t k = 0; k < m; ++k) {
            cplx akp = A[k * m + p];
            cplx akq = A[k * m + q];
            A[k * m + p] = c * akp - s_q * akq;
            A[k * m + q] = s_p * akp + c * akq;
        }
        for (size_t k = 0; k < m; ++k) {
            cplx apk = A[p * m + k];
            cplx aqk = A[q * m + k];
            A[p * m + k] = c * apk - conj(s_q) * aqk;
            A[q * m + k] = conj(s_p) * apk + c * aqk;
        }
        for (size_t k = 0; k < m; ++k) {
            cplx vkp = evec[k * m + p];
            cplx vkq = evec[k * m + q];
            evec[k * m + p] = c * vkp - s_q * vkq;
            evec[k * m + q] = s_p * vkp + c * vkq;
        }
    }

    for (size_t i = 0; i < m; ++i)
        eval[i] = creal(A[i * m + i]);
}

/* Sort eigenpairs by descending eigenvalue. */
static void sort_eig_desc(double *eval, cplx *evec, size_t m) {
    for (size_t i = 0; i < m; ++i) {
        size_t best = i;
        for (size_t j = i + 1; j < m; ++j)
            if (eval[j] > eval[best]) best = j;
        if (best != i) {
            double t = eval[i]; eval[i] = eval[best]; eval[best] = t;
            for (size_t k = 0; k < m; ++k) {
                cplx tv = evec[k * m + i];
                evec[k * m + i] = evec[k * m + best];
                evec[k * m + best] = tv;
            }
        }
    }
}

/* ===================================================================
 * Beamforming
 * =================================================================== */

/* Map an angle-grid index to an angle in [-pi/2, +pi/2]. */
static double grid_angle(size_t i, size_t nangles) {
    return -M_PI / 2.0 + M_PI * (double)i / (double)(nangles - 1);
}

void dsp_beamform_conventional(const cplx *R, size_t nsensors,
                               double spacing,
                               double *power, size_t nangles) {
    cplx *a = malloc(nsensors * sizeof(cplx));
    if (!a) return;

    for (size_t ang = 0; ang < nangles; ++ang) {
        double theta = grid_angle(ang, nangles);
        dsp_array_steering(nsensors, spacing, theta, a);

        /* Response = a^H R a (real, since R is Hermitian). */
        cplx acc = 0.0;
        for (size_t i = 0; i < nsensors; ++i) {
            cplx Ra = 0.0;
            for (size_t j = 0; j < nsensors; ++j)
                Ra += R[i * nsensors + j] * a[j];
            acc += conj(a[i]) * Ra;
        }
        power[ang] = creal(acc);
    }
    free(a);
}

int dsp_beamform_mvdr(const cplx *R, size_t nsensors, double spacing,
                      double *power, size_t nangles) {
    cplx *a = malloc(nsensors * sizeof(cplx));
    cplx *w = malloc(nsensors * sizeof(cplx));
    if (!a || !w) { free(a); free(w); return -1; }

    for (size_t ang = 0; ang < nangles; ++ang) {
        double theta = grid_angle(ang, nangles);
        dsp_array_steering(nsensors, spacing, theta, a);

        /* Solve R w = a, then the MVDR power is 1 / (a^H w). */
        if (hermitian_solve(R, nsensors, a, w) != 0) {
            free(a); free(w);
            return -1;
        }
        cplx denom = 0.0;
        for (size_t i = 0; i < nsensors; ++i)
            denom += conj(a[i]) * w[i];

        double d = creal(denom);
        power[ang] = (fabs(d) > 1e-300) ? 1.0 / d : 1e300;
    }
    free(a); free(w);
    return 0;
}

/* ===================================================================
 * Spatial MUSIC
 * =================================================================== */

int dsp_doa_music(const cplx *R, size_t nsensors, double spacing,
                  size_t nsources, double *pseudo, size_t nangles) {
    if (nsources == 0 || nsources >= nsensors)
        return -1;

    cplx   *A    = malloc(nsensors * nsensors * sizeof(cplx));
    double *eval = malloc(nsensors * sizeof(double));
    cplx   *evec = malloc(nsensors * nsensors * sizeof(cplx));
    cplx   *a    = malloc(nsensors * sizeof(cplx));
    if (!A || !eval || !evec || !a) {
        free(A); free(eval); free(evec); free(a);
        return -1;
    }
    memcpy(A, R, nsensors * nsensors * sizeof(cplx));

    hermitian_eig(A, nsensors, eval, evec);
    sort_eig_desc(eval, evec, nsensors);

    /* The top `nsources` eigenvectors span the signal subspace; the
     * rest span the noise subspace. MUSIC peaks where a steering
     * vector is orthogonal to the noise subspace. */
    for (size_t ang = 0; ang < nangles; ++ang) {
        double theta = grid_angle(ang, nangles);
        dsp_array_steering(nsensors, spacing, theta, a);

        double denom = 0.0;
        for (size_t v = nsources; v < nsensors; ++v) {
            /* Projection of a onto noise eigenvector v: a^H v. */
            cplx dot = 0.0;
            for (size_t k = 0; k < nsensors; ++k)
                dot += conj(a[k]) * evec[k * nsensors + v];
            denom += creal(dot * conj(dot));
        }
        pseudo[ang] = (denom > 1e-300) ? 1.0 / denom : 1e300;
    }

    free(A); free(eval); free(evec); free(a);
    return 0;
}

/* ===================================================================
 * Spatial ESPRIT
 *
 * Splits the array into two maximally-overlapping subarrays (rows
 * 0..M-2 and rows 1..M-1 of the signal-subspace eigenvectors). They
 * are related by a diagonal phase rotation whose eigenvalues are
 * exp(j*2*pi*d*sin(theta)) - one per source. The angle is read off
 * each eigenvalue's phase.
 * =================================================================== */

int dsp_doa_esprit(const cplx *R, size_t nsensors, double spacing,
                   size_t nsources, double *angles) {
    if (nsources == 0 || nsources >= nsensors)
        return -1;

    cplx   *A    = malloc(nsensors * nsensors * sizeof(cplx));
    double *eval = malloc(nsensors * sizeof(double));
    cplx   *evec = malloc(nsensors * nsensors * sizeof(cplx));
    if (!A || !eval || !evec) {
        free(A); free(eval); free(evec);
        return -1;
    }
    memcpy(A, R, nsensors * nsensors * sizeof(cplx));

    hermitian_eig(A, nsensors, eval, evec);
    sort_eig_desc(eval, evec, nsensors);

    /* Signal subspace Es: the top nsources eigenvectors (columns).
     * Form the two subarrays Es1 (rows 0..M-2) and Es2 (rows 1..M-1).
     * Solve Es1 * Psi = Es2 in the least-squares sense; Psi is
     * nsources x nsources and its eigenvalues carry the angles. */
    size_t d   = nsources;
    size_t sub = nsensors - 1;

    /* Build the d x d normal-equation matrices:
     *   G = Es1^H Es1,  C = Es1^H Es2. */
    cplx *G = calloc(d * d, sizeof(cplx));
    cplx *C = calloc(d * d, sizeof(cplx));
    cplx *Psi = calloc(d * d, sizeof(cplx));
    if (!G || !C || !Psi) {
        free(A); free(eval); free(evec);
        free(G); free(C); free(Psi);
        return -1;
    }
    for (size_t i = 0; i < sub; ++i) {
        for (size_t r = 0; r < d; ++r) {
            cplx e1r = evec[i * nsensors + r];        /* Es1[i][r] */
            for (size_t c = 0; c < d; ++c) {
                cplx e1c = evec[i * nsensors + c];
                cplx e2c = evec[(i + 1) * nsensors + c];  /* Es2 */
                G[r * d + c] += conj(e1r) * e1c;
                C[r * d + c] += conj(e1r) * e2c;
            }
        }
    }

    /* Solve G * Psi = C, one column of Psi at a time. */
    for (size_t col = 0; col < d; ++col) {
        cplx *rhs = malloc(d * sizeof(cplx));
        cplx *sol = malloc(d * sizeof(cplx));
        if (!rhs || !sol) { free(rhs); free(sol); break; }
        for (size_t r = 0; r < d; ++r)
            rhs[r] = C[r * d + col];
        if (hermitian_solve(G, d, rhs, sol) == 0)
            for (size_t r = 0; r < d; ++r)
                Psi[r * d + col] = sol[r];
        free(rhs); free(sol);
    }

    /* The eigenvalues of Psi are exp(j*2*pi*d*sin(theta)) - they lie
     * on the unit circle, all of magnitude one. Plain QR converges
     * very slowly on such a matrix, so a Wilkinson-style shift is
     * used: subtract the trailing diagonal entry, run a QR step, add
     * it back. This deflates one eigenvalue per few steps. */
    cplx *H  = malloc(d * d * sizeof(cplx));
    cplx *Q  = malloc(d * d * sizeof(cplx));
    cplx *Rm = malloc(d * d * sizeof(cplx));
    cplx *lam = malloc(d * sizeof(cplx));
    if (!H || !Q || !Rm || !lam) {
        free(A); free(eval); free(evec);
        free(G); free(C); free(Psi);
        free(H); free(Q); free(Rm); free(lam);
        return -1;
    }
    memcpy(H, Psi, d * d * sizeof(cplx));

    /* Shifted QR with deflation: shrink the active block from size m
     * down to 1, recording one eigenvalue each time it deflates. */
    size_t m = d;
    int found = 0;
    while (m > 1 && found < (int)nsources) {
        for (int iter = 0; iter < 200; ++iter) {
            /* Wilkinson shift = current trailing diagonal entry. */
            cplx mu = H[(m - 1) * d + (m - 1)];
            for (size_t i = 0; i < m; ++i)
                H[i * d + i] -= mu;

            /* Complex Gram-Schmidt QR of the active m x m block. */
            for (size_t i = 0; i < d * d; ++i) Q[i] = 0.0;
            for (size_t c = 0; c < m; ++c)
                for (size_t r = 0; r < m; ++r)
                    Rm[r * d + c] = H[r * d + c];
            for (size_t c = 0; c < m; ++c) {
                for (size_t k = 0; k < c; ++k) {
                    cplx dot = 0.0;
                    for (size_t r = 0; r < m; ++r)
                        dot += conj(Q[r * d + k]) * H[r * d + c];
                    for (size_t r = 0; r < m; ++r)
                        Rm[r * d + c] -= dot * Q[r * d + k];
                }
                double nrm = 0.0;
                for (size_t r = 0; r < m; ++r)
                    nrm += creal(Rm[r * d + c] * conj(Rm[r * d + c]));
                nrm = sqrt(nrm);
                if (nrm < 1e-300) nrm = 1e-300;
                for (size_t r = 0; r < m; ++r)
                    Q[r * d + c] = Rm[r * d + c] / nrm;
            }
            /* H <- R Q  with  R = Q^H H, on the active block. */
            for (size_t r = 0; r < m; ++r)
                for (size_t c = 0; c < m; ++c) {
                    cplx s = 0.0;
                    for (size_t k = 0; k < m; ++k)
                        s += conj(Q[k * d + r]) * H[k * d + c];
                    Rm[r * d + c] = s;
                }
            for (size_t r = 0; r < m; ++r)
                for (size_t c = 0; c < m; ++c) {
                    cplx s = 0.0;
                    for (size_t k = 0; k < m; ++k)
                        s += Rm[r * d + k] * Q[k * d + c];
                    H[r * d + c] = s;
                }
            /* Undo the shift. */
            for (size_t i = 0; i < m; ++i)
                H[i * d + i] += mu;

            /* Converged once the sub-diagonal entry is negligible. */
            if (cabs(H[(m - 1) * d + (m - 2)]) < 1e-12)
                break;
        }
        /* Deflate: the trailing diagonal entry is an eigenvalue. */
        lam[found++] = H[(m - 1) * d + (m - 1)];
        --m;
    }
    if (m == 1 && found < (int)nsources)
        lam[found++] = H[0];

    /* Convert each eigenvalue exp(j*2*pi*d*sin(theta)) to an angle. */
    int got = 0;
    for (int i = 0; i < found && got < (int)nsources; ++i) {
        double ang = carg(lam[i]);            /* 2*pi*d*sin(theta) */
        double sintheta = ang / (2.0 * M_PI * spacing);
        if (sintheta >  1.0) sintheta =  1.0;
        if (sintheta < -1.0) sintheta = -1.0;
        angles[got++] = asin(sintheta);
    }

    free(A); free(eval); free(evec);
    free(G); free(C); free(Psi);
    free(H); free(Q); free(Rm); free(lam);
    return got;
}
