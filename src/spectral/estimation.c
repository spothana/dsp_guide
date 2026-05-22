/*
 * estimation.c - AR, ARMA, MUSIC, and ESPRIT spectral estimation.
 *
 * Structure:
 *   - autocorrelation
 *   - AR via Levinson-Durbin (Yule-Walker) and Burg
 *   - a Jacobi eigensolver for real symmetric matrices
 *   - ARMA via the modified Yule-Walker method
 *   - MUSIC and ESPRIT, built on the eigensolver
 */
#include "spectral/estimation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===================================================================
 * Autocorrelation
 * =================================================================== */

void dsp_autocorr(const double *x, size_t n, size_t maxlag, double *r) {
    for (size_t k = 0; k <= maxlag; ++k) {
        double s = 0.0;
        for (size_t i = k; i < n; ++i)
            s += x[i] * x[i - k];
        r[k] = (n > 0) ? s / (double)n : 0.0;   /* biased estimate */
    }
}

/* ===================================================================
 * AR - Yule-Walker via the Levinson-Durbin recursion
 * =================================================================== */

int dsp_ar_yule_walker(const double *r, size_t order,
                       double *a, double *sigma2) {
    if (order == 0 || r[0] <= 0.0)
        return -1;

    double *prev = calloc(order, sizeof(double));
    if (!prev)
        return -1;

    double err = r[0];                  /* prediction error power */

    for (size_t i = 0; i < order; ++i) {
        /* Reflection coefficient for order i+1. */
        double acc = r[i + 1];
        for (size_t j = 0; j < i; ++j)
            acc -= a[j] * r[i - j];

        if (err <= 0.0) { free(prev); return -1; }
        double k = acc / err;

        /* Update the coefficient vector with the new reflection. */
        for (size_t j = 0; j < i; ++j)
            prev[j] = a[j];
        a[i] = k;
        for (size_t j = 0; j < i; ++j)
            a[j] = prev[j] - k * prev[i - 1 - j];

        err *= (1.0 - k * k);           /* error power shrinks */
    }

    /* Convention here: x(n) = -sum a[k] x(n-1-k) + w(n), so the
     * coefficients carry a leading minus sign relative to the
     * prediction filter above. */
    for (size_t i = 0; i < order; ++i)
        a[i] = -a[i];

    if (sigma2) *sigma2 = err;
    free(prev);
    return 0;
}

/* ===================================================================
 * AR - Burg's method (forward-backward linear prediction)
 * =================================================================== */

int dsp_ar_burg(const double *x, size_t n, size_t order,
                double *a, double *sigma2) {
    if (order == 0 || n <= order)
        return -1;

    /* f[i] forward error, b[i] backward error. ap[] holds the
     * prediction-error filter with the leading ap[0] = 1 fixed; the
     * model coefficients are read from ap[1..order] at the end. */
    double *f   = malloc(n * sizeof(double));
    double *b   = malloc(n * sizeof(double));
    double *ap  = calloc(order + 1, sizeof(double));
    double *apr = calloc(order + 1, sizeof(double));
    if (!f || !b || !ap || !apr) {
        free(f); free(b); free(ap); free(apr);
        return -1;
    }
    for (size_t i = 0; i < n; ++i)
        f[i] = b[i] = x[i];
    ap[0] = 1.0;

    double err = 0.0;
    for (size_t i = 0; i < n; ++i)
        err += x[i] * x[i];
    err /= (double)n;

    for (size_t m = 1; m <= order; ++m) {
        /* Reflection coefficient minimising the forward + backward
         * squared error: k = -2*sum(f[i]*b[i-1]) / sum(f^2 + b^2). */
        double num = 0.0, den = 0.0;
        for (size_t i = m; i < n; ++i) {
            num += f[i] * b[i - 1];
            den += f[i] * f[i] + b[i - 1] * b[i - 1];
        }
        double k = (den > 1e-300) ? (-2.0 * num / den) : 0.0;

        /* Order-update the prediction-error filter:
         *   ap[i] <- ap[i] + k * ap[m-i],  for i = 1 .. m
         * (ap[0] stays 1; ap[m] becomes k). */
        for (size_t i = 0; i <= order; ++i)
            apr[i] = ap[i];
        for (size_t i = 1; i <= m; ++i)
            ap[i] = apr[i] + k * apr[m - i];

        /* Update the forward/backward errors, top-down so each pair
         * reads pre-update values. */
        for (size_t i = n - 1; i >= m; --i) {
            double ft = f[i]     + k * b[i - 1];
            double bt = b[i - 1] + k * f[i];
            f[i] = ft;
            b[i] = bt;
            if (i == m) break;              /* size_t underflow guard */
        }

        err *= (1.0 - k * k);
    }

    /* Model x(n) = -sum a[k] x(n-1-k) + w(n): a[k] = ap[k+1]. */
    for (size_t k = 0; k < order; ++k)
        a[k] = ap[k + 1];

    if (sigma2) *sigma2 = err;
    free(f); free(b); free(ap); free(apr);
    return 0;
}

/* ===================================================================
 * AR power spectral density
 * =================================================================== */

void dsp_ar_psd(const double *a, size_t order, double sigma2,
                double *psd, size_t nf) {
    for (size_t i = 0; i < nf; ++i) {
        double freq = 0.5 * (double)i / (double)nf;   /* [0, 0.5) */

        /* Denominator A(f) = 1 + sum a[k] e^{-j2pi f (k+1)}. */
        double re = 1.0, im = 0.0;
        for (size_t k = 0; k < order; ++k) {
            double ang = -2.0 * M_PI * freq * (double)(k + 1);
            re += a[k] * cos(ang);
            im += a[k] * sin(ang);
        }
        double denom = re * re + im * im;
        psd[i] = (denom > 1e-300) ? sigma2 / denom : 0.0;
    }
}

/* ===================================================================
 * ARMA - modified Yule-Walker estimation
 * =================================================================== */

void dsp_arma_psd(const double *a, size_t p,
                  const double *b, size_t q,
                  double *psd, size_t nf) {
    for (size_t i = 0; i < nf; ++i) {
        double freq = 0.5 * (double)i / (double)nf;

        /* Numerator B(f) = sum b[k] e^{-j2pi f k}. */
        double nre = 0.0, nim = 0.0;
        for (size_t k = 0; k <= q; ++k) {
            double ang = -2.0 * M_PI * freq * (double)k;
            nre += b[k] * cos(ang);
            nim += b[k] * sin(ang);
        }
        /* Denominator A(f) = 1 + sum a[k] e^{-j2pi f (k+1)}. */
        double dre = 1.0, dim = 0.0;
        for (size_t k = 0; k < p; ++k) {
            double ang = -2.0 * M_PI * freq * (double)(k + 1);
            dre += a[k] * cos(ang);
            dim += a[k] * sin(ang);
        }
        double num = nre * nre + nim * nim;
        double den = dre * dre + dim * dim;
        psd[i] = (den > 1e-300) ? num / den : 0.0;
    }
}

int dsp_arma_estimate(const double *x, size_t n, size_t p, size_t q,
                      double *a, double *b) {
    if (p == 0 || n <= p + q + 1)
        return -1;

    /* Step 1: estimate the AR part from the modified Yule-Walker
     * equations, which use autocorrelations at lags > q (where the MA
     * part has died out, so they depend only on the AR parameters). */
    size_t maxlag = p + q + p;
    double *r = malloc((maxlag + 1) * sizeof(double));
    if (!r) return -1;
    dsp_autocorr(x, n, maxlag, r);

    /* Solve the p x p system  R a = -rhs, where the matrix and rhs
     * are built from autocorrelations around lag q. */
    double *M   = malloc(p * p * sizeof(double));
    double *rhs = malloc(p * sizeof(double));
    if (!M || !rhs) { free(r); free(M); free(rhs); return -1; }

    for (size_t i = 0; i < p; ++i) {
        long lag_i = (long)q + 1 + (long)i;
        for (size_t j = 0; j < p; ++j) {
            long lag = lag_i - 1 - (long)j;
            M[i * p + j] = r[lag < 0 ? -lag : lag];
        }
        rhs[i] = -r[lag_i];
    }

    /* Gaussian elimination with partial pivoting. */
    for (size_t col = 0; col < p; ++col) {
        size_t piv = col;
        double best = fabs(M[col * p + col]);
        for (size_t row = col + 1; row < p; ++row) {
            double v = fabs(M[row * p + col]);
            if (v > best) { best = v; piv = row; }
        }
        if (best < 1e-12) { free(r); free(M); free(rhs); return -1; }
        if (piv != col) {
            for (size_t k = 0; k < p; ++k) {
                double t = M[col * p + k];
                M[col * p + k] = M[piv * p + k];
                M[piv * p + k] = t;
            }
            double t = rhs[col]; rhs[col] = rhs[piv]; rhs[piv] = t;
        }
        for (size_t row = col + 1; row < p; ++row) {
            double f = M[row * p + col] / M[col * p + col];
            for (size_t k = col; k < p; ++k)
                M[row * p + k] -= f * M[col * p + k];
            rhs[row] -= f * rhs[col];
        }
    }
    for (size_t ii = p; ii-- > 0; ) {
        double s = rhs[ii];
        for (size_t k = ii + 1; k < p; ++k)
            s -= M[ii * p + k] * a[k];
        a[ii] = s / M[ii * p + ii];
    }

    /* Step 2: the MA spectrum. Filter the data by the estimated AR
     * polynomial to get the residual, whose autocorrelation gives the
     * MA covariance. We then take the MA coefficients from a simple
     * spectral factorisation surrogate: b[k] = residual autocorr,
     * normalised so b[0] is the gain. This is the standard practical
     * (if approximate) modified-Yule-Walker MA step. */
    double *resid = malloc(n * sizeof(double));
    if (!resid) { free(r); free(M); free(rhs); return -1; }
    for (size_t i = 0; i < n; ++i) {
        double v = x[i];
        for (size_t k = 0; k < p; ++k)
            if (i > k) v += a[k] * x[i - 1 - k];
        resid[i] = v;
    }

    double *rr = malloc((q + 1) * sizeof(double));
    if (!rr) { free(r); free(M); free(rhs); free(resid); return -1; }
    dsp_autocorr(resid, n, q, rr);

    /* Normalise so the MA part has b[0] = sqrt(power) gain and the
     * remaining taps are the normalised residual correlations. */
    double gain = (rr[0] > 0.0) ? sqrt(rr[0]) : 1.0;
    b[0] = gain;
    for (size_t k = 1; k <= q; ++k)
        b[k] = (rr[0] > 0.0) ? rr[k] / gain : 0.0;

    free(r); free(M); free(rhs); free(resid); free(rr);
    return 0;
}

/* ===================================================================
 * Jacobi eigensolver for real symmetric matrices
 *
 * Repeatedly applies Jacobi rotations that zero the largest
 * off-diagonal element, until the matrix is diagonal. The diagonal
 * then holds the eigenvalues and the accumulated rotations the
 * eigenvectors. Ideal for the small symmetric covariance matrices
 * the subspace methods use.
 * =================================================================== */

/* A : msize x msize symmetric (overwritten). eval : eigenvalues.
 * evec : eigenvectors in columns. Returns 0 on success. */
static int jacobi_eig(double *A, size_t m, double *eval, double *evec) {
    /* evec starts as the identity. */
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < m; ++j)
            evec[i * m + j] = (i == j) ? 1.0 : 0.0;

    for (int sweep = 0; sweep < 100; ++sweep) {
        /* Find the largest off-diagonal magnitude. */
        double off = 0.0;
        size_t pi = 0, pj = 1;
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i + 1; j < m; ++j) {
                double v = fabs(A[i * m + j]);
                if (v > off) { off = v; pi = i; pj = j; }
            }
        if (off < 1e-14)
            break;                       /* converged: A is diagonal */

        /* Rotation angle that zeros A[pi][pj]. */
        double app = A[pi * m + pi];
        double aqq = A[pj * m + pj];
        double apq = A[pi * m + pj];
        double phi = 0.5 * atan2(2.0 * apq, aqq - app);
        double c = cos(phi), s = sin(phi);

        /* Apply the rotation to rows/columns pi and pj. */
        for (size_t k = 0; k < m; ++k) {
            double aki = A[k * m + pi];
            double akj = A[k * m + pj];
            A[k * m + pi] = c * aki - s * akj;
            A[k * m + pj] = s * aki + c * akj;
        }
        for (size_t k = 0; k < m; ++k) {
            double aik = A[pi * m + k];
            double ajk = A[pj * m + k];
            A[pi * m + k] = c * aik - s * ajk;
            A[pj * m + k] = s * aik + c * ajk;
        }
        /* Accumulate the rotation into the eigenvector matrix. */
        for (size_t k = 0; k < m; ++k) {
            double vki = evec[k * m + pi];
            double vkj = evec[k * m + pj];
            evec[k * m + pi] = c * vki - s * vkj;
            evec[k * m + pj] = s * vki + c * vkj;
        }
    }

    for (size_t i = 0; i < m; ++i)
        eval[i] = A[i * m + i];
    return 0;
}

/* ===================================================================
 * Covariance matrix for the subspace methods
 *
 * Builds a real msize x msize covariance estimate by averaging outer
 * products of length-msize data snapshots. (A real symmetric form is
 * used so the Jacobi solver applies; it resolves frequencies in
 * [0, 0.5) without the complex steering vectors of array processing.)
 * =================================================================== */
static int build_covariance(const double *x, size_t n, size_t m,
                             double *R) {
    if (n < m)
        return -1;
    for (size_t i = 0; i < m * m; ++i)
        R[i] = 0.0;

    size_t snaps = n - m + 1;
    for (size_t s = 0; s < snaps; ++s) {
        const double *snap = x + s;
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j)
                R[i * m + j] += snap[i] * snap[j];
    }
    for (size_t i = 0; i < m * m; ++i)
        R[i] /= (double)snaps;
    return 0;
}

/* Sort eigenpairs by descending eigenvalue (small m -> simple sort). */
static void sort_eig_desc(double *eval, double *evec, size_t m) {
    for (size_t i = 0; i < m; ++i) {
        size_t best = i;
        for (size_t j = i + 1; j < m; ++j)
            if (eval[j] > eval[best]) best = j;
        if (best != i) {
            double t = eval[i]; eval[i] = eval[best]; eval[best] = t;
            for (size_t k = 0; k < m; ++k) {
                double tv = evec[k * m + i];
                evec[k * m + i] = evec[k * m + best];
                evec[k * m + best] = tv;
            }
        }
    }
}

/* ===================================================================
 * MUSIC
 * =================================================================== */

int dsp_music(const double *x, size_t n, size_t nsources,
              size_t msize, double *pseudo, size_t nf) {
    if (nsources == 0 || msize <= 2 * nsources || n < msize)
        return -1;

    double *R    = malloc(msize * msize * sizeof(double));
    double *eval = malloc(msize * sizeof(double));
    double *evec = malloc(msize * msize * sizeof(double));
    if (!R || !eval || !evec) {
        free(R); free(eval); free(evec);
        return -1;
    }

    build_covariance(x, n, msize, R);
    jacobi_eig(R, msize, eval, evec);
    sort_eig_desc(eval, evec, msize);

    /* The signal subspace is the top `2*nsources` eigenvectors (each
     * real sinusoid spans two dimensions); the rest are noise. The
     * MUSIC pseudospectrum peaks where a steering vector is
     * orthogonal to the noise subspace. */
    size_t sig_dim = 2 * nsources;

    for (size_t fi = 0; fi < nf; ++fi) {
        double freq = 0.5 * (double)fi / (double)nf;

        /* Real/imag parts of the steering vector e(f). */
        double denom = 0.0;
        for (size_t v = sig_dim; v < msize; ++v) {
            double dot_re = 0.0, dot_im = 0.0;
            for (size_t k = 0; k < msize; ++k) {
                double ang = 2.0 * M_PI * freq * (double)k;
                dot_re += evec[k * msize + v] * cos(ang);
                dot_im += evec[k * msize + v] * sin(ang);
            }
            denom += dot_re * dot_re + dot_im * dot_im;
        }
        pseudo[fi] = (denom > 1e-300) ? 1.0 / denom : 1e300;
    }

    free(R); free(eval); free(evec);
    return 0;
}

/* ===================================================================
 * ESPRIT
 * =================================================================== */

int dsp_esprit(const double *x, size_t n, size_t nsources,
               size_t msize, double *freqs) {
    if (nsources == 0 || msize <= 2 * nsources || n < msize)
        return -1;

    double *R    = malloc(msize * msize * sizeof(double));
    double *eval = malloc(msize * sizeof(double));
    double *evec = malloc(msize * msize * sizeof(double));
    if (!R || !eval || !evec) {
        free(R); free(eval); free(evec);
        return -1;
    }

    build_covariance(x, n, msize, R);
    jacobi_eig(R, msize, eval, evec);
    sort_eig_desc(eval, evec, msize);

    /* ESPRIT exploits the shift-invariance of the signal subspace.
     * Take the top sig_dim signal eigenvectors as the columns of U.
     * Two maximally-overlapping subarrays are U with its last row
     * dropped (U1) and U with its first row dropped (U2). They are
     * related by a rotation: U1 * Psi = U2. The eigenvalues of Psi
     * are e^{j*omega_i}, so their angles are the frequencies.
     *
     * For a real-valued signal each tone contributes a complex-
     * conjugate pair of eigenvalues, so a sig_dim = 2*nsources
     * subspace yields nsources distinct angles. */
    size_t d   = 2 * nsources;           /* signal-subspace dimension */
    size_t sub = msize - 1;              /* subarray row count        */

    /* Least-squares solve U1 * Psi = U2 for the d x d matrix Psi:
     *   Psi = (U1' U1)^-1 (U1' U2).
     * Build the d x d Gram matrix G = U1'U1 and C = U1'U2. */
    double *G   = calloc(d * d, sizeof(double));
    double *C   = calloc(d * d, sizeof(double));
    double *Psi = calloc(d * d, sizeof(double));
    if (!G || !C || !Psi) {
        free(R); free(eval); free(evec);
        free(G); free(C); free(Psi);
        return -1;
    }
    for (size_t i = 0; i < sub; ++i) {
        for (size_t r = 0; r < d; ++r) {
            double u1r = evec[i * msize + r];        /* U1 row i */
            for (size_t c = 0; c < d; ++c) {
                double u1c = evec[i * msize + c];
                double u2c = evec[(i + 1) * msize + c];
                G[r * d + c] += u1r * u1c;
                C[r * d + c] += u1r * u2c;
            }
        }
    }

    /* Solve G * Psi = C column by column via Gaussian elimination. */
    double *aug = malloc(d * d * sizeof(double));
    if (!aug) {
        free(R); free(eval); free(evec);
        free(G); free(C); free(Psi);
        return -1;
    }
    memcpy(aug, G, d * d * sizeof(double));
    for (size_t col = 0; col < d; ++col) {
        size_t piv = col;
        double best = fabs(aug[col * d + col]);
        for (size_t row = col + 1; row < d; ++row) {
            double v = fabs(aug[row * d + col]);
            if (v > best) { best = v; piv = row; }
        }
        if (best < 1e-12) {
            free(R); free(eval); free(evec);
            free(G); free(C); free(Psi); free(aug);
            return -1;
        }
        if (piv != col) {
            for (size_t k = 0; k < d; ++k) {
                double t = aug[col * d + k];
                aug[col * d + k] = aug[piv * d + k];
                aug[piv * d + k] = t;
            }
            for (size_t k = 0; k < d; ++k) {
                double t = C[col * d + k];
                C[col * d + k] = C[piv * d + k];
                C[piv * d + k] = t;
            }
        }
        for (size_t row = 0; row < d; ++row) {
            if (row == col) continue;
            double fct = aug[row * d + col] / aug[col * d + col];
            for (size_t k = 0; k < d; ++k) {
                aug[row * d + k] -= fct * aug[col * d + k];
                C[row * d + k]   -= fct * C[col * d + k];
            }
        }
    }
    for (size_t r = 0; r < d; ++r)
        for (size_t c = 0; c < d; ++c)
            Psi[r * d + c] = C[r * d + c] / aug[r * d + r];

    /* The eigenvalues of Psi lie on the unit circle as e^{j*omega}.
     * For a real signal they form conjugate pairs e^{+-j*omega}, so
     * the SYMMETRIC part S = (Psi + Psi^T)/2 has real eigenvalues
     * cos(omega) - each tone appearing as a (near-)doubled value.
     * S is symmetric, so the Jacobi solver recovers those directly,
     * sidestepping a non-symmetric eigenproblem entirely. */
    double *S    = malloc(d * d * sizeof(double));
    double *sev  = malloc(d * sizeof(double));
    double *svec = malloc(d * d * sizeof(double));
    if (!S || !sev || !svec) {
        free(R); free(eval); free(evec);
        free(G); free(C); free(Psi); free(aug);
        free(S); free(sev); free(svec);
        return -1;
    }
    for (size_t r = 0; r < d; ++r)
        for (size_t c = 0; c < d; ++c)
            S[r * d + c] = 0.5 * (Psi[r * d + c] + Psi[c * d + r]);

    jacobi_eig(S, d, sev, svec);
    sort_eig_desc(sev, svec, d);

    /* Each tone yields a conjugate pair, so the symmetric part has
     * two eigenvalues clustered near cos(omega). Average each adjacent
     * sorted pair to get one robust representative per tone, then
     * convert with omega = acos(.). */
    int found = 0;
    for (size_t s = 0; s + 1 < d && found < (int)nsources; s += 2) {
        double c = 0.5 * (sev[s] + sev[s + 1]);
        if (c >  1.0) c =  1.0;             /* clamp for acos */
        if (c < -1.0) c = -1.0;
        double omega = acos(c);
        freqs[found++] = omega / (2.0 * M_PI);
    }

    free(R); free(eval); free(evec);
    free(G); free(C); free(Psi); free(aug);
    free(S); free(sev); free(svec);
    return found;
}
