/*
 * kalman.c - Linear Kalman filter, EKF, and tracking helpers.
 *
 * Structure:
 *   - small dense matrix helpers (multiply, transpose, add, invert)
 *   - the linear Kalman predict / update cycle
 *   - the Extended Kalman Filter (nonlinear, via caller Jacobians)
 *   - a constant-velocity tracker set-up and a sensor-fusion helper
 */
#include "adaptive/kalman.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ===================================================================
 * Small dense matrix helpers (row-major)
 * =================================================================== */

/* C = A (r x k) * B (k x c). */
static void mat_mul(const double *A, const double *B, double *C,
                    size_t r, size_t k, size_t c) {
    for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < c; ++j) {
            double s = 0.0;
            for (size_t t = 0; t < k; ++t)
                s += A[i * k + t] * B[t * c + j];
            C[i * c + j] = s;
        }
}

/* C = A (r x k) * B^T, where B is (c x k). */
static void mat_mul_bt(const double *A, const double *B, double *C,
                       size_t r, size_t k, size_t c) {
    for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < c; ++j) {
            double s = 0.0;
            for (size_t t = 0; t < k; ++t)
                s += A[i * k + t] * B[j * k + t];
            C[i * c + j] = s;
        }
}

/* y = A (r x c) * v (length c). */
static void mat_vec(const double *A, const double *v, double *y,
                    size_t r, size_t c) {
    for (size_t i = 0; i < r; ++i) {
        double s = 0.0;
        for (size_t j = 0; j < c; ++j)
            s += A[i * c + j] * v[j];
        y[i] = s;
    }
}

/* C = A + B, element-wise, n elements. */
static void mat_add(const double *A, const double *B, double *C,
                    size_t n) {
    for (size_t i = 0; i < n; ++i)
        C[i] = A[i] + B[i];
}

/*
 * Invert an m x m matrix by Gauss-Jordan with partial pivoting.
 * Returns 0 on success, -1 if singular. `inv` may not alias `A`.
 */
static int mat_invert(const double *A, double *inv, size_t m) {
    double *aug = malloc(m * m * sizeof(double));
    if (!aug)
        return -1;
    memcpy(aug, A, m * m * sizeof(double));

    /* inv starts as the identity. */
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < m; ++j)
            inv[i * m + j] = (i == j) ? 1.0 : 0.0;

    for (size_t col = 0; col < m; ++col) {
        /* Partial pivot. */
        size_t piv = col;
        double best = fabs(aug[col * m + col]);
        for (size_t row = col + 1; row < m; ++row) {
            double v = fabs(aug[row * m + col]);
            if (v > best) { best = v; piv = row; }
        }
        if (best < 1e-14) { free(aug); return -1; }
        if (piv != col) {
            for (size_t k = 0; k < m; ++k) {
                double t = aug[col * m + k];
                aug[col * m + k] = aug[piv * m + k];
                aug[piv * m + k] = t;
                t = inv[col * m + k];
                inv[col * m + k] = inv[piv * m + k];
                inv[piv * m + k] = t;
            }
        }
        /* Scale the pivot row to a unit diagonal. */
        double d = aug[col * m + col];
        for (size_t k = 0; k < m; ++k) {
            aug[col * m + k] /= d;
            inv[col * m + k] /= d;
        }
        /* Eliminate the column from every other row. */
        for (size_t row = 0; row < m; ++row) {
            if (row == col) continue;
            double f = aug[row * m + col];
            for (size_t k = 0; k < m; ++k) {
                aug[row * m + k] -= f * aug[col * m + k];
                inv[row * m + k] -= f * inv[col * m + k];
            }
        }
    }
    free(aug);
    return 0;
}

/* ===================================================================
 * Linear Kalman filter
 * =================================================================== */

int dsp_kalman_init(dsp_kalman *kf, size_t nstate, size_t nmeas) {
    if (nstate == 0 || nmeas == 0)
        return -1;
    kf->nstate = nstate;
    kf->nmeas  = nmeas;
    kf->x = calloc(nstate, sizeof(double));
    kf->P = calloc(nstate * nstate, sizeof(double));
    kf->F = calloc(nstate * nstate, sizeof(double));
    kf->Q = calloc(nstate * nstate, sizeof(double));
    kf->H = calloc(nmeas  * nstate, sizeof(double));
    kf->R = calloc(nmeas  * nmeas,  sizeof(double));
    if (!kf->x || !kf->P || !kf->F || !kf->Q || !kf->H || !kf->R) {
        dsp_kalman_free(kf);
        return -1;
    }
    return 0;
}

void dsp_kalman_free(dsp_kalman *kf) {
    free(kf->x); free(kf->P); free(kf->F);
    free(kf->Q); free(kf->H); free(kf->R);
    kf->x = kf->P = kf->F = kf->Q = kf->H = kf->R = NULL;
    kf->nstate = kf->nmeas = 0;
}

void dsp_kalman_predict(dsp_kalman *kf) {
    size_t n = kf->nstate;

    /* x <- F x */
    double *xnew = malloc(n * sizeof(double));
    if (!xnew) return;
    mat_vec(kf->F, kf->x, xnew, n, n);
    memcpy(kf->x, xnew, n * sizeof(double));
    free(xnew);

    /* P <- F P F^T + Q */
    double *FP   = malloc(n * n * sizeof(double));
    double *FPFt = malloc(n * n * sizeof(double));
    if (!FP || !FPFt) { free(FP); free(FPFt); return; }
    mat_mul(kf->F, kf->P, FP, n, n, n);
    mat_mul_bt(FP, kf->F, FPFt, n, n, n);
    mat_add(FPFt, kf->Q, kf->P, n * n);
    free(FP); free(FPFt);
}

/*
 * Shared measurement update. The state has already been predicted;
 * `H` is the measurement matrix (or the EKF's measurement Jacobian),
 * and `hx` is the predicted measurement h(x) - which for the linear
 * filter is simply H x. Corrects x and P in place.
 */
static int kalman_update_core(double *x, double *P,
                              const double *H, const double *R,
                              const double *z, const double *hx,
                              size_t n, size_t m) {
    /* Innovation  y = z - h(x). */
    double *y = malloc(m * sizeof(double));
    /* Innovation covariance  S = H P H^T + R. */
    double *HP   = malloc(m * n * sizeof(double));
    double *HPHt = malloc(m * m * sizeof(double));
    double *S    = malloc(m * m * sizeof(double));
    double *Sinv = malloc(m * m * sizeof(double));
    /* Kalman gain  K = P H^T S^-1. */
    double *PHt  = malloc(n * m * sizeof(double));
    double *K    = malloc(n * m * sizeof(double));
    double *Ky   = malloc(n * sizeof(double));
    double *KH   = malloc(n * n * sizeof(double));
    double *KHP  = malloc(n * n * sizeof(double));
    if (!y || !HP || !HPHt || !S || !Sinv || !PHt || !K || !Ky
        || !KH || !KHP) {
        free(y); free(HP); free(HPHt); free(S); free(Sinv);
        free(PHt); free(K); free(Ky); free(KH); free(KHP);
        return -1;
    }

    for (size_t i = 0; i < m; ++i)
        y[i] = z[i] - hx[i];

    mat_mul(H, P, HP, m, n, n);              /* H P            */
    mat_mul_bt(HP, H, HPHt, m, n, m);        /* H P H^T        */
    mat_add(HPHt, R, S, m * m);              /* S = H P H^T + R*/

    if (mat_invert(S, Sinv, m) != 0) {
        free(y); free(HP); free(HPHt); free(S); free(Sinv);
        free(PHt); free(K); free(Ky); free(KH); free(KHP);
        return -1;
    }

    mat_mul_bt(P, H, PHt, n, n, m);          /* P H^T          */
    mat_mul(PHt, Sinv, K, n, m, m);          /* K = P H^T S^-1 */

    /* State correction  x <- x + K y. */
    mat_vec(K, y, Ky, n, m);
    for (size_t i = 0; i < n; ++i)
        x[i] += Ky[i];

    /* Covariance update  P <- (I - K H) P. */
    mat_mul(K, H, KH, n, m, n);              /* K H            */
    mat_mul(KH, P, KHP, n, n, n);            /* K H P          */
    for (size_t i = 0; i < n * n; ++i)
        P[i] -= KHP[i];

    free(y); free(HP); free(HPHt); free(S); free(Sinv);
    free(PHt); free(K); free(Ky); free(KH); free(KHP);
    return 0;
}

int dsp_kalman_update(dsp_kalman *kf, const double *z) {
    size_t n = kf->nstate, m = kf->nmeas;

    /* For the linear filter the predicted measurement is h(x) = H x. */
    double *hx = malloc(m * sizeof(double));
    if (!hx) return -1;
    mat_vec(kf->H, kf->x, hx, m, n);

    int rc = kalman_update_core(kf->x, kf->P, kf->H, kf->R,
                                z, hx, n, m);
    free(hx);
    return rc;
}

/* ===================================================================
 * Extended Kalman Filter
 * =================================================================== */

int dsp_ekf_init(dsp_ekf *ekf, size_t nstate, size_t nmeas) {
    /* The EKF reuses the linear filter's storage; F and H hold the
     * Jacobians, refreshed each step from the caller's callbacks. */
    if (dsp_kalman_init(&ekf->base, nstate, nmeas) != 0)
        return -1;
    ekf->f = NULL; ekf->h = NULL;
    ekf->Fjac = NULL; ekf->Hjac = NULL;
    ekf->user = NULL;
    return 0;
}

void dsp_ekf_free(dsp_ekf *ekf) {
    dsp_kalman_free(&ekf->base);
}

void dsp_ekf_predict(dsp_ekf *ekf) {
    dsp_kalman *kf = &ekf->base;
    size_t n = kf->nstate;

    /* Refresh the transition Jacobian at the current state, then
     * propagate the covariance through it:  P <- Fjac P Fjac^T + Q. */
    ekf->Fjac(kf->x, n, kf->F, ekf->user);

    double *FP   = malloc(n * n * sizeof(double));
    double *FPFt = malloc(n * n * sizeof(double));
    double *xnew = malloc(n * sizeof(double));
    if (!FP || !FPFt || !xnew) {
        free(FP); free(FPFt); free(xnew);
        return;
    }
    mat_mul(kf->F, kf->P, FP, n, n, n);
    mat_mul_bt(FP, kf->F, FPFt, n, n, n);
    mat_add(FPFt, kf->Q, kf->P, n * n);

    /* Propagate the state through the true nonlinear transition f. */
    ekf->f(kf->x, n, xnew, ekf->user);
    memcpy(kf->x, xnew, n * sizeof(double));

    free(FP); free(FPFt); free(xnew);
}

int dsp_ekf_update(dsp_ekf *ekf, const double *z) {
    dsp_kalman *kf = &ekf->base;
    size_t n = kf->nstate, m = kf->nmeas;

    /* Refresh the measurement Jacobian and the predicted measurement
     * from the nonlinear measurement model. */
    ekf->Hjac(kf->x, n, kf->H, m, ekf->user);

    double *hx = malloc(m * sizeof(double));
    if (!hx) return -1;
    ekf->h(kf->x, n, hx, m, ekf->user);

    int rc = kalman_update_core(kf->x, kf->P, kf->H, kf->R,
                                z, hx, n, m);
    free(hx);
    return rc;
}

/* ===================================================================
 * Tracking helpers
 * =================================================================== */

int dsp_kalman_tracker_init(dsp_kalman *kf, size_t ndim, double dt,
                            double proc_std, double meas_std) {
    if (ndim == 0 || ndim > 3)
        return -1;

    /* State = [pos_0, vel_0, pos_1, vel_1, ...], 2 entries per axis;
     * measurement = the position of each axis. */
    size_t n = 2 * ndim;
    size_t m = ndim;
    if (dsp_kalman_init(kf, n, m) != 0)
        return -1;

    /* State transition: a constant-velocity model per axis,
     *   pos <- pos + vel*dt,   vel <- vel.
     * F is block-diagonal with [[1, dt], [0, 1]] blocks. */
    for (size_t d = 0; d < ndim; ++d) {
        size_t p = 2 * d;            /* position index */
        size_t v = 2 * d + 1;        /* velocity index */
        kf->F[p * n + p] = 1.0;
        kf->F[p * n + v] = dt;
        kf->F[v * n + v] = 1.0;
    }

    /* Process noise: a small acceleration disturbance, the standard
     * constant-velocity Q built from proc_std. Per axis:
     *   [[dt^4/4, dt^3/2], [dt^3/2, dt^2]] * proc_std^2. */
    double q = proc_std * proc_std;
    double dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt3 * dt;
    for (size_t d = 0; d < ndim; ++d) {
        size_t p = 2 * d, v = 2 * d + 1;
        kf->Q[p * n + p] = q * dt4 / 4.0;
        kf->Q[p * n + v] = q * dt3 / 2.0;
        kf->Q[v * n + p] = q * dt3 / 2.0;
        kf->Q[v * n + v] = q * dt2;
    }

    /* Measurement matrix: each axis's position is observed directly. */
    for (size_t d = 0; d < ndim; ++d)
        kf->H[d * n + 2 * d] = 1.0;

    /* Measurement-noise covariance: diagonal, from meas_std. */
    double r = meas_std * meas_std;
    for (size_t d = 0; d < ndim; ++d)
        kf->R[d * m + d] = r;

    return 0;
}

double dsp_kalman_fuse(const double *meas, const double *var,
                       size_t nsensors, double *fused_var) {
    /* Inverse-variance ("precision") weighting: the maximum-likelihood
     * fusion of independent Gaussian measurements. Each sensor
     * contributes in proportion to 1/variance, so an accurate sensor
     * dominates a noisy one, and the fused variance is smaller than
     * any individual variance. This is the static, single-sample form
     * of the Kalman update. */
    double sum_prec = 0.0;
    double sum_wx   = 0.0;
    for (size_t i = 0; i < nsensors; ++i) {
        double v = (var[i] > 1e-300) ? var[i] : 1e-300;
        double prec = 1.0 / v;
        sum_prec += prec;
        sum_wx   += prec * meas[i];
    }
    if (sum_prec <= 0.0) {
        if (fused_var) *fused_var = 0.0;
        return 0.0;
    }
    if (fused_var)
        *fused_var = 1.0 / sum_prec;
    return sum_wx / sum_prec;
}
