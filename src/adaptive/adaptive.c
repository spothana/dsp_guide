/*
 * adaptive.c - LMS, NLMS, and RLS adaptive filters.
 *
 * All three share the same structure: a tapped delay line of recent
 * inputs, a weight vector, and a per-sample update rule. They differ
 * only in how the weights are updated from the error.
 */
#include "adaptive/adaptive.h"
#include <stdlib.h>
#include <string.h>

/* Shift a new sample into a delay line; newest at index 0. */
static void push_line(double *line, size_t n, double sample) {
    for (size_t i = n - 1; i > 0; --i)
        line[i] = line[i - 1];
    line[0] = sample;
}

/* Dot product of weights and delay line: the filter output. */
static double dot(const double *a, const double *b, size_t n) {
    double s = 0.0;
    for (size_t i = 0; i < n; ++i)
        s += a[i] * b[i];
    return s;
}

/* ===================================================================
 * LMS - Least Mean Squares
 * =================================================================== */

int dsp_lms_init(dsp_lms *f, size_t ntaps, double mu) {
    if (ntaps == 0) return -1;
    f->weights = calloc(ntaps, sizeof(double));
    f->line    = calloc(ntaps, sizeof(double));
    if (!f->weights || !f->line) {
        free(f->weights); free(f->line);
        f->weights = f->line = NULL;
        return -1;
    }
    f->ntaps = ntaps;
    f->mu    = mu;
    return 0;
}

void dsp_lms_free(dsp_lms *f) {
    free(f->weights);
    free(f->line);
    f->weights = f->line = NULL;
    f->ntaps = 0;
}

double dsp_lms_update(dsp_lms *f, double input, double desired,
                      double *err_out) {
    push_line(f->line, f->ntaps, input);

    double y   = dot(f->weights, f->line, f->ntaps);
    double err = desired - y;

    /* Stochastic-gradient step: w += mu * error * x. */
    for (size_t i = 0; i < f->ntaps; ++i)
        f->weights[i] += f->mu * err * f->line[i];

    if (err_out) *err_out = err;
    return y;
}

double dsp_lms_train(dsp_lms *f, const double *input,
                     const double *desired, size_t n) {
    double mse = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double err;
        dsp_lms_update(f, input[i], desired[i], &err);
        mse += err * err;
    }
    return (n > 0) ? mse / (double)n : 0.0;
}

/* ===================================================================
 * NLMS - Normalized Least Mean Squares
 * =================================================================== */

int dsp_nlms_init(dsp_nlms *f, size_t ntaps, double mu, double eps) {
    if (ntaps == 0) return -1;
    f->weights = calloc(ntaps, sizeof(double));
    f->line    = calloc(ntaps, sizeof(double));
    if (!f->weights || !f->line) {
        free(f->weights); free(f->line);
        f->weights = f->line = NULL;
        return -1;
    }
    f->ntaps = ntaps;
    f->mu    = mu;
    f->eps   = (eps > 0.0) ? eps : 1e-6;
    return 0;
}

void dsp_nlms_free(dsp_nlms *f) {
    free(f->weights);
    free(f->line);
    f->weights = f->line = NULL;
    f->ntaps = 0;
}

double dsp_nlms_update(dsp_nlms *f, double input, double desired,
                       double *err_out) {
    push_line(f->line, f->ntaps, input);

    double y   = dot(f->weights, f->line, f->ntaps);
    double err = desired - y;

    /* Normalise the step by the input power ||x||^2. This makes the
     * effective learning rate independent of the input signal level,
     * the one weakness of plain LMS. */
    double power = dot(f->line, f->line, f->ntaps);
    double step  = f->mu / (power + f->eps);

    for (size_t i = 0; i < f->ntaps; ++i)
        f->weights[i] += step * err * f->line[i];

    if (err_out) *err_out = err;
    return y;
}

double dsp_nlms_train(dsp_nlms *f, const double *input,
                      const double *desired, size_t n) {
    double mse = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double err;
        dsp_nlms_update(f, input[i], desired[i], &err);
        mse += err * err;
    }
    return (n > 0) ? mse / (double)n : 0.0;
}

/* ===================================================================
 * RLS - Recursive Least Squares
 *
 * RLS keeps P, an estimate of the inverse input-correlation matrix.
 * Each sample it forms the Kalman gain k = P x / (lambda + x' P x),
 * updates the weights by k * error, and updates P. The result is a
 * least-squares fit over all past samples (exponentially weighted by
 * lambda), which converges far faster than gradient methods at the
 * cost of the O(L^2) matrix work.
 * =================================================================== */

int dsp_rls_init(dsp_rls *f, size_t ntaps, double lambda,
                 double delta) {
    if (ntaps == 0 || lambda <= 0.0 || lambda > 1.0)
        return -1;

    f->weights = calloc(ntaps, sizeof(double));
    f->line    = calloc(ntaps, sizeof(double));
    f->P       = calloc(ntaps * ntaps, sizeof(double));
    f->gain    = calloc(ntaps, sizeof(double));
    f->tmp     = calloc(ntaps, sizeof(double));
    if (!f->weights || !f->line || !f->P || !f->gain || !f->tmp) {
        free(f->weights); free(f->line); free(f->P);
        free(f->gain); free(f->tmp);
        f->weights = f->line = f->P = f->gain = f->tmp = NULL;
        return -1;
    }

    f->ntaps  = ntaps;
    f->lambda = lambda;
    f->delta  = (delta > 0.0) ? delta : 100.0;

    /* Initialise P = delta * I. A large delta gives a gentle,
     * numerically stable start. */
    for (size_t i = 0; i < ntaps; ++i)
        f->P[i * ntaps + i] = f->delta;

    return 0;
}

void dsp_rls_free(dsp_rls *f) {
    free(f->weights);
    free(f->line);
    free(f->P);
    free(f->gain);
    free(f->tmp);
    f->weights = f->line = f->P = f->gain = f->tmp = NULL;
    f->ntaps = 0;
}

double dsp_rls_update(dsp_rls *f, double input, double desired,
                      double *err_out) {
    size_t L = f->ntaps;
    push_line(f->line, L, input);

    /* a priori output and error, using the weights from last step. */
    double y   = dot(f->weights, f->line, L);
    double err = desired - y;

    /* tmp = P * x   (matrix-vector product, O(L^2)). */
    for (size_t i = 0; i < L; ++i) {
        double s = 0.0;
        for (size_t j = 0; j < L; ++j)
            s += f->P[i * L + j] * f->line[j];
        f->tmp[i] = s;
    }

    /* denominator = lambda + x' * P * x. */
    double denom = f->lambda + dot(f->line, f->tmp, L);
    if (denom < 1e-12) denom = 1e-12;

    /* Kalman gain k = (P x) / denom. */
    for (size_t i = 0; i < L; ++i)
        f->gain[i] = f->tmp[i] / denom;

    /* Weight update: w += k * error. */
    for (size_t i = 0; i < L; ++i)
        f->weights[i] += f->gain[i] * err;

    /* P update: P = (P - k (P x)') / lambda. tmp still holds P x, so
     * the outer product k * tmp' is the rank-1 correction. */
    double inv_lambda = 1.0 / f->lambda;
    for (size_t i = 0; i < L; ++i)
        for (size_t j = 0; j < L; ++j)
            f->P[i * L + j] =
                (f->P[i * L + j] - f->gain[i] * f->tmp[j]) * inv_lambda;

    if (err_out) *err_out = err;
    return y;
}

double dsp_rls_train(dsp_rls *f, const double *input,
                     const double *desired, size_t n) {
    double mse = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double err;
        dsp_rls_update(f, input[i], desired[i], &err);
        mse += err * err;
    }
    return (n > 0) ? mse / (double)n : 0.0;
}
