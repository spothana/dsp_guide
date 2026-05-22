/*
 * equalize.c - LMS adaptive FIR equalizer.
 *
 * Each call to dsp_lms_update does three things:
 *   1. shift the new sample into the delay line,
 *   2. compute the filter output (a dot product),
 *   3. nudge every tap by mu * error * (the sample it multiplies).
 * Over a training sequence the taps settle toward the channel's
 * inverse, flattening the overall response.
 */
#include "coding/equalize.h"
#include <stdlib.h>
#include <string.h>

int dsp_lms_init(dsp_lms *eq, size_t ntaps, double mu) {
    if (ntaps == 0) return -1;
    eq->weights = calloc(ntaps, sizeof(double));
    eq->line    = calloc(ntaps, sizeof(double));
    if (!eq->weights || !eq->line) {
        free(eq->weights);
        free(eq->line);
        eq->weights = eq->line = NULL;
        return -1;
    }
    eq->ntaps = ntaps;
    eq->mu    = mu;
    /* A common initialisation: a single unit "centre spike" tap, so
     * the equalizer starts as a pass-through and adapts from there. */
    eq->weights[ntaps / 2] = 1.0;
    return 0;
}

void dsp_lms_free(dsp_lms *eq) {
    free(eq->weights);
    free(eq->line);
    eq->weights = eq->line = NULL;
    eq->ntaps = 0;
}

double dsp_lms_update(dsp_lms *eq, double input, double desired,
                      double *err_out) {
    /* 1. Shift the delay line; newest sample at index 0. */
    for (size_t i = eq->ntaps - 1; i > 0; --i)
        eq->line[i] = eq->line[i - 1];
    eq->line[0] = input;

    /* 2. Filter output = weights . line. */
    double y = 0.0;
    for (size_t i = 0; i < eq->ntaps; ++i)
        y += eq->weights[i] * eq->line[i];

    /* 3. LMS tap update along the negative gradient of squared error. */
    double err = desired - y;
    for (size_t i = 0; i < eq->ntaps; ++i)
        eq->weights[i] += eq->mu * err * eq->line[i];

    if (err_out)
        *err_out = err;
    return y;
}

double dsp_lms_train(dsp_lms *eq, const double *input,
                     const double *desired, size_t n) {
    double mse = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double err;
        dsp_lms_update(eq, input[i], desired[i], &err);
        mse += err * err;
    }
    return (n > 0) ? mse / (double)n : 0.0;
}
