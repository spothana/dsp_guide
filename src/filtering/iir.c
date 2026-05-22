/*
 * iir.c - Direct-form-I biquad: design, apply, and stability test.
 */
#include "filtering/iir.h"

void dsp_iir_design_lowpass(dsp_biquad *bq, double fc, double q) {
    /* RBJ cookbook low-pass biquad. */
    double w0    = 2.0 * M_PI * fc;
    double cosw0 = cos(w0);
    double sinw0 = sin(w0);
    double alpha = sinw0 / (2.0 * q);

    /* Unnormalised coefficients. */
    double b0 = (1.0 - cosw0) / 2.0;
    double b1 =  1.0 - cosw0;
    double b2 = (1.0 - cosw0) / 2.0;
    double a0 =  1.0 + alpha;
    double a1 = -2.0 * cosw0;
    double a2 =  1.0 - alpha;

    /* Normalise so a0 == 1. */
    bq->b0 = b0 / a0;
    bq->b1 = b1 / a0;
    bq->b2 = b2 / a0;
    bq->a1 = a1 / a0;
    bq->a2 = a2 / a0;
}

void dsp_iir_apply(const dsp_biquad *bq,
                   const double *x, double *y, size_t n) {
    /* Direct form I delay registers, all starting at zero. */
    double x1 = 0.0, x2 = 0.0;   /* x[n-1], x[n-2] */
    double y1 = 0.0, y2 = 0.0;   /* y[n-1], y[n-2] */

    for (size_t i = 0; i < n; ++i) {
        double xn = x[i];
        double yn = bq->b0 * xn + bq->b1 * x1 + bq->b2 * x2
                                 - bq->a1 * y1 - bq->a2 * y2;
        y[i] = yn;

        /* Shift the delay line. */
        x2 = x1; x1 = xn;
        y2 = y1; y1 = yn;
    }
}

int dsp_iir_is_stable(const dsp_biquad *bq) {
    /* Schur-Cohn / triangle stability region for a 2nd-order section. */
    double a1 = bq->a1, a2 = bq->a2;
    if (fabs(a2) >= 1.0) return 0;
    if (fabs(a1) >= 1.0 + a2) return 0;
    return 1;
}
