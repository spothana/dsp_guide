/*
 * kalman.h - Kalman filtering and state estimation
 *
 * PROBLEM SOLVED
 *   The adaptive filters in adaptive.h LEARN an unknown system. The
 *   Kalman filter solves a different problem: it ESTIMATES the hidden
 *   STATE of a KNOWN dynamic system from a stream of noisy
 *   measurements. "State" is whatever you are tracking - a position
 *   and velocity, an orientation, a temperature.
 *
 *   The filter keeps two things: a state estimate and a covariance
 *   matrix saying how uncertain that estimate is. Each step it runs a
 *   two-phase cycle:
 *
 *     PREDICT - push the state forward by the system's motion model;
 *               uncertainty grows.
 *     UPDATE  - fold in a new measurement; the Kalman gain weighs
 *               model vs measurement by their relative uncertainty,
 *               and the covariance shrinks.
 *
 *   It is the optimal (minimum-mean-square-error) estimator when the
 *   system is linear and the noise is Gaussian.
 *
 * SENSOR FUSION
 *   The headline application: several noisy sensors measuring the
 *   same quantity. The update step naturally COMBINES them - each
 *   measurement pulls the estimate in proportion to its precision -
 *   so the fused estimate beats any single sensor. A GPS + IMU
 *   navigation filter is the classic example.
 *
 * RELATION TO adaptive.h
 *   RLS is, mathematically, a Kalman filter for a particular state
 *   model (the "state" being the filter weights, assumed constant).
 *   The two share the recursive predict/correct structure; the Kalman
 *   filter is the general form.
 *
 * THIS MODULE PROVIDES
 *   - a linear Kalman filter in general matrix form;
 *   - the Extended Kalman Filter (EKF) for nonlinear systems, which
 *     linearises the model at each step via caller-supplied Jacobians;
 *   - ready-made constant-velocity tracking helpers.
 *
 * Matrices are plain row-major double arrays; dimensions are passed
 * explicitly. The filters are written for the small state sizes
 * (typically 2-9) that real tracking and fusion problems use.
 */
#ifndef DSP_KALMAN_H
#define DSP_KALMAN_H

#include <stddef.h>

/* ===================================================================
 * Linear Kalman filter
 * =================================================================== */

/*
 * A linear Kalman filter.
 *   nstate : state-vector dimension
 *   nmeas  : measurement-vector dimension
 *   x      : current state estimate, length nstate
 *   P      : state covariance, nstate x nstate, row-major
 *   F      : state-transition matrix, nstate x nstate
 *   Q      : process-noise covariance, nstate x nstate
 *   H      : measurement matrix, nmeas x nstate
 *   R      : measurement-noise covariance, nmeas x nmeas
 * F, Q, H, R describe the system and are set by the caller (directly
 * or via dsp_kalman_init); x and P evolve as the filter runs.
 */
typedef struct {
    size_t  nstate;
    size_t  nmeas;
    double *x;
    double *P;
    double *F;
    double *Q;
    double *H;
    double *R;
} dsp_kalman;

/*
 * Allocate a Kalman filter for the given dimensions. All matrices are
 * zero-initialised; fill F, Q, H, R and the initial x, P before use.
 * Returns 0 on success, -1 on allocation failure. Pair with
 * dsp_kalman_free.
 */
int dsp_kalman_init(dsp_kalman *kf, size_t nstate, size_t nmeas);

/* Release a Kalman filter's matrices. */
void dsp_kalman_free(dsp_kalman *kf);

/*
 * Predict step: advance the state and grow the covariance by the
 * motion model.
 *   x <- F x
 *   P <- F P F^T + Q
 */
void dsp_kalman_predict(dsp_kalman *kf);

/*
 * Update step: fold in a measurement vector `z` (length nmeas).
 * Computes the Kalman gain, corrects the state, and shrinks the
 * covariance. Returns 0 on success, -1 if the innovation covariance
 * is singular.
 */
int dsp_kalman_update(dsp_kalman *kf, const double *z);

/* ===================================================================
 * Extended Kalman Filter (nonlinear systems)
 * =================================================================== */

/*
 * Caller-supplied model functions for the EKF.
 *   state transition  f : x_next = f(x)
 *   measurement model h : z = h(x)
 * and their Jacobians (partial-derivative matrices), evaluated at the
 * current state:
 *   Fjac : d f / d x, nstate x nstate
 *   Hjac : d h / d x, nmeas  x nstate
 * `user` is an opaque pointer passed straight through, for parameters.
 */
typedef void (*dsp_ekf_transition)(const double *x, size_t nstate,
                                   double *x_next, void *user);
typedef void (*dsp_ekf_measurement)(const double *x, size_t nstate,
                                    double *z, size_t nmeas,
                                    void *user);
typedef void (*dsp_ekf_jacobian_f)(const double *x, size_t nstate,
                                   double *Fjac, void *user);
typedef void (*dsp_ekf_jacobian_h)(const double *x, size_t nstate,
                                   double *Hjac, size_t nmeas,
                                   void *user);

/*
 * An Extended Kalman Filter. It reuses the linear filter's storage for
 * x, P, Q, R (F and H hold the Jacobians, refreshed each step); the
 * nonlinear behaviour comes from the four callbacks.
 */
typedef struct {
    dsp_kalman          base;     /* x, P, Q, R, and the Jacobian stores */
    dsp_ekf_transition  f;
    dsp_ekf_measurement h;
    dsp_ekf_jacobian_f  Fjac;
    dsp_ekf_jacobian_h  Hjac;
    void               *user;     /* opaque, passed to the callbacks */
} dsp_ekf;

/*
 * Allocate an EKF. Set Q, R and the initial x, P after this call, and
 * assign the four model callbacks (and optional `user` pointer).
 * Returns 0 on success, -1 on allocation failure. Pair with
 * dsp_ekf_free.
 */
int dsp_ekf_init(dsp_ekf *ekf, size_t nstate, size_t nmeas);

/* Release an EKF's storage. */
void dsp_ekf_free(dsp_ekf *ekf);

/*
 * EKF predict: propagate the state through the nonlinear transition f,
 * and the covariance through its Jacobian:  P <- Fjac P Fjac^T + Q.
 */
void dsp_ekf_predict(dsp_ekf *ekf);

/*
 * EKF update: fold in measurement `z` (length nmeas). The innovation
 * uses the nonlinear h; the gain and covariance use the Jacobian Hjac.
 * Returns 0 on success, -1 on a singular innovation covariance.
 */
int dsp_ekf_update(dsp_ekf *ekf, const double *z);

/* ===================================================================
 * Ready-made tracking helpers
 *
 * The most common Kalman set-up: a constant-velocity motion model.
 * The state is position and velocity per spatial axis; the filter is
 * fed noisy position measurements and recovers smoothed position AND
 * the unmeasured velocity.
 * =================================================================== */

/*
 * Configure `kf` as a constant-velocity tracker.
 *   ndim       : number of spatial dimensions (1, 2, or 3)
 *   dt         : time step between updates
 *   proc_std   : process-noise standard deviation (how much the
 *                velocity is allowed to wander - model "looseness")
 *   meas_std   : measurement-noise standard deviation of the position
 *                sensor
 * The filter is allocated inside this call with state size 2*ndim
 * (position and velocity per axis) and measurement size ndim. The
 * caller still sets the initial state x and covariance P, then runs
 * the usual predict/update cycle. Returns 0 on success, -1 on error.
 * Pair with dsp_kalman_free.
 */
int dsp_kalman_tracker_init(dsp_kalman *kf, size_t ndim, double dt,
                            double proc_std, double meas_std);

/*
 * Fuse measurements of one quantity from several sensors of differing
 * accuracy into a single best estimate - a one-line sensor-fusion
 * helper, the static (single-sample) form of the Kalman update.
 *   meas      : the sensor readings, length nsensors
 *   var       : each sensor's noise variance, length nsensors
 *   nsensors  : number of sensors
 *   fused_var : if non-NULL, receives the variance of the fused
 *               estimate (smaller than any individual sensor's)
 * Returns the inverse-variance-weighted (maximum-likelihood) estimate.
 */
double dsp_kalman_fuse(const double *meas, const double *var,
                       size_t nsensors, double *fused_var);

#endif /* DSP_KALMAN_H */
