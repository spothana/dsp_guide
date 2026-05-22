/*
 * array.h - Array signal processing: beamforming and DOA estimation
 *
 * PROBLEM SOLVED
 *   A single sensor (microphone, antenna) measures a signal but not
 *   the DIRECTION it came from. An ARRAY of sensors can: a wave
 *   arriving at an angle reaches the sensors at slightly different
 *   times, and that pattern of delays encodes the direction.
 *
 *   Two tasks build on this:
 *     BEAMFORMING - combine the sensor signals so the array "listens"
 *       preferentially in one direction, amplifying a wave from there
 *       and suppressing others. A spatial filter.
 *     DIRECTION OF ARRIVAL (DOA) - estimate the angle(s) the
 *       wave(s) came from.
 *
 * THE UNIFORM LINEAR ARRAY (ULA)
 *   This module models a ULA: M sensors equally spaced by d (in
 *   wavelengths). A plane wave from angle theta reaches sensor k with
 *   a phase advance of 2*pi*d*k*sin(theta) relative to sensor 0. The
 *   vector of those phases is the STEERING VECTOR a(theta) - the
 *   array's complex response to a wave from that angle.
 *
 * THE COVARIANCE MATRIX
 *   With several time snapshots of the M sensors, the M x M sample
 *   covariance matrix R captures the spatial structure. Every method
 *   here works from R:
 *     delay-and-sum - steer a(theta), measure a' R a
 *     MVDR/Capon    - the minimum-variance distortionless beamformer
 *     MUSIC/ESPRIT  - eigen-split R into signal and noise subspaces
 *
 * RELATION TO spectral/estimation.h
 *   MUSIC and ESPRIT also appear there as TEMPORAL frequency
 *   estimators. They are the same algorithms: a temporal frequency
 *   and a spatial direction both reduce to estimating the phase slope
 *   of a complex exponential. The spatial versions here use a complex
 *   Hermitian covariance and complex steering vectors, so they need a
 *   complex eigensolver rather than the real one used there.
 */
#ifndef DSP_ARRAY_H
#define DSP_ARRAY_H

#include "common.h"

/* ===================================================================
 * Array geometry and signal synthesis
 * =================================================================== */

/*
 * Steering vector of a uniform linear array.
 *   nsensors : number of array elements M
 *   spacing  : element spacing d, in wavelengths (0.5 is standard)
 *   theta    : arrival angle in radians, measured from broadside
 *   a        : output steering vector, length nsensors; element k is
 *              exp(j * 2*pi * spacing * k * sin(theta))
 */
void dsp_array_steering(size_t nsensors, double spacing,
                        double theta, cplx *a);

/*
 * Synthesise array snapshots for a set of plane-wave sources plus
 * sensor noise - a test-signal generator for the methods below.
 *   nsensors   : array element count M
 *   spacing    : element spacing in wavelengths
 *   nsnapshots : number of time snapshots to generate
 *   angles     : source arrival angles in radians, length nsources
 *   nsources   : number of sources
 *   noise_std  : standard deviation of the complex sensor noise
 *   seed       : RNG seed
 *   snapshots  : output, nsnapshots * nsensors complex values,
 *                row-major; row t is the M sensor outputs at time t
 */
void dsp_array_synthesize(size_t nsensors, double spacing,
                          size_t nsnapshots,
                          const double *angles, size_t nsources,
                          double noise_std, unsigned seed,
                          cplx *snapshots);

/*
 * Sample spatial covariance matrix from array snapshots.
 *   snapshots  : nsnapshots * nsensors complex values, row-major
 *   nsensors   : array element count M
 *   nsnapshots : number of snapshots
 *   R          : output M x M Hermitian covariance, row-major;
 *                R = (1/nsnapshots) sum_t x(t) x(t)^H
 */
void dsp_array_covariance(const cplx *snapshots, size_t nsensors,
                          size_t nsnapshots, cplx *R);

/* ===================================================================
 * Beamforming
 * =================================================================== */

/*
 * Conventional (delay-and-sum) beamformer power response.
 *   R        : M x M spatial covariance matrix
 *   nsensors : array element count M
 *   spacing  : element spacing in wavelengths
 *   power    : output power vs angle, length nangles, over the angle
 *              grid theta = -pi/2 .. +pi/2
 *   nangles  : number of angle grid points
 * For each angle the response is a(theta)^H R a(theta): steering the
 * array there and measuring the received power. Its peaks indicate
 * source directions, with a resolution limited by the array aperture.
 */
void dsp_beamform_conventional(const cplx *R, size_t nsensors,
                               double spacing,
                               double *power, size_t nangles);

/*
 * MVDR / Capon (minimum-variance distortionless) beamformer power.
 *   R        : M x M spatial covariance matrix
 *   nsensors : array element count M
 *   spacing  : element spacing in wavelengths
 *   power    : output power vs angle, length nangles
 *   nangles  : number of angle grid points
 * The MVDR response is 1 / (a(theta)^H R^-1 a(theta)): it keeps unit
 * gain on the look direction while minimising total output power, so
 * it places sharp nulls on interferers and resolves sources far
 * better than the conventional beamformer.
 * Returns 0 on success, -1 if R is singular.
 */
int dsp_beamform_mvdr(const cplx *R, size_t nsensors, double spacing,
                      double *power, size_t nangles);

/* ===================================================================
 * Direction-of-arrival estimation
 * =================================================================== */

/*
 * Spatial MUSIC pseudospectrum.
 *   R        : M x M spatial covariance matrix
 *   nsensors : array element count M
 *   spacing  : element spacing in wavelengths
 *   nsources : number of sources (signal-subspace dimension)
 *   pseudo   : output pseudospectrum, length nangles, over
 *              theta = -pi/2 .. +pi/2; its PEAKS are the arrival
 *              angles
 *   nangles  : number of angle grid points
 * Eigen-decomposes R, then for each angle measures how orthogonal the
 * steering vector is to the noise subspace. Super-resolution: it
 * separates sources closer than the conventional beamwidth.
 * Returns 0 on success, -1 on bad parameters.
 */
int dsp_doa_music(const cplx *R, size_t nsensors, double spacing,
                  size_t nsources, double *pseudo, size_t nangles);

/*
 * Spatial ESPRIT - direct DOA estimation, no angle search.
 *   R        : M x M spatial covariance matrix
 *   nsensors : array element count M (>= nsources + 1)
 *   spacing  : element spacing in wavelengths
 *   nsources : number of sources
 *   angles   : output arrival angles in radians, length nsources;
 *              not sorted in any particular order
 * Exploits the shift-invariance of the ULA: two overlapping
 * subarrays differ only by a phase rotation, whose eigenvalues give
 * the angles directly.
 * Returns the number of angles found (<= nsources), or -1 on error.
 */
int dsp_doa_esprit(const cplx *R, size_t nsensors, double spacing,
                   size_t nsources, double *angles);

#endif /* DSP_ARRAY_H */
