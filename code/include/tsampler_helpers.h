#ifndef TSAMPLER_HELPERS_H
#define TSAMPLER_HELPERS_H

/*************************************************
* Internal TSampler helper functions shared
* between tsampler_sampling.c and tsampler_sign.c.
*
* These implement base-b arithmetic, Gaussian
* coset sampling, and polynomial inversion that
* are specific to the TSampler's dual-gadget
* approach and not present in the AC project.
**************************************************/

#include "params.h"
#include "tsampler_params.h"
#include "arith.h"

/* Base-b modular reduction (coefficient-wise, result in [0, b-1]) */
void tsampler_poly_q_mod_bl(poly_q res, const poly_q arg);
void tsampler_poly_q_mod_bh(poly_q res, const poly_q arg);

/* Exact division by base-b (coefficient-wise) */
void tsampler_poly_q_div_bl(poly_q res, const poly_q arg);
void tsampler_poly_q_div_bh(poly_q res, const poly_q arg);

/* Discrete Gaussian coset sampling: sample from D_{b*Z + center, sigma} */
void tsampler_poly_q_gaussian_coset_sl(poly_q res, const poly_q arg);
void tsampler_poly_q_gaussian_coset_sh(poly_q res, const poly_q arg);

/* Invert polynomial modulo bH in Z[x]/(bH, x^n+1) */
void tsampler_poly_q_invert_mod_bh(poly_q res, const poly_q arg);

/* Sample scalar polynomial from D_{Z, 0, s4} */
void tsampler_poly_q_sample_gaussian_s4(poly_q res);

/* Sample d-vector from discrete Gaussian perturbations */
void tsampler_poly_q_vec_d_gaussian_sqrt_s2sq_sLsq(poly_q_vec_d res);
void tsampler_poly_q_vec_d_gaussian_sqrt_s4sq_sHsq(poly_q_vec_d res);

#endif /* TSAMPLER_HELPERS_H */
