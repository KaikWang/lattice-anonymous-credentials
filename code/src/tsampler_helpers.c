/*************************************************
* TSampler helper functions — base-b arithmetic,
* Gaussian coset sampling, and polynomial inversion.
*
* Ported from tag-friendly_sampler's poly_q.c.
**************************************************/

#include "tsampler_helpers.h"
#include "random.h"

#include <flint/nmod_poly.h>

/*************************************************
* Name:        tsampler_poly_q_mod_bl
*
* Description: Reduce each coefficient modulo PARAM_BL_TS.
*              Result is in [0, BL_TS-1].
**************************************************/
void tsampler_poly_q_mod_bl(poly_q res, const poly_q arg) {
  coeff_q c;
  for (size_t j = 0; j < PARAM_N; ++j) {
    c = poly_q_get_coeff(arg, j);
    c %= PARAM_BL_TS;
    c += (c >> (sizeof(coeff_q)*8-1)) & PARAM_BL_TS;
    poly_q_set_coeff(res, j, c);
  }
}

/*************************************************
* Name:        tsampler_poly_q_mod_bh
*
* Description: Reduce each coefficient modulo PARAM_BH_TS.
*              Result is in [0, BH_TS-1].
**************************************************/
void tsampler_poly_q_mod_bh(poly_q res, const poly_q arg) {
  coeff_q c;
  for (size_t j = 0; j < PARAM_N; ++j) {
    c = poly_q_get_coeff(arg, j);
    c %= PARAM_BH_TS;
    c += (c >> (sizeof(coeff_q)*8-1)) & PARAM_BH_TS;
    poly_q_set_coeff(res, j, c);
  }
}

/*************************************************
* Name:        tsampler_poly_q_div_bl
*
* Description: Exact integer division of each
*              coefficient by PARAM_BL_TS.
**************************************************/
void tsampler_poly_q_div_bl(poly_q res, const poly_q arg) {
  coeff_q c;
  for (size_t i = 0; i < PARAM_N; i++) {
    c = poly_q_get_coeff(arg, i);
    c /= PARAM_BL_TS;
    poly_q_set_coeff(res, i, c);
  }
}

/*************************************************
* Name:        tsampler_poly_q_div_bh
*
* Description: Exact integer division of each
*              coefficient by PARAM_BH_TS.
**************************************************/
void tsampler_poly_q_div_bh(poly_q res, const poly_q arg) {
  coeff_q c;
  for (size_t i = 0; i < PARAM_N; i++) {
    c = poly_q_get_coeff(arg, i);
    c /= PARAM_BH_TS;
    poly_q_set_coeff(res, i, c);
  }
}

/*************************************************
* Name:        tsampler_poly_q_gaussian_coset_sl
*
* Description: Sample from the discrete Gaussian
*   over bL*Z + center, with width sL/bL.
*   output = center + bL * z
*   where z ~ D_{Z, sigma, -center/bL}
**************************************************/
void tsampler_poly_q_gaussian_coset_sl(poly_q res, const poly_q arg) {
  coeff_q c;
  int64_t z;
  double neg_c_div_bl;
  for (size_t i = 0; i < PARAM_N; i++) {
    c = poly_q_get_coeff(arg, i);
    neg_c_div_bl = -((double)c) / PARAM_BL_TS;
    z = SampleZ(neg_c_div_bl, PARAM_SL_DIV_BL_TS);
    c += (PARAM_BL_TS * z);
    poly_q_set_coeff(res, i, c);
  }
}

/*************************************************
* Name:        tsampler_poly_q_gaussian_coset_sh
*
* Description: Sample from the discrete Gaussian
*   over bH*Z + center, with width sH/bH.
**************************************************/
void tsampler_poly_q_gaussian_coset_sh(poly_q res, const poly_q arg) {
  coeff_q c;
  int64_t z;
  double neg_c_div_bh;
  for (size_t i = 0; i < PARAM_N; i++) {
    c = poly_q_get_coeff(arg, i);
    neg_c_div_bh = -((double)c) / PARAM_BH_TS;
    z = SampleZ(neg_c_div_bh, PARAM_SH_DIV_BH_TS);
    c += (PARAM_BH_TS * z);
    poly_q_set_coeff(res, i, c);
  }
}

/*************************************************
* Name:        tsampler_poly_q_invert_mod_bh
*
* Description: Invert a polynomial modulo PARAM_BH_TS
*   in Z[x]/(bH, x^n+1). Uses FLINT's nmod_poly_xgcd.
*   The tag polynomial must be invertible modulo bH.
**************************************************/
void tsampler_poly_q_invert_mod_bh(poly_q res, const poly_q arg) {
  nmod_poly_t arg_bh, res_bh, f_mod_bh, gcd_poly, aux_poly;
  coeff_q c;

  nmod_poly_init(arg_bh, PARAM_BH_TS);
  nmod_poly_init(res_bh, PARAM_BH_TS);
  nmod_poly_init(f_mod_bh, PARAM_BH_TS);
  nmod_poly_init(gcd_poly, PARAM_BH_TS);
  nmod_poly_init(aux_poly, PARAM_BH_TS);

  /* f_mod_bh = x^n + 1  (mod bH) */
  nmod_poly_zero(f_mod_bh);
  nmod_poly_set_coeff_ui(f_mod_bh, 0, 1);
  nmod_poly_set_coeff_ui(f_mod_bh, PARAM_N, 1);

  /* arg_bh = arg mod bH */
  for (size_t j = 0; j < PARAM_N; j++) {
    c = poly_q_get_coeff(arg, j);
    c %= PARAM_BH_TS;
    c += (c >> (sizeof(coeff_q)*8-1)) & PARAM_BH_TS;
    nmod_poly_set_coeff_ui(arg_bh, j, (ulong)c);
  }

  /* Extended GCD: f_mod_bh * aux + arg_bh * res_bh = gcd */
  nmod_poly_xgcd(gcd_poly, aux_poly, res_bh, f_mod_bh, arg_bh);

  if (!nmod_poly_is_one(gcd_poly)) {
    /* tag is not invertible — set res to zero as signal */
    poly_q_zero(res);
  } else {
    /* convert res_bh back to poly_q */
    for (size_t j = 0; j < PARAM_N; j++) {
      c = (coeff_q)nmod_poly_get_coeff_ui(res_bh, j);
      poly_q_set_coeff(res, j, c);
    }
  }

  nmod_poly_clear(arg_bh);
  nmod_poly_clear(res_bh);
  nmod_poly_clear(f_mod_bh);
  nmod_poly_clear(gcd_poly);
  nmod_poly_clear(aux_poly);
}

/*************************************************
* Name:        tsampler_poly_q_sample_gaussian_s4
*
* Description: Sample a scalar polynomial from
*   the discrete Gaussian D_{Z, 0, s4}.
**************************************************/
void tsampler_poly_q_sample_gaussian_s4(poly_q res) {
  coeff_q cj;
  for (size_t j = 0; j < PARAM_N; j++) {
    cj = SampleZ(0, PARAM_S4_TS);
    poly_q_set_coeff(res, j, cj);
  }
}

/*************************************************
* Name:        tsampler_poly_q_vec_d_gaussian_sqrt_s2sq_sLsq
*
* Description: Sample each entry of a d-vector
*   from D_{Z, 0, sqrt(s2² - sL²)}.
**************************************************/
void tsampler_poly_q_vec_d_gaussian_sqrt_s2sq_sLsq(poly_q_vec_d res) {
  coeff_q cj;
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      cj = SampleZ(0, PARAM_SQRT_S2SQ_SLSQ_TS);
      poly_q_set_coeff(res->entries[i], j, cj);
    }
  }
}

/*************************************************
* Name:        tsampler_poly_q_vec_d_gaussian_sqrt_s4sq_sHsq
*
* Description: Sample each entry of a d-vector
*   from D_{Z, 0, sqrt(s4² - sH²)}.
**************************************************/
void tsampler_poly_q_vec_d_gaussian_sqrt_s4sq_sHsq(poly_q_vec_d res) {
  coeff_q cj;
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      cj = SampleZ(0, PARAM_SQRT_S4SQ_SHSQ_TS);
      poly_q_set_coeff(res->entries[i], j, cj);
    }
  }
}
