/*************************************************
* SEP signature with TSampler-style gadget sampler.
*
* This is a variant of sep.c where the internal
* Klein/ZF gadget sampler is replaced with the
* TSampler's base-decomposition approach.
*
* Reuses: sep_keygen(), sep_verify(),
*         _sep_verify_from_commitment()
* Replaces: sample_klein() → tsampler_gadget_sample_klein_replacement()
**************************************************/

#include "sep_tsampler.h"
#include "tsampler_sampling.h"
#include "arith.h"
#include "randombytes.h"
#include "poly_q_sampling.h"
#include "fips202.h"
#include "covariance.h"

/* Forward declaration from sep.c (internal) */
void _sep_sign_commitment(
    sep_sig_t *sig, uint8_t state[STATE_BYTES],
    const sep_sk_t *sk, const sep_pk_t *pk,
    const poly_q_vec_d cmt);

/* ================================================
 * Modified preimage sampler using TSampler gadget
 * ================================================ */

/*************************************************
* Name:        poly_q_vec_2d_dk_sample_pre_tsampler
*
* Description: Same as poly_q_vec_2d_dk_sample_pre
*   but replaces the Klein sampler (sample_klein)
*   with base-decomposition gadget sampling.
**************************************************/
static void poly_q_vec_2d_dk_sample_pre_tsampler(
    poly_q_vec_d v1[2],
    poly_q_vec_d v2[PARAM_K],
    const poly_q_mat_d_d R[2][PARAM_K],
    const poly_q_mat_d_d A,
    const poly_q_mat_d_d B[PARAM_K],
    const poly_q_vec_d u,
    const poly_q tag,
    const poly_real_mat_2d_2d S) {
  size_t i;
  poly_q_vec_d p1[2];
  poly_q_vec_d p2[PARAM_K];
  poly_q taginv;
  poly_q_vec_d w, tmp, y[PARAM_K];
  int64_t bexpi;

  /* init */
  poly_q_init(taginv);
  poly_q_vec_d_init(w);
  poly_q_vec_d_init(tmp);
  poly_q_vec_d_init(p1[0]);
  poly_q_vec_d_init(p1[1]);
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_init(p2[i]);
    poly_q_vec_d_init(y[i]);
    poly_q_vec_d_zero(y[i]);
  }

  /* Step 1: perturbation sampler (same as original SEP) */
  {
    /* poly_q_vec_2d_dk_sample_perturb is static in poly_q_sampling.c,
     * so we inline the equivalent logic here. */
    size_t jj;
    poly_real tmp_sub, tmp_mul;
    poly_real_vec_2d c;

    poly_real_init(tmp_sub);
    poly_real_init(tmp_mul);
    poly_real_vec_2d_init(c);

    /* sample p2 from Gaussian */
    for (i = 0; i < PARAM_K; i++) {
      poly_q_vec_d_gaussian_sqrt_s2sq_sGsq(p2[i]);
    }

    /* p1 = R * p2 */
    poly_q_mat_d_d_mul_vec_d(p1[0], R[0][0], p2[0]);
    poly_q_mat_d_d_mul_vec_d(p1[1], R[1][0], p2[0]);
    for (i = 1; i < PARAM_K; i++) {
      poly_q_mat_d_d_muladd_vec_d(p1[0], R[0][i], p2[i]);
      poly_q_mat_d_d_muladd_vec_d(p1[1], R[1][i], p2[i]);
    }

    /* elliptical center computation */
    for (i = 0; i < PARAM_D; i++) {
      poly_real_from_poly_q(c->entries[i          ], p1[0]->entries[i]);
      poly_real_from_poly_q(c->entries[i + PARAM_D], p1[1]->entries[i]);
      poly_real_mul_scalar(c->entries[i],
          c->entries[i], PARAM_NEGSGSQ_DIV_S2SQ_SGSQ);
      poly_real_mul_scalar(c->entries[i + PARAM_D],
          c->entries[i + PARAM_D], PARAM_NEGSGSQ_DIV_S2SQ_SGSQ);
    }

    for (i = 2*PARAM_D; i > 0; i--) {
      poly_q_samplefz(p1[(i-1)/PARAM_D]->entries[(i-1)%PARAM_D],
          S->rows[i-1]->entries[i-1], c->entries[i-1]);
      poly_real_sub_poly_real_poly_q(tmp_sub,
          p1[(i-1)/PARAM_D]->entries[(i-1)%PARAM_D], c->entries[i-1]);
      for (jj = 0; jj < i-1; jj++) {
        poly_real_mul(tmp_mul, tmp_sub, S->rows[jj]->entries[i-1]);
        poly_real_add(c->entries[jj], c->entries[jj], tmp_mul);
      }
    }

    poly_real_vec_2d_clear(c);
    poly_real_clear(tmp_mul);
    poly_real_clear(tmp_sub);
  }

  /* Step 2: compute w (same as original SEP) */
  poly_q_invert(taginv, tag);
  poly_q_vec_d_sub(w, u, p1[0]);
  poly_q_mat_d_d_mulsub_vec_d(w, A, p1[1]);
  for (i = 0; i < PARAM_K; i++) {
    poly_q_mat_d_d_muladd_vec_d(w, B[i], p2[i]);
  }
  poly_q_vec_d_mul_poly(w, w, taginv);
  poly_q_vec_d_sub(w, w, p2[0]);
  bexpi = 1;
  for (i = 1; i < PARAM_K; i++) {
    bexpi *= PARAM_B;
    poly_q_vec_d_mul_scalar(tmp, p2[i], bexpi);
    poly_q_vec_d_sub(w, w, tmp);
  }

  /* Step 3: gadget sampling — TSampler base decomposition
   * instead of sample_klein(y, w).
   * The base decomposition directly gives z such that G*z = w,
   * unlike Klein which gives a lattice point y that needs
   * y[0] += w to get z. */
  tsampler_gadget_sample_klein_replacement(y, w);
  /* NO y[0] += w here — base decomp already includes it */

  /* Step 4: combine (same as original SEP) */
  poly_q_mat_d_d_mul_vec_d(v1[0], R[0][0], y[0]);
  poly_q_mat_d_d_mul_vec_d(v1[1], R[1][0], y[0]);
  poly_q_vec_d_add(v1[0], v1[0], p1[0]);
  poly_q_vec_d_add(v1[1], v1[1], p1[1]);
  for (i = 1; i < PARAM_K; i++) {
    poly_q_mat_d_d_muladd_vec_d(v1[0], R[0][i], y[i]);
    poly_q_mat_d_d_muladd_vec_d(v1[1], R[1][i], y[i]);
  }

  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_add(v2[i], p2[i], y[i]);
  }

  /* clean up */
  poly_q_clear(taginv);
  poly_q_vec_d_clear(w);
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(p1[0]);
  poly_q_vec_d_clear(p1[1]);
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_clear(p2[i]);
    poly_q_vec_d_clear(y[i]);
  }
}

/* ================================================
 * Modified sign commitment using TSampler gadget
 * ================================================ */

static void _sep_sign_commitment_tsampler(
    sep_sig_t *sig, uint8_t state[STATE_BYTES],
    const sep_sk_t *sk, const sep_pk_t *pk,
    const poly_q_vec_d cmt) {
  size_t i;
  poly_q_mat_d_k A3;
  poly_q_mat_d_d A;
  poly_q_vec_d u, tmp;
  poly_q_vec_d v1[2];
  uint64_t norm2sq_v1, norm2sq_v2, norm2sq_v3;

  poly_q_mat_d_d_init(A);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(u);
  poly_q_vec_d_init(tmp);
  poly_q_vec_d_init(v1[0]);
  poly_q_vec_d_init(v1[1]);

  /* compute tag from state */
  poly_q_binary_fixed_weight(sig->tag, state);

  /* increment state */
  uint64_t stateinc = ((uint64_t)*(uint32_t*)state) + 1;
  *(uint32_t*)state = (uint32_t)stateinc;
  uint64_t carry = stateinc >> 32;
  for (i = 1; i < STATE_BYTES/4; i++) {
    stateinc = ((uint64_t)*(uint32_t*)&state[i*4]) + carry;
    *(uint32_t*)&state[i*4] = (uint32_t)stateinc;
    carry = stateinc >> 32;
  }

  /* expand uniform matrices */
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);

  int reject_cnt = 0;
reject_signature:
  if (++reject_cnt > 100) {
    fprintf(stderr, "sep_tsampler: >100 rejections, last: v1=%lu/%lu v2=%lu/%lu v3=%lu/%lu\n",
            (unsigned long)norm2sq_v1, (unsigned long)PARAM_B1SQ,
            (unsigned long)norm2sq_v2, (unsigned long)PARAM_B2SQ,
            (unsigned long)norm2sq_v3, (unsigned long)PARAM_B3SQ);
    goto done;
  }
  /* sample v3 */
  poly_q_vec_k_sample_gaussian_s2(sig->v3);

  /* compute u + cmt - A3.v3 */
  poly_q_mat_d_k_mul_vec_k(tmp, A3, sig->v3);
  poly_q_vec_d_sub(tmp, cmt, tmp);
  poly_q_vec_d_add(tmp, u, tmp);

  /* call TSampler-style preimage sampler */
  poly_q_vec_2d_dk_sample_pre_tsampler(
      v1, sig->v2, sk->R, A, pk->B, tmp, sig->tag, sk->S);

  /* norm check */
  norm2sq_v1 = poly_q_vec_d_norm2(v1[0]);
  norm2sq_v1 += poly_q_vec_d_norm2(v1[1]);
  norm2sq_v2 = poly_q_vec_d_norm2(sig->v2[0]);
  for (i = 1; i < PARAM_K; i++) {
    norm2sq_v2 += poly_q_vec_d_norm2(sig->v2[i]);
  }
  norm2sq_v3 = poly_q_vec_k_norm2(sig->v3);
  if ((norm2sq_v1 > PARAM_B1SQ) ||
      (norm2sq_v2 > PARAM_B2SQ) ||
      (norm2sq_v3 > PARAM_B3SQ)) {
    goto reject_signature;
  }

done:
  poly_q_vec_d_set(sig->v12, v1[1]);

  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(v1[0]);
  poly_q_vec_d_clear(v1[1]);
}

/* ================================================
 * Public interface
 * ================================================ */

void sep_sign_tsampler(
    sep_sig_t *sig, uint8_t state[STATE_BYTES],
    const sep_sk_t *sk, const sep_pk_t *pk,
    const uint8_t msg[PARAM_M * PARAM_N / 8]) {
  size_t i;
  poly_q_vec_d cmt;
  poly_q_mat_d_m D;
  poly_q_vec_m m;

  poly_q_vec_d_init(cmt);
  poly_q_mat_d_m_init(D);
  poly_q_vec_m_init(m);

  /* commitment cmt = D.m */
  for (i = 0; i < PARAM_M; i++) {
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  }
  poly_q_mat_d_m_uniform(D, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_m_mul_vec_m(cmt, D, m);

  _sep_sign_commitment_tsampler(sig, state, sk, pk, cmt);

  poly_q_vec_d_clear(cmt);
  poly_q_mat_d_m_clear(D);
  poly_q_vec_m_clear(m);
}
