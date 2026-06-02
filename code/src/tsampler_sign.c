/*************************************************
* TSampler-based signature scheme for ablation.
*
* Implements keygen, sign, verify using the
* tag-friendly dual-gadget sampler.
*
* Reuses the AC project's compute_covariance and
* sk_sq_spectral_norm by zero-padding R to K=5.
**************************************************/

#include "tsampler_sign.h"
#include "tsampler_sampling.h"
#include "tsampler_params.h"
#include "tsampler_helpers.h"
#include "arith.h"
#include "randombytes.h"
#include "poly_q_sampling.h"
#include "fips202.h"
#include "covariance.h"

/* ================================================
 * Init / clear helpers
 * ================================================ */

void tsampler_keys_init(tsampler_pk_t *pk, tsampler_sk_t *sk) {
  for (size_t i = 0; i < PARAM_KH_TS; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(sk->R[j][i]);
    }
    poly_q_mat_d_d_init(pk->B[i]);
  }
  poly_real_mat_2d_2d_init(sk->S);
}

void tsampler_keys_clear(tsampler_pk_t *pk, tsampler_sk_t *sk) {
  for (size_t i = 0; i < PARAM_KH_TS; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(sk->R[j][i]);
    }
    poly_q_mat_d_d_clear(pk->B[i]);
  }
  poly_real_mat_2d_2d_clear(sk->S);
}

void tsampler_sig_init(tsampler_sig_t *sig) {
  poly_q_init(sig->tag);
  for (size_t i = 0; i < PARAM_KH_TS; i++) {
    poly_q_vec_d_init(sig->v2[i]);
  }
  poly_q_vec_d_init(sig->v12);
  poly_q_init(sig->v3);
}

void tsampler_sig_clear(tsampler_sig_t *sig) {
  poly_q_clear(sig->tag);
  for (size_t i = 0; i < PARAM_KH_TS; i++) {
    poly_q_vec_d_clear(sig->v2[i]);
  }
  poly_q_vec_d_clear(sig->v12);
  poly_q_clear(sig->v3);
}

/* ================================================
 * Key generation
 * ================================================ */

/*************************************************
* Name:        tsampler_keygen
*
* Description: Generate TSampler keypair.
*   Zero-pads R to K=5 dimensions internally
*   so that the existing sk_sq_spectral_norm
*   and compute_covariance functions can be
*   reused without modification.
*
*   Only the first KH_TS=3 columns of B are
*   stored in the public key.
**************************************************/
void tsampler_keygen(tsampler_pk_t *pk, tsampler_sk_t *sk) {
  uint8_t root_seed[SEED_BYTES], seeds[SEED_BYTES*2];
  uint8_t *public_seed = seeds, *secret_seed = &seeds[SEED_BYTES];
  poly_q_mat_d_d A;
  poly_q_mat_d_d RRstar[2][2];
  uint32_t kappa;
  size_t i, j;

  /* We create R_padded[2][PARAM_K] (= R_padded[2][5]) so that
   * sk_sq_spectral_norm(RRstar, R_padded) type-checks.
   * Only the first PARAM_KH_TS=3 columns are non-zero. */
  poly_q_mat_d_d R_padded[2][PARAM_K];
  poly_q_mat_d_d B_padded[PARAM_K];

  /* generate random seeds */
  randombytes(root_seed, SEED_BYTES);
  sha3_512(seeds, root_seed, SEED_BYTES);
#if SEED_BYTES != 32
#error "SEED_BYTES must be 32."
#endif

  /* init matrices */
  poly_q_mat_d_d_init(A);
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(RRstar[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_q_mat_d_d_init(R_padded[i][j]);
    }
  }
  for (j = 0; j < PARAM_K; j++) {
    poly_q_mat_d_d_init(B_padded[j]);
  }

  /* expand uniform A' from seed */
  poly_q_mat_d_d_uniform(A, public_seed, TS_DOMSEP_A, 0);

  /* sample R from centered binomial B_1,
   * reject if spectral norm too large */
  kappa = 0;
  {
    int iter = 0;
    const int max_iter = 100;
    double norm_val;
    do {
      for (i = 0; i < PARAM_KH_TS; i++) {
        for (j = 0; j < 2; j++) {
          poly_q_mat_d_d_binomial(R_padded[j][i], secret_seed,
              kappa++, TS_DOMSEP_R);
        }
      }
      /* zero out the padding columns (KH_TS..K-1) */
      for (i = PARAM_KH_TS; i < PARAM_K; i++) {
        for (j = 0; j < 2; j++) {
          poly_q_mat_d_d_zero(R_padded[j][i]);
        }
      }
      norm_val = sk_sq_spectral_norm(RRstar, R_padded);
      if (++iter >= max_iter) {
        fprintf(stderr, "TSampler keygen: spectral norm "
                "%.2f exceeded bound %.2f after %d iters\n",
                norm_val, PARAM_R_MAX_SQ_SPECTRAL_NORM_TS, max_iter);
        break;
      }
    } while (norm_val > PARAM_R_MAX_SQ_SPECTRAL_NORM_TS);
  }

  /* copy the non-padded R columns to the secret key */
  for (i = 0; i < PARAM_KH_TS; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_set(sk->R[j][i], R_padded[j][i]);
    }
  }

  /* compute B = (I | A') * R (full padded version) */
  for (i = 0; i < PARAM_K; i++) {
    poly_q_mat_d_d_mul_mat_d_d(B_padded[i], A, R_padded[1][i]);
    poly_q_mat_d_d_add(B_padded[i], B_padded[i], R_padded[0][i]);
  }

  /* store only the first KH_TS columns in the public key */
  for (i = 0; i < PARAM_KH_TS; i++) {
    poly_q_mat_d_d_set(pk->B[i], B_padded[i]);
  }

  /* copy public seed */
  for (i = 0; i < SEED_BYTES; i++) {
    pk->seed[i] = public_seed[i];
  }

  /* compute Schur complement S for perturbation sampling
   * (reuses existing compute_covariance — RRstar was computed
   *  by sk_sq_spectral_norm above) */
  compute_covariance(sk->S, RRstar);

  /* clean up */
  poly_q_mat_d_d_clear(A);
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(RRstar[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_q_mat_d_d_clear(R_padded[i][j]);
    }
  }
  for (j = 0; j < PARAM_K; j++) {
    poly_q_mat_d_d_clear(B_padded[j]);
  }
}

/* ================================================
 * Tag generation (same as SEP / original TSampler)
 * ================================================ */

static void _tsampler_tag_gen(poly_q tag, uint8_t state[STATE_BYTES_TS]) {
  /* compute tag from state (fixed-weight binary polynomial) */
  poly_q_binary_fixed_weight(tag, state);

  /* increment state */
#if (STATE_BYTES_TS % 4) != 0
#error "STATE_BYTES_TS must be multiple of 4"
#endif
  uint64_t stateinc = ((uint64_t)*(uint32_t*)state) + 1;
  *(uint32_t*)state = (uint32_t)stateinc;
  uint64_t carry = stateinc >> 32;
  for (size_t i = 1; i < STATE_BYTES_TS/4; i++) {
    stateinc = ((uint64_t)*(uint32_t*)&state[i*4]) + carry;
    *(uint32_t*)&state[i*4] = (uint32_t)stateinc;
    carry = stateinc >> 32;
  }
}

/* ================================================
 * Sign
 * ================================================ */

/*************************************************
* Name:        tsampler_sign
*
* Description: Sign a message with TSampler.
*
*   Follows the original TSampler signing flow:
*   1. Expand A', a3, u, d from public seed
*   2. Compute u = u + d * m  (message embedding)
*   3. Generate tag from state
*   4. Compute tag^{-1} mod bH
*   5. Sample scalar v3 from Gaussian D_{s4}
*   6. Rejection loop:
*      a. Sample v3
*      b. u' = u + d*m - a3*v3
*      c. Call tsampler_sample_pre
*      d. Norm check → goto reject if too large
*   7. Output sig = (tag, v12, v2[0..KH_TS-1], v3)
**************************************************/
void tsampler_sign(
    tsampler_sig_t      *sig,
    uint8_t              state[STATE_BYTES_TS],
    const tsampler_sk_t *sk,
    const tsampler_pk_t *pk,
    const uint8_t        msg[PARAM_N/8]) {
  size_t i;
  poly_q_mat_d_d A;
  poly_q_vec_d u, a3, d, tmp;
  poly_q_vec_d v11;
  poly_q m;
  uint64_t norm2sq_v11, norm2sq_v12, norm2sq_v2, norm2sq_v3;

  /* init */
  poly_q_mat_d_d_init(A);
  poly_q_init(m);
  poly_q_vec_d_init(u);
  poly_q_vec_d_init(a3);
  poly_q_vec_d_init(d);
  poly_q_vec_d_init(tmp);
  poly_q_vec_d_init(v11);

  /* expand uniform matrices from public seed */
  poly_q_mat_d_d_uniform(A,  pk->seed, TS_DOMSEP_A,  0);
  poly_q_vec_d_uniform(   u,  pk->seed, TS_DOMSEP_U);
  poly_q_vec_d_uniform(   a3, pk->seed, TS_DOMSEP_A3);
  poly_q_vec_d_uniform(   d,  pk->seed, TS_DOMSEP_D);

  /* message embedding: u = u + d*m.
   * Use tmp to avoid in-place nmod_poly_mul (FLINT restriction). */
  poly_q_from_bits(m, msg);
  poly_q_vec_d_mul_poly(tmp, d, m);
  poly_q_vec_d_add(u, u, tmp);
  /* tmp now holds d_original * m.  Save for later: d = tmp */
  poly_q_vec_d_set(d, tmp);

  /* generate tag + taginv.  Retry if tag is not invertible mod bH. */
  {
    int tag_tries = 0;
    int reject_loop_count;
  try_new_tag:
    if (++tag_tries > 100) {
      fprintf(stderr, "TSampler sign: failed to find invertible tag\n");
      goto sign_cleanup;
    }
    _tsampler_tag_gen(sig->tag, state);
    tsampler_poly_q_invert_mod_bh(m, sig->tag);

    /* check if tag is invertible: result poly should be non-zero */
    {
      int is_zero = 1;
      for (size_t ci = 0; ci < PARAM_N; ci++) {
        if (poly_q_get_coeff(m, ci) != 0) { is_zero = 0; break; }
      }
      if (is_zero) goto try_new_tag;
    }

    reject_loop_count = 0;
  reject_signature:
    if (++reject_loop_count > 100) {
      goto try_new_tag;  /* try a new tag after many rejections */
    }

    /* sample scalar v3 from discrete Gaussian */
    tsampler_poly_q_sample_gaussian_s4(sig->v3);

    /* compute u + d*m - a3*v3  (d reused as temp) */
    poly_q_vec_d_mul_poly(d, a3, sig->v3);
    poly_q_vec_d_sub(d, u, d);

    /* call TSampler preimage sampler */
    tsampler_sample_pre(v11, sig->v12, sig->v2,
        sk->R, A, pk->B, d, sig->tag, m, sk->S);
    /* note: m holds taginv at this point */

    /* norm check: reject if any vector too long */
    norm2sq_v11 = poly_q_vec_d_norm2(v11);
    norm2sq_v12 = poly_q_vec_d_norm2(sig->v12);
    norm2sq_v2  = poly_q_vec_d_norm2(sig->v2[0]);
    for (i = 1; i < PARAM_KH_TS; i++) {
      norm2sq_v2 += poly_q_vec_d_norm2(sig->v2[i]);
    }
    norm2sq_v3 = poly_q_sq_norm2(sig->v3);

    if ((norm2sq_v11 > PARAM_B11SQ_TS) ||
        (norm2sq_v12 > PARAM_B12SQ_TS) ||
        (norm2sq_v2  > PARAM_B2SQ_TS)  ||
        (norm2sq_v3  > PARAM_B3SQ_TS)) {
      goto reject_signature;
    }

  }

sign_cleanup:

  /* clean up */
  poly_q_mat_d_d_clear(A);
  poly_q_clear(m);
  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(a3);
  poly_q_vec_d_clear(d);
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(v11);
}

/* ================================================
 * Verify
 * ================================================ */

/*************************************************
* Name:        tsampler_verify
*
* Description: Verify a TSampler signature.
*
*   Reconstructs v11 from the verification equation:
*     v11 = u + d*m - A'*v12 - a3*v3
*           + Σ(B_i - tag * qL * bH^i * I) * v2_i
*
*   Checks:
*     ||v11||² ≤ B11SQ
*     ||v12||² ≤ B12SQ
*     ||v2||²  ≤ B2SQ
*     ||v3||²  ≤ B3SQ
*     tag has Hamming weight W
**************************************************/
int tsampler_verify(
    const tsampler_sig_t *sig,
    const uint8_t         msg[PARAM_N/8],
    const tsampler_pk_t  *pk) {
  size_t i, j;
  int64_t bexpi;
  poly_q m;
  poly_q_vec_d u, a3, v11;
  poly_q_mat_d_d A;
  poly_q_mat_d_d Btmp;
  poly_q tag_times_bexpi;
  uint64_t norm2sq_v11, norm2sq_v12, norm2sq_v2, norm2sq_v3;
  int64_t tag_weight;

  /* init */
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(Btmp);
  poly_q_vec_d_init(u);
  poly_q_vec_d_init(a3);
  poly_q_vec_d_init(v11);
  poly_q_init(m);
  poly_q_init(tag_times_bexpi);

  /* expand uniform matrices from public seed */
  poly_q_mat_d_d_uniform(A,  pk->seed, TS_DOMSEP_A,  0);
  poly_q_vec_d_uniform(   u,  pk->seed, TS_DOMSEP_U);
  poly_q_vec_d_uniform(   a3, pk->seed, TS_DOMSEP_A3);
  poly_q_vec_d_uniform(   v11, pk->seed, TS_DOMSEP_D); /* v11 as temp for d */

  /* v11 = u_orig + d_orig*m.  Avoid in-place poly mul. */
  {
    poly_q_vec_d d_times_m;
    poly_q_vec_d_init(d_times_m);
    poly_q_from_bits(m, msg);
    poly_q_vec_d_mul_poly(d_times_m, v11, m);  /* d_times_m = d_orig * m */
    poly_q_vec_d_add(v11, u, d_times_m);        /* v11 = u_orig + d_orig*m */
    poly_q_vec_d_clear(d_times_m);
  }

  /* subtract A'*v12 and a3*v3 */
  {
    poly_q_vec_d sub;
    poly_q_vec_d_init(sub);
    poly_q_vec_d_mul_poly(sub, a3, sig->v3);
    poly_q_mat_d_d_muladd_vec_d(sub, A, sig->v12);
    poly_q_vec_d_sub(v11, v11, sub);
    poly_q_vec_d_clear(sub);
  }

  /* add (B_i - tag * qL * bH^i * I) * v2_i */
  bexpi = PARAM_QL_TS;
  for (i = 0; i < PARAM_KH_TS; i++) {
    poly_q_mul_scalar(tag_times_bexpi, sig->tag, bexpi);
    poly_q_mat_d_d_set(Btmp, pk->B[i]);
    for (j = 0; j < PARAM_D; j++) {
      poly_q_sub(Btmp->rows[j]->entries[j],
          Btmp->rows[j]->entries[j], tag_times_bexpi);
    }
    poly_q_mat_d_d_muladd_vec_d(v11, Btmp, sig->v2[i]);
    bexpi *= PARAM_BH_TS;
  }

  /* compute norms */
  norm2sq_v11 = poly_q_vec_d_norm2(v11);
  norm2sq_v12 = poly_q_vec_d_norm2(sig->v12);
  norm2sq_v2  = poly_q_vec_d_norm2(sig->v2[0]);
  for (i = 1; i < PARAM_KH_TS; i++) {
    norm2sq_v2 += poly_q_vec_d_norm2(sig->v2[i]);
  }
  norm2sq_v3  = poly_q_sq_norm2(sig->v3);
  tag_weight  = poly_q_weight(sig->tag);

  /* clean up */
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(Btmp);
  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(a3);
  poly_q_vec_d_clear(v11);
  poly_q_clear(m);
  poly_q_clear(tag_times_bexpi);

  /* verify all checks pass */
  {
    int ok = (norm2sq_v11 <= PARAM_B11SQ_TS) &&
             (norm2sq_v12 <= PARAM_B12SQ_TS) &&
             (norm2sq_v2  <= PARAM_B2SQ_TS)  &&
             (norm2sq_v3  <= PARAM_B3SQ_TS)  &&
             (tag_weight  == PARAM_W_TS);
    if (!ok) {
      static int dbg = 0;
      if (dbg++ < 3) {
        fprintf(stderr, "  [VERIFY] norms: v11=%lu/%lu v12=%lu/%lu "
                "v2=%lu/%lu v3=%lu/%lu tag_w=%ld/%d\n",
                (unsigned long)norm2sq_v11, (unsigned long)PARAM_B11SQ_TS,
                (unsigned long)norm2sq_v12, (unsigned long)PARAM_B12SQ_TS,
                (unsigned long)norm2sq_v2,  (unsigned long)PARAM_B2SQ_TS,
                (unsigned long)norm2sq_v3,  (unsigned long)PARAM_B3SQ_TS,
                (long)tag_weight, PARAM_W_TS);
      }
    }
    return ok;
  }
}
