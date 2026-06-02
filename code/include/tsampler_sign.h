#ifndef TSAMPLER_SIGN_H
#define TSAMPLER_SIGN_H

#include <stdint.h>
#include "params.h"           /* for PARAM_N, PARAM_D, SEED_BYTES */
#include "tsampler_params.h"  /* for PARAM_KH_TS etc. */
#include "arith.h"

/*************************************************
* TSampler-based signature scheme for ablation.
*
* This is a SEPARATE signature scheme from the
* original SEP. It uses the tag-friendly dual-gadget
* structure with base-decomposition sampling.
*
* Comparison with SEP:
*   SEP:     single gadget G=[1|b|...|b^4], K=5, Klein
*   TSampler: dual gadget [G_L | qL*G_{H,t}],
*             KL=2, KH=3, base decomposition
**************************************************/

/* ================================================
 * Key and signature types
 * ================================================ */

typedef struct {
  poly_q_mat_d_d    R[2][PARAM_KH_TS];  /* secret trapdoor */
  poly_real_mat_2d_2d S;                /* Schur complement */
} tsampler_sk_t;

typedef struct {
  poly_q_mat_d_d B[PARAM_KH_TS];        /* public matrices */
  uint8_t        seed[SEED_BYTES];      /* seed for A', a3, u, d */
} tsampler_pk_t;

typedef struct {
  poly_q       tag;                     /* tag polynomial */
  poly_q_vec_d v12;                     /* middle preimage */
  poly_q_vec_d v2[PARAM_KH_TS];         /* bottom preimage */
  poly_q       v3;                      /* scalar Gaussian (not vector!) */
} tsampler_sig_t;

/* ================================================
 * Init / clear
 * ================================================ */

void tsampler_keys_init(tsampler_pk_t *pk, tsampler_sk_t *sk);
void tsampler_keys_clear(tsampler_pk_t *pk, tsampler_sk_t *sk);
void tsampler_sig_init(tsampler_sig_t *sig);
void tsampler_sig_clear(tsampler_sig_t *sig);

/* ================================================
 * Core operations
 * ================================================ */

/*************************************************
* Name:        tsampler_keygen
*
* Description: Generate TSampler public/private key.
*   - Samples R from centered binomial distribution
*   - Checks spectral norm (rejection)
*   - Computes B = (I | A') * R
*   - Computes Schur complement S for sampling
**************************************************/
void tsampler_keygen(tsampler_pk_t *pk, tsampler_sk_t *sk);

/*************************************************
* Name:        tsampler_sign
*
* Description: Sign a message using TSampler.
*   - Expands public matrices from seed
*   - Generates tag from state (fixed-weight binary)
*   - Samples v3 from discrete Gaussian (scalar)
*   - Calls tsampler_sample_pre for short preimage
*   - Norm check + rejection sampling
*
*   Signature: sig = (tag, v12, v2[0..KH_TS-1], v3)
**************************************************/
void tsampler_sign(
    tsampler_sig_t *sig,
    uint8_t         state[STATE_BYTES_TS],
    const tsampler_sk_t *sk,
    const tsampler_pk_t *pk,
    const uint8_t   msg[PARAM_N/8]);

/*************************************************
* Name:        tsampler_verify
*
* Description: Verify a TSampler signature.
*   - Expands public matrices from seed
*   - Reconstructs v11 from the verification equation
*   - Checks norm bounds and tag weight
*
*   Returns 1 if valid, 0 otherwise.
**************************************************/
int tsampler_verify(
    const tsampler_sig_t *sig,
    const uint8_t         msg[PARAM_N/8],
    const tsampler_pk_t  *pk);

#endif /* TSAMPLER_SIGN_H */
