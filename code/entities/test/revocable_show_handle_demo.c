#include <stdio.h>

#include "arith.h"
#include "sep.h"
#include "osig.h"
#include "random.h"
#include "randombytes.h"
#include "show.h"

static void init_show_lang(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                           poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                           poly_qshow_mat_k_k A3_embed[PARAM_D],
                           poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                           poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                           poly_qshow_mat_k_k Dx_embed[PARAM_D],
                           poly_qshow_vec_k u_embed[PARAM_D]) {
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(u_embed[i]);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_embed[i][j]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_D * PARAM_K; j++) {
      poly_qshow_mat_k_k_init(B_embed[i][j]);
    }
    poly_qshow_mat_k_k_init(A3_embed[i]);
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qshow_mat_k_k_init(D_embed[i][j]);
    }
    poly_qshow_mat_k_k_init(Dx_embed[i]);
  }
}

static void clear_show_lang(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                            poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                            poly_qshow_mat_k_k A3_embed[PARAM_D],
                            poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                            poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                            poly_qshow_mat_k_k Dx_embed[PARAM_D],
                            poly_qshow_vec_k u_embed[PARAM_D]) {
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(u_embed[i]);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_embed[i][j]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_D * PARAM_K; j++) {
      poly_qshow_mat_k_k_clear(B_embed[i][j]);
    }
    poly_qshow_mat_k_k_clear(A3_embed[i]);
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qshow_mat_k_k_clear(D_embed[i][j]);
    }
    poly_qshow_mat_k_k_clear(Dx_embed[i]);
  }
}

int main(void) {
  sep_pk_t pk;
  sep_sk_t sk;
  user_pk_t upk;
  user_sk_t usk;
  sep_sig_t sig;
  show_proof_t proof;
  poly_q_vec_d r[2], cmt;
  poly_qshow_vec_m1 s1;
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K];
  poly_qshow_mat_k_k A3_embed[PARAM_D];
  poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M];
  poly_qshow_mat_k_k Dx_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];
  uint8_t msg[PARAM_M * PARAM_N / 8];
  uint8_t state[STATE_BYTES];
  uint8_t crs_seed[CRS_SEED_BYTES];
  uint32_t handle = 0xa5c3017bu;
  int ok;

  arith_setup();
  random_init();

  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  sep_sig_init(&sig);
  show_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  poly_qshow_vec_m1_init(s1);
  init_show_lang(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);

  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(state, sizeof(state));
  randombytes(crs_seed, sizeof(crs_seed));
  randombytes(msg, sizeof(msg));

  osig_user_commit_handle(r, cmt, msg, &upk, handle);
  osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
  osig_user_sig_complete(&sig, r);

  if (!osig_user_verify_handle(&sig, &pk, &upk, msg, handle)) {
    printf("handle-bound signature did not verify before show.\n");
    ok = 0;
    goto cleanup;
  }

  show_user_embed_handle(A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                         Dx_embed, u_embed, s1, &upk, &usk, &pk, &sig,
                         msg, handle);
  show_user_prove_handle(&proof, A_embed, B_embed, A3_embed, Ds_embed,
                         D_embed, Dx_embed, s1, crs_seed, upk.seed);
  ok = show_verify_handle(&proof, A_embed, B_embed, A3_embed, Ds_embed,
                          D_embed, Dx_embed, u_embed, crs_seed, upk.seed);

cleanup:
  clear_show_lang(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  poly_qshow_vec_m1_clear(s1);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  show_proof_clear(&proof);
  sep_sig_clear(&sig);
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  arith_teardown();

  if (ok) {
    printf("revocable show handle proof passed.\n");
    return 0;
  }
  printf("revocable show handle proof FAILED.\n");
  return 1;
}
