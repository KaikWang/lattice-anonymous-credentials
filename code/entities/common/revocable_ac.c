#include "revocable_ac.h"

#include "fips202.h"
#include "poly_q_sampling.h"
#include <stdlib.h>
#include <string.h>

static void revac_absorb_u64(keccak_state *state, uint64_t x) {
  uint8_t buf[8];

  for (size_t i = 0; i < sizeof(buf); i++) {
    buf[i] = (uint8_t)(x >> (8 * i));
  }
  shake256_absorb(state, buf, sizeof(buf));
}

static void revac_absorb_poly_q(keccak_state *state, const poly_q p) {
  for (size_t i = 0; i < PARAM_N; i++) {
    revac_absorb_u64(state, (uint64_t)poly_q_get_coeff(p, i));
  }
}

uint32_t revac_handle_from_identity(const user_pk_t *pk, uint64_t id) {
  keccak_state state;
  uint8_t out[4];

  shake256_init(&state);
  revac_absorb_u64(&state, id);
  shake256_absorb(&state, pk->seed, sizeof(pk->seed));
  for (size_t i = 0; i < PARAM_D; i++) {
    revac_absorb_poly_q(&state, pk->t->entries[i]);
  }
  shake256_finalize(&state);
  shake256_squeeze(out, sizeof(out), &state);
  return ((uint32_t)out[0]) |
         ((uint32_t)out[1] << 8) |
         ((uint32_t)out[2] << 16) |
         ((uint32_t)out[3] << 24);
}

void revac_show_proof_init(revac_show_proof_t *proof) {
  show_proof_init(&proof->sig_show);
  memset(proof->nonce, 0, sizeof(proof->nonce));
  proof->acc_zk_kind = REVAC_ACC_ZK_NONE;
  proof->acc_zk_proof = NULL;
  proof->acc_zk_proof_len = 0;
}

void revac_show_proof_clear(revac_show_proof_t *proof) {
  free(proof->acc_zk_proof);
  proof->acc_zk_proof = NULL;
  proof->acc_zk_proof_len = 0;
  proof->acc_zk_kind = REVAC_ACC_ZK_NONE;
  show_proof_clear(&proof->sig_show);
}

int revac_show_proof_set_acc_zk(revac_show_proof_t *proof,
                                revac_acc_zk_kind_t kind,
                                const uint8_t *bytes,
                                size_t len) {
  uint8_t *copy = NULL;

  if ((bytes == NULL && len != 0) ||
      (kind == REVAC_ACC_ZK_NONE && len != 0)) {
    return 0;
  }
  if (len != 0) {
    copy = malloc(len);
    if (copy == NULL) {
      return 0;
    }
    memcpy(copy, bytes, len);
  }
  free(proof->acc_zk_proof);
  proof->acc_zk_proof = copy;
  proof->acc_zk_proof_len = len;
  proof->acc_zk_kind = len == 0 ? REVAC_ACC_ZK_NONE : kind;
  return 1;
}

void revac_show_context_init(revac_show_context_t *ctx) {
  memset(ctx->user_seed, 0, sizeof(ctx->user_seed));
  ctx->acc_seq = 0;
  poly_q_vec_d_init(ctx->acc);
}

void revac_show_context_clear(revac_show_context_t *ctx) {
  poly_q_vec_d_clear(ctx->acc);
}

void revac_show_context_from_ta(revac_show_context_t *ctx,
                                const revac_ta_t *ta,
                                const uint8_t user_seed[SEED_BYTES]) {
  memcpy(ctx->user_seed, user_seed, SEED_BYTES);
  ctx->acc_seq = ta->acc_st.ctr;
  poly_q_vec_d_set(ctx->acc, ta->acc);
}

void revac_show_lang_init(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
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

void revac_show_lang_clear(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
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

void revac_show_public_embed(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                             poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                             poly_qshow_mat_k_k A3_embed[PARAM_D],
                             poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                             poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                             poly_qshow_mat_k_k Dx_embed[PARAM_D],
                             poly_qshow_vec_k u_embed[PARAM_D],
                             const sep_pk_t *ta_pk,
                             const uint8_t user_seed[SEED_BYTES]) {
  poly_q_mat_d_d A, Ds[2], Dx;
  poly_q_mat_d_m D;
  poly_q_vec_d a3, u;

  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(Ds[0]);
  poly_q_mat_d_d_init(Ds[1]);
  poly_q_mat_d_d_init(Dx);
  poly_q_mat_d_m_init(D);
  poly_q_vec_d_init(a3);
  poly_q_vec_d_init(u);

  poly_q_mat_d_d_uniform(A, ta_pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_vec_d_uniform(a3, ta_pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_mat_d_d_uniform(Ds[0], ta_pk->seed, DOMAIN_SEPARATOR_DS, 0);
  poly_q_mat_d_d_uniform(Ds[1], ta_pk->seed, DOMAIN_SEPARATOR_DS, PARAM_D);
  poly_q_mat_d_m_uniform(D, ta_pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_d_uniform(Dx, user_seed, DOMAIN_SEPARATOR_DX, 0);
  poly_q_vec_d_uniform(u, ta_pk->seed, DOMAIN_SEPARATOR_U);

  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(u_embed[i], u->entries[i], PARAM_Q1_SHOW);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_Q1_SHOW);
      poly_qshow_subring_embed_mat_k_k(Ds_embed[i][j], Ds[0]->rows[i]->entries[j], PARAM_Q1_SHOW);
      poly_qshow_subring_embed_mat_k_k(Ds_embed[i][j + PARAM_D], Ds[1]->rows[i]->entries[j], PARAM_Q1_SHOW);
    }
    for (size_t j = 0; j < PARAM_D * PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(B_embed[i][j], ta_pk->B[j / PARAM_D]->rows[i]->entries[j % PARAM_D], PARAM_Q1_SHOW);
    }
    poly_qshow_subring_embed_mat_k_k(A3_embed[i], a3->entries[i], PARAM_Q1_SHOW);
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qshow_subring_embed_mat_k_k(D_embed[i][j], D->rows[i]->entries[j], PARAM_Q1_SHOW);
    }
    poly_qshow_subring_embed_mat_k_k(Dx_embed[i], Dx->rows[i]->entries[0], PARAM_Q1_SHOW);
  }

  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(a3);
  poly_q_mat_d_m_clear(D);
  poly_q_mat_d_d_clear(Dx);
  poly_q_mat_d_d_clear(Ds[1]);
  poly_q_mat_d_d_clear(Ds[0]);
  poly_q_mat_d_d_clear(A);
}

void revac_show_derive_crs_seed(uint8_t crs_seed[CRS_SEED_BYTES],
                                const sep_pk_t *ta_pk,
                                const uint8_t user_seed[SEED_BYTES],
                                const uint8_t nonce[REVAC_SHOW_NONCE_BYTES],
                                uint32_t acc_seq,
                                const poly_q_vec_d acc) {
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, ta_pk->seed, SEED_BYTES);
  shake256_absorb(&state, user_seed, SEED_BYTES);
  shake256_absorb(&state, nonce, REVAC_SHOW_NONCE_BYTES);
  shake256_absorb(&state, (const uint8_t *)&acc_seq, sizeof(acc_seq));
  for (size_t i = 0; i < PARAM_D; i++) {
    revac_absorb_poly_q(&state, acc->entries[i]);
  }
  shake256_finalize(&state);
  shake256_squeeze(crs_seed, CRS_SEED_BYTES, &state);
}
