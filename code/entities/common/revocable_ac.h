#ifndef REVOCABLE_AC_H
#define REVOCABLE_AC_H

#include <stdint.h>
#include <stddef.h>

#include "acc_genpp.h"
#include "osig.h"
#include "show.h"

#define REVAC_MAX_UAVS ACC_GENPP_MAX_DICT_ENTRIES
#define REVAC_GCS_CERT_VALID 1
#define REVAC_SHOW_NONCE_BYTES 32

typedef enum {
  REVAC_ACC_ZK_NONE = 0,
  REVAC_ACC_ZK_LABRADOR_MEMVER = 1,
  REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE = 2
} revac_acc_zk_kind_t;

typedef struct {
  uint64_t id;
  uint32_t handle;
  user_pk_t pk;
  int used;
  int issued;
  int revoked;
} revac_uav_record_t;

typedef struct {
  sep_pk_t ta_pk;
  sep_sk_t ta_sk;
  acc_genpp_pp_t acc_pp;
  acc_genpp_td_t acc_td;
  acc_genpp_state_t acc_st;
  poly_q_vec_d acc;
  revac_uav_record_t uavs[REVAC_MAX_UAVS];
} revac_ta_t;

typedef struct {
  uint64_t id;
  user_pk_t pk;
  user_sk_t sk;
  uint32_t handle;
  sep_sig_t sig;
  acc_genpp_wit_t witness;
  uint32_t last_seq;
  int registered;
  int issued;
} revac_uav_t;

typedef struct {
  uint64_t id;
  user_pk_t pk;
  user_sk_t sk;
  uint64_t timestamp;
  int registered;
  int cert_valid;
} revac_gcs_t;

typedef struct {
  show_proof_t sig_show;
  uint8_t nonce[REVAC_SHOW_NONCE_BYTES];
  revac_acc_zk_kind_t acc_zk_kind;
  uint8_t *acc_zk_proof;
  size_t acc_zk_proof_len;
} revac_show_proof_t;

typedef struct {
  uint8_t user_seed[SEED_BYTES];
  uint32_t acc_seq;
  poly_q_vec_d acc;
} revac_show_context_t;

uint32_t revac_handle_from_identity(const user_pk_t *pk, uint64_t id);
void revac_show_proof_init(revac_show_proof_t *proof);
void revac_show_proof_clear(revac_show_proof_t *proof);
int revac_show_proof_set_acc_zk(revac_show_proof_t *proof,
                                revac_acc_zk_kind_t kind,
                                const uint8_t *bytes,
                                size_t len);
void revac_show_context_init(revac_show_context_t *ctx);
void revac_show_context_clear(revac_show_context_t *ctx);
void revac_show_context_from_ta(revac_show_context_t *ctx,
                                const revac_ta_t *ta,
                                const uint8_t user_seed[SEED_BYTES]);

void revac_show_lang_init(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                          poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                          poly_qshow_mat_k_k A3_embed[PARAM_D],
                          poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                          poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                          poly_qshow_mat_k_k Dx_embed[PARAM_D],
                          poly_qshow_vec_k u_embed[PARAM_D]);

void revac_show_lang_clear(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                           poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                           poly_qshow_mat_k_k A3_embed[PARAM_D],
                           poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                           poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                           poly_qshow_mat_k_k Dx_embed[PARAM_D],
                           poly_qshow_vec_k u_embed[PARAM_D]);

void revac_show_public_embed(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                             poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                             poly_qshow_mat_k_k A3_embed[PARAM_D],
                             poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                             poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                             poly_qshow_mat_k_k Dx_embed[PARAM_D],
                             poly_qshow_vec_k u_embed[PARAM_D],
                             const sep_pk_t *ta_pk,
                             const uint8_t user_seed[SEED_BYTES]);

void revac_show_derive_crs_seed(uint8_t crs_seed[CRS_SEED_BYTES],
                                const sep_pk_t *ta_pk,
                                const uint8_t user_seed[SEED_BYTES],
                                const uint8_t nonce[REVAC_SHOW_NONCE_BYTES],
                                uint32_t acc_seq,
                                const poly_q_vec_d acc);

#endif
