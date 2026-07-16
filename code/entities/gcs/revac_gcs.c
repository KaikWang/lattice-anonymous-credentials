#include "revac_gcs.h"

#include <string.h>

#include "revac_labrador_bridge.h"

void revac_gcs_init(revac_gcs_t *gcs) {
  memset(gcs, 0, sizeof(*gcs));
  user_keys_init(&gcs->pk, &gcs->sk);
}

void revac_gcs_clear(revac_gcs_t *gcs) {
  user_keys_clear(&gcs->pk, &gcs->sk);
}

void revac_gcs_keygen(revac_gcs_t *gcs, uint64_t id, const uint8_t seed[SEED_BYTES], uint64_t timestamp) {
  gcs->id = id;
  gcs->timestamp = timestamp;
  osig_user_keygen(&gcs->pk, &gcs->sk, seed);
}

int revac_gcs_verify_uav(const revac_ta_t *ta,
                         const revac_gcs_t *gcs,
                         const revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8]) {
  if (!gcs->registered || gcs->cert_valid != REVAC_GCS_CERT_VALID || !uav->issued) {
    return 0;
  }
  if (!osig_user_verify_handle(&uav->sig, &ta->ta_pk, &uav->pk, msg, uav->handle)) {
    return 0;
  }
  return acc_genpp_memver(&ta->acc_pp, ta->acc, uav->handle, &uav->witness);
}

int revac_gcs_verify_show_signature(const revac_ta_t *ta,
                                    const revac_gcs_t *gcs,
                                    const revac_show_context_t *ctx,
                                    const revac_show_proof_t *proof) {
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K];
  poly_qshow_mat_k_k A3_embed[PARAM_D];
  poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M];
  poly_qshow_mat_k_k Dx_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];
  uint8_t expected_crs[CRS_SEED_BYTES];
  int ok;

  if (!gcs->registered || gcs->cert_valid != REVAC_GCS_CERT_VALID) {
    return 0;
  }
  if (ctx->acc_seq != ta->acc_st.ctr || !poly_q_vec_d_equal(ctx->acc, ta->acc)) {
    return 0;
  }
  revac_show_derive_crs_seed(expected_crs, &ta->ta_pk, ctx->user_seed,
                             proof->nonce, ctx->acc_seq, ctx->acc);

  revac_show_lang_init(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  revac_show_public_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                          Dx_embed, u_embed, &ta->ta_pk, ctx->user_seed);
  ok = show_verify_handle(&proof->sig_show, A_embed, B_embed, A3_embed,
                          Ds_embed, D_embed, Dx_embed, u_embed,
                          expected_crs, ctx->user_seed);
  revac_show_lang_clear(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  return ok;
}

int revac_gcs_verify_show(const revac_ta_t *ta,
                          const revac_gcs_t *gcs,
                          const revac_show_context_t *ctx,
                          const revac_show_proof_t *proof) {
  if (!revac_gcs_verify_show_signature(ta, gcs, ctx, proof)) {
    return 0;
  }
  if (proof->acc_zk_kind != REVAC_ACC_ZK_LABRADOR_MEMVER ||
      proof->acc_zk_proof == NULL ||
      proof->acc_zk_proof_len == 0) {
    return 0;
  }
  return revac_labrador_memver_verify(ta, ctx, proof);
}
