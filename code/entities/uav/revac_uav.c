#include "revac_uav.h"

#include <string.h>

#include "randombytes.h"
#include "revac_labrador_bridge.h"

void revac_uav_init(revac_uav_t *uav) {
  memset(uav, 0, sizeof(*uav));
  user_keys_init(&uav->pk, &uav->sk);
  sep_sig_init(&uav->sig);
  acc_genpp_wit_init(&uav->witness);
}

void revac_uav_clear(revac_uav_t *uav) {
  acc_genpp_wit_clear(&uav->witness);
  sep_sig_clear(&uav->sig);
  user_keys_clear(&uav->pk, &uav->sk);
}

void revac_uav_keygen(revac_uav_t *uav, uint64_t id, const uint8_t seed[SEED_BYTES]) {
  uav->id = id;
  osig_user_keygen(&uav->pk, &uav->sk, seed);
}

void revac_uav_apply_update(revac_uav_t *uav,
                            const revac_ta_t *ta,
                            const acc_genpp_update_t *upmsg) {
  acc_genpp_memwitup(&uav->witness, &ta->acc_pp, uav->handle, upmsg);
  if (upmsg->seq > uav->last_seq) {
    uav->last_seq = upmsg->seq;
  }
}

int revac_uav_recover_updates(revac_uav_t *uav,
                              const revac_ta_t *ta) {
  acc_genpp_update_t missing[ACC_GENPP_BML_CAPACITY];
  size_t missing_len = 0;

  if (acc_genpp_bml_get_missing(&ta->acc_st, uav->last_seq, missing,
                                ACC_GENPP_BML_CAPACITY, &missing_len)) {
    for (size_t i = 0; i < missing_len; i++) {
      revac_uav_apply_update(uav, ta, &missing[i]);
    }
    return 1;
  }

  if (!acc_genpp_memwitsync(&uav->witness, &ta->acc_pp, &ta->acc_td,
                            &ta->acc_st, ta->acc, uav->handle)) {
    return 0;
  }
  uav->last_seq = ta->acc_st.ctr;
  return 1;
}

int revac_uav_show_prove(revac_show_proof_t *proof,
                         const revac_ta_t *ta,
                         const revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8],
                         const uint8_t nonce[REVAC_SHOW_NONCE_BYTES]) {
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K];
  poly_qshow_mat_k_k A3_embed[PARAM_D];
  poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M];
  poly_qshow_mat_k_k Dx_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];
  poly_qshow_vec_m1 s1;
  uint8_t crs_seed[CRS_SEED_BYTES];
  int ok = 0;

  if (!uav->issued) {
    return 0;
  }
  if (!acc_genpp_memver(&ta->acc_pp, ta->acc, uav->handle, &uav->witness)) {
    return 0;
  }

  memcpy(proof->nonce, nonce, REVAC_SHOW_NONCE_BYTES);
  if (!revac_show_proof_set_acc_zk(proof, REVAC_ACC_ZK_NONE, NULL, 0)) {
    return 0;
  }
  revac_show_derive_crs_seed(crs_seed, &ta->ta_pk, uav->pk.seed,
                             proof->nonce, ta->acc_st.ctr, ta->acc);

  poly_qshow_vec_m1_init(s1);
  revac_show_lang_init(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);

  show_user_embed_handle(A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                         Dx_embed, u_embed, s1, &uav->pk, &uav->sk,
                         &ta->ta_pk, &uav->sig, msg, uav->handle);
  show_user_prove_handle(&proof->sig_show, A_embed, B_embed, A3_embed,
                         Ds_embed, D_embed, Dx_embed, s1,
                         crs_seed, uav->pk.seed);
  if (revac_labrador_online_enabled() &&
      !revac_labrador_memver_prove(proof, ta, uav)) {
    goto cleanup;
  }
  ok = 1;

cleanup:
  revac_show_lang_clear(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  poly_qshow_vec_m1_clear(s1);
  return ok;
}
