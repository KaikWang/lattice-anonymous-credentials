#include "revac_ta.h"

#include <string.h>

#include "randombytes.h"

static void revac_copy_user_pk(user_pk_t *dst, const user_pk_t *src) {
  poly_q_vec_d_set(dst->t, src->t);
  memcpy(dst->seed, src->seed, sizeof(dst->seed));
}

static int revac_find_uav_by_handle(const revac_ta_t *ta, uint32_t handle) {
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    if (ta->uavs[i].used && ta->uavs[i].handle == handle) {
      return (int)i;
    }
  }
  return -1;
}

void revac_ta_init(revac_ta_t *ta) {
  sep_keys_init(&ta->ta_pk, &ta->ta_sk);
  acc_genpp_pp_init(&ta->acc_pp);
  acc_genpp_td_init(&ta->acc_td);
  acc_genpp_state_init(&ta->acc_st);
  poly_q_vec_d_init(ta->acc);
  memset(ta->uavs, 0, sizeof(ta->uavs));
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    poly_q_vec_d_init(ta->uavs[i].pk.t);
    memset(ta->uavs[i].pk.seed, 0, sizeof(ta->uavs[i].pk.seed));
  }
}

void revac_ta_clear(revac_ta_t *ta) {
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    poly_q_vec_d_clear(ta->uavs[i].pk.t);
  }
  poly_q_vec_d_clear(ta->acc);
  acc_genpp_state_clear(&ta->acc_st);
  acc_genpp_td_clear(&ta->acc_td);
  acc_genpp_pp_clear(&ta->acc_pp);
  sep_keys_clear(&ta->ta_pk, &ta->ta_sk);
}

void revac_ta_setup(revac_ta_t *ta) {
  sep_keygen(&ta->ta_pk, &ta->ta_sk);
  acc_genpp_gen(&ta->acc_pp, &ta->acc_td, &ta->acc_st, ta->acc);
}

int revac_ta_register_uav(revac_ta_t *ta, revac_uav_t *uav) {
  uint32_t handle = revac_handle_from_identity(&uav->pk, uav->id);

  if (revac_find_uav_by_handle(ta, handle) >= 0) {
    return 0;
  }
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    if (!ta->uavs[i].used) {
      ta->uavs[i].used = 1;
      ta->uavs[i].id = uav->id;
      ta->uavs[i].handle = handle;
      revac_copy_user_pk(&ta->uavs[i].pk, &uav->pk);
      uav->handle = handle;
      uav->registered = 1;
      return 1;
    }
  }
  return 0;
}

int revac_ta_register_gcs(revac_ta_t *ta, revac_gcs_t *gcs) {
  (void)ta;
  gcs->registered = 1;
  gcs->cert_valid = REVAC_GCS_CERT_VALID;
  return 1;
}

int revac_ta_issue_uav(revac_ta_t *ta,
                       revac_uav_t *uav,
                       const uint8_t msg[PARAM_M * PARAM_N / 8]) {
  uint8_t state[STATE_BYTES];
  poly_q_vec_d r[2], cmt, next_acc;
  acc_genpp_update_t upmsg;
  int db_pos;
  int ok = 0;

  if (!uav->registered) {
    return 0;
  }
  db_pos = revac_find_uav_by_handle(ta, uav->handle);
  if (db_pos < 0 || ta->uavs[db_pos].revoked) {
    return 0;
  }

  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(next_acc);

  randombytes(state, sizeof(state));
  osig_user_commit_handle(r, cmt, msg, &uav->pk, uav->handle);
  osig_signer_sign_commitment(&uav->sig, state, &ta->ta_sk, &ta->ta_pk, cmt);
  osig_user_sig_complete(&uav->sig, r);
  if (!osig_user_verify_handle(&uav->sig, &ta->ta_pk, &uav->pk, msg, uav->handle)) {
    goto cleanup;
  }
  if (!acc_genpp_add(next_acc, &uav->witness, &upmsg, &ta->acc_st,
                     &ta->acc_pp, &ta->acc_td, ta->acc, uav->handle)) {
    goto cleanup;
  }

  poly_q_vec_d_set(ta->acc, next_acc);
  uav->last_seq = ta->acc_st.ctr;
  uav->issued = 1;
  ta->uavs[db_pos].issued = 1;
  ok = 1;

cleanup:
  poly_q_vec_d_clear(next_acc);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  return ok;
}

int revac_ta_revoke_uav(revac_ta_t *ta,
                        uint32_t handle,
                        acc_genpp_update_t *upmsg) {
  acc_genpp_wit_t revoke_wit;
  poly_q_vec_d next_acc;
  int db_pos = revac_find_uav_by_handle(ta, handle);
  int ok = 0;

  if (db_pos < 0 || !ta->uavs[db_pos].issued || ta->uavs[db_pos].revoked) {
    return 0;
  }

  acc_genpp_wit_init(&revoke_wit);
  poly_q_vec_d_init(next_acc);

  if (!acc_genpp_memwitsync(&revoke_wit, &ta->acc_pp, &ta->acc_td, &ta->acc_st, ta->acc, handle)) {
    goto cleanup;
  }
  if (!acc_genpp_delete(next_acc, upmsg, &ta->acc_st, &ta->acc_pp, &ta->acc_td, ta->acc, handle, &revoke_wit)) {
    goto cleanup;
  }
  poly_q_vec_d_set(ta->acc, next_acc);
  ta->uavs[db_pos].revoked = 1;
  ok = 1;

cleanup:
  poly_q_vec_d_clear(next_acc);
  acc_genpp_wit_clear(&revoke_wit);
  return ok;
}
