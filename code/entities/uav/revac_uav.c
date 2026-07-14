#include "revac_uav.h"

#include <string.h>

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
