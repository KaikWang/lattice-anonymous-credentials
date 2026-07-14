#include "revac_gcs.h"

#include <string.h>

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
