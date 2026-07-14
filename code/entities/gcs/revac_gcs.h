#ifndef REVAC_GCS_H
#define REVAC_GCS_H

#include "../common/revocable_ac.h"

void revac_gcs_init(revac_gcs_t *gcs);
void revac_gcs_clear(revac_gcs_t *gcs);
void revac_gcs_keygen(revac_gcs_t *gcs, uint64_t id, const uint8_t seed[SEED_BYTES], uint64_t timestamp);

int revac_gcs_verify_uav(const revac_ta_t *ta,
                         const revac_gcs_t *gcs,
                         const revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8]);

#endif
