#ifndef REVAC_UAV_H
#define REVAC_UAV_H

#include "../common/revocable_ac.h"

void revac_uav_init(revac_uav_t *uav);
void revac_uav_clear(revac_uav_t *uav);
void revac_uav_keygen(revac_uav_t *uav, uint64_t id, const uint8_t seed[SEED_BYTES]);

void revac_uav_apply_update(revac_uav_t *uav,
                            const revac_ta_t *ta,
                            const acc_genpp_update_t *upmsg);

int revac_uav_recover_updates(revac_uav_t *uav,
                              const revac_ta_t *ta);

int revac_uav_show_prove(revac_show_proof_t *proof,
                         const revac_ta_t *ta,
                         const revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8],
                         const uint8_t nonce[REVAC_SHOW_NONCE_BYTES]);

#endif
