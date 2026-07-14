#ifndef REVAC_TA_H
#define REVAC_TA_H

#include "../common/revocable_ac.h"

void revac_ta_init(revac_ta_t *ta);
void revac_ta_clear(revac_ta_t *ta);
void revac_ta_setup(revac_ta_t *ta);

int revac_ta_register_uav(revac_ta_t *ta, revac_uav_t *uav);
int revac_ta_register_gcs(revac_ta_t *ta, revac_gcs_t *gcs);

int revac_ta_issue_uav(revac_ta_t *ta,
                       revac_uav_t *uav,
                       const uint8_t msg[PARAM_M * PARAM_N / 8]);

int revac_ta_revoke_uav(revac_ta_t *ta,
                        uint32_t handle,
                        acc_genpp_update_t *upmsg);

#endif
