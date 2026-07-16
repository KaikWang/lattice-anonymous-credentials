#ifndef REVAC_LABRADOR_BRIDGE_H
#define REVAC_LABRADOR_BRIDGE_H

#include "revocable_ac.h"

int revac_labrador_online_enabled(void);

int revac_labrador_memver_prove(revac_show_proof_t *proof,
                                const revac_ta_t *ta,
                                const revac_uav_t *uav);

int revac_labrador_memver_verify(const revac_ta_t *ta,
                                 const revac_show_context_t *ctx,
                                 const revac_show_proof_t *proof);

#endif
