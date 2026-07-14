#ifndef REVOCABLE_AC_H
#define REVOCABLE_AC_H

#include <stdint.h>
#include <stddef.h>

#include "acc_genpp.h"
#include "osig.h"

#define REVAC_MAX_UAVS ACC_GENPP_MAX_DICT_ENTRIES
#define REVAC_GCS_CERT_VALID 1

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

uint32_t revac_handle_from_identity(const user_pk_t *pk, uint64_t id);

#endif
