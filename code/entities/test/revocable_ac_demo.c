#include <stdio.h>
#include <string.h>

#include "arith.h"
#include "random.h"
#include "randombytes.h"
#include "revac_gcs.h"
#include "revac_labrador_bridge.h"
#include "revac_ta.h"
#include "revac_uav.h"

static int expect_true(int cond, const char *msg) {
  if (!cond) {
    printf("%s\n", msg);
    return 0;
  }
  return 1;
}

static int expect_false(int cond, const char *msg) {
  return expect_true(!cond, msg);
}

static int revac_plan_process_demo(void) {
  int ok = 1;
  revac_ta_t ta;
  revac_uav_t uav1, uav2, uav3;
  revac_gcs_t gcs;
  revac_show_proof_t show1;
  revac_show_context_t show1_ctx;
  acc_genpp_update_t revoke_uav1, revoke_uav2;
  uint8_t nonce1[REVAC_SHOW_NONCE_BYTES];
  uint8_t msg1[PARAM_M * PARAM_N / 8];
  uint8_t msg2[PARAM_M * PARAM_N / 8];
  uint8_t msg3[PARAM_M * PARAM_N / 8];

  printf("\nrevocable_ac_plan_process_demo\n");

  revac_ta_init(&ta);
  revac_uav_init(&uav1);
  revac_uav_init(&uav2);
  revac_uav_init(&uav3);
  revac_gcs_init(&gcs);
  revac_show_proof_init(&show1);
  revac_show_context_init(&show1_ctx);

  printf("[phase 1] TA setup, anonymous credential keygen, AccGen, BML init\n");
  revac_ta_setup(&ta);
  if (!expect_true(ta.acc_st.ctr == 0, "TA accumulator counter was not initialized.") ||
      !expect_true(ta.acc_st.bml_len == 0, "TA BML was not initialized empty.") ||
      !expect_true(gcs.cert_valid == 0, "GCS cert should be invalid before registration.")) {
    ok = 0;
    goto cleanup;
  }

  printf("[phase 2] UAV/GCS keygen and registration\n");
  revac_uav_keygen(&uav1, UINT64_C(0x1000000000000001), ta.ta_pk.seed);
  revac_uav_keygen(&uav2, UINT64_C(0x1000000000000002), ta.ta_pk.seed);
  revac_uav_keygen(&uav3, UINT64_C(0x1000000000000003), ta.ta_pk.seed);
  revac_gcs_keygen(&gcs, UINT64_C(0x2000000000000001),
                   ta.ta_pk.seed, UINT64_C(2026071401));
  randombytes(msg1, sizeof(msg1));
  randombytes(msg2, sizeof(msg2));
  randombytes(msg3, sizeof(msg3));
  randombytes(nonce1, sizeof(nonce1));

  if (!expect_true(revac_ta_register_uav(&ta, &uav1), "UAV1 registration failed.") ||
      !expect_true(revac_ta_register_uav(&ta, &uav2), "UAV2 registration failed.") ||
      !expect_true(revac_ta_register_uav(&ta, &uav3), "UAV3 registration failed.") ||
      !expect_true(revac_ta_register_gcs(&ta, &gcs), "GCS registration failed.") ||
      !expect_true(uav1.handle == revac_handle_from_identity(&uav1.pk, uav1.id),
                   "UAV1 handle is not H(PK_UAV, ID_UAV).") ||
      !expect_true(uav2.handle == revac_handle_from_identity(&uav2.pk, uav2.id),
                   "UAV2 handle is not H(PK_UAV, ID_UAV).") ||
      !expect_true(uav3.handle == revac_handle_from_identity(&uav3.pk, uav3.id),
                   "UAV3 handle is not H(PK_UAV, ID_UAV).") ||
      !expect_true(gcs.registered && gcs.cert_valid == REVAC_GCS_CERT_VALID,
                   "GCS certificate state was not installed by TA.")) {
    ok = 0;
    goto cleanup;
  }
  printf(": registration and handle derivation passed\n");

  printf("[phase 3] anonymous credential issuance and accumulator Add\n");
  if (!expect_true(revac_ta_issue_uav(&ta, &uav1, msg1), "UAV1 issuance failed.") ||
      !expect_true(revac_ta_issue_uav(&ta, &uav2, msg2), "UAV2 issuance failed.") ||
      !expect_true(revac_ta_issue_uav(&ta, &uav3, msg3), "UAV3 issuance failed.") ||
      !expect_true(revac_gcs_verify_uav(&ta, &gcs, &uav1, msg1),
                   "GCS verification failed for UAV1 before revocation.") ||
      !expect_true(revac_gcs_verify_uav(&ta, &gcs, &uav2, msg2),
                   "GCS verification failed for UAV2 before revocation.") ||
      !expect_true(revac_gcs_verify_uav(&ta, &gcs, &uav3, msg3),
                   "GCS verification failed for UAV3 before revocation.")) {
    ok = 0;
    goto cleanup;
  }
  gcs.cert_valid = 0;
  if (!expect_false(revac_gcs_verify_uav(&ta, &gcs, &uav1, msg1),
                    "GCS accepted a UAV credential with an invalid GCS certificate.")) {
    ok = 0;
    goto cleanup;
  }
  gcs.cert_valid = REVAC_GCS_CERT_VALID;
  printf(": issuance, signature verification, MemVer, and GCS cert gate passed\n");

  printf("[phase 3b] anonymous show proof API without clear handle/witness at GCS\n");
  revac_show_context_from_ta(&show1_ctx, &ta, uav1.pk.seed);
  if (!expect_true(revac_uav_show_prove(&show1, &ta, &uav1, msg1, nonce1),
                   "UAV1 failed to produce handle-bound show proof.") ||
      !expect_true(revac_gcs_verify_show_signature(&ta, &gcs, &show1_ctx, &show1),
                   "GCS failed to verify handle-bound signature show proof.")) {
    ok = 0;
    goto cleanup;
  }
  if (revac_labrador_online_enabled()) {
    if (!expect_true(revac_gcs_verify_show(&ta, &gcs, &show1_ctx, &show1),
                     "GCS failed to verify full Labrador-backed revocable show proof.")) {
      ok = 0;
      goto cleanup;
    }
    printf(": signature show proof and Labrador MemVer proof passed\n");
  } else {
    if (!expect_false(revac_gcs_verify_show(&ta, &gcs, &show1_ctx, &show1),
                      "GCS accepted a full revocable show proof without AccGenpp ZK proof.")) {
      ok = 0;
      goto cleanup;
    }
    printf(": signature show proof passed; full revocation-ZK gate correctly requires AccGenpp proof\n");
  }

  printf("[phase 4] Delete, broadcast update, MemWitUp, BML recovery, MemWitSync fallback\n");
  if (!expect_true(revac_ta_revoke_uav(&ta, uav1.handle, &revoke_uav1),
                   "UAV1 revocation failed.") ||
      !expect_true(revoke_uav1.seq == 1 && revoke_uav1.handle == uav1.handle,
                   "UAV1 revocation broadcast did not carry the expected seq/handle.") ||
      !expect_true(ta.acc_st.bml_len == 1, "BML did not store the latest revocation broadcast.")) {
    ok = 0;
    goto cleanup;
  }

  revac_uav_apply_update(&uav2, &ta, &revoke_uav1);
  if (!expect_false(revac_gcs_verify_uav(&ta, &gcs, &uav1, msg1),
                    "Revoked UAV1 still verified.") ||
      !expect_true(revac_gcs_verify_uav(&ta, &gcs, &uav2, msg2),
                   "Online UAV2 failed after direct MemWitUp.") ||
      !expect_false(revac_gcs_verify_uav(&ta, &gcs, &uav3, msg3),
                    "Offline UAV3 with stale witness unexpectedly verified.")) {
    ok = 0;
    goto cleanup;
  }

  if (!expect_true(revac_uav_recover_updates(&uav3, &ta),
                   "BML recovery failed for offline UAV3.") ||
      !expect_true(revac_gcs_verify_uav(&ta, &gcs, &uav3, msg3),
                   "GCS verification failed after UAV3 BML recovery.") ||
      !expect_true(uav3.last_seq == revoke_uav1.seq,
                   "UAV3 last_seq was not updated after BML recovery.")) {
    ok = 0;
    goto cleanup;
  }
  printf(": direct broadcast update and BML recovery passed\n");

  if (!expect_true(revac_ta_revoke_uav(&ta, uav2.handle, &revoke_uav2),
                   "UAV2 revocation failed.") ||
      !expect_false(revac_gcs_verify_uav(&ta, &gcs, &uav2, msg2),
                    "Revoked UAV2 still verified.") ||
      !expect_false(revac_gcs_verify_uav(&ta, &gcs, &uav3, msg3),
                    "UAV3 stale witness unexpectedly verified after second revocation.")) {
    ok = 0;
    goto cleanup;
  }

  ta.acc_st.bml_len = 0;
  ta.acc_st.bml_start = 0;
  memset(ta.acc_st.bml, 0, sizeof(ta.acc_st.bml));
  if (!expect_true(revac_uav_recover_updates(&uav3, &ta),
                   "MemWitSync fallback failed when BML was unavailable.") ||
      !expect_true(revac_gcs_verify_uav(&ta, &gcs, &uav3, msg3),
                   "GCS verification failed after MemWitSync fallback.") ||
      !expect_true(uav3.last_seq == revoke_uav2.seq,
                   "UAV3 last_seq was not updated after MemWitSync fallback.")) {
    ok = 0;
    goto cleanup;
  }

  if (!expect_false(revac_ta_register_uav(&ta, &uav1),
                    "TA accepted a duplicate/revoked UAV handle.")) {
    ok = 0;
    goto cleanup;
  }
  printf(": revocation, stale rejection, and MemWitSync fallback passed\n");

cleanup:
  revac_show_proof_clear(&show1);
  revac_show_context_clear(&show1_ctx);
  revac_gcs_clear(&gcs);
  revac_uav_clear(&uav3);
  revac_uav_clear(&uav2);
  revac_uav_clear(&uav1);
  revac_ta_clear(&ta);
  return ok;
}

int main(void) {
  int ok;

  arith_setup();
  random_init();
  ok = revac_plan_process_demo();
  arith_teardown();

  if (ok) {
    printf("\npassed.\n");
    return 0;
  }
  printf("\nFAILED!\n");
  return 1;
}
