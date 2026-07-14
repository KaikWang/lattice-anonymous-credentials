#include <stdio.h>

#include "arith.h"
#include "random.h"
#include "randombytes.h"
#include "revac_gcs.h"
#include "revac_ta.h"
#include "revac_uav.h"

static int revac_lifecycle_test(void) {
  int ok = 1;
  revac_ta_t ta;
  revac_uav_t uav1, uav2;
  revac_gcs_t gcs;
  acc_genpp_update_t upmsg;
  uint8_t msg1[PARAM_M * PARAM_N / 8], msg2[PARAM_M * PARAM_N / 8];

  printf("\nrevocable_ac_lifecycle_test\n");

  revac_ta_init(&ta);
  revac_uav_init(&uav1);
  revac_uav_init(&uav2);
  revac_gcs_init(&gcs);

  revac_ta_setup(&ta);
  revac_uav_keygen(&uav1, UINT64_C(0x1000000000000001), ta.ta_pk.seed);
  revac_uav_keygen(&uav2, UINT64_C(0x1000000000000002), ta.ta_pk.seed);
  revac_gcs_keygen(&gcs, UINT64_C(0x2000000000000001), ta.ta_pk.seed, UINT64_C(2026071401));
  randombytes(msg1, sizeof(msg1));
  randombytes(msg2, sizeof(msg2));

  if (!revac_ta_register_uav(&ta, &uav1) ||
      !revac_ta_register_uav(&ta, &uav2) ||
      !revac_ta_register_gcs(&ta, &gcs)) {
    printf("registration failed.\n");
    ok = 0;
    goto cleanup;
  }
  printf(":");
  fflush(stdout);

  if (!revac_ta_issue_uav(&ta, &uav1, msg1) ||
      !revac_ta_issue_uav(&ta, &uav2, msg2)) {
    printf("issuance failed.\n");
    ok = 0;
    goto cleanup;
  }
  if (!revac_gcs_verify_uav(&ta, &gcs, &uav1, msg1) ||
      !revac_gcs_verify_uav(&ta, &gcs, &uav2, msg2)) {
    printf("GCS verification failed before revocation.\n");
    ok = 0;
    goto cleanup;
  }
  printf(":");
  fflush(stdout);

  if (!revac_ta_revoke_uav(&ta, uav1.handle, &upmsg)) {
    printf("revocation failed.\n");
    ok = 0;
    goto cleanup;
  }
  if (revac_gcs_verify_uav(&ta, &gcs, &uav1, msg1)) {
    printf("revoked UAV still verified.\n");
    ok = 0;
    goto cleanup;
  }
  if (revac_gcs_verify_uav(&ta, &gcs, &uav2, msg2)) {
    printf("UAV with stale witness unexpectedly verified after revocation.\n");
    ok = 0;
    goto cleanup;
  }
  printf(":");
  fflush(stdout);

  if (!revac_uav_recover_updates(&uav2, &ta)) {
    printf("BML/witness recovery failed.\n");
    ok = 0;
    goto cleanup;
  }
  if (!revac_gcs_verify_uav(&ta, &gcs, &uav2, msg2)) {
    printf("GCS verification failed after BML recovery.\n");
    ok = 0;
    goto cleanup;
  }
  if (uav2.last_seq != upmsg.seq) {
    printf("UAV last_seq was not updated.\n");
    ok = 0;
    goto cleanup;
  }
  printf(":");
  fflush(stdout);

cleanup:
  revac_gcs_clear(&gcs);
  revac_uav_clear(&uav2);
  revac_uav_clear(&uav1);
  revac_ta_clear(&ta);
  return ok;
}

int main(void) {
  int ok;

  arith_setup();
  random_init();
  ok = revac_lifecycle_test();
  arith_teardown();

  if (ok) {
    printf("\npassed.\n");
    return 0;
  }
  printf("\nFAILED!\n");
  return 1;
}
