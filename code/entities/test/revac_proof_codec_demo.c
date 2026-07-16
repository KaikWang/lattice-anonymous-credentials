#include <stdio.h>
#include <stdlib.h>

#include "arith.h"
#include "random.h"
#include "randombytes.h"
#include "revac_gcs.h"
#include "revac_proof_codec.h"
#include "revac_ta.h"
#include "revac_uav.h"

static int expect_true(int cond, const char *msg) {
  if (!cond) {
    fprintf(stderr, "FAIL: %s\n", msg);
    return 0;
  }
  return 1;
}

static void print_size(const char *label, size_t bytes) {
  printf("%s: %zu bytes (%.2f KiB)\n", label, bytes, bytes / 1024.0);
}

int main(void) {
  revac_ta_t ta;
  revac_uav_t uav;
  revac_gcs_t gcs;
  revac_show_proof_t proof;
  revac_show_proof_t unpacked;
  revac_show_context_t ctx;
  uint8_t msg[PARAM_M * PARAM_N / 8];
  uint8_t nonce[REVAC_SHOW_NONCE_BYTES];
  uint8_t *wire = NULL;
  size_t wire_len = 0;
  int online;
  int ok = 1;

  arith_setup();
  random_init();

  revac_ta_init(&ta);
  revac_uav_init(&uav);
  revac_gcs_init(&gcs);
  revac_show_proof_init(&proof);
  revac_show_proof_init(&unpacked);
  revac_show_context_init(&ctx);

  revac_ta_setup(&ta);
  revac_uav_keygen(&uav, UINT64_C(0x1000000000000123), ta.ta_pk.seed);
  revac_gcs_keygen(&gcs, UINT64_C(0x2000000000000123),
                   ta.ta_pk.seed, UINT64_C(123456));
  randombytes(msg, sizeof(msg));
  randombytes(nonce, sizeof(nonce));

  ok = ok && expect_true(revac_ta_register_uav(&ta, &uav),
                         "UAV registration failed");
  ok = ok && expect_true(revac_ta_register_gcs(&ta, &gcs),
                         "GCS registration failed");
  ok = ok && expect_true(revac_ta_issue_uav(&ta, &uav, msg),
                         "UAV issuance failed");
  if (!ok) {
    goto cleanup;
  }

  revac_show_context_from_ta(&ctx, &ta, uav.pk.seed);
  ok = expect_true(revac_uav_show_prove(&proof, &ta, &uav, msg, nonce),
                   "show proof generation failed");
  if (!ok) {
    goto cleanup;
  }

  print_size("signature show proof raw packed estimate",
             revac_show_signature_raw_size_bytes());
  print_size("signature show proof compact bound",
             revac_show_signature_compact_bound_bytes());
  print_size("in-memory/wire-equivalent proof before codec",
             revac_show_proof_wire_size_bytes(&proof));
  print_size("accumulator zk proof payload before codec",
             proof.acc_zk_proof_len);

  ok = expect_true(revac_show_proof_pack_compact(&wire, &wire_len, &proof, 1),
                   "compact proof packing failed");
  if (!ok) {
    goto cleanup;
  }
  print_size("compact revocable show proof wire bytes", wire_len);

  ok = expect_true(revac_show_proof_unpack_compact(&unpacked, wire, wire_len),
                   "compact proof unpacking failed");
  if (!ok) {
    goto cleanup;
  }
  print_size("accumulator zk proof payload after codec",
             unpacked.acc_zk_proof_len);

  ok = expect_true(revac_gcs_verify_show_signature(&ta, &gcs, &ctx, &unpacked),
                   "unpacked signature show proof did not verify");
  if (!ok) {
    goto cleanup;
  }

  online = unpacked.acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER ||
           unpacked.acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE;
  if (online) {
    ok = expect_true(revac_gcs_verify_show(&ta, &gcs, &ctx, &unpacked),
                     "unpacked revocable show proof did not verify");
    if (!ok) {
      goto cleanup;
    }
  }

  printf("codec demo passed%s.\n",
         online ? " with accumulator proof verification" :
                  " with signature proof verification");

cleanup:
  free(wire);
  revac_show_context_clear(&ctx);
  revac_show_proof_clear(&unpacked);
  revac_show_proof_clear(&proof);
  revac_gcs_clear(&gcs);
  revac_uav_clear(&uav);
  revac_ta_clear(&ta);
  arith_teardown();
  return ok ? 0 : 1;
}
