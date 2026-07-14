#include <stdio.h>

#include "acc_genpp.h"
#include "arith.h"
#include "random.h"

static int acc_genpp_lifecycle_test(void) {
  int ok = 1;
  acc_genpp_pp_t pp;
  acc_genpp_td_t td;
  acc_genpp_state_t st;
  acc_genpp_wit_t sx1, sx2, sx1_sync;
  acc_genpp_update_t msg;
  acc_genpp_update_t missing[ACC_GENPP_BML_CAPACITY];
  size_t missing_len = 0;
  poly_q_vec_d A0, A1, A2;
  const uint32_t x1 = 0x10203040u;
  const uint32_t x2 = 0x55667788u;

  printf("\nacc_genpp_lifecycle_test\n");
  printf("acc size: %.2f KB, compact witness size: %.2f KB, full witness struct payload: %.2f KB\n",
         acc_genpp_acc_size_bytes() / 1024.0,
         acc_genpp_wit_compact_size_bytes() / 1024.0,
         acc_genpp_wit_full_size_bytes() / 1024.0);

  acc_genpp_pp_init(&pp);
  acc_genpp_td_init(&td);
  acc_genpp_state_init(&st);
  acc_genpp_wit_init(&sx1);
  acc_genpp_wit_init(&sx2);
  acc_genpp_wit_init(&sx1_sync);
  poly_q_vec_d_init(A0);
  poly_q_vec_d_init(A1);
  poly_q_vec_d_init(A2);

  acc_genpp_gen(&pp, &td, &st, A0);

  if (!acc_genpp_add(A1, &sx1, &msg, &st, &pp, &td, A0, x1)) {
    printf("Add(x1) failed.\n");
    ok = 0;
    goto cleanup;
  }
  if (!acc_genpp_memver(&pp, A1, x1, &sx1)) {
    printf("MemVer(x1) failed after Add.\n");
    ok = 0;
    goto cleanup;
  }
  if (acc_genpp_add(A1, &sx1, &msg, &st, &pp, &td, A1, x1)) {
    printf("Add(x1) unexpectedly accepted duplicate handle.\n");
    ok = 0;
    goto cleanup;
  }

  if (!acc_genpp_add(A1, &sx2, &msg, &st, &pp, &td, A1, x2)) {
    printf("Add(x2) failed.\n");
    ok = 0;
    goto cleanup;
  }
  if (!acc_genpp_memver(&pp, A1, x2, &sx2)) {
    printf("MemVer(x2) failed after Add.\n");
    ok = 0;
    goto cleanup;
  }

  if (!acc_genpp_delete(A2, &msg, &st, &pp, &td, A1, x1, &sx1)) {
    printf("Delete(x1) failed.\n");
    ok = 0;
    goto cleanup;
  }
  if (msg.seq != 1 || msg.op != ACC_GENPP_OP_DEL || msg.handle != x1) {
    printf("Delete(x1) produced an unexpected upmsg.\n");
    ok = 0;
    goto cleanup;
  }
  if (acc_genpp_memver(&pp, A2, x1, &sx1)) {
    printf("MemVer(x1) unexpectedly accepted a deleted member without update.\n");
    ok = 0;
    goto cleanup;
  }

  if (!acc_genpp_bml_get_missing(&st, 0, missing, ACC_GENPP_BML_CAPACITY, &missing_len) ||
      missing_len != 1 || missing[0].seq != msg.seq || missing[0].handle != msg.handle) {
    printf("BML failed to return the missed Delete(x1) broadcast.\n");
    ok = 0;
    goto cleanup;
  }
  acc_genpp_memwitup(&sx2, &pp, x2, &missing[0]);
  if (!acc_genpp_memver(&pp, A2, x2, &sx2)) {
    printf("MemWitUp/MemVer(x2) failed after Delete(x1).\n");
    ok = 0;
    goto cleanup;
  }
  if (!acc_genpp_bml_get_missing(&st, msg.seq, missing, ACC_GENPP_BML_CAPACITY, &missing_len) ||
      missing_len != 0) {
    printf("BML unexpectedly returned broadcasts for an up-to-date witness.\n");
    ok = 0;
    goto cleanup;
  }

  if (!acc_genpp_memwitsync(&sx1_sync, &pp, &td, &st, A2, x2)) {
    printf("MemWitSync(x2) failed.\n");
    ok = 0;
    goto cleanup;
  }
  if (!acc_genpp_memver(&pp, A2, x2, &sx1_sync)) {
    printf("MemVer synced witness failed.\n");
    ok = 0;
    goto cleanup;
  }

  printf(":");

cleanup:
  poly_q_vec_d_clear(A2);
  poly_q_vec_d_clear(A1);
  poly_q_vec_d_clear(A0);
  acc_genpp_wit_clear(&sx1_sync);
  acc_genpp_wit_clear(&sx2);
  acc_genpp_wit_clear(&sx1);
  acc_genpp_state_clear(&st);
  acc_genpp_td_clear(&td);
  acc_genpp_pp_clear(&pp);
  return ok;
}

int main(void) {
  int ok;

  arith_setup();
  random_init();
  ok = acc_genpp_lifecycle_test();
  arith_teardown();

  if (ok) {
    printf("\npassed.\n");
    return 0;
  }
  printf("\nFAILED!\n");
  return 1;
}
