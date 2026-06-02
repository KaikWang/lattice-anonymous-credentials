/*************************************************
* Ablation test: SEP (Klein) vs SEP (TSampler)
*
* Both use identical keygen, verify, and parameters.
* The only difference is the gadget sampling method:
*   SEP:        Klein/ZF lattice sampler
*   SEP-TSamp:  Base-decomposition coset sampler
**************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "params.h"
#include "arith.h"
#include "sep.h"
#include "sep_tsampler.h"
#include "randombytes.h"
#include "random.h"

static double get_time_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

#define N_TRIALS 3

int main(void) {
  double t_sep_k = 0, t_sep_s = 0, t_sep_v = 0;
  double t_ts_k  = 0, t_ts_s  = 0, t_ts_v  = 0;
  int sep_ok = 1, ts_ok = 1;

  printf("=== Ablation: SEP(Klein) vs SEP(TSampler) ===\n");
  printf("N=%d D=%d Q=%d K=%d B=%d\n\n",
         PARAM_N, PARAM_D, (int)PARAM_Q, PARAM_K, PARAM_B);

  arith_setup();
  random_init();

  /* ── Original SEP (Klein) ── */
  printf("--- SEP (Klein sampler) ---\n");
  for (int t = 0; t < N_TRIALS; t++) {
    sep_pk_t pk; sep_sk_t sk; sep_sig_t sig;
    uint8_t msg[PARAM_M * PARAM_N / 8];
    uint8_t st[STATE_BYTES] = {0};
    sep_keys_init(&pk, &sk); sep_sig_init(&sig);
    randombytes(msg, sizeof(msg));

    double t0 = get_time_ms();
    sep_keygen(&pk, &sk);
    double t1 = get_time_ms();
    sep_sign(&sig, st, &sk, &pk, msg);
    double t2 = get_time_ms();
    int ok = sep_verify(&sig, msg, &pk);
    double t3 = get_time_ms();
    if (!ok) sep_ok = 0;

    t_sep_k += (t1-t0); t_sep_s += (t2-t1); t_sep_v += (t3-t2);
    printf("  trial %d: keygen=%.0f sign=%.0f verify=%.2fms %s\n",
           t+1, t1-t0, t2-t1, t3-t2, ok?"OK":"FAIL");
    fflush(stdout);
    sep_keys_clear(&pk, &sk); sep_sig_clear(&sig);
  }
  printf("  Avg: keygen=%.0f sign=%.0f verify=%.2fms [%s]\n\n",
         t_sep_k/N_TRIALS, t_sep_s/N_TRIALS, t_sep_v/N_TRIALS,
         sep_ok?"PASS":"FAIL");

  /* ── SEP with TSampler gadget ── */
  printf("--- SEP (TSampler base-decomp) ---\n");
  for (int t = 0; t < N_TRIALS; t++) {
    sep_pk_t pk; sep_sk_t sk; sep_sig_t sig;
    uint8_t msg[PARAM_M * PARAM_N / 8];
    uint8_t st[STATE_BYTES] = {0};
    sep_keys_init(&pk, &sk); sep_sig_init(&sig);
    randombytes(msg, sizeof(msg));

    double t0 = get_time_ms();
    sep_keygen(&pk, &sk);
    double t1 = get_time_ms();
    sep_sign_tsampler(&sig, st, &sk, &pk, msg);
    double t2 = get_time_ms();
    int ok = sep_verify(&sig, msg, &pk);
    double t3 = get_time_ms();
    if (!ok) ts_ok = 0;

    t_ts_k += (t1-t0); t_ts_s += (t2-t1); t_ts_v += (t3-t2);
    printf("  trial %d: keygen=%.0f sign=%.0f verify=%.2fms %s\n",
           t+1, t1-t0, t2-t1, t3-t2, ok?"OK":"FAIL");
    fflush(stdout);
    sep_keys_clear(&pk, &sk); sep_sig_clear(&sig);
  }
  printf("  Avg: keygen=%.0f sign=%.0f verify=%.2fms [%s]\n\n",
         t_ts_k/N_TRIALS, t_ts_s/N_TRIALS, t_ts_v/N_TRIALS,
         ts_ok?"PASS":"FAIL");

  printf("=== Summary: SEP=%s  TSampler=%s ===\n",
         sep_ok?"PASS":"FAIL", ts_ok?"PASS":"FAIL");

  arith_teardown();
  return (sep_ok && ts_ok) ? 0 : 1;
}
