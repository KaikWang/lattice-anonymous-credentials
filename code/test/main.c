#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "sep.h"
#include "osig.h"
#include "show.h"
#include "randombytes.h"
#include "random.h"

/* ── Toggle: uncomment to print proof/witness sizes ── */
//#define BENCHMARK_PROOF_SIZES// 注释掉后，不打印零知识证明大小
/* ── Toggle: uncomment to dump proof/witness as hex ── */
//#define DUMP_PROOF_WITNESS_HEX

#define NTESTS 1
#define NSUBTESTS 1

static int sep_test(void)
{
  int rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  randombytes(state, STATE_BYTES);

  printf("\nsep_test\n");

  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      sep_sign(&sig, state, &sk, &pk, msg);
      if (!sep_verify(&sig, msg, &pk))
      {
        printf("sep_verify returned zero for a valid signature.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      msg[0] ^= 1;
      if (sep_verify(&sig, msg, &pk))
      {
        printf("sep_verify returned non-zero for a valid signature on the wrong message.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

  sep_test_cleanup:
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return rval;
}

static int osig_signing_test(void)
{
  int i,j,rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  user_sk_t usk;
  user_pk_t upk;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  randombytes(state, STATE_BYTES);

  printf("\nosig_signing_test\n");

  // init
  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);
  user_keys_init(&upk, &usk);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);

  for (i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    osig_user_keygen(&upk, &usk, pk.seed);
    for (j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      osig_user_commit(r, cmt, msg, &upk);
      osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
      osig_user_sig_complete(&sig, r);
      if (!osig_user_verify(&sig, &pk, &upk, msg))
      {
        printf("osig_user_verify returned zero for a valid signature.\n");
        rval = 0;
        goto osig_signing_test_cleanup;
      }
      msg[0] ^= 1;
      if (osig_user_verify(&sig, &pk, &upk, msg))
      {
        printf("osig_user_verify returned non-zero for a valid signature on the wrong message.\n");
        rval = 0;
        goto osig_signing_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

osig_signing_test_cleanup:
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return rval;
}

static int osig_proof_test(void)
{
  int i,j,rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  osig_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  coeff_qiss coeff;
  poly_qiss_vec_k u[2*PARAM_D], s1[PARAM_M1_K_ISS];
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  printf("\nosig_proof_test\n");

  // init
  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  osig_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  for (i = 0; i < 2*PARAM_D; i++)
  {
    poly_qiss_vec_k_init(u[i]);
  }
  for (i = 0; i < PARAM_M1_K_ISS; i++)
  {
    poly_qiss_vec_k_init(s1[i]);
  }
  for (i = 0; i < PARAM_D; i++)
  {
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qiss_mat_k_k_init(D_embed[i][j]);
    }
  }

  for (i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    osig_user_keygen(&upk, &usk, pk.seed);
    randombytes(crs_seed, CRS_SEED_BYTES);
    for (j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      osig_user_commit(r, cmt, msg, &upk);
      osig_user_embed(A_embed, Ds_embed, D_embed, u, s1, &upk, &usk, cmt, r, msg);
      osig_user_prove(&proof, A_embed, Ds_embed, D_embed, u, s1, crs_seed, upk.seed);
      if (!osig_signer_verify(&proof, A_embed, Ds_embed, D_embed, u, crs_seed, upk.seed))
      {
        printf("osig_signer_verify returned zero for a valid proof.\n");
        rval = 0;
        goto osig_proof_test_cleanup;
      }
      coeff = poly_qiss_get_coeff(u[0]->entries[0], 0) + 1;
      poly_qiss_set_coeff(u[0]->entries[0], 0, coeff);
      if (osig_signer_verify(&proof, A_embed, Ds_embed, D_embed, u, crs_seed, upk.seed))
      {
        printf("osig_signer_verify returned non-zero for a valid proof but tampered commitment.\n");
        rval = 0;
        goto osig_proof_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

osig_proof_test_cleanup:
  for (i = 0; i < 2*PARAM_D; i++)
  {
    poly_qiss_vec_k_clear(u[i]);
  }
  for (i = 0; i < PARAM_M1_K_ISS; i++)
  {
    poly_qiss_vec_k_clear(s1[i]);
  }
  for (i = 0; i < PARAM_D; i++)
  {
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qiss_mat_k_k_clear(D_embed[i][j]);
    }
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  osig_proof_clear(&proof);
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  return rval;
}

static int show_proof_test(void)
{
  int i,j,rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  sep_sig_t sig;
  show_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  coeff_qshow coeff;
  poly_qshow_vec_m1 s1;
  poly_qshow_vec_k u_embed[PARAM_D];
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_D*PARAM_K]; poly_qshow_vec_k A3_embed[PARAM_D];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  printf("\nshow_proof_test\n");

  // init
  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  sep_sig_init(&sig);
  show_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  poly_qshow_vec_m1_init(s1);
  for (i = 0; i < PARAM_D; i++)
  {
    poly_qshow_vec_k_init(u_embed[i]);
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qshow_mat_k_k_init(A_embed[i][j]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_D*PARAM_K; j++)
    {
      poly_qshow_mat_k_k_init(B_embed[i][j]);
    }
    poly_qshow_vec_k_init(A3_embed[i]);
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qshow_mat_k_k_init(D_embed[i][j]);
    }
  }

  for (i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    osig_user_keygen(&upk, &usk, pk.seed);
    randombytes(crs_seed, CRS_SEED_BYTES);
    for (j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      osig_user_commit(r, cmt, msg, &upk);
      osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
      osig_user_sig_complete(&sig, r);
      if (!osig_user_verify(&sig, &pk, &upk, msg))
      {
        printf("osig_user_verify returned zero for a valid signature.\n");
        rval = 0;
        goto show_proof_test_cleanup;
      }
      show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1, &upk, &usk, &pk, &sig, msg);
      show_user_prove(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, s1, crs_seed, upk.seed);
      if (!show_verify(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, crs_seed, upk.seed))
      {
        printf("show_verify returned zero for a valid proof.\n");
        rval = 0;
        goto show_proof_test_cleanup;
      }
      coeff = poly_qshow_get_coeff(u_embed[0]->entries[0], 0) + 1;
      poly_qshow_set_coeff(u_embed[0]->entries[0], 0, coeff);
      if (show_verify(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, crs_seed, upk.seed))
      {
        printf("show_verify returned non-zero for a valid proof but wrong statement.\n");
        rval = 0;
        goto show_proof_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

show_proof_test_cleanup:
  sep_keys_clear(&pk, &sk);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  show_proof_clear(&proof);
  poly_q_vec_d_clear(r[0]);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(cmt);
  poly_qshow_vec_m1_clear(s1);
  for (i = 0; i < PARAM_D; i++)
  {
    poly_qshow_vec_k_clear(u_embed[i]);
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qshow_mat_k_k_clear(A_embed[i][j]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_D*PARAM_K; j++)
    {
      poly_qshow_mat_k_k_clear(B_embed[i][j]);
    }
    poly_qshow_vec_k_clear(A3_embed[i]);
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qshow_mat_k_k_clear(D_embed[i][j]);
    }
  }
  return rval;
}

#ifdef DUMP_PROOF_WITNESS_HEX
/* ── Hex dump helper (writes to FILE*) ── */

static void hex_dump(FILE *fp, const char *label, const uint8_t *buf, size_t len)
{
  fprintf(fp, "\n── %s (%zu bytes) ──\n", label, len);
  for (size_t off = 0; off < len; off += 16) {
    fprintf(fp, "%08zx  ", off);
    for (size_t i = 0; i < 16; i++) {
      if (off + i < len)
        fprintf(fp, "%02x ", buf[off + i]);
      else
        fprintf(fp, "   ");
      if (i == 7) fprintf(fp, " ");
    }
    fprintf(fp, " |");
    for (size_t i = 0; i < 16 && off + i < len; i++) {
      unsigned char c = buf[off + i];
      fputc(isprint(c) ? c : '.', fp);
    }
    fprintf(fp, "|\n");
  }
}

#define DUMP_OUTPUT_FILE "/home/wkk/pq_ac/lattice-anonymous-credentials/code/test/proof_witness_hex.txt"

/* ── Dump osig (issuance) proof and witness as hex ── */

static void dump_osig_proof_and_witness(void)
{
  int i, j;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  osig_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  poly_qiss_vec_k u[2 * PARAM_D], s1[PARAM_M1_K_ISS];
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2 * PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M * PARAM_N / 8], crs_seed[CRS_SEED_BYTES];
  uint8_t *packed_buf = NULL;
  size_t buf_max = 0, pos;

  FILE *fp = fopen(DUMP_OUTPUT_FILE, "w");
  if (!fp) return;

  printf("\nDumping osig proof & witness to %s ...\n", DUMP_OUTPUT_FILE);

  fprintf(fp, "========================================\n");
  fprintf(fp, "=== OSIG (Issuance) Proof & Witness ===\n");
  fprintf(fp, "========================================\n");

  /* compute max buffer size needed */
  {
    size_t s = POLYQISS_VECD_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = POLYQISS_VEC256L_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = PARAM_ARP_ISS * COEFFQISS_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = POLYQISS_VECL_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = POLYQISS_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = PARAM_M1_K_ISS * POLYQISS_VECK_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = PARAM_M2_ISS * POLYQISS_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    /* witness s1 */
    s = PARAM_M1_K_ISS * POLYQISS_VECK_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
  }
  packed_buf = (uint8_t *)malloc(buf_max + 1);
  if (!packed_buf) { fclose(fp); return; }

  /* init */
  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  osig_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  for (i = 0; i < 2 * PARAM_D; i++) poly_qiss_vec_k_init(u[i]);
  for (i = 0; i < PARAM_M1_K_ISS; i++) poly_qiss_vec_k_init(s1[i]);
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + 0]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_M; j++) poly_qiss_mat_k_k_init(D_embed[i][j]);
  }

  /* run one issuance */
  sep_keygen(&pk, &sk);
  randombytes(state, STATE_BYTES);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(crs_seed, CRS_SEED_BYTES);
  randombytes(msg, PARAM_M * PARAM_N / 8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_user_embed(A_embed, Ds_embed, D_embed, u, s1, &upk, &usk, cmt, r, msg);
  osig_user_prove(&proof, A_embed, Ds_embed, D_embed, u, s1, crs_seed, upk.seed);

  fprintf(fp, "\n--- osig_proof_t fields ---\n");

  /* tA */
  poly_qiss_vec_d_pack(packed_buf, proof.tA);
  hex_dump(fp, "proof.tA (poly_qiss_vec_d)", packed_buf, POLYQISS_VECD_PACKEDBYTES);

  /* tB */
  poly_qiss_vec_256_l_pack(packed_buf, proof.tB);
  hex_dump(fp, "proof.tB (poly_qiss_vec_256_l)", packed_buf, POLYQISS_VEC256L_PACKEDBYTES);

  /* z3[PARAM_ARP_ISS] */
  pos = 0;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    coeff_qiss_pack(&packed_buf[pos], proof.z3[i]);
    pos += COEFFQISS_PACKEDBYTES;
  }
  hex_dump(fp, "proof.z3 (coeff_qiss[256])", packed_buf, pos);

  /* h */
  poly_qiss_vec_l_pack(packed_buf, proof.h);
  hex_dump(fp, "proof.h (poly_qiss_vec_l)", packed_buf, POLYQISS_VECL_PACKEDBYTES);

  /* t1 */
  poly_qiss_pack(packed_buf, proof.t1);
  hex_dump(fp, "proof.t1 (poly_qiss)", packed_buf, POLYQISS_PACKEDBYTES);

  /* c */
  poly_qiss_pack(packed_buf, proof.c);
  hex_dump(fp, "proof.c (poly_qiss)", packed_buf, POLYQISS_PACKEDBYTES);

  /* ctr_c */
  hex_dump(fp, "proof.ctr_c (uint32_t)", (const uint8_t *)&proof.ctr_c, sizeof(uint32_t));

  /* z1[PARAM_M1_K_ISS] */
  pos = 0;
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_pack(&packed_buf[pos], proof.z1[i]);
    pos += POLYQISS_VECK_PACKEDBYTES;
  }
  hex_dump(fp, "proof.z1 (poly_qiss_vec_k[26])", packed_buf, pos);

  /* z2 (poly_qiss_vec_m2 via get_poly + pack) */
  {
    poly_qiss tmp;
    poly_qiss_init(tmp);
    pos = 0;
    for (i = 0; i < PARAM_M2_ISS; i++) {
      poly_qiss_vec_m2_get_poly(tmp, proof.z2, i);
      poly_qiss_pack(&packed_buf[pos], tmp);
      pos += POLYQISS_PACKEDBYTES;
    }
    hex_dump(fp, "proof.z2 (poly_qiss_vec_m2)", packed_buf, pos);
    poly_qiss_clear(tmp);
  }

  fprintf(fp, "\n--- osig witness s1[PARAM_M1_K_ISS] ---\n");
  pos = 0;
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_pack(&packed_buf[pos], s1[i]);
    pos += POLYQISS_VECK_PACKEDBYTES;
  }
  hex_dump(fp, "witness.s1 (poly_qiss_vec_k[26])", packed_buf, pos);

  /* cleanup */
  sep_keys_clear(&pk, &sk);
  user_keys_clear(&upk, &usk);
  osig_proof_clear(&proof);
  poly_q_vec_d_clear(r[0]);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(cmt);
  for (i = 0; i < 2 * PARAM_D; i++) poly_qiss_vec_k_clear(u[i]);
  for (i = 0; i < PARAM_M1_K_ISS; i++) poly_qiss_vec_k_clear(s1[i]);
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + 0]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_M; j++) poly_qiss_mat_k_k_clear(D_embed[i][j]);
  }
  free(packed_buf);
  fclose(fp);
}

/* ── Dump show proof and witness as hex ── */

static void dump_show_proof_and_witness(void)
{
  int i, j;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  sep_sig_t sig;
  show_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  poly_qshow_vec_m1 s1;
  poly_qshow_vec_k u_embed[PARAM_D];
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_D*PARAM_K]; poly_qshow_vec_k A3_embed[PARAM_D];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2 * PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M * PARAM_N / 8], crs_seed[CRS_SEED_BYTES];
  uint8_t *packed_buf = NULL;
  size_t buf_max = 0, pos;

  FILE *fp = fopen(DUMP_OUTPUT_FILE, "a");
  if (!fp) return;

  printf("Dumping show proof & witness to %s ...\n", DUMP_OUTPUT_FILE);

  fprintf(fp, "\n========================================\n");
  fprintf(fp, "=== SHOW Proof & Witness ===\n");
  fprintf(fp, "========================================\n");

  /* compute max buffer size needed */
  {
    size_t s = POLYQSHOW_VECD_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = POLYQSHOW_VEC256L_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = PARAM_ARP_SHOW * COEFFQSHOW_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = POLYQSHOW_VECL_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = POLYQSHOW_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = PARAM_M1_SHOW * POLYQSHOW_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
    s = PARAM_M2_SHOW * POLYQSHOW_PACKEDBYTES;
    if (s > buf_max) buf_max = s;
  }
  packed_buf = (uint8_t *)malloc(buf_max + 1);
  if (!packed_buf) { fclose(fp); return; }

  /* init */
  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  sep_sig_init(&sig);
  show_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  poly_qshow_vec_m1_init(s1);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(u_embed[i]);
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_embed[i][j]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + 0]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_D * PARAM_K; j++)
      poly_qshow_mat_k_k_init(B_embed[i][j]);
    for (j = 0; j < PARAM_K; j++)
      poly_qshow_vec_k_init(A3_embed[i]); /* was mat_k_k_init */
    for (j = 0; j < PARAM_M; j++)
      poly_qshow_mat_k_k_init(D_embed[i][j]);
  }

  /* run one full issuance + show */
  sep_keygen(&pk, &sk);
  randombytes(state, STATE_BYTES);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(crs_seed, CRS_SEED_BYTES);
  randombytes(msg, PARAM_M * PARAM_N / 8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
  osig_user_sig_complete(&sig, r);
  show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1,
                  &upk, &usk, &pk, &sig, msg);
  show_user_prove(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                  s1, crs_seed, upk.seed);

  fprintf(fp, "\n--- show_proof_t fields ---\n");

  /* tA */
  poly_qshow_vec_d_pack(packed_buf, proof.tA);
  hex_dump(fp, "proof.tA (poly_qshow_vec_d)", packed_buf, POLYQSHOW_VECD_PACKEDBYTES);

  /* tB */
  poly_qshow_vec_256_l_pack(packed_buf, proof.tB);
  hex_dump(fp, "proof.tB (poly_qshow_vec_256_l)", packed_buf, POLYQSHOW_VEC256L_PACKEDBYTES);

  /* z3[PARAM_ARP_SHOW] */
  pos = 0;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow_pack(&packed_buf[pos], proof.z3[i]);
    pos += COEFFQSHOW_PACKEDBYTES;
  }
  hex_dump(fp, "proof.z3 (coeff_qshow[256])", packed_buf, pos);

  /* h */
  poly_qshow_vec_l_pack(packed_buf, proof.h);
  hex_dump(fp, "proof.h (poly_qshow_vec_l)", packed_buf, POLYQSHOW_VECL_PACKEDBYTES);

  /* t1 */
  poly_qshow_pack(packed_buf, proof.t1);
  hex_dump(fp, "proof.t1 (poly_qshow)", packed_buf, POLYQSHOW_PACKEDBYTES);

  /* c */
  poly_qshow_pack(packed_buf, proof.c);
  hex_dump(fp, "proof.c (poly_qshow)", packed_buf, POLYQSHOW_PACKEDBYTES);

  /* ctr_c */
  hex_dump(fp, "proof.ctr_c (uint32_t)", (const uint8_t *)&proof.ctr_c, sizeof(uint32_t));

  /* z1 (poly_qshow_vec_m1 = PARAM_M1_SHOW poly_qshow entries) */
  {
    poly_qshow tmp;
    poly_qshow_init(tmp);
    pos = 0;
    for (i = 0; i < PARAM_M1_SHOW; i++) {
      poly_qshow_vec_m1_get_poly(tmp, proof.z1, i);
      poly_qshow_pack(&packed_buf[pos], tmp);
      pos += POLYQSHOW_PACKEDBYTES;
    }
    hex_dump(fp, "proof.z1 (poly_qshow_vec_m1)", packed_buf, pos);
    poly_qshow_clear(tmp);
  }

  /* z2 (poly_qshow_vec_m2 = PARAM_M2_SHOW poly_qshow entries) */
  {
    poly_qshow tmp;
    poly_qshow_init(tmp);
    pos = 0;
    for (i = 0; i < PARAM_M2_SHOW; i++) {
      poly_qshow_vec_m2_get_poly(tmp, proof.z2, i);
      poly_qshow_pack(&packed_buf[pos], tmp);
      pos += POLYQSHOW_PACKEDBYTES;
    }
    hex_dump(fp, "proof.z2 (poly_qshow_vec_m2)", packed_buf, pos);
    poly_qshow_clear(tmp);
  }

  fprintf(fp, "\n--- show witness s1 (poly_qshow_vec_m1) ---\n");
  {
    poly_qshow tmp;
    poly_qshow_init(tmp);
    pos = 0;
    for (i = 0; i < PARAM_M1_SHOW; i++) {
      poly_qshow_vec_m1_get_poly(tmp, s1, i);
      poly_qshow_pack(&packed_buf[pos], tmp);
      pos += POLYQSHOW_PACKEDBYTES;
    }
    hex_dump(fp, "witness.s1 (poly_qshow_vec_m1)", packed_buf, pos);
    poly_qshow_clear(tmp);
  }

  /* cleanup */
  sep_keys_clear(&pk, &sk);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  show_proof_clear(&proof);
  poly_q_vec_d_clear(r[0]);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(cmt);
  poly_qshow_vec_m1_clear(s1);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(u_embed[i]);
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_embed[i][j]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + 0]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_D * PARAM_K; j++)
      poly_qshow_mat_k_k_clear(B_embed[i][j]);
    for (j = 0; j < PARAM_K; j++)
      poly_qshow_vec_k_clear(A3_embed[i]); /* was mat_k_k_clear */
    for (j = 0; j < PARAM_M; j++)
      poly_qshow_mat_k_k_clear(D_embed[i][j]);
  }
  free(packed_buf);
  fclose(fp);
}

#endif /* DUMP_PROOF_WITNESS_HEX */

#ifdef BENCHMARK_PROOF_SIZES
/* ── Proof / witness size computation ── */

/* osig_proof_t  (issuance)  —  all fields in packed bytes */
#define OSIG_PROOF_SIZE_BYTES \
    (POLYQISS_VECD_PACKEDBYTES                       /* tA  */ \
   + POLYQISS_VEC256L_PACKEDBYTES                    /* tB  */ \
   + PARAM_ARP_ISS * COEFFQISS_PACKEDBYTES           /* z3  */ \
   + POLYQISS_VECL_PACKEDBYTES                       /* h   */ \
   + POLYQISS_PACKEDBYTES                            /* t1  */ \
   + POLYQISS_PACKEDBYTES                            /* c   */ \
   + sizeof(uint32_t)                                /* ctr_c */ \
   + PARAM_M1_K_ISS * POLYQISS_VECK_PACKEDBYTES      /* z1  */ \
   + PARAM_M2_ISS   * POLYQISS_PACKEDBYTES)          /* z2  */

/* osig witness s1:  PARAM_M1_K_ISS × vec_k  →  total polynomials × packed */
#define OSIG_WITNESS_SIZE_BYTES \
    (PARAM_M1_K_ISS * PARAM_K_ISS * POLYQISS_PACKEDBYTES)

/* show_proof_t  (show)  —  all fields in packed bytes */
#define SHOW_PROOF_SIZE_BYTES \
    (POLYQSHOW_VECD_PACKEDBYTES                      /* tA  */ \
   + POLYQSHOW_VEC256L_PACKEDBYTES                   /* tB  */ \
   + PARAM_ARP_SHOW * COEFFQSHOW_PACKEDBYTES         /* z3  */ \
   + POLYQSHOW_VECL_PACKEDBYTES                      /* h   */ \
   + POLYQSHOW_PACKEDBYTES                           /* t1  */ \
   + POLYQSHOW_PACKEDBYTES                           /* c   */ \
   + sizeof(uint32_t)                                /* ctr_c */ \
   + PARAM_M1_SHOW * POLYQSHOW_PACKEDBYTES           /* z1  */ \
   + PARAM_M2_SHOW * POLYQSHOW_PACKEDBYTES)          /* z2  */

/* show witness s1:  PARAM_M1_SHOW polynomials */
#define SHOW_WITNESS_SIZE_BYTES \
    (PARAM_M1_SHOW * POLYQSHOW_PACKEDBYTES)

static void print_proof_sizes(void)
{
    printf("\n========== Proof / Witness Sizes ==========\n");
    printf("  osig proof  : %zu bytes  (%.2f KB)\n",
           (size_t)OSIG_PROOF_SIZE_BYTES,
           OSIG_PROOF_SIZE_BYTES / 1024.0);
    printf("  osig witness: %zu bytes  (%.2f KB)\n",
           (size_t)OSIG_WITNESS_SIZE_BYTES,
           OSIG_WITNESS_SIZE_BYTES / 1024.0);
    printf("  show proof  : %zu bytes  (%.2f KB)\n",
           (size_t)SHOW_PROOF_SIZE_BYTES,
           SHOW_PROOF_SIZE_BYTES / 1024.0);
    printf("  show witness: %zu bytes  (%.2f KB)\n",
           (size_t)SHOW_WITNESS_SIZE_BYTES,
           SHOW_WITNESS_SIZE_BYTES / 1024.0);
    printf("============================================\n\n");
}
#endif /* BENCHMARK_PROOF_SIZES */

int main(void) {
  int pass = 1;
  arith_setup();
  random_init();
  printf("wkk Hello from the unit tests.\n");

#ifdef DUMP_PROOF_WITNESS_HEX
  dump_osig_proof_and_witness();
  dump_show_proof_and_witness();
#endif

#ifdef BENCHMARK_PROOF_SIZES
  print_proof_sizes();
#endif

  for (int i = 0; i < NTESTS; i++)
  {
    pass &= sep_test();
    pass &= osig_signing_test();
    pass &= osig_proof_test();
    pass &= show_proof_test();
    if (!pass)
    {
      printf("FAILED!\n");
      break;
    } else {
      printf(".");
    }
  }
  if (pass)
  {
    printf("passed.\n");
  }
  arith_teardown();
  return 0;
}
