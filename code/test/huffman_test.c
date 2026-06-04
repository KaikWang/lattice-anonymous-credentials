/*************************************************
 * Full Show Proof Compression Test
 *
 * Compresses EVERY field of show_proof_t:
 *   Uniform (58b/coeff):  tA, tB, h, t1
 *   Challenge (5b/coeff): c
 *   Raw:                   ctr_c (4B)
 *   Huffman (Gaussian):    z1, z2, z3
 *
 * Verifies: roundtrip + proof re-verification.
 *************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "sep.h"
#include "osig.h"
#include "show.h"
#include "huffman_compress.h"
#include "randombytes.h"
#include "random.h"

/* ── Extract coefficients (non-negative, for uniform elements) ── */
static void poly_extract(int64_t *c, const poly_qshow p) {
    for (int j = 0; j < PARAM_N_SHOW; j++) c[j] = poly_qshow_get_coeff(p, j);
}
static void poly_writeback(poly_qshow p, const int64_t *c) {
    for (int j = 0; j < PARAM_N_SHOW; j++) poly_qshow_set_coeff(p, j, c[j]);
}

/* ── Extract centered coefficients (for challenge c: [-8,8]) ── */
static void poly_extract_centered(int64_t *c, const poly_qshow p) {
    for (int j = 0; j < PARAM_N_SHOW; j++) c[j] = poly_qshow_get_coeff_centered(p, j);
}

/* ── vec_d (uniform) ── */
static void vd_extract(int64_t *c, const poly_qshow_vec_d v) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_D_SHOW; i++) {
        poly_qshow_vec_d_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) *c++ = poly_qshow_get_coeff(t, j);
    }
    poly_qshow_clear(t);
}
static void vd_writeback(poly_qshow_vec_d v, const int64_t *c) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_D_SHOW; i++) {
        poly_qshow_vec_d_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) poly_qshow_set_coeff(t, j, *c++);
        poly_qshow_vec_d_set_poly(v, t, i);
    }
    poly_qshow_clear(t);
}

/* ── vec_256_l (uniform) ── */
static void vl256_extract(int64_t *c, const poly_qshow_vec_256_l v) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
        poly_qshow_vec_256_l_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) *c++ = poly_qshow_get_coeff(t, j);
    }
    poly_qshow_clear(t);
}
static void vl256_writeback(poly_qshow_vec_256_l v, const int64_t *c) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
        poly_qshow_vec_256_l_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) poly_qshow_set_coeff(t, j, *c++);
        poly_qshow_vec_256_l_set_poly(v, t, i);
    }
    poly_qshow_clear(t);
}

/* ── vec_l (uniform) ── */
static void vl_extract(int64_t *c, const poly_qshow_vec_l v) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_vec_l_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) *c++ = poly_qshow_get_coeff(t, j);
    }
    poly_qshow_clear(t);
}
static void vl_writeback(poly_qshow_vec_l v, const int64_t *c) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_vec_l_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) poly_qshow_set_coeff(t, j, *c++);
        poly_qshow_vec_l_set_poly(v, t, i);
    }
    poly_qshow_clear(t);
}

/* ── vec_m1 / vec_m2 (Gaussian → centered) ── */
static void vm1_extract(int64_t *c, const poly_qshow_vec_m1 v) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_M1_SHOW; i++) {
        poly_qshow_vec_m1_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) *c++ = poly_qshow_get_coeff_centered(t, j);
    }
    poly_qshow_clear(t);
}
static void vm1_writeback(poly_qshow_vec_m1 v, const int64_t *c) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_M1_SHOW; i++) {
        poly_qshow_vec_m1_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) poly_qshow_set_coeff(t, j, *c++);
        poly_qshow_vec_m1_set_poly(v, t, i);
    }
    poly_qshow_clear(t);
}
static void vm2_extract(int64_t *c, const poly_qshow_vec_m2 v) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_M2_SHOW; i++) {
        poly_qshow_vec_m2_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) *c++ = poly_qshow_get_coeff_centered(t, j);
    }
    poly_qshow_clear(t);
}
static void vm2_writeback(poly_qshow_vec_m2 v, const int64_t *c) {
    poly_qshow t; poly_qshow_init(t);
    for (int i = 0; i < PARAM_M2_SHOW; i++) {
        poly_qshow_vec_m2_get_poly(t, v, i);
        for (int j = 0; j < PARAM_N_SHOW; j++) poly_qshow_set_coeff(t, j, *c++);
        poly_qshow_vec_m2_set_poly(v, t, i);
    }
    poly_qshow_clear(t);
}

/* ── Main test ── */

int main(void) {
    arith_setup();
    random_init();

    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║    FULL Show Proof Compression Test                  ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    /* ── Run full Show proof protocol ── */
    sep_sk_t   sk;     sep_pk_t   pk;
    user_sk_t  usk;    user_pk_t  upk;
    sep_sig_t  sig;    show_proof_t proof;

    poly_q_vec_d r[2], cmt;
    poly_qshow_vec_m1 s1;
    poly_qshow_vec_k  u_embed[PARAM_D];
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_D*PARAM_K];
    poly_qshow_vec_k A3_embed[PARAM_D];
    poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
    uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];

    sep_keys_init(&pk, &sk); user_keys_init(&upk, &usk);
    sep_sig_init(&sig); show_proof_init(&proof);
    poly_q_vec_d_init(r[0]); poly_q_vec_d_init(r[1]); poly_q_vec_d_init(cmt);
    poly_qshow_vec_m1_init(s1);
    for (int i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_init(u_embed[i]);
        for (int j = 0; j < PARAM_D; j++) {
            poly_qshow_mat_k_k_init(A_embed[i][j]);
            poly_qshow_mat_k_k_init(Ds_embed[i][j]);
            poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
        }
        for (int j = 0; j < PARAM_D*PARAM_K; j++) poly_qshow_mat_k_k_init(B_embed[i][j]);
        poly_qshow_vec_k_init(A3_embed[i]);
        for (int j = 0; j < PARAM_M; j++)       poly_qshow_mat_k_k_init(D_embed[i][j]);
    }

    sep_keygen(&pk, &sk);
    osig_user_keygen(&upk, &usk, pk.seed);
    randombytes(state, STATE_BYTES);
    randombytes(crs_seed, CRS_SEED_BYTES);
    randombytes(msg, PARAM_M*PARAM_N/8);
    osig_user_commit(r, cmt, msg, &upk);
    osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
    osig_user_sig_complete(&sig, r);
    if (!osig_user_verify(&sig, &pk, &upk, msg)) {
        printf("ERROR: credential signature verification failed!\n"); goto cleanup;
    }
    show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1,
                    &upk, &usk, &pk, &sig, msg);
    show_user_prove(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                    s1, crs_seed, upk.seed);

    printf("Show proof generated.\n\n");

    /* ── Field definitions ── */
    #define N_TA  (PARAM_D_SHOW * PARAM_N_SHOW)             /* 1472 */
    #define N_TB  (PARAM_ARP_DIV_N_L_SHOW * PARAM_N_SHOW)  /* 704  */
    #define N_H   (PARAM_L_SHOW * PARAM_N_SHOW)             /* 448  */
    #define N_T1  PARAM_N_SHOW                              /* 64   */
    #define N_C   PARAM_N_SHOW                              /* 64   */
    #define N_Z1  (PARAM_M1_SHOW * PARAM_N_SHOW)            /* 13504 */
    #define N_Z2  (PARAM_M2_SHOW * PARAM_N_SHOW)            /* 4736  */
    #define N_Z3  PARAM_ARP_SHOW                            /* 256   */

    /* allocate buffers */
    enum { F_TA=0, F_TB, F_H, F_T1, F_C, F_CTRC, F_Z3, F_Z1, F_Z2, N_FIELDS };
    #define NF N_FIELDS
    int64_t *orig[NF], *rest[NF];
    uint8_t *comp[NF];
    size_t   counts[NF], naive_sizes[NF], comp_sizes[NF], max_sizes[NF];
    const char *names[NF];

    counts[F_TA]=N_TA; counts[F_TB]=N_TB; counts[F_H]=N_H; counts[F_T1]=N_T1;
    counts[F_C]=N_C; counts[F_CTRC]=1; counts[F_Z3]=N_Z3; counts[F_Z1]=N_Z1; counts[F_Z2]=N_Z2;

    names[F_TA]="tA"; names[F_TB]="tB"; names[F_H]="h"; names[F_T1]="t1";
    names[F_C]="c"; names[F_CTRC]="ctr_c"; names[F_Z3]="z3"; names[F_Z1]="z1"; names[F_Z2]="z2";

    for (int f = 0; f < N_FIELDS; f++) {
        orig[f] = (int64_t *)malloc((counts[f] + 1) * sizeof(int64_t));
        rest[f] = (int64_t *)calloc(counts[f] + 1, sizeof(int64_t));
        naive_sizes[f] = counts[f] * sizeof(int64_t);
    }
    /* special: ctr_c */
    naive_sizes[F_CTRC] = sizeof(uint32_t);

    /* ── Extract all fields from the proof ── */
    vd_extract(orig[F_TA], proof.tA);
    vl256_extract(orig[F_TB], proof.tB);
    vl_extract(orig[F_H], proof.h);
    poly_extract(orig[F_T1], proof.t1);
    poly_extract_centered(orig[F_C], proof.c);
    orig[F_CTRC][0] = (int64_t)(uint64_t)proof.ctr_c;
    vm1_extract(orig[F_Z1], proof.z1);
    vm2_extract(orig[F_Z2], proof.z2);
    for (size_t i = 0; i < N_Z3; i++) orig[F_Z3][i] = proof.z3[i];

    /* ── Compute max compressed buffer sizes ── */
    max_sizes[F_TA] = COMPACT_UNIFORM_BYTES(N_TA) + 8;
    max_sizes[F_TB] = COMPACT_UNIFORM_BYTES(N_TB) + 8;
    max_sizes[F_H]  = COMPACT_UNIFORM_BYTES(N_H) + 8;
    max_sizes[F_T1] = COMPACT_UNIFORM_BYTES(N_T1) + 8;
    max_sizes[F_C]  = COMPACT_CHALLENGE_BYTES(N_C) + 8;
    max_sizes[F_CTRC] = 8;
    max_sizes[F_Z3] = HUFFMAN_Z3_COMPRESSED_BYTES_MAX;
    max_sizes[F_Z1] = HUFFMAN_Z1_COMPRESSED_BYTES_MAX;
    max_sizes[F_Z2] = HUFFMAN_Z2_COMPRESSED_BYTES_MAX;
    for (int f = 0; f < N_FIELDS; f++) comp[f] = (uint8_t *)malloc(max_sizes[f]);

    /* ── Compress each field ── */
    comp_sizes[F_TA]   = compact_pack_uniform(comp[F_TA], orig[F_TA], N_TA);
    comp_sizes[F_TB]   = compact_pack_uniform(comp[F_TB], orig[F_TB], N_TB);
    comp_sizes[F_H]    = compact_pack_uniform(comp[F_H],  orig[F_H],  N_H);
    comp_sizes[F_T1]   = compact_pack_uniform(comp[F_T1], orig[F_T1], N_T1);
    comp_sizes[F_C]    = compact_pack_challenge(comp[F_C], orig[F_C], N_C);
    comp_sizes[F_CTRC] = 4; memcpy(comp[F_CTRC], &proof.ctr_c, 4);
    comp_sizes[F_Z3]   = huffman_pack_z3(comp[F_Z3], orig[F_Z3], N_Z3);
    comp_sizes[F_Z1]   = huffman_pack_z1(comp[F_Z1], orig[F_Z1], N_Z1);
    comp_sizes[F_Z2]   = huffman_pack_z2(comp[F_Z2], orig[F_Z2], N_Z2);

    /* ── Decompress each field ── */
    compact_unpack_uniform(rest[F_TA], comp[F_TA], comp_sizes[F_TA], N_TA);
    compact_unpack_uniform(rest[F_TB], comp[F_TB], comp_sizes[F_TB], N_TB);
    compact_unpack_uniform(rest[F_H],  comp[F_H],  comp_sizes[F_H],  N_H);
    compact_unpack_uniform(rest[F_T1], comp[F_T1], comp_sizes[F_T1], N_T1);
    compact_unpack_challenge(rest[F_C], comp[F_C], comp_sizes[F_C], N_C);
    { uint32_t tmp; memcpy(&tmp, comp[F_CTRC], 4); rest[F_CTRC][0] = (int64_t)(uint64_t)tmp; }
    huffman_unpack_z3(rest[F_Z3], comp[F_Z3], comp_sizes[F_Z3], N_Z3);
    huffman_unpack_z1(rest[F_Z1], comp[F_Z1], comp_sizes[F_Z1], N_Z1);
    huffman_unpack_z2(rest[F_Z2], comp[F_Z2], comp_sizes[F_Z2], N_Z2);

    /* ── Verify roundtrip ── */
    int all_ok = 1;
    for (int f = 0; f < N_FIELDS; f++) {
        int errs = 0;
        for (size_t i = 0; i < counts[f]; i++) {
            if (orig[f][i] != rest[f][i]) {
                if (errs < 3) printf("  MISMATCH %s[%zu]: %" PRId64 " != %" PRId64 "\n",
                                      names[f], i, orig[f][i], rest[f][i]);
                errs++;
            }
        }
        if (errs) { printf("  *** %s: %d errors ***\n", names[f], errs); all_ok = 0; }
    }

    /* ── Write back and re-verify proof ── */
    if (all_ok) {
        vd_writeback(proof.tA, rest[F_TA]);
        vl256_writeback(proof.tB, rest[F_TB]);
        vl_writeback(proof.h, rest[F_H]);
        poly_writeback(proof.t1, rest[F_T1]);
        poly_writeback(proof.c, rest[F_C]);
        proof.ctr_c = (uint32_t)(uint64_t)rest[F_CTRC][0];
        vm1_writeback(proof.z1, rest[F_Z1]);
        vm2_writeback(proof.z2, rest[F_Z2]);
        for (size_t i = 0; i < N_Z3; i++) proof.z3[i] = rest[F_Z3][i];

        int ver = show_verify(&proof, A_embed, B_embed, A3_embed,
                              Ds_embed, D_embed, u_embed, crs_seed, upk.seed);
        printf("\nProof re-verification after full roundtrip: %s\n",
               ver ? "PASS" : "FAIL");
    }

    /* ── Print per-field sizes ── */
    printf("\n┌───────┬──────────┬──────────┬──────────┬──────────┬──────────┐\n");
    printf("│ Field │  coeffs  │  Naive   │Compressed│  Ratio   │  Method  │\n");
    printf("├───────┼──────────┼──────────┼──────────┼──────────┼──────────┤\n");

    size_t tot_naive = 0, tot_comp = 0;
    const char *method[] = {"58b","58b","58b","58b","5b","raw","Huff","Huff","Huff"};

    for (int f = 0; f < N_FIELDS; f++) {
        printf("│ %-5s │ %8zu │ %7zu B │ %7zu B │ %7.2fx │ %-8s │\n",
               names[f], counts[f], naive_sizes[f], comp_sizes[f],
               (double)naive_sizes[f] / (double)comp_sizes[f], method[f]);
        tot_naive += naive_sizes[f];
        tot_comp  += comp_sizes[f];
    }

    printf("├───────┼──────────┼──────────┼──────────┼──────────┼──────────┤\n");
    printf("│ TOTAL │ %8zu │ %7zu B │ %7zu B │ %7.2fx │          │\n",
           counts[F_TA]+counts[F_TB]+counts[F_H]+counts[F_T1]+counts[F_C]
           +counts[F_CTRC]+counts[F_Z3]+counts[F_Z1]+counts[F_Z2],
           tot_naive, tot_comp, (double)tot_naive/(double)tot_comp);
    printf("└───────┴──────────┴──────────┴──────────┴──────────┴──────────┘\n");

    printf("\n  Total naive:      %7zu B = %.2f KB\n", tot_naive, tot_naive/1024.0);
    printf("  Total compressed:  %7zu B = %.2f KB\n", tot_comp,  tot_comp/1024.0);
    printf("  Savings:           %7zu B = %.2f KB (%.1f%%)\n",
           tot_naive - tot_comp, (tot_naive-tot_comp)/1024.0,
           100.0*(tot_naive-tot_comp)/tot_naive);

    /* ── Theoretical minimum breakdown ── */
    printf("\n═══════════════════════════════════════════\n");
    printf("  THEORETICAL MINIMUM (per field)\n");
    printf("═══════════════════════════════════════════\n");
    printf("  Uniform (58b/coeff):\n");
    printf("    tA: %5d B  tB: %5d B  h: %5d B  t1: %4d B\n",
           COMPACT_TA_BYTES, COMPACT_TB_BYTES, COMPACT_H_BYTES, COMPACT_T1_BYTES);
    printf("  Challenge (5b/coeff):\n");
    printf("    c: %5d B\n", COMPACT_C_BYTES);
    printf("  Raw:  ctr_c: %5d B\n", COMPACT_CTRC_BYTES);
    printf("  Huffman (Gaussian):\n");
    printf("    z3: ~926 B  z1: ~52818 B  z2: ~12091 B\n");
    printf("  ───────────────────────────────\n");
    { int theory_total = COMPACT_TA_BYTES + COMPACT_TB_BYTES + COMPACT_H_BYTES +
                         COMPACT_T1_BYTES + COMPACT_C_BYTES + COMPACT_CTRC_BYTES +
                         926 + 52818 + 12091;
      printf("  Theoretical total minimum: ~%d B = %.2f KB\n",
             theory_total, theory_total / 1024.0); }
    printf("\n");

    /* cleanup */
    for (int f = 0; f < N_FIELDS; f++) { free(orig[f]); free(rest[f]); free(comp[f]); }

cleanup:
    poly_qshow_vec_m1_clear(s1);
    poly_q_vec_d_clear(cmt); poly_q_vec_d_clear(r[1]); poly_q_vec_d_clear(r[0]);
    show_proof_clear(&proof); sep_sig_clear(&sig);
    user_keys_clear(&upk, &usk); sep_keys_clear(&pk, &sk);
    for (int i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_clear(u_embed[i]);
        for (int j = 0; j < PARAM_D; j++) {
            poly_qshow_mat_k_k_clear(A_embed[i][j]);
            poly_qshow_mat_k_k_clear(Ds_embed[i][j]);
            poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
        }
        for (int j = 0; j < PARAM_D*PARAM_K; j++) poly_qshow_mat_k_k_clear(B_embed[i][j]);
        poly_qshow_vec_k_clear(A3_embed[i]);
        for (int j = 0; j < PARAM_M; j++)       poly_qshow_mat_k_k_clear(D_embed[i][j]);
    }
    arith_teardown();
    return 0;
}
