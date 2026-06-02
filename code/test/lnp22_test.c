/**
 * LNP22 Zero-Knowledge Proof for Linear Relation with Short Witness
 *
 * Proves knowledge of a vector s such that:
 *   P · s = u  (mod qR)
 *   ||s|| < B_s
 *
 * where P (matrix), u (vector), B_s (scalar bound) are public.
 *
 * The protocol follows the 5-round structure from the show proof,
 * but simplified to a single linear constraint and a single norm bound.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "arith.h"
#include "randombytes.h"
#include "random.h"
#include "poly_qshow_sampling.h"
#include "sep.h"
#include "show.h"
#include "macros.h"
#include "fips202.h"
#include "four_squares.h"

/* ────────────────────────────────────────────
 * Domain separators for LNP22 proof
 * ──────────────────────────────────────────── */
#define DOMAIN_SEPARATOR_A1_LNP22    128
#define DOMAIN_SEPARATOR_A2_LNP22    129
#define DOMAIN_SEPARATOR_BYG_LNP22   130
#define DOMAIN_SEPARATOR_B_LNP22     131
#define DOMAIN_SEPARATOR_CHAL1_LNP22 132
#define DOMAIN_SEPARATOR_CHAL2_LNP22 133
#define DOMAIN_SEPARATOR_CHAL3_LNP22 134
#define DOMAIN_SEPARATOR_CHAL4_LNP22 135
#define DOMAIN_SEPARATOR_RAND_S2_LNP22 136
#define DOMAIN_SEPARATOR_RAND_G_LNP22  137

/* ────────────────────────────────────────────
 * LNP22 proof parameters
 *
 * The witness s has LNP22_M entries (polynomials in R_qshow).
 * The slack polynomial a is appended to make ||s||² + ||a||² = B_s².
 * Total witness s1 has LNP22_M + 1 entries; we use the existing
 * poly_qshow_vec_m1 (PARAM_M1_SHOW = 211 entries), filling only
 * the first LNP22_M+1, and zeroing the rest.
 *
 * Linear constraint:  P is D × LNP22_M, u is D vector.
 * (We use PARAM_D = 4 rows to match the μ-challenge count.)
 * ──────────────────────────────────────────── */
#define LNP22_M      8       /* witness dimension (without slack) */
#define LNP22_BS     10      /* norm bound B_s on the witness s   */
#define LNP22_BSSQ   ((uint64_t)(LNP22_BS * LNP22_BS))

/* Gaussian rejection sampling parameters (test-appropriate values) */
#define LNP22_REJ1   2
#define LNP22_REJ2   2
#define LNP22_REJ3   2

/* Squared verification bounds for z1, z2, z3.
 * Use the same bounds as the show proof since we reuse the same Gaussian widths. */
#define LNP22_BOUND_Z1_LO  PARAM_B1SQ_SHOW_LOW64
#define LNP22_BOUND_Z1_HI  PARAM_B1SQ_SHOW_HIGH64
#define LNP22_BOUND_Z2     PARAM_B2SQ_SHOW
#define LNP22_BOUND_Z3     PARAM_B3SQ_SHOW

/* Number of quadratic constraints (just 1: the norm bound).
 * NOTE: vec_qshow_uniform always writes PARAM_ARP_SHOW + 6 entries,
 * so chal_2 arrays must be sized accordingly, even though we only
 * use the first PARAM_ARP_SHOW + 1 entries. */
#define LNP22_NUM_QUAD      1
#define LNP22_CHAL2_WIDTH   (PARAM_ARP_SHOW + 6)  /* must match vec_qshow_uniform */

/* Byte-length computations */
#define LNP22_CHAL_BASE_BYTES (1 + CRS_SEED_BYTES + SEED_BYTES)
#define LNP22_CHAL1_BYTES \
    (LNP22_CHAL_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES + POLYQSHOW_VEC256L_PACKEDBYTES)
#define LNP22_CHAL2_BYTES \
    (LNP22_CHAL1_BYTES + PARAM_ARP_SHOW*COEFFQSHOW_PACKEDBYTES)
#define LNP22_CHAL3_BYTES \
    (LNP22_CHAL2_BYTES + POLYQSHOW_VECL_PACKEDBYTES)
#define LNP22_CHAL4_BYTES \
    (LNP22_CHAL3_BYTES + 2*POLYQSHOW_PACKEDBYTES)

/* ─────────────────────────────────────────────────────
 * LNP22 proof structure (mirrors show_proof_t)
 * ───────────────────────────────────────────────────── */
typedef struct {
    poly_qshow_vec_d      tA;
    poly_qshow_vec_256_l  tB;
    coeff_qshow           z3[PARAM_ARP_SHOW];
    poly_qshow_vec_l      h;
    poly_qshow            t1;
    poly_qshow            c;
    uint32_t              ctr_c;
    poly_qshow_vec_m1     z1;
    poly_qshow_vec_m2     z2;
} lnp22_proof_t;

/* ────────────────────────────────────────────
 * Rejection sampling (same as show_prove.c)
 * ──────────────────────────────────────────── */
static int _lnp22_reject_exp(double threshold) {
    uint64_t u = 0;
    double udbl = 0;
    randombytes((uint8_t *)&u, 6);
    udbl = (double)u / (double)(1ul << 48);
    if (udbl > threshold) return 1;
    return 0;
}

/* ────────────────────────────────────────────
 * Initialize / clear proof structure
 * ──────────────────────────────────────────── */
static void lnp22_proof_init(lnp22_proof_t *proof) {
    poly_qshow_vec_m1_init(proof->z1);
    poly_qshow_vec_d_init(proof->tA);
    poly_qshow_vec_256_l_init(proof->tB);
    poly_qshow_vec_l_init(proof->h);
    poly_qshow_init(proof->t1);
    poly_qshow_init(proof->c);
    poly_qshow_vec_m2_init(proof->z2);
}

static void lnp22_proof_clear(lnp22_proof_t *proof) {
    poly_qshow_vec_m1_clear(proof->z1);
    poly_qshow_vec_d_clear(proof->tA);
    poly_qshow_vec_256_l_clear(proof->tB);
    poly_qshow_vec_l_clear(proof->h);
    poly_qshow_clear(proof->t1);
    poly_qshow_clear(proof->c);
    poly_qshow_vec_m2_clear(proof->z2);
}

/* ─────────────────────────────────────────────────────────
 * lnp22_embed:  Build witness s1 from secret s.
 *
 *   s1 = [ s[0..M-1] || slack || 0...0 ]
 *
 * where slack is the four-square decomposition of B_s² - ||s||².
 *
 * Also initialises the statement matrices P_embed and u_embed
 * (the public linear-constraint parameters).
 * ───────────────────────────────────────────────────────── */
static void lnp22_embed(
    poly_qshow_vec_m1          s1,
    poly_qshow                 P_embed[PARAM_D][LNP22_M],
    poly_qshow                 u_embed[PARAM_D],
    const poly_qshow           s[LNP22_M],
    const poly_qshow           P[PARAM_D][LNP22_M],
    const poly_qshow           u[PARAM_D],
    uint64_t                   b_s_sq)
{
    size_t i, j;
    uint64_t norm2_s, gap;
    uint64_t four_squares_res[4];
    poly_qshow tmp;

    poly_qshow_init(tmp);

    /* ── Compute ||s||² ── */
    norm2_s = 0;
    for (i = 0; i < LNP22_M; i++) {
        norm2_s += (uint64_t)poly_qshow_sq_norm2(s[i]);
    }

    /* ── Compute slack via four-square decomposition ── */
    assert(norm2_s <= b_s_sq);
    gap = b_s_sq - norm2_s;
    four_squares(four_squares_res, gap);

    /* ── Build s1: [ s[0..M-1] | slack | 0...0 ] ── */
    for (i = 0; i < LNP22_M; i++) {
        poly_qshow_set(s1->entries[i], s[i]);
    }
    /* slack polynomial: set first 4 coefficients to the 4-square roots */
    poly_qshow_zero(s1->entries[LNP22_M]);
    for (i = 0; i < 4; i++) {
        poly_qshow_set_coeff(s1->entries[LNP22_M], i, (coeff_qshow)four_squares_res[i]);
    }
    /* remaining entries stay zero */
    for (i = LNP22_M + 1; i < PARAM_M1_SHOW; i++) {
        poly_qshow_zero(s1->entries[i]);
    }

    /* ── Copy P and u into the embed arrays ── */
    for (i = 0; i < PARAM_D; i++) {
        for (j = 0; j < LNP22_M; j++) {
            poly_qshow_set(P_embed[i][j], P[i][j]);
        }
        poly_qshow_set(u_embed[i], u[i]);
    }

    poly_qshow_clear(tmp);
}

/* ═══════════════════════════════════════════════════════
 * ROUND 1: Commitments and first challenge
 * ═══════════════════════════════════════════════════════ */
static void lnp22_prove_round1(
    lnp22_proof_t              *proof,
    poly_qshow_vec_m1          chal_1[PARAM_ARP_SHOW],
    uint8_t                    buf[LNP22_CHAL1_BYTES],
    poly_qshow_vec_m2          s2,
    poly_qshow_vec_m1          y1,
    poly_qshow_vec_m2          y2,
    poly_qshow_vec_256_l       y3_g,
    const poly_qshow_vec_m1    s1,
    const poly_qshow_mat_d_m1  A1,
    const poly_qshow_mat_d_m2  A2,
    const poly_qshow_mat_256l_m2 Byg,
    const uint8_t              randomness_seed[SEED_BYTES],
    const uint32_t             kappa)
{
    size_t i, j;
    uint32_t kpp = kappa;
    coeff_qshow tmp_coeff;
    poly_qshow_vec_d w, tmp_vec_d;
    uint8_t binomial_seed[SEED_BYTES];

    poly_qshow_vec_d_init(w);
    poly_qshow_vec_d_init(tmp_vec_d);

    /* Sample ABDLOP commitment randomness */
    poly_qshow_vec_m2_binomial(s2, randomness_seed, kpp++,
                               DOMAIN_SEPARATOR_RAND_S2_LNP22);

    /* Sample Gaussian masks for c·s1, c·s2, R·s1 */
    poly_qshow_vec_m1_sample_gaussian_s1(y1);
    poly_qshow_vec_m2_sample_gaussian_s2(y2);
    for (i = 0; i < PARAM_ARP_DIV_N_SHOW; i++) {
        for (j = 0; j < PARAM_N_SHOW; j++) {
            tmp_coeff = SampleZ(0, PARAM_S3_SHOW);
            poly_qshow_set_coeff(y3_g->entries[i], j, tmp_coeff);
        }
    }
    /* Uniform masks for automorphism equations */
    for (i = PARAM_ARP_DIV_N_SHOW; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
        poly_qshow_uniform_but_zero(y3_g->entries[i], randomness_seed, kpp++,
                                    DOMAIN_SEPARATOR_RAND_G_LNP22);
    }

    /* tA = A1·s1 + A2·s2,  w = A1·y1 + A2·y2,  tB = Byg·s2 + y3_g */
    poly_qshow_mat_d_m1_mul_vec_m1(proof->tA, A1, s1);
    poly_qshow_mat_d_m2_mul_vec_m2(tmp_vec_d, A2, s2);
    poly_qshow_vec_d_add(proof->tA, proof->tA, tmp_vec_d);

    poly_qshow_mat_d_m1_mul_vec_m1(w, A1, y1);
    poly_qshow_mat_d_m2_mul_vec_m2(tmp_vec_d, A2, y2);
    poly_qshow_vec_d_add(w, w, tmp_vec_d);

    poly_qshow_mat_256l_m2_mul_vec_m2(proof->tB, Byg, s2);
    poly_qshow_vec_256_l_add(proof->tB, proof->tB, y3_g);

    /* First challenge */
    buf[0] = 1;
    poly_qshow_vec_d_pack(
        buf + LNP22_CHAL_BASE_BYTES, proof->tA);
    poly_qshow_vec_d_pack(
        buf + LNP22_CHAL_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w);
    poly_qshow_vec_256_l_pack(
        buf + LNP22_CHAL_BASE_BYTES + 2 * POLYQSHOW_VECD_PACKEDBYTES,
        proof->tB);
    shake256(binomial_seed, SEED_BYTES, buf, LNP22_CHAL1_BYTES);
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        poly_qshow_vec_m1_binomial(chal_1[i], binomial_seed,
                                   DOMAIN_SEPARATOR_CHAL1_LNP22, i, SEED_BYTES);
    }

    poly_qshow_vec_d_clear(w);
    poly_qshow_vec_d_clear(tmp_vec_d);
}

/* ═══════════════════════════════════════════════════════
 * ROUND 2: Approximate Range Proof + second challenge
 * ═══════════════════════════════════════════════════════ */
static int lnp22_prove_round2(
    lnp22_proof_t              *proof,
    coeff_qshow                chal_2[PARAM_L_SHOW][LNP22_CHAL2_WIDTH],
    uint8_t                    buf[LNP22_CHAL2_BYTES],
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],
    const poly_qshow_vec_256_l y3_g)
{
    size_t i, j, k;
    int64_t tmp;
    uint64_t sq_norm_y3 = 0, sq_norm_z3 = 0;
    uint8_t challenge_seed[SEED_BYTES];
    coeff_qshow tmp_coeff;

    /* z3 = R·s1 + y3  (in Z, not in the ring) */
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        tmp = 0;
        for (j = 0; j < PARAM_M1_SHOW; j++) {
            for (k = 0; k < PARAM_N_SHOW; k++) {
                tmp += poly_qshow_get_coeff_centered(chal_1[i]->entries[j], k) *
                       poly_qshow_get_coeff_centered(s1->entries[j], k);
            }
        }
        tmp_coeff = poly_qshow_get_coeff_centered(
            y3_g->entries[i / PARAM_N_SHOW], i % PARAM_N_SHOW);
        tmp += tmp_coeff;
        CHK_UI_OVF_ADDITION(sq_norm_y3, tmp_coeff * tmp_coeff);
        CHK_UI_OVF_ADDITION(sq_norm_z3, tmp * tmp);
        proof->z3[i] = tmp;
    }

    /* Rejection sampling */
    if (_lnp22_reject_exp(
            exp(M_PI * (double)(sq_norm_y3 - sq_norm_z3) /
                (double)PARAM_S3SQ_SHOW) /
            (double)LNP22_REJ3)) {
        return 0;
    }

    /* Second challenge */
    buf[0] = 2;
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        coeff_qshow_pack(
            buf + LNP22_CHAL1_BYTES + i * COEFFQSHOW_PACKEDBYTES,
            proof->z3[i]);
    }
    shake256(challenge_seed, SEED_BYTES, buf, LNP22_CHAL2_BYTES);
    for (i = 0; i < PARAM_L_SHOW; i++) {
        vec_qshow_uniform(chal_2[i], challenge_seed,
                          DOMAIN_SEPARATOR_CHAL2_LNP22, i, SEED_BYTES);
    }
    return 1;
}

/* ═══════════════════════════════════════════════════════
 * ROUND 3: h-vector + third challenge
 * ═══════════════════════════════════════════════════════ */
static void lnp22_prove_round3(
    lnp22_proof_t              *proof,
    poly_qshow_vec_l           chal_3_l,
    poly_qshow_vec_k           chal_3_dk[PARAM_D],
    uint8_t                    buf[LNP22_CHAL3_BYTES],
    poly_qshow_vec_256         sum_gamma_e_star[PARAM_L_SHOW],
    poly_qshow_vec_m1          sum_gamma_r_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],
    const coeff_qshow          chal_2[PARAM_L_SHOW][LNP22_CHAL2_WIDTH],
    const poly_qshow_vec_256_l y3_g,
    const poly_qshow           quadratic_precomp)
{
    size_t i, j, k;
    poly_qshow tmp_poly;
    poly_qshow_vec_m1 tmp_vec_m1, rj_star;
    uint8_t challenge_seed[SEED_BYTES];

    poly_qshow_init(tmp_poly);
    poly_qshow_vec_m1_init(tmp_vec_m1);
    poly_qshow_vec_m1_init(rj_star);

    /* Σ γ_{ij} · r_j* */
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
        poly_qshow_vec_m1_conjugate(rj_star, chal_1[j]);
        for (i = 0; i < PARAM_L_SHOW; i++) {
            if (j == 0) {
                poly_qshow_vec_m1_mul_scalar(
                    sum_gamma_r_star[i], rj_star, chal_2[i][j]);
            } else {
                poly_qshow_vec_m1_mul_scalar(
                    tmp_vec_m1, rj_star, chal_2[i][j]);
                poly_qshow_vec_m1_add(
                    sum_gamma_r_star[i], sum_gamma_r_star[i], tmp_vec_m1);
            }
        }
    }

    /* Σ γ_{ij} · e_j* = conjugate(τ⁻¹([γ_{i,0} | ... | γ_{i,ARP}])) */
    for (i = 0; i < PARAM_L_SHOW; i++) {
        for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
            poly_qshow_set_coeff(
                sum_gamma_e_star[i]->entries[j], 0,
                chal_2[i][j * PARAM_N_SHOW]);
            for (k = 1; k < PARAM_N_SHOW; k++) {
                poly_qshow_set_coeff(
                    sum_gamma_e_star[i]->entries[j], k,
                    -chal_2[i][(j + 1) * PARAM_N_SHOW - k]);
            }
        }
    }

    /*
     * h_i =  g_i
     *       - Σ_j γ_{ij}·z3_j
     *       + Σ_j γ_{ij}·e_j*·y3
     *       + Σ_j γ_{ij}·r_j*·s1
     *       + γ_{i,ARP} · (⟨s1*,s1⟩ − B_s²)
     */
    for (i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_set(tmp_poly,
                       y3_g->entries[PARAM_ARP_DIV_N_SHOW + i]);

        /* − Σ γ_{ij}·z3_j */
        for (j = 0; j < PARAM_ARP_SHOW; j++) {
            poly_qshow_muladd_constant(tmp_poly, chal_2[i][j],
                                       -proof->z3[j]);
        }

        poly_qshow_set(proof->h->entries[i], tmp_poly);

        /* + Σ γ_{ij}·e_j*·y3 */
        for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
            poly_qshow_mul(tmp_poly, sum_gamma_e_star[i]->entries[j],
                           y3_g->entries[j]);
            poly_qshow_add(proof->h->entries[i],
                           proof->h->entries[i], tmp_poly);
        }

        /* + Σ γ_{ij}·r_j*·s1 */
        poly_qshow_vec_m1_mul_inner(tmp_poly,
                                    sum_gamma_r_star[i], s1);
        poly_qshow_add(proof->h->entries[i],
                       proof->h->entries[i], tmp_poly);

        /* + γ_{i,ARP} · quadratic_precomp */
        poly_qshow_mul_scalar(tmp_poly, quadratic_precomp,
                              chal_2[i][PARAM_ARP_SHOW]);
        poly_qshow_add(proof->h->entries[i],
                       proof->h->entries[i], tmp_poly);
    }

    /* Third challenge */
    buf[0] = 3;
    poly_qshow_vec_l_pack(buf + LNP22_CHAL2_BYTES, proof->h);
    poly_qshow_vec_l_uniform(chal_3_l, buf, DOMAIN_SEPARATOR_CHAL3_LNP22,
                             LNP22_CHAL3_BYTES);
    shake256(challenge_seed, SEED_BYTES, buf, LNP22_CHAL3_BYTES);
    for (i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_uniform(chal_3_dk[i], challenge_seed,
                                 DOMAIN_SEPARATOR_CHAL3_LNP22,
                                 i + PARAM_L_SHOW, SEED_BYTES);
    }

    poly_qshow_clear(tmp_poly);
    poly_qshow_vec_m1_clear(tmp_vec_m1);
    poly_qshow_vec_m1_clear(rj_star);
}

/* ═════════════════════════════════════════════════════════
 * ROUND 4: Garbage terms t0, t1 + fourth challenge c
 *
 *   e0 = Σ μ_i·γ_{i,ARP} · ⟨y1*, y1⟩
 *
 *   e1 = Σ μ_i·γ_{i,ARP} · (⟨y1*, s1⟩ + ⟨s1*, y1⟩)
 *        + Σ_{i<d} μ_{ℓ+i} · (P · y1)_i
 *
 *   t0 = ⟨b, y2⟩ + e0
 *   t1 = ⟨b, s2⟩ + e1
 * ═══════════════════════════════════════════════════════ */
static void lnp22_prove_round4(
    lnp22_proof_t              *proof,
    uint8_t                    buf[LNP22_CHAL4_BYTES],
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m2    s2,
    const poly_qshow_mat_256l_m2 Byg,
    const poly_qshow_vec_m2    b,
    const poly_qshow           P_embed[PARAM_D][LNP22_M],
    const poly_qshow_vec_256   sum_gamma_e_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1    sum_gamma_r_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1    y1,
    const poly_qshow_vec_m2    y2,
    const coeff_qshow          chal_2[PARAM_L_SHOW][LNP22_CHAL2_WIDTH],
    const poly_qshow_vec_l     chal_3_l,
    const poly_qshow_vec_k     chal_3_dk[PARAM_D])
{
    size_t i, j;
    uint32_t kappa_c;
    poly_qshow sum_mu_gamma, tmp_poly, e0, e1, t0;
    poly_qshow_vec_256_l tmp_vec_256_l;
    uint8_t challenge_seed[SEED_BYTES];

    poly_qshow_init(sum_mu_gamma);
    poly_qshow_zero(sum_mu_gamma);
    poly_qshow_init(tmp_poly);
    poly_qshow_init(e0);
    poly_qshow_zero(e0);
    poly_qshow_init(e1);
    poly_qshow_zero(e1);
    poly_qshow_init(t0);
    poly_qshow_vec_256_l_init(tmp_vec_256_l);

    /* ── Σ μ_i · γ_{i,ARP} ── */
    for (i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i],
                              chal_2[i][PARAM_ARP_SHOW]);
        poly_qshow_add(sum_mu_gamma, sum_mu_gamma, tmp_poly);
    }

    /* ── Quadratic terms:  ⟨y1*, y1⟩  and  ⟨y1*, s1⟩ + ⟨s1*, y1⟩ ── */
    {
        poly_qshow y1i_star, y1s_y1, y1s_s1;
        poly_qshow_init(y1i_star);
        poly_qshow_init(y1s_y1);
        poly_qshow_init(y1s_s1);
        poly_qshow_zero(y1s_y1);
        poly_qshow_zero(y1s_s1);

        for (i = 0; i < PARAM_M1_SHOW; i++) {
            poly_qshow_conjugate(y1i_star, y1->entries[i]);

            poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
            poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);

            poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
            poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
        }

        /* e0 = sum_mu_gamma · ⟨y1*, y1⟩ */
        poly_qshow_mul(e0, sum_mu_gamma, y1s_y1);

        /* e1 = sum_mu_gamma · (⟨y1*, s1⟩ + ⟨s1*, y1⟩) */
        poly_qshow_conjugate(tmp_poly, y1s_s1);
        poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
        poly_qshow_mul(e1, sum_mu_gamma, y1s_s1);

        poly_qshow_clear(y1i_star);
        poly_qshow_clear(y1s_y1);
        poly_qshow_clear(y1s_s1);
    }

    /* ── Linear terms from sum_gamma_r_star / sum_gamma_e_star ── */
    poly_qshow_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, y2);
    for (i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow y1s_y1;
        poly_qshow_init(y1s_y1);

        poly_qshow_vec_m1_mul_inner(y1s_y1, sum_gamma_r_star[i], y1);
        poly_qshow_sub(y1s_y1, y1s_y1,
                       tmp_vec_256_l->entries[PARAM_ARP_DIV_N_SHOW + i]);
        for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
            poly_qshow_mul(tmp_poly, sum_gamma_e_star[i]->entries[j],
                           tmp_vec_256_l->entries[j]);
            poly_qshow_sub(y1s_y1, y1s_y1, tmp_poly);
        }
        poly_qshow_mul(y1s_y1, y1s_y1, chal_3_l->entries[i]);
        poly_qshow_add(e1, e1, y1s_y1);

        poly_qshow_clear(y1s_y1);
    }

    /* ── Linear constraint:  Σ_{i<d} μ_{ℓ+i} · (P · y1)_i ── */
    for (i = 0; i < PARAM_D; i++) {
        poly_qshow y1s_y1;
        poly_qshow_init(y1s_y1);
        poly_qshow_zero(y1s_y1);

        for (j = 0; j < LNP22_M; j++) {
            poly_qshow_mul(tmp_poly, P_embed[i][j], y1->entries[j]);
            poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
        }
        poly_qshow_mul(y1s_y1, y1s_y1,
                       chal_3_dk[i]->entries[0]);
        poly_qshow_add(e1, e1, y1s_y1);

        poly_qshow_clear(y1s_y1);
    }

    /* ── Garbage terms ── */
    poly_qshow_vec_m2_mul_inner(t0, b, y2);
    poly_qshow_add(t0, t0, e0);

    poly_qshow_vec_m2_mul_inner(tmp_poly, b, s2);
    poly_qshow_add(proof->t1, tmp_poly, e1);

    /* ── Fourth challenge c ── */
    kappa_c = 0;
    buf[0] = 4;
    poly_qshow_pack(buf + LNP22_CHAL3_BYTES, t0);
    poly_qshow_pack(buf + LNP22_CHAL3_BYTES + POLYQSHOW_PACKEDBYTES,
                    proof->t1);
    shake256(challenge_seed, SEED_BYTES, buf, LNP22_CHAL4_BYTES);
    do {
        poly_qshow_sample_challenge(proof->c, challenge_seed,
                                    DOMAIN_SEPARATOR_CHAL4_LNP22,
                                    kappa_c++, SEED_BYTES);
    } while (challenge_size_show(proof->c) > PARAM_ETA_SHOW);
    proof->ctr_c = kappa_c - 1;

    poly_qshow_clear(sum_mu_gamma);
    poly_qshow_clear(tmp_poly);
    poly_qshow_clear(e0);
    poly_qshow_clear(e1);
    poly_qshow_clear(t0);
    poly_qshow_vec_256_l_clear(tmp_vec_256_l);
}

/* ═══════════════════════════════════════════════════════
 * ROUND 5: Response + rejection sampling
 * ═══════════════════════════════════════════════════════ */
static int lnp22_prove_round5(
    lnp22_proof_t            *proof,
    const poly_qshow_vec_m1  s1,
    const poly_qshow_vec_m2  s2,
    const poly_qshow_vec_m1  y1,
    const poly_qshow_vec_m2  y2)
{
    uint64_t sq_norm_y1, sq_norm_z1, sq_norm_y2, sq_norm_z2;

    /* z1 = y1 + c·s1 */
    poly_qshow_vec_m1_mul_poly_qshow(proof->z1, s1, proof->c);
    poly_qshow_vec_m1_add(proof->z1, proof->z1, y1);
    sq_norm_y1 = (uint64_t)poly_qshow_vec_m1_norm2(y1);
    sq_norm_z1 = (uint64_t)poly_qshow_vec_m1_norm2(proof->z1);

    /* z2 = y2 + c·s2 */
    poly_qshow_vec_m2_mul_poly_qshow(proof->z2, s2, proof->c);
    poly_qshow_vec_m2_add(proof->z2, proof->z2, y2);
    sq_norm_y2 = poly_qshow_vec_m2_norm2(y2);
    sq_norm_z2 = poly_qshow_vec_m2_norm2(proof->z2);

    /* Rejection sampling */
    if (_lnp22_reject_exp(
            exp(M_PI * (double)(sq_norm_y1 - sq_norm_z1) /
                (double)PARAM_S1SQ_SHOW) /
            (double)LNP22_REJ1) ||
        _lnp22_reject_exp(
            exp(M_PI * (double)(sq_norm_y2 - sq_norm_z2) /
                (double)PARAM_S2SQ_SHOW) /
            (double)LNP22_REJ2)) {
        return 0;
    }
    return 1;
}

/* ═══════════════════════════════════════════════════════
 * lnp22_prove:  Generate the zero-knowledge proof.
 *
 * Parameters:
 *   proof     – output proof structure (initialised by caller)
 *   P_embed   – constraint matrix  (D_SHOW × LNP22_M)
 *   u_embed   – constraint target  (D_SHOW vector)
 *   s1        – witness vector     (M1_SHOW entries)
 *   crs_seed  – CRS randomness seed
 *   seed      – public parameter seed
 * ═══════════════════════════════════════════════════════ */
static void lnp22_prove(
    lnp22_proof_t             *proof,
    const poly_qshow          P_embed[PARAM_D][LNP22_M],
    const poly_qshow          u_embed[PARAM_D],
    const poly_qshow_vec_m1   s1,
    const uint8_t             crs_seed[CRS_SEED_BYTES],
    const uint8_t             seed[SEED_BYTES])
{
    size_t i;
    uint8_t randomness_seed[SEED_BYTES];
    (void)u_embed;  /* u_embed only needed in verify, not prove */
    uint8_t buf[LNP22_CHAL4_BYTES];
    uint32_t kappa;
    poly_qshow tmp_poly, s1i_star, one, quadratic_precomp;
    poly_qshow_mat_d_m1 A1;
    poly_qshow_mat_d_m2 A2;
    poly_qshow_mat_256l_m2 Byg;
    poly_qshow_vec_m1 y1;
    poly_qshow_vec_m1 sum_gamma_r_star[PARAM_L_SHOW];
    poly_qshow_vec_256 sum_gamma_e_star[PARAM_L_SHOW];
    poly_qshow_vec_m2 s2, y2, b;
    poly_qshow_vec_256_l y3_g;

    /* Challenges */
    poly_qshow_vec_m1 chal_1[PARAM_ARP_SHOW];
    coeff_qshow chal_2[PARAM_L_SHOW][LNP22_CHAL2_WIDTH];
    poly_qshow_vec_l chal_3_l;
    poly_qshow_vec_k chal_3_dk[PARAM_D];

    /* ── Init ── */
    poly_qshow_init(tmp_poly);
    poly_qshow_init(s1i_star);
    poly_qshow_init(one);
    poly_qshow_init(quadratic_precomp);
    poly_qshow_zero(quadratic_precomp);

    poly_qshow_mat_d_m1_init(A1);
    poly_qshow_vec_m1_init(y1);
    for (i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_vec_m1_init(sum_gamma_r_star[i]);
        poly_qshow_vec_256_init(sum_gamma_e_star[i]);
    }
    poly_qshow_mat_d_m2_init(A2);
    poly_qshow_mat_256l_m2_init(Byg);
    poly_qshow_vec_m2_init(s2);
    poly_qshow_vec_m2_init(y2);
    poly_qshow_vec_m2_init(b);
    poly_qshow_vec_256_l_init(y3_g);
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        poly_qshow_vec_m1_init(chal_1[i]);
    }
    poly_qshow_vec_l_init(chal_3_l);
    for (i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_init(chal_3_dk[i]);
    }

    /* Random secret seed */
    randombytes(randomness_seed, SEED_BYTES);

    /* Expand CRS */
    poly_qshow_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_LNP22);
    poly_qshow_mat_d_m2_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_LNP22);
    poly_qshow_mat_256l_m2_uniform(Byg, crs_seed,
                                   DOMAIN_SEPARATOR_BYG_LNP22);
    poly_qshow_vec_m2_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_LNP22);

    /* Byte-pack the statement */
    memset(buf, 0, LNP22_CHAL4_BYTES);
    for (i = 0; i < CRS_SEED_BYTES; i++) {
        buf[1 + i] = crs_seed[i];
    }
    for (i = 0; i < SEED_BYTES; i++) {
        buf[1 + CRS_SEED_BYTES + i] = seed[i];
    }

    /* ── Precompute quadratic constraint: ⟨s1*, s1⟩ − B_s² ── */
    for (i = 0; i < PARAM_N_SHOW; i++) {
        poly_qshow_set_coeff(one, i, 1);
    }
    for (i = 0; i < PARAM_M1_SHOW; i++) {
        poly_qshow_conjugate(s1i_star, s1->entries[i]);
        poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
        poly_qshow_add(quadratic_precomp, quadratic_precomp, tmp_poly);
    }
    /* subtract B_s² from constant term */
    poly_qshow_muladd_constant(quadratic_precomp,
                               -(coeff_qshow)LNP22_BSSQ, 1);

    /* ── 5-round protocol with rejection sampling loop ── */
    kappa = 0;
reject:
    lnp22_prove_round1(proof, chal_1, buf, s2, y1, y2, y3_g,
                       s1, A1, A2, Byg, randomness_seed, kappa);
    kappa += PARAM_ARP_DIV_N_L_SHOW - PARAM_ARP_DIV_N_SHOW + 1;

    if (!lnp22_prove_round2(proof, chal_2, buf, s1, chal_1, y3_g)) {
        goto reject;
    }

    lnp22_prove_round3(proof, chal_3_l, chal_3_dk, buf,
                       sum_gamma_e_star, sum_gamma_r_star,
                       s1, chal_1, chal_2, y3_g, quadratic_precomp);

    lnp22_prove_round4(proof, buf, s1, s2, Byg, b,
                       P_embed, sum_gamma_e_star, sum_gamma_r_star,
                       y1, y2, chal_2, chal_3_l, chal_3_dk);

    if (!lnp22_prove_round5(proof, s1, s2, y1, y2)) {
        goto reject;
    }

    /* ── Clean up ── */
    poly_qshow_clear(s1i_star);
    poly_qshow_clear(tmp_poly);
    poly_qshow_clear(one);
    poly_qshow_clear(quadratic_precomp);

    poly_qshow_mat_d_m1_clear(A1);
    poly_qshow_vec_m1_clear(y1);
    for (i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_vec_m1_clear(sum_gamma_r_star[i]);
        poly_qshow_vec_256_clear(sum_gamma_e_star[i]);
    }
    poly_qshow_mat_d_m2_clear(A2);
    poly_qshow_mat_256l_m2_clear(Byg);
    poly_qshow_vec_m2_clear(s2);
    poly_qshow_vec_m2_clear(y2);
    poly_qshow_vec_m2_clear(b);
    poly_qshow_vec_256_l_clear(y3_g);
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        poly_qshow_vec_m1_clear(chal_1[i]);
    }
    poly_qshow_vec_l_clear(chal_3_l);
    for (i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_clear(chal_3_dk[i]);
    }
}

/* ═══════════════════════════════════════════════════════
 * lnp22_verify:  Verify the zero-knowledge proof.
 *
 * Returns 1 if valid, 0 otherwise.
 * ═══════════════════════════════════════════════════════ */
static int lnp22_verify(
    const lnp22_proof_t       *proof,
    const poly_qshow          P_embed[PARAM_D][LNP22_M],
    const poly_qshow          u_embed[PARAM_D],
    const uint8_t             crs_seed[CRS_SEED_BYTES],
    const uint8_t             seed[SEED_BYTES])
{
    size_t i, j;
    int is_valid = 1;
    uint8_t buf[LNP22_CHAL4_BYTES];
    uint8_t challenge_seed[SEED_BYTES];
    uint64_t sq_norm_z3, sq_norm_z2;
    uint128 sq_norm_z1;
    coeff_qshow chal_2[PARAM_L_SHOW][LNP22_CHAL2_WIDTH];
    coeff_qshow tmp_coeff;
    poly_qshow tmp_poly, z1s_z1, z1s_c_one, sum_mu_gamma, t0, c, c_one;
    poly_qshow_mat_d_m1 A1;
    poly_qshow_mat_d_m2 A2;
    poly_qshow_mat_256l_m2 Byg;
    poly_qshow_vec_d tmp_vec_d, w;
    poly_qshow_vec_m1 sum_gamma_r_star_i, tmp_vec_m1;
    poly_qshow sum_gamma_e_star_ij;
    poly_qshow_vec_m1 chal_1[PARAM_ARP_SHOW];
    poly_qshow_vec_m2 b;
    poly_qshow_vec_256_l tmp_vec_256_l, Z;
    poly_qshow_vec_l chal_3_l;
    poly_qshow_vec_k chal_3_dk[PARAM_D];

    /* ── Init ── */
    poly_qshow_init(tmp_poly);
    poly_qshow_init(z1s_z1);
    poly_qshow_init(z1s_c_one);
    poly_qshow_init(c);
    poly_qshow_init(c_one);
    poly_qshow_init(t0);
    poly_qshow_init(sum_mu_gamma);
    poly_qshow_zero(sum_mu_gamma);

    poly_qshow_mat_d_m1_init(A1);
    poly_qshow_mat_d_m2_init(A2);
    poly_qshow_mat_256l_m2_init(Byg);
    poly_qshow_vec_m1_init(tmp_vec_m1);
    poly_qshow_vec_m1_init(sum_gamma_r_star_i);
    poly_qshow_vec_d_init(tmp_vec_d);
    poly_qshow_vec_d_init(w);
    poly_qshow_init(sum_gamma_e_star_ij);
    poly_qshow_vec_m2_init(b);
    poly_qshow_vec_256_l_init(tmp_vec_256_l);
    poly_qshow_vec_256_l_init(Z);
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        poly_qshow_vec_m1_init(chal_1[i]);
    }
    poly_qshow_vec_l_init(chal_3_l);
    for (i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_init(chal_3_dk[i]);
    }

    /* ── Expand CRS ── */
    poly_qshow_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_LNP22);
    poly_qshow_mat_d_m2_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_LNP22);
    poly_qshow_mat_256l_m2_uniform(Byg, crs_seed,
                                   DOMAIN_SEPARATOR_BYG_LNP22);
    poly_qshow_vec_m2_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_LNP22);

    /* Byte-pack the statement */
    memset(buf, 0, LNP22_CHAL4_BYTES);
    for (i = 0; i < CRS_SEED_BYTES; i++) {
        buf[1 + i] = crs_seed[i];
    }
    for (i = 0; i < SEED_BYTES; i++) {
        buf[1 + CRS_SEED_BYTES + i] = seed[i];
    }

    /* c_one = polynomial with all coefficients = 1 */
    for (i = 0; i < PARAM_N_SHOW; i++) {
        poly_qshow_set_coeff(c_one, i, 1);
    }
    poly_qshow_mul(c_one, c_one, proof->c);

    /* ═══════════════ Reconstruct first challenge ═══════════════ */
    /* w = A1·z1 + A2·z2 − c·tA */
    poly_qshow_vec_d_mul_poly_qshow(tmp_vec_d, proof->tA, proof->c);
    poly_qshow_mat_d_m2_mul_vec_m2(w, A2, proof->z2);
    poly_qshow_vec_d_sub(w, w, tmp_vec_d);
    poly_qshow_mat_d_m1_mul_vec_m1(tmp_vec_d, A1, proof->z1);
    poly_qshow_vec_d_add(w, w, tmp_vec_d);

    buf[0] = 1;
    poly_qshow_vec_d_pack(
        buf + LNP22_CHAL_BASE_BYTES, proof->tA);
    poly_qshow_vec_d_pack(
        buf + LNP22_CHAL_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w);
    poly_qshow_vec_256_l_pack(
        buf + LNP22_CHAL_BASE_BYTES + 2 * POLYQSHOW_VECD_PACKEDBYTES,
        proof->tB);
    shake256(challenge_seed, SEED_BYTES, buf, LNP22_CHAL1_BYTES);
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        poly_qshow_vec_m1_binomial(chal_1[i], challenge_seed,
                                   DOMAIN_SEPARATOR_CHAL1_LNP22, i, SEED_BYTES);
    }

    /* ═══════════════ Reconstruct second challenge ═══════════════ */
    buf[0] = 2;
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        coeff_qshow_pack(
            buf + LNP22_CHAL1_BYTES + i * COEFFQSHOW_PACKEDBYTES,
            proof->z3[i]);
    }
    shake256(challenge_seed, SEED_BYTES, buf, LNP22_CHAL2_BYTES);
    for (i = 0; i < PARAM_L_SHOW; i++) {
        vec_qshow_uniform(chal_2[i], challenge_seed,
                          DOMAIN_SEPARATOR_CHAL2_LNP22, i, SEED_BYTES);
    }

    /* ═══════════════ Reconstruct third challenge ═══════════════ */
    buf[0] = 3;
    poly_qshow_vec_l_pack(buf + LNP22_CHAL2_BYTES, proof->h);
    poly_qshow_vec_l_uniform(chal_3_l, buf, DOMAIN_SEPARATOR_CHAL3_LNP22,
                             LNP22_CHAL3_BYTES);
    shake256(challenge_seed, SEED_BYTES, buf, LNP22_CHAL3_BYTES);
    for (i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_uniform(chal_3_dk[i], challenge_seed,
                                 DOMAIN_SEPARATOR_CHAL3_LNP22,
                                 i + PARAM_L_SHOW, SEED_BYTES);
    }

    /* ═══════════════ Reconstruct fourth challenge ═══════════════ */
    /*
     * Recompute t0:
     *
     *   t0 = ⟨b, z2⟩ − c·t1                                    [garbage]
     *      + Σ μ_i·γ_{i,ARP} · (⟨z1*,z1⟩ − c²·B_s²)           [norm bound]
     *      + Σ_i μ_{ℓ+i} · (P·z1 − c·u)_i                      [linear constraint]
     *      + Σ μ_i·(Σ γ_{ij}·z3_j + Σ γ_{ij}·e_j*·Z + Z_{arp+i})
     */

    /* ── ⟨b, z2⟩ − c·t1 ── */
    poly_qshow_vec_m2_mul_inner(t0, b, proof->z2);
    poly_qshow_mul(tmp_poly, proof->c, proof->t1);
    poly_qshow_sub(t0, t0, tmp_poly);

    /* c holds proof->c² for now */
    poly_qshow_mul(c, proof->c, proof->c);

    /* ── Σ μ_i·γ_{i,ARP} ── */
    for (i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i],
                              chal_2[i][PARAM_ARP_SHOW]);
        poly_qshow_add(sum_mu_gamma, sum_mu_gamma, tmp_poly);
    }

    /* ── Z = c·tB − Byg·z2 ── */
    poly_qshow_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, proof->z2);
    poly_qshow_vec_256_l_mul_poly_qshow(Z, proof->tB, proof->c);
    poly_qshow_vec_256_l_sub(Z, Z, tmp_vec_256_l);

    /* ── Quadratic norm term:  Σμγ·(⟨z1*,z1⟩ − c²·B_s²) ── */
    poly_qshow_zero(z1s_z1);
    for (i = 0; i < PARAM_M1_SHOW; i++) {
        poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
        poly_qshow_mul(tmp_poly, tmp_poly, proof->z1->entries[i]);
        poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
    }
    /* subtract c² · B_s² */
    poly_qshow_mul_scalar(tmp_poly, c, LNP22_BSSQ);
    poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
    poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma);
    poly_qshow_add(t0, t0, z1s_z1);

    /* ── Linear constraint:  Σ_i μ_{ℓ+i} · (P·z1 − c·u)_i ── */
    /* Conjugate all chal_1 for the r-star / e-star sums */
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
        poly_qshow_vec_m1_conjugate(tmp_vec_m1, chal_1[j]);
        poly_qshow_vec_m1_set(chal_1[j], tmp_vec_m1);
    }

    poly_qshow_zero(z1s_c_one);
    for (i = 0; i < PARAM_L_SHOW; i++) {
        poly_qshow_set(tmp_poly, proof->h->entries[i]);
        /* h_i must have zero constant coefficient */
        is_valid = is_valid && (poly_qshow_get_coeff(tmp_poly, 0) == 0);
        if (!is_valid) goto lnp22_verify_cleanup;

        for (j = 0; j < PARAM_ARP_SHOW; j++) {
            poly_qshow_muladd_constant(tmp_poly, chal_2[i][j],
                                       proof->z3[j]);
        }
        poly_qshow_mul(tmp_poly, tmp_poly, proof->c);
        poly_qshow_sub(tmp_poly,
                       Z->entries[PARAM_ARP_DIV_N_SHOW + i], tmp_poly);

        /* Σ γ_{ij}·r_j*·z1 */
        for (j = 0; j < PARAM_ARP_SHOW; j++) {
            if (j == 0) {
                poly_qshow_vec_m1_mul_scalar(
                    sum_gamma_r_star_i, chal_1[j], chal_2[i][j]);
            } else {
                poly_qshow_vec_m1_mul_scalar(
                    tmp_vec_m1, chal_1[j], chal_2[i][j]);
                poly_qshow_vec_m1_add(
                    sum_gamma_r_star_i, sum_gamma_r_star_i, tmp_vec_m1);
            }
        }
        poly_qshow_vec_m1_mul_inner(z1s_z1, sum_gamma_r_star_i, proof->z1);
        poly_qshow_add(tmp_poly, tmp_poly, z1s_z1);

        /* Σ γ_{ij}·e_j*·Z */
        for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
            poly_qshow_set_coeff(
                sum_gamma_e_star_ij, 0,
                chal_2[i][j * PARAM_N_SHOW]);
            for (size_t k = 1; k < PARAM_N_SHOW; k++) {
                poly_qshow_set_coeff(
                    sum_gamma_e_star_ij, k,
                    -chal_2[i][(j + 1) * PARAM_N_SHOW - k]);
            }
            poly_qshow_mul(sum_gamma_e_star_ij,
                           sum_gamma_e_star_ij, Z->entries[j]);
            poly_qshow_add(tmp_poly, tmp_poly, sum_gamma_e_star_ij);
        }
        poly_qshow_mul(tmp_poly, tmp_poly, chal_3_l->entries[i]);
        poly_qshow_add(z1s_c_one, z1s_c_one, tmp_poly);
    }
    poly_qshow_mul(z1s_c_one, z1s_c_one, proof->c);
    poly_qshow_add(t0, t0, z1s_c_one);

    /* ── Linear constraint into z1s_c_one, added to t0 ── */
    poly_qshow_zero(z1s_c_one);
    for (i = 0; i < PARAM_D; i++) {
        /* z1s_z1 used as temp: P_i · z1[0..M-1] */
        poly_qshow_zero(z1s_z1);
        for (j = 0; j < LNP22_M; j++) {
            poly_qshow_mul(tmp_poly, P_embed[i][j],
                           proof->z1->entries[j]);
            poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
        }
        /* subtract c · u_i */
        poly_qshow_mul(tmp_poly, proof->c, u_embed[i]);
        poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
        /* weight by μ_{ℓ+i} */
        poly_qshow_mul(z1s_z1, z1s_z1,
                       chal_3_dk[i]->entries[0]);
        poly_qshow_add(z1s_c_one, z1s_c_one, z1s_z1);
    }
    poly_qshow_mul(z1s_c_one, z1s_c_one, proof->c);
    poly_qshow_add(t0, t0, z1s_c_one);

    /* ── Compute fourth challenge c' ── */
    poly_qshow_zero(c);
    buf[0] = 4;
    poly_qshow_pack(buf + LNP22_CHAL3_BYTES, t0);
    poly_qshow_pack(buf + LNP22_CHAL3_BYTES + POLYQSHOW_PACKEDBYTES,
                    proof->t1);
    shake256(challenge_seed, SEED_BYTES, buf, LNP22_CHAL4_BYTES);
    poly_qshow_sample_challenge(tmp_poly, challenge_seed,
                                DOMAIN_SEPARATOR_CHAL4_LNP22,
                                proof->ctr_c, SEED_BYTES);

    /* ═══════════════ Checks ═══════════════ */
    /* Norm bounds */
    sq_norm_z1 = poly_qshow_vec_m1_norm2(proof->z1);
    sq_norm_z2 = poly_qshow_vec_m2_norm2(proof->z2);
    sq_norm_z3 = 0;
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        tmp_coeff = proof->z3[i];
        CHK_UI_OVF_ADDITION(sq_norm_z3, (uint64_t)(tmp_coeff * tmp_coeff));
    }
    is_valid = is_valid &&
               (sq_norm_z1 <= ((uint128)LNP22_BOUND_Z1_LO +
                               (((uint128)LNP22_BOUND_Z1_HI) << 64)));
    is_valid = is_valid && (sq_norm_z2 <= LNP22_BOUND_Z2);
    is_valid = is_valid && (sq_norm_z3 <= LNP22_BOUND_Z3);
    is_valid = is_valid && poly_qshow_equal(proof->c, tmp_poly);

lnp22_verify_cleanup:
    /* Clean up */
    poly_qshow_clear(tmp_poly);
    poly_qshow_clear(z1s_z1);
    poly_qshow_clear(z1s_c_one);
    poly_qshow_clear(c);
    poly_qshow_clear(c_one);
    poly_qshow_clear(t0);
    poly_qshow_clear(sum_mu_gamma);
    poly_qshow_clear(sum_gamma_e_star_ij);

    poly_qshow_mat_d_m1_clear(A1);
    poly_qshow_mat_d_m2_clear(A2);
    poly_qshow_mat_256l_m2_clear(Byg);
    poly_qshow_vec_m1_clear(tmp_vec_m1);
    poly_qshow_vec_m1_clear(sum_gamma_r_star_i);
    poly_qshow_vec_d_clear(tmp_vec_d);
    poly_qshow_vec_d_clear(w);
    poly_qshow_vec_m2_clear(b);
    poly_qshow_vec_256_l_clear(tmp_vec_256_l);
    poly_qshow_vec_256_l_clear(Z);
    for (i = 0; i < PARAM_ARP_SHOW; i++) {
        poly_qshow_vec_m1_clear(chal_1[i]);
    }
    poly_qshow_vec_l_clear(chal_3_l);
    for (i = 0; i < PARAM_D; i++) {
        poly_qshow_vec_k_clear(chal_3_dk[i]);
    }

    return is_valid;
}

/* ═══════════════════════════════════════════════════════
 * Helper: generate a random small-norm witness
 *
 * Fills s[0..M-1] with small-coefficient polynomials such
 * that ||s||² < b_s_sq. Uses centred coefficients in {−1,0,1}.
 * ═══════════════════════════════════════════════════════ */
static void lnp22_random_small_witness(
    poly_qshow s[LNP22_M],
    uint64_t b_s_sq)
{
    size_t i, j;
    int64_t c;
    uint64_t total_norm = 0;

    for (i = 0; i < LNP22_M; i++) {
        poly_qshow_zero(s[i]);
        for (j = 0; j < PARAM_N_SHOW; j++) {
            /* Ternary coefficients */
            uint64_t r;
            randombytes((uint8_t *)&r, sizeof(r));
            r = r % 3;
            c = (r == 0) ? -1 : ((r == 1) ? 0 : 1);
            poly_qshow_set_coeff(s[i], j, c);
            total_norm += (c * c);
        }
    }

    /* If too large, zero out coefficients until within bound */
    i = 0;
    j = 0;
    while (total_norm > b_s_sq) {
        c = poly_qshow_get_coeff(s[i], j);
        if (c != 0) {
            poly_qshow_set_coeff(s[i], j, 0);
            total_norm -= (c * c);
        }
        j++;
        if (j >= PARAM_N_SHOW) { j = 0; i++; }
        if (i >= LNP22_M) break; /* shouldn't happen, LNP22_BSSQ is large enough */
    }
}

/* ═══════════════════════════════════════════════════════
 * Helper: generate a random matrix P and compute u = P·s
 * ═══════════════════════════════════════════════════════ */
static void lnp22_random_statement(
    poly_qshow P[PARAM_D][LNP22_M],
    poly_qshow u[PARAM_D],
    const poly_qshow s[LNP22_M])
{
    size_t i, j;
    uint64_t r;
    coeff_qshow c;
    poly_qshow tmp;

    poly_qshow_init(tmp);

    for (i = 0; i < PARAM_D; i++) {
        /* Random matrix row */
        for (j = 0; j < LNP22_M; j++) {
            poly_qshow_zero(P[i][j]);
            for (size_t k = 0; k < PARAM_N_SHOW; k++) {
                randombytes((uint8_t *)&r, sizeof(r));
                c = (coeff_qshow)(r % 17); /* small-ish coefficients */
                poly_qshow_set_coeff(P[i][j], k, c);
            }
        }

        /* u_i = Σ_j P[i][j] · s[j] */
        poly_qshow_zero(u[i]);
        for (j = 0; j < LNP22_M; j++) {
            poly_qshow_mul(tmp, P[i][j], s[j]);
            poly_qshow_add(u[i], u[i], tmp);
        }
    }

    poly_qshow_clear(tmp);
}

/* ═══════════════════════════════════════════════════════
 * Test: main entry point
 * ═══════════════════════════════════════════════════════ */
static int lnp22_test(void) {
    int rval = 1;
    size_t i, j;
    coeff_qshow coeff;

    /* Statement */
    poly_qshow P[PARAM_D][LNP22_M];
    poly_qshow u[PARAM_D];
    poly_qshow P_embed[PARAM_D][LNP22_M];
    poly_qshow u_embed[PARAM_D];

    /* Witness */
    poly_qshow s[LNP22_M];
    poly_qshow_vec_m1 s1;
    lnp22_proof_t proof;

    uint8_t crs_seed[CRS_SEED_BYTES];
    uint8_t seed[SEED_BYTES];

    printf("\nlnp22_test\n");

    /* Init */
    lnp22_proof_init(&proof);
    poly_qshow_vec_m1_init(s1);
    for (i = 0; i < LNP22_M; i++) {
        poly_qshow_init(s[i]);
    }
    for (i = 0; i < PARAM_D; i++) {
        for (j = 0; j < LNP22_M; j++) {
            poly_qshow_init(P[i][j]);
            poly_qshow_init(P_embed[i][j]);
        }
        poly_qshow_init(u[i]);
        poly_qshow_init(u_embed[i]);
    }

    randombytes(crs_seed, CRS_SEED_BYTES);
    randombytes(seed, SEED_BYTES);

    for (int sub = 0; sub < 5; sub++) {
        /* Step 1: generate random short witness */
        lnp22_random_small_witness(s, LNP22_BSSQ);

        /* Step 2: generate random P and compute u = P·s */
        lnp22_random_statement(P, u, s);

        /* Step 3: embed witness with four-square slack */
        lnp22_embed(s1, P_embed, u_embed, s, P, u, LNP22_BSSQ);

        /* Step 4: generate proof */
        lnp22_prove(&proof, P_embed, u_embed, s1, crs_seed, seed);

        /* Step 5: verify valid proof */
        if (!lnp22_verify(&proof, P_embed, u_embed, crs_seed, seed)) {
            printf("lnp22_verify returned zero for a valid proof.\n");
            rval = 0;
            goto lnp22_test_cleanup;
        }

        /* Step 6: tamper with u → verification should fail */
        coeff = poly_qshow_get_coeff(u_embed[0], 0) + 1;
        poly_qshow_set_coeff(u_embed[0], 0, coeff);
        if (lnp22_verify(&proof, P_embed, u_embed, crs_seed, seed)) {
            printf("lnp22_verify returned non-zero for a wrong u.\n");
            rval = 0;
            goto lnp22_test_cleanup;
        }
        /* restore u */
        poly_qshow_set_coeff(u_embed[0], 0,
                             poly_qshow_get_coeff(u_embed[0], 0) - 1);

        printf(":");
        fflush(stdout);
    }

lnp22_test_cleanup:
    lnp22_proof_clear(&proof);
    poly_qshow_vec_m1_clear(s1);
    for (i = 0; i < LNP22_M; i++) {
        poly_qshow_clear(s[i]);
    }
    for (i = 0; i < PARAM_D; i++) {
        for (j = 0; j < LNP22_M; j++) {
            poly_qshow_clear(P[i][j]);
            poly_qshow_clear(P_embed[i][j]);
        }
        poly_qshow_clear(u[i]);
        poly_qshow_clear(u_embed[i]);
    }
    return rval;
}

/* ═══════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════ */

/* Proof size in bytes (sum of all packed fields of lnp22_proof_t) */
#define LNP22_PROOF_SIZE_BYTES \
    (POLYQSHOW_VECD_PACKEDBYTES       /* tA: 23*512 = 11776 */ \
   + POLYQSHOW_VEC256L_PACKEDBYTES    /* tB: 11*512 = 5632  */ \
   + PARAM_ARP_SHOW * COEFFQSHOW_PACKEDBYTES  /* z3: 256*8 = 2048 */ \
   + POLYQSHOW_VECL_PACKEDBYTES       /* h:  7*512 = 3584  */ \
   + POLYQSHOW_PACKEDBYTES            /* t1: 64*8 = 512    */ \
   + POLYQSHOW_PACKEDBYTES            /* c:  64*8 = 512    */ \
   + sizeof(uint32_t)                 /* ctr_c: 4          */ \
   + PARAM_M1_SHOW * POLYQSHOW_PACKEDBYTES  /* z1: 211*512 = 108032 */ \
   + PARAM_M2_SHOW * POLYQSHOW_PACKEDBYTES) /* z2: 74*512 = 37888  */

/* Witness s1 size in bytes (packed) */
#define LNP22_WITNESS_SIZE_BYTES \
    (PARAM_M1_SHOW * POLYQSHOW_PACKEDBYTES)  /* 211*512 = 108032 */

static double _time_diff_sec(struct timespec *start, struct timespec *end) {
    return (end->tv_sec - start->tv_sec)
         + (end->tv_nsec - start->tv_nsec) / 1e9;
}

int main(void) {
    struct timespec t0, t1, t2;
    double prove_time, verify_time;

    arith_setup();
    random_init();
    printf("lnp22 Hello from the LNP22 unit test.\n");

    /* ── Benchmark a single run ── */
    {
        size_t i, j;
        poly_qshow P[PARAM_D][LNP22_M], u[PARAM_D];
        poly_qshow P_embed[PARAM_D][LNP22_M], u_embed[PARAM_D];
        poly_qshow s[LNP22_M];
        poly_qshow_vec_m1 s1;
        lnp22_proof_t proof;
        uint8_t crs_seed[CRS_SEED_BYTES], seed[SEED_BYTES];

        lnp22_proof_init(&proof);
        poly_qshow_vec_m1_init(s1);
        for (i = 0; i < LNP22_M; i++) poly_qshow_init(s[i]);
        for (i = 0; i < PARAM_D; i++) {
            for (j = 0; j < LNP22_M; j++) {
                poly_qshow_init(P[i][j]);
                poly_qshow_init(P_embed[i][j]);
            }
            poly_qshow_init(u[i]);
            poly_qshow_init(u_embed[i]);
        }
        randombytes(crs_seed, CRS_SEED_BYTES);
        randombytes(seed, SEED_BYTES);

        /* Generate one random instance */
        lnp22_random_small_witness(s, LNP22_BSSQ);
        lnp22_random_statement(P, u, s);
        lnp22_embed(s1, P_embed, u_embed, s, P, u, LNP22_BSSQ);

        /* Time prove */
        clock_gettime(CLOCK_MONOTONIC, &t0);
        lnp22_prove(&proof, P_embed, u_embed, s1, crs_seed, seed);//调用零知识证明生成函数
        clock_gettime(CLOCK_MONOTONIC, &t1);

        /* Time verify */
        lnp22_verify(&proof, P_embed, u_embed, crs_seed, seed);//调用零知识证明验证函数
        clock_gettime(CLOCK_MONOTONIC, &t2);

        prove_time  = _time_diff_sec(&t0, &t1);
        verify_time = _time_diff_sec(&t1, &t2);

        printf("\n");
        printf("  Proof size : %.2f KB (%zu bytes)\n",
               LNP22_PROOF_SIZE_BYTES / 1024.0, (size_t)LNP22_PROOF_SIZE_BYTES);
        printf("  Witness s1 : %.2f KB (%zu bytes packed)\n",
               LNP22_WITNESS_SIZE_BYTES / 1024.0, (size_t)LNP22_WITNESS_SIZE_BYTES);
        printf("  Prove time : %.3f ms\n", prove_time * 1000.0);
        printf("  Verify time: %.3f ms\n", verify_time * 1000.0);

        /* Clean up benchmark */
        lnp22_proof_clear(&proof);
        poly_qshow_vec_m1_clear(s1);
        for (i = 0; i < LNP22_M; i++) poly_qshow_clear(s[i]);
        for (i = 0; i < PARAM_D; i++) {
            for (j = 0; j < LNP22_M; j++) {
                poly_qshow_clear(P[i][j]);
                poly_qshow_clear(P_embed[i][j]);
            }
            poly_qshow_clear(u[i]);
            poly_qshow_clear(u_embed[i]);
        }
    }

    /* ── Run correctness tests ── */
    int pass = lnp22_test();

    if (pass) {
        printf("\nlnp22: passed.\n");
    } else {
        printf("\nlnp22: FAILED!\n");
    }

    arith_teardown();
    return pass ? 0 : 1;
}
