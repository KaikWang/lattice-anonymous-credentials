#include "acc_genpp.h"

#include <string.h>

#include "fips202.h"
#include "poly_q_sampling.h"
#include "randombytes.h"
#include "random.h"

#define ACC_TRUNC_S2SQ 201183.59824153673253022134
#define ACC_TRUNC_SQRT_S4SQ_SGSQ 63.52563493692179008576
#define ACC_TRUNC_SGINVSQ_S4INVSQ -2991.20481842988692733343
#define ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ -0.49560117302052814070
#define ACC_TRUNC_NEG_SGSQ -2000.00165310704187504598

typedef struct {
  poly_q_mat_d_d blocks[ACC_GENPP_GADGET_DIM][ACC_GENPP_GADGET_DIM];
} acc_genpp_ginv_mat_t;

static uint64_t acc_genpp_gadget_power(size_t j) {
  return ((uint64_t)1) << j;
}

static int acc_genpp_handle_bit(uint32_t x, size_t pos) {
  return (int)((x >> pos) & 1u);
}

static void acc_genpp_gadget_block(poly_q_mat_d_d out, size_t j) {
  poly_q_mat_d_d_zero(out);
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_q_set_coeff(out->rows[i]->entries[i], 0, (coeff_q)acc_genpp_gadget_power(j));
  }
}

static void acc_genpp_Di(acc_genpp_gmat_t *out, const acc_genpp_Ba_t *B, size_t i, int bit) {
  poly_q_mat_d_d tmp;
  poly_q_mat_d_d_init(tmp);

  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    if (bit) {
      poly_q_mat_d_d_set(out->blocks[j], B->bit[i].blocks[j]);
    } else {
      acc_genpp_gadget_block(tmp, j);
      poly_q_mat_d_d_sub(out->blocks[j], tmp, B->bit[i].blocks[j]);
    }
  }

  poly_q_mat_d_d_clear(tmp);
}

static void acc_genpp_ginv_mat_init(acc_genpp_ginv_mat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_init(m->blocks[i][j]);
    }
  }
}

static void acc_genpp_ginv_mat_clear(acc_genpp_ginv_mat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_clear(m->blocks[i][j]);
    }
  }
}

static void acc_genpp_ginv_vec(poly_q_vec_d out[ACC_GENPP_GADGET_DIM], const poly_q_vec_d in) {
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    poly_q_vec_d_zero(out[j]);
  }

  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      coeff_q c = poly_q_get_coeff(in->entries[row], coeff);
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        if ((c >> j) & 1) {
          poly_q_set_coeff(out[j]->entries[row], coeff, 1);
        }
      }
    }
  }
}

static void acc_genpp_ginv_mat(acc_genpp_ginv_mat_t *out, const acc_genpp_gmat_t *in) {
  for (size_t r = 0; r < ACC_GENPP_GADGET_DIM; r++) {
    for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
      poly_q_mat_d_d_zero(out->blocks[r][c]);
    }
  }

  for (size_t cblock = 0; cblock < ACC_GENPP_GADGET_DIM; cblock++) {
    for (size_t r = 0; r < PARAM_D; r++) {
      for (size_t c = 0; c < PARAM_D; c++) {
        for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
          coeff_q v = poly_q_get_coeff(in->blocks[cblock]->rows[r]->entries[c], coeff);
          for (size_t bit = 0; bit < ACC_GENPP_GADGET_DIM; bit++) {
            if ((v >> bit) & 1) {
              poly_q_set_coeff(out->blocks[bit][cblock]->rows[r]->entries[c], coeff, 1);
            }
          }
        }
      }
    }
  }
}

static void acc_genpp_gmat_mul_ginv(acc_genpp_gmat_t *out,
                                    const acc_genpp_gmat_t *lhs,
                                    const acc_genpp_ginv_mat_t *rhs) {
  poly_q_mat_d_d acc, tmp;
  poly_q_mat_d_d_init(acc);
  poly_q_mat_d_d_init(tmp);

  for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
    poly_q_mat_d_d_zero(acc);
    for (size_t k = 0; k < ACC_GENPP_GADGET_DIM; k++) {
      poly_q_mat_d_d_mul_mat_d_d(tmp, lhs->blocks[k], rhs->blocks[k][c]);
      poly_q_mat_d_d_add(acc, acc, tmp);
    }
    poly_q_mat_d_d_set(out->blocks[c], acc);
  }

  poly_q_mat_d_d_clear(tmp);
  poly_q_mat_d_d_clear(acc);
}

static void acc_genpp_gmat_mul_gvec(poly_q_vec_d out,
                                    const acc_genpp_gmat_t *lhs,
                                    const poly_q_vec_d rhs[ACC_GENPP_GADGET_DIM]) {
  poly_q_vec_d_zero(out);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_mat_d_d_muladd_vec_d(out, lhs->blocks[i], rhs[i]);
  }
}

static void acc_genpp_sample_small_vec_d(poly_q_vec_d out,
                                         const uint8_t seed[ACC_GENPP_SEED_BYTES],
                                         uint32_t x,
                                         size_t bit,
                                         size_t gadget) {
  keccak_state state;
  uint8_t buf[SHAKE256_RATE];
  size_t pos = sizeof(buf);

  poly_q_vec_d_zero(out);
  shake256_init(&state);
  shake256_absorb(&state, seed, ACC_GENPP_SEED_BYTES);
  shake256_absorb(&state, (const uint8_t *)&x, sizeof(x));
  shake256_absorb(&state, (const uint8_t *)&bit, sizeof(bit));
  shake256_absorb(&state, (const uint8_t *)&gadget, sizeof(gadget));
  shake256_finalize(&state);

  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      uint8_t sample;
      do {
        if (pos >= sizeof(buf)) {
          shake256_squeeze(buf, sizeof(buf), &state);
          pos = 0;
        }
        sample = buf[pos++] & 3u;
      } while (sample == 3u);
      poly_q_set_coeff(out->entries[i], j, (coeff_q)sample - 1);
    }
  }
}

static void acc_genpp_mul_Mx_right(poly_q_vec_d out,
                                   const acc_genpp_pp_t *pp,
                                   uint32_t x,
                                   const poly_q_vec_d right[ACC_GENPP_HANDLE_BITS][ACC_GENPP_GADGET_DIM]) {
  poly_q_mat_d_d Gblock, Mblock;

  poly_q_vec_d_zero(out);
  poly_q_mat_d_d_init(Gblock);
  poly_q_mat_d_d_init(Mblock);

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    int xbit = acc_genpp_handle_bit(x, i);
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_set(Mblock, pp->Ba.bit[i].blocks[j]);
      if (xbit) {
        acc_genpp_gadget_block(Gblock, j);
        poly_q_mat_d_d_sub(Mblock, Mblock, Gblock);
      }
      poly_q_mat_d_d_muladd_vec_d(out, Mblock, right[i][j]);
    }
  }

  poly_q_mat_d_d_clear(Mblock);
  poly_q_mat_d_d_clear(Gblock);
}

static void acc_trunc_vec_gaussian_sqrt_s4sq_sgsq(poly_q_vec_d out) {
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      poly_q_set_coeff(out->entries[i], j, SampleZ(0, ACC_TRUNC_SQRT_S4SQ_SGSQ));
    }
  }
}

static void acc_trunc_compute_rrstar(poly_q_mat_d_d RRstar[2][2],
                                     const poly_q_mat_d_d R[2][ACC_GENPP_TRUNC_K_L]) {
  poly_q_mat_d_d Rstar[ACC_GENPP_TRUNC_K_L][2], tmp;

  poly_q_mat_d_d_init(tmp);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(Rstar[i][j]);
      poly_q_mat_d_d_conjugate(Rstar[i][j], R[j][i]);
    }
  }

  for (size_t i = 0; i < 2; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_mul_mat_d_d(RRstar[i][j], R[i][0], Rstar[0][j]);
      for (size_t k = 1; k < ACC_GENPP_TRUNC_K_L; k++) {
        poly_q_mat_d_d_mul_mat_d_d(tmp, R[i][k], Rstar[k][j]);
        poly_q_mat_d_d_add(RRstar[i][j], RRstar[i][j], tmp);
      }
    }
  }

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(Rstar[i][j]);
    }
  }
  poly_q_mat_d_d_clear(tmp);
}

static void acc_trunc_compute_covariance(poly_real_mat_2d_2d S,
                                         const poly_q_mat_d_d RRstar[2][2]) {
  poly_real finv, tmp;

  poly_real_init(finv);
  poly_real_init(tmp);

  for (size_t i = 0; i < 2 * PARAM_D; i++) {
    for (size_t j = 0; j < 2 * PARAM_D; j++) {
      poly_real_from_poly_q(S->rows[i]->entries[j],
          RRstar[i / PARAM_D][j / PARAM_D]->rows[i % PARAM_D]->entries[j % PARAM_D]);
      poly_real_mul_scalar(S->rows[i]->entries[j], S->rows[i]->entries[j], ACC_TRUNC_SGINVSQ_S4INVSQ);
    }
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_real_add_constant(S->rows[i]->entries[i], PARAM_S1SQ);
  }
  for (size_t i = PARAM_D; i < 2 * PARAM_D; i++) {
    poly_real_add_constant(S->rows[i]->entries[i], PARAM_S3SQ);
  }

  for (size_t i = 2 * PARAM_D - 1; i > PARAM_D - 1; i--) {
    poly_real_invert(finv, S->rows[i]->entries[i]);
    for (size_t j = 0; j < i; j++) {
      poly_real_mul(S->rows[j]->entries[i], S->rows[j]->entries[i], finv);
    }
    for (size_t j = 0; j < i; j++) {
      for (size_t k = 0; k < i; k++) {
        poly_real_mul(tmp, S->rows[j]->entries[i], S->rows[i]->entries[k]);
        poly_real_sub(S->rows[j]->entries[k], S->rows[j]->entries[k], tmp);
      }
    }
  }

  poly_real_clear(tmp);
  poly_real_clear(finv);
}

static void acc_trunc_sample_perturb(
    poly_q_vec_d pL[ACC_GENPP_TRUNC_L],
    poly_q_vec_d p12,
    poly_q_vec_d p2[ACC_GENPP_TRUNC_K_L],
    const poly_real m_sG_t_tstar,
    const poly_q_mat_d_d R[2][ACC_GENPP_TRUNC_K_L],
    const poly_real_mat_2d_2d S) {
  poly_real tmp_sub, tmp_mul, tmp_sample;
  poly_real_vec_2d c;
  poly_real S_d[PARAM_D][PARAM_D];

  poly_real_init(tmp_sub);
  poly_real_init(tmp_mul);
  poly_real_init(tmp_sample);
  poly_real_vec_2d_init(c);
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_real_init(S_d[i][j]);
    }
  }

  arb_poly_set(tmp_sub, m_sG_t_tstar);
  poly_real_add_constant(tmp_sub, (coeff_real)ACC_TRUNC_S2SQ);
  arb_poly_zero(tmp_mul);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_L; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_q_samplefz(pL[i]->entries[j], tmp_sub, tmp_mul);
    }
  }

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    acc_trunc_vec_gaussian_sqrt_s4sq_sgsq(p2[i]);
  }

  poly_q_mat_d_d_mul_vec_d(pL[0], R[0][0], p2[0]);
  poly_q_mat_d_d_mul_vec_d(p12, R[1][0], p2[0]);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(pL[0], R[0][i], p2[i]);
    poly_q_mat_d_d_muladd_vec_d(p12, R[1][i], p2[i]);
  }

  for (size_t i = 0; i < PARAM_D; i++) {
    poly_real_from_poly_q(c->entries[i], pL[0]->entries[i]);
    poly_real_from_poly_q(c->entries[i + PARAM_D], p12->entries[i]);
    poly_real_mul_scalar(c->entries[i], c->entries[i], ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ);
    poly_real_mul_scalar(c->entries[i + PARAM_D], c->entries[i + PARAM_D], ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ);
  }

  for (size_t i = 2 * PARAM_D; i > PARAM_D; i--) {
    poly_q_samplefz(p12->entries[i - 1 - PARAM_D], S->rows[i - 1]->entries[i - 1], c->entries[i - 1]);
    poly_real_from_poly_q(tmp_sample, p12->entries[i - 1 - PARAM_D]);
    poly_real_sub(tmp_sub, tmp_sample, c->entries[i - 1]);
    for (size_t j = 0; j < i - 1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S->rows[j]->entries[i - 1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      arb_poly_set(S_d[i][j], S->rows[i]->entries[j]);
    }
    poly_real_add(S_d[i][i], S_d[i][i], m_sG_t_tstar);
  }

  for (size_t i = PARAM_D; i > 0; i--) {
    poly_q_samplefz(pL[0]->entries[i - 1], S_d[i - 1][i - 1], c->entries[i - 1]);
    poly_real_invert(tmp_sub, S_d[i - 1][i - 1]);
    for (size_t j = 0; j < i - 1; j++) {
      poly_real_mul(S_d[j][i - 1], S_d[j][i - 1], tmp_sub);
    }
    for (size_t j = 0; j < i - 1; j++) {
      for (size_t k = 0; k < i - 1; k++) {
        poly_real_mul(tmp_mul, S_d[j][i - 1], S_d[i - 1][k]);
        poly_real_sub(S_d[j][k], S_d[j][k], tmp_mul);
      }
    }
    poly_real_from_poly_q(tmp_sample, pL[0]->entries[i - 1]);
    poly_real_sub(tmp_sub, tmp_sample, c->entries[i - 1]);
    for (size_t j = 0; j < i - 1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S_d[j][i - 1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_real_clear(S_d[i][j]);
    }
  }
  poly_real_vec_2d_clear(c);
  poly_real_clear(tmp_sample);
  poly_real_clear(tmp_mul);
  poly_real_clear(tmp_sub);
}

static const double acc_trunc_klein_widths[ACC_GENPP_TRUNC_K] = {
  3.42997312037747725810,
  3.44004613446786633446,
  3.44010564823841225035,
  3.44010600038776326315,
  3.45610035816311800261
};

static const double acc_trunc_neg_scaled_gadget_gso[ACC_GENPP_TRUNC_K][ACC_GENPP_TRUNC_K] = {
  {-0.07647058823529412352, -0.00588214820229020913, -0.00045247284545076348, -0.00003480560345379704, -0.00000269779563119137},
  { 0.00588235294117647051, -0.07646792662977272559, -0.00588214699085992491, -0.00045247284489936151, -0.00003507134320546393},
  { 0.00000000000000000000,  0.00591695381295464801, -0.07646791088117903257, -0.00588214698369169970, -0.00045592746167107091},
  { 0.00000000000000000000,  0.00000000000000000000,  0.00591715854467268331, -0.07646791078799210217, -0.00592705700172388959},
  { 0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000,  0.00591715975610271599, -0.07705174102241059420}
};

static const int64_t acc_trunc_gadget_basis[ACC_GENPP_TRUNC_K][ACC_GENPP_TRUNC_K] = {
  {14,  0,  0,  0,  7},
  {-1, 14,  0,  0,  3},
  { 0, -1, 14,  0, 12},
  { 0,  0, -1, 14,  7},
  { 0,  0,  0, -1, 10}
};

static void acc_trunc_sample_klein(poly_q_vec_d v[ACC_GENPP_TRUNC_K], const poly_q_vec_d w) {
  for (int64_t i = PARAM_N * PARAM_D * ACC_GENPP_TRUNC_K; i > 0; i--) {
    int64_t i_1 = (i - 1) / (PARAM_N * PARAM_D);
    int64_t i_2 = (i - 1) % (PARAM_N * PARAM_D);
    int64_t i_2_1 = i_2 / PARAM_N;
    int64_t i_2_2 = i_2 % PARAM_N;
    coeff_q w_coeff = poly_q_get_coeff(w->entries[i_2_1], (size_t)i_2_2);
    double di;

    if (i == PARAM_N * PARAM_D * ACC_GENPP_TRUNC_K) {
      di = acc_trunc_neg_scaled_gadget_gso[0][i_1] * ((double)w_coeff);
    } else {
      coeff_q v_coeff = poly_q_get_coeff_centered(v[0]->entries[i_2_1], (size_t)i_2_2);
      di = acc_trunc_neg_scaled_gadget_gso[0][i_1] * ((double)(w_coeff + v_coeff));
      for (size_t j = 1; j < ACC_GENPP_TRUNC_K; j++) {
        v_coeff = poly_q_get_coeff_centered(v[j]->entries[i_2_1], (size_t)i_2_2);
        di += acc_trunc_neg_scaled_gadget_gso[j][i_1] * ((double)v_coeff);
      }
    }

    int64_t zi = SampleZ(di, acc_trunc_klein_widths[i_1]);
    for (size_t j = 0; j < ACC_GENPP_TRUNC_K; j++) {
      coeff_q v_coeff = poly_q_get_coeff(v[j]->entries[i_2_1], (size_t)i_2_2);
      v_coeff += zi * acc_trunc_gadget_basis[j][i_1];
      poly_q_set_coeff(v[j]->entries[i_2_1], (size_t)i_2_2, v_coeff);
    }
  }
}

static void acc_trunc_sampler(
    poly_q_vec_d v11,
    poly_q_vec_d v12,
    poly_q_vec_d v2[ACC_GENPP_TRUNC_K_L],
    const poly_q_mat_d_d R[2][ACC_GENPP_TRUNC_K_L],
    const poly_q_mat_d_d A,
    const poly_q_mat_d_d B[ACC_GENPP_TRUNC_K_L],
    const poly_q_vec_d u,
    const poly_q tag,
    const poly_real_mat_2d_2d S) {
  poly_q_vec_d pL[ACC_GENPP_TRUNC_L], p12, p2[ACC_GENPP_TRUNC_K_L];
  poly_q_vec_d w, tmp, y[ACC_GENPP_TRUNC_K];
  poly_q taginv;
  poly_real m_sG_t_tstar;
  int64_t bexpi;

  poly_q_init(taginv);
  poly_real_init(m_sG_t_tstar);
  poly_q_vec_d_init(w);
  poly_q_vec_d_init(tmp);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_init(pL[i]);
    poly_q_vec_d_init(y[i]);
  }
  poly_q_vec_d_init(p12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_init(p2[i]);
    poly_q_vec_d_init(y[ACC_GENPP_TRUNC_L + i]);
  }

  poly_q_conjugate(taginv, tag);
  poly_q_mul(taginv, taginv, tag);
  poly_real_from_poly_q(m_sG_t_tstar, taginv);
  poly_real_mul_scalar(m_sG_t_tstar, m_sG_t_tstar, (coeff_real)ACC_TRUNC_NEG_SGSQ);

  acc_trunc_sample_perturb(pL, p12, p2, m_sG_t_tstar, R, S);

  bexpi = ACC_GENPP_TRUNC_B;
  for (size_t i = 1; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_mul_scalar(tmp, pL[i], bexpi);
    poly_q_vec_d_add(pL[0], pL[0], tmp);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  poly_q_invert(taginv, tag);
  poly_q_vec_d_sub(w, u, pL[0]);
  poly_q_mat_d_d_mulsub_vec_d(w, A, p12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(w, B[i], p2[i]);
  }
  poly_q_vec_d_mul_poly(w, w, taginv);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_mul_scalar(tmp, p2[i], bexpi);
    poly_q_vec_d_sub(w, w, tmp);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  acc_trunc_sample_klein(y, w);
  poly_q_vec_d_add(y[0], y[0], w);

  poly_q_vec_d_mul_poly(y[0], y[0], tag);
  poly_q_vec_d_add(pL[0], pL[0], y[0]);
  bexpi = ACC_GENPP_TRUNC_B;
  for (size_t i = 1; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_mul_poly(y[i], y[i], tag);
    poly_q_vec_d_mul_scalar(y[i], y[i], bexpi);
    poly_q_vec_d_add(pL[0], pL[0], y[i]);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  poly_q_mat_d_d_mul_vec_d(v11, R[0][0], y[ACC_GENPP_TRUNC_L]);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(v11, R[0][i], y[ACC_GENPP_TRUNC_L + i]);
  }
  poly_q_vec_d_add(v11, v11, pL[0]);

  poly_q_mat_d_d_mul_vec_d(v12, R[1][0], y[ACC_GENPP_TRUNC_L]);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(v12, R[1][i], y[ACC_GENPP_TRUNC_L + i]);
  }
  poly_q_vec_d_add(v12, v12, p12);

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_add(v2[i], p2[i], y[ACC_GENPP_TRUNC_L + i]);
  }

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_clear(p2[i]);
    poly_q_vec_d_clear(y[ACC_GENPP_TRUNC_L + i]);
  }
  poly_q_vec_d_clear(p12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_clear(pL[i]);
    poly_q_vec_d_clear(y[i]);
  }
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(w);
  poly_real_clear(m_sG_t_tstar);
  poly_q_clear(taginv);
}

static void acc_genpp_sample_Ba(acc_genpp_Ba_t *B, const uint8_t seed[SEED_BYTES]) {
  uint8_t block_seed[SEED_BYTES];
  keccak_state state;

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      shake256_init(&state);
      shake256_absorb(&state, seed, SEED_BYTES);
      shake256_absorb(&state, (const uint8_t *)&i, sizeof(i));
      shake256_absorb(&state, (const uint8_t *)&j, sizeof(j));
      shake256_finalize(&state);
      shake256_squeeze(block_seed, SEED_BYTES, &state);
      poly_q_mat_d_d_uniform(B->bit[i].blocks[j], block_seed, DOMAIN_SEPARATOR_DS, 0);
    }
  }
}

static int acc_genpp_dict_find(const acc_genpp_state_t *st, uint32_t x) {
  for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
    if (st->dict[i].used && st->dict[i].handle == x) {
      return (int)i;
    }
  }
  return -1;
}

static int acc_genpp_dict_put(acc_genpp_state_t *st, uint32_t x, acc_genpp_op_t op, uint32_t index) {
  int pos = acc_genpp_dict_find(st, x);
  if (pos < 0) {
    for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
      if (!st->dict[i].used) {
        pos = (int)i;
        st->dict[i].used = 1;
        st->dict[i].handle = x;
        break;
      }
    }
  }
  if (pos < 0) {
    return 0;
  }
  st->dict[pos].op = op;
  st->dict[pos].index = index;
  return 1;
}

static void acc_genpp_bml_push(acc_genpp_state_t *st, const acc_genpp_update_t *upmsg) {
  size_t pos;

  if (st->bml_len < ACC_GENPP_BML_CAPACITY) {
    pos = (st->bml_start + st->bml_len) % ACC_GENPP_BML_CAPACITY;
    st->bml_len++;
  } else {
    pos = st->bml_start;
    st->bml_start = (st->bml_start + 1) % ACC_GENPP_BML_CAPACITY;
  }
  st->bml[pos] = *upmsg;
}

void acc_genpp_gmat_init(acc_genpp_gmat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_mat_d_d_init(m->blocks[i]);
  }
}

void acc_genpp_gmat_clear(acc_genpp_gmat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_mat_d_d_clear(m->blocks[i]);
  }
}

void acc_genpp_Ba_init(acc_genpp_Ba_t *B) {
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    acc_genpp_gmat_init(&B->bit[i]);
  }
}

void acc_genpp_Ba_clear(acc_genpp_Ba_t *B) {
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    acc_genpp_gmat_clear(&B->bit[i]);
  }
}

void acc_genpp_state_init(acc_genpp_state_t *st) {
  memset(st->ka, 0, sizeof(st->ka));
  st->ctr = 0;
  poly_q_vec_d_init(st->acc0);
  memset(st->dict, 0, sizeof(st->dict));
  memset(st->bml, 0, sizeof(st->bml));
  st->bml_start = 0;
  st->bml_len = 0;
}

void acc_genpp_state_clear(acc_genpp_state_t *st) {
  poly_q_vec_d_clear(st->acc0);
}

void acc_genpp_pp_init(acc_genpp_pp_t *pp) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_mat_d_d_init(pp->A.B[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_init(pp->trunc_B[i]);
  }
  memset(pp->A.seed, 0, sizeof(pp->A.seed));
  acc_genpp_Ba_init(&pp->Ba);
  poly_q_vec_d_init(pp->ua);
}

void acc_genpp_pp_clear(acc_genpp_pp_t *pp) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_mat_d_d_clear(pp->A.B[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_clear(pp->trunc_B[i]);
  }
  acc_genpp_Ba_clear(&pp->Ba);
  poly_q_vec_d_clear(pp->ua);
}

void acc_genpp_td_init(acc_genpp_td_t *td) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(td->RA.R[j][i]);
    }
  }
  poly_real_mat_2d_2d_init(td->RA.S);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(td->trunc_R[j][i]);
    }
  }
  poly_real_mat_2d_2d_init(td->trunc_S);
}

void acc_genpp_td_clear(acc_genpp_td_t *td) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(td->RA.R[j][i]);
    }
  }
  poly_real_mat_2d_2d_clear(td->RA.S);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(td->trunc_R[j][i]);
    }
  }
  poly_real_mat_2d_2d_clear(td->trunc_S);
}

void acc_genpp_wit_init(acc_genpp_wit_t *wit) {
  poly_q_vec_d_init(wit->left);
  poly_q_vec_d_init(wit->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_init(wit->v2[i]);
  }
  poly_q_init(wit->tag);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_init(wit->right[i][j]);
    }
  }
}

void acc_genpp_wit_clear(acc_genpp_wit_t *wit) {
  poly_q_vec_d_clear(wit->left);
  poly_q_vec_d_clear(wit->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_clear(wit->v2[i]);
  }
  poly_q_clear(wit->tag);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_clear(wit->right[i][j]);
    }
  }
}

size_t acc_genpp_acc_size_bytes(void) {
  return (PARAM_D * PARAM_N * PARAM_Q_BITLEN + 7) / 8;
}

size_t acc_genpp_wit_compact_size_bytes(void) {
  const size_t vec_bits = PARAM_D * PARAM_N * PARAM_Q_BITLEN;
  const size_t tag_bits = PARAM_N;
  const size_t vec_count = 2 + ACC_GENPP_TRUNC_K_L;
  return (vec_count * vec_bits + tag_bits + 7) / 8;
}

size_t acc_genpp_wit_full_size_bytes(void) {
  const size_t vec_bits = PARAM_D * PARAM_N * PARAM_Q_BITLEN;
  const size_t tag_bits = PARAM_N;
  const size_t vec_count = 2 + ACC_GENPP_TRUNC_K_L + ACC_GENPP_HANDLE_BITS * ACC_GENPP_GADGET_DIM;
  return (vec_count * vec_bits + tag_bits + 7) / 8;
}

uint64_t acc_genpp_wit_norm2(const acc_genpp_wit_t *wit) {
  uint64_t norm = 0;

  norm += poly_q_vec_d_norm2(wit->left);
  norm += poly_q_vec_d_norm2(wit->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    norm += poly_q_vec_d_norm2(wit->v2[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      norm += poly_q_vec_d_norm2(wit->right[i][j]);
    }
  }
  return norm;
}

void acc_genpp_Fk(uint8_t r[ACC_GENPP_SEED_BYTES],
                  const uint8_t ka[ACC_GENPP_SEED_BYTES],
                  uint32_t x,
                  const poly_q_vec_d acc) {
  keccak_state state;
  coeff_q coeff;

  shake256_init(&state);
  shake256_absorb(&state, ka, ACC_GENPP_SEED_BYTES);
  shake256_absorb(&state, (const uint8_t *)&x, sizeof(x));
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      coeff = poly_q_get_coeff(acc->entries[i], j);
      shake256_absorb(&state, (const uint8_t *)&coeff, sizeof(coeff));
    }
  }
  shake256_finalize(&state);
  shake256_squeeze(r, ACC_GENPP_SEED_BYTES, &state);
}

void acc_genpp_evalf(acc_genpp_gmat_t *out,
                     const acc_genpp_Ba_t *B,
                     uint32_t selector) {
  acc_genpp_gmat_t cur, D;
  acc_genpp_ginv_mat_t inv;

  acc_genpp_gmat_init(&cur);
  acc_genpp_gmat_init(&D);
  acc_genpp_ginv_mat_init(&inv);

  acc_genpp_Di(&cur, B, ACC_GENPP_HANDLE_BITS - 1,
               acc_genpp_handle_bit(selector, ACC_GENPP_HANDLE_BITS - 1));

  for (size_t off = 1; off < ACC_GENPP_HANDLE_BITS; off++) {
    size_t i = ACC_GENPP_HANDLE_BITS - 1 - off;
    acc_genpp_ginv_mat(&inv, &cur);
    acc_genpp_Di(&D, B, i, acc_genpp_handle_bit(selector, i));
    acc_genpp_gmat_mul_ginv(&cur, &D, &inv);
  }

  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    poly_q_mat_d_d_set(out->blocks[j], cur.blocks[j]);
  }

  acc_genpp_ginv_mat_clear(&inv);
  acc_genpp_gmat_clear(&D);
  acc_genpp_gmat_clear(&cur);
}

void acc_genpp_evalfx_apply(acc_genpp_wit_t *delta,
                            const acc_genpp_Ba_t *B,
                            uint32_t selector,
                            uint32_t point,
                            const poly_q_vec_d u) {
  acc_genpp_gmat_t tail[ACC_GENPP_HANDLE_BITS + 1], D;
  poly_q_vec_d ginv_u[ACC_GENPP_GADGET_DIM];
  poly_q_vec_d tmp_in[ACC_GENPP_GADGET_DIM], tmp_out[ACC_GENPP_GADGET_DIM];
  acc_genpp_ginv_mat_t inv;
  int prefix_equal = 1;

  for (size_t i = 0; i <= ACC_GENPP_HANDLE_BITS; i++) {
    acc_genpp_gmat_init(&tail[i]);
  }
  acc_genpp_gmat_init(&D);
  acc_genpp_ginv_mat_init(&inv);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_init(ginv_u[i]);
    poly_q_vec_d_init(tmp_in[i]);
    poly_q_vec_d_init(tmp_out[i]);
  }

  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    acc_genpp_gadget_block(tail[ACC_GENPP_HANDLE_BITS].blocks[j], j);
  }

  for (size_t off = 0; off < ACC_GENPP_HANDLE_BITS; off++) {
    size_t i = ACC_GENPP_HANDLE_BITS - 1 - off;
    acc_genpp_ginv_mat(&inv, &tail[i + 1]);
    acc_genpp_Di(&D, B, i, acc_genpp_handle_bit(selector, i));
    acc_genpp_gmat_mul_ginv(&tail[i], &D, &inv);
  }

  acc_genpp_ginv_vec(ginv_u, u);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_zero(delta->right[i][j]);
    }
  }
  poly_q_vec_d_zero(delta->left);

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    if (!prefix_equal) {
      break;
    }

    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_set(tmp_in[j], ginv_u[j]);
    }

    if (i + 1 < ACC_GENPP_HANDLE_BITS) {
      acc_genpp_ginv_mat(&inv, &tail[i + 1]);
      for (size_t r = 0; r < ACC_GENPP_GADGET_DIM; r++) {
        poly_q_vec_d_zero(tmp_out[r]);
        for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
          poly_q_mat_d_d_muladd_vec_d(tmp_out[r], inv.blocks[r][c], tmp_in[c]);
        }
      }
    } else {
      for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
        poly_q_vec_d_set(tmp_out[c], tmp_in[c]);
      }
    }

    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (acc_genpp_handle_bit(selector, i)) {
        poly_q_vec_d_set(delta->right[i][j], tmp_out[j]);
      } else {
        poly_q_vec_d_neg(delta->right[i][j], tmp_out[j]);
      }
    }

    if (acc_genpp_handle_bit(selector, i) != acc_genpp_handle_bit(point, i)) {
      prefix_equal = 0;
    }
  }

  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_clear(tmp_out[i]);
    poly_q_vec_d_clear(tmp_in[i]);
    poly_q_vec_d_clear(ginv_u[i]);
  }
  acc_genpp_ginv_mat_clear(&inv);
  acc_genpp_gmat_clear(&D);
  for (size_t i = 0; i <= ACC_GENPP_HANDLE_BITS; i++) {
    acc_genpp_gmat_clear(&tail[i]);
  }
}

void acc_genpp_sampleleft(acc_genpp_wit_t *e,
                          const acc_genpp_pp_t *pp,
                          const acc_genpp_td_t *td,
                          const poly_q_vec_d u,
                          uint32_t x,
                          const uint8_t ra[ACC_GENPP_SEED_BYTES]) {
  poly_q_mat_d_d A;
  poly_q_vec_d right_contrib, target;

  poly_q_mat_d_d_init(A);
  poly_q_vec_d_init(right_contrib);
  poly_q_vec_d_init(target);

  poly_q_zero(e->tag);
  poly_q_set_coeff(e->tag, 0, 1);
  poly_q_mat_d_d_uniform(A, pp->A.seed, DOMAIN_SEPARATOR_A, 0);

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      acc_genpp_sample_small_vec_d(e->right[i][j], ra, x, i, j);
    }
  }

  acc_genpp_mul_Mx_right(right_contrib, pp, x, e->right);
  poly_q_vec_d_sub(target, u, right_contrib);
  acc_trunc_sampler(e->left, e->v12, e->v2, td->trunc_R, A, pp->trunc_B, target, e->tag, td->trunc_S);

  poly_q_vec_d_clear(target);
  poly_q_vec_d_clear(right_contrib);
  poly_q_mat_d_d_clear(A);
}

void acc_genpp_gen(acc_genpp_pp_t *pp,
                   acc_genpp_td_t *td,
                   acc_genpp_state_t *st,
                   poly_q_vec_d A0) {
  uint8_t seed[SEED_BYTES], rseed[SEED_BYTES];
  poly_q_mat_d_d A, RRstar[2][2];

  sep_keygen(&pp->A, &td->RA);
  randombytes(seed, sizeof(seed));
  acc_genpp_sample_Ba(&pp->Ba, seed);
  randombytes(rseed, sizeof(rseed));

  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(RRstar[0][0]);
  poly_q_mat_d_d_init(RRstar[0][1]);
  poly_q_mat_d_d_init(RRstar[1][0]);
  poly_q_mat_d_d_init(RRstar[1][1]);

  poly_q_mat_d_d_uniform(A, pp->A.seed, DOMAIN_SEPARATOR_A, 0);
  uint32_t kappa = 0;
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_binomial(td->trunc_R[j][i], rseed, kappa++, DOMAIN_SEPARATOR_R);
    }
    poly_q_mat_d_d_mul_mat_d_d(pp->trunc_B[i], A, td->trunc_R[1][i]);
    poly_q_mat_d_d_add(pp->trunc_B[i], pp->trunc_B[i], td->trunc_R[0][i]);
  }
  acc_trunc_compute_rrstar(RRstar, td->trunc_R);
  acc_trunc_compute_covariance(td->trunc_S, RRstar);

  randombytes(st->ka, sizeof(st->ka));
  poly_q_vec_d_uniform(pp->ua, seed, DOMAIN_SEPARATOR_U);
  poly_q_vec_d_uniform(A0, seed, DOMAIN_SEPARATOR_DS);
  poly_q_vec_d_set(st->acc0, A0);
  st->ctr = 0;
  memset(st->dict, 0, sizeof(st->dict));

  poly_q_mat_d_d_clear(RRstar[1][1]);
  poly_q_mat_d_d_clear(RRstar[1][0]);
  poly_q_mat_d_d_clear(RRstar[0][1]);
  poly_q_mat_d_d_clear(RRstar[0][0]);
  poly_q_mat_d_d_clear(A);
}

int acc_genpp_add(poly_q_vec_d A_out,
                  acc_genpp_wit_t *sx,
                  acc_genpp_update_t *upmsg,
                  acc_genpp_state_t *st,
                  const acc_genpp_pp_t *pp,
                  const acc_genpp_td_t *td,
                  const poly_q_vec_d A,
                  uint32_t x) {
  uint8_t ra[ACC_GENPP_SEED_BYTES];
  int pos = acc_genpp_dict_find(st, x);

  if (pos >= 0 && st->dict[pos].op != ACC_GENPP_OP_EMPTY) {
    return 0;
  }

  acc_genpp_Fk(ra, st->ka, x, A);
  acc_genpp_sampleleft(sx, pp, td, A, x, ra);
  if (!acc_genpp_dict_put(st, x, ACC_GENPP_OP_ADD, st->ctr)) {
    return 0;
  }

  poly_q_vec_d_set(A_out, A);
  upmsg->seq = st->ctr;
  upmsg->op = ACC_GENPP_OP_ADD;
  upmsg->handle = 0;
  return 1;
}

int acc_genpp_delete(poly_q_vec_d A_out,
                     acc_genpp_update_t *upmsg,
                     acc_genpp_state_t *st,
                     const acc_genpp_pp_t *pp,
                     const acc_genpp_td_t *td,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx) {
  acc_genpp_gmat_t B1x;
  poly_q_vec_d ginv_ua[ACC_GENPP_GADGET_DIM], delta;
  int pos;

  (void)td;
  if (st->ctr >= ACC_GENPP_MAX_DELETIONS) {
    return 0;
  }
  pos = acc_genpp_dict_find(st, x);
  if (pos < 0 || st->dict[pos].op != ACC_GENPP_OP_ADD) {
    return 0;
  }
  if (!acc_genpp_memver(pp, A, x, sx)) {
    return 0;
  }

  acc_genpp_gmat_init(&B1x);
  poly_q_vec_d_init(delta);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_init(ginv_ua[i]);
  }

  acc_genpp_evalf(&B1x, &pp->Ba, x);
  acc_genpp_ginv_vec(ginv_ua, pp->ua);
  acc_genpp_gmat_mul_gvec(delta, &B1x, ginv_ua);
  poly_q_vec_d_add(A_out, A, delta);

  st->dict[pos].op = ACC_GENPP_OP_DEL;
  st->dict[pos].index = st->ctr + 1;
  st->ctr++;
  upmsg->seq = st->ctr;
  upmsg->op = ACC_GENPP_OP_DEL;
  upmsg->handle = x;
  acc_genpp_bml_push(st, upmsg);

  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_clear(ginv_ua[i]);
  }
  poly_q_vec_d_clear(delta);
  acc_genpp_gmat_clear(&B1x);
  return 1;
}

void acc_genpp_memwitup(acc_genpp_wit_t *sx,
                        const acc_genpp_pp_t *pp,
                        uint32_t x,
                        const acc_genpp_update_t *upmsg) {
  acc_genpp_wit_t delta;

  if (upmsg->op != ACC_GENPP_OP_DEL) {
    return;
  }

  acc_genpp_wit_init(&delta);
  acc_genpp_evalfx_apply(&delta, &pp->Ba, upmsg->handle, x, pp->ua);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_add(sx->right[i][j], sx->right[i][j], delta.right[i][j]);
    }
  }
  acc_genpp_wit_clear(&delta);
}

int acc_genpp_bml_get_missing(const acc_genpp_state_t *st,
                              uint32_t last_seq,
                              acc_genpp_update_t *out,
                              size_t out_cap,
                              size_t *out_len) {
  size_t count = 0;
  uint32_t oldest_seq;

  if (out_len == NULL || last_seq > st->ctr) {
    return 0;
  }
  *out_len = 0;
  if (last_seq == st->ctr) {
    return 1;
  }
  if (st->bml_len == 0) {
    return 0;
  }

  oldest_seq = st->ctr - (uint32_t)st->bml_len + 1;
  if (last_seq + 1 < oldest_seq) {
    return 0;
  }

  for (size_t i = 0; i < st->bml_len; i++) {
    size_t pos = (st->bml_start + i) % ACC_GENPP_BML_CAPACITY;
    const acc_genpp_update_t *msg = &st->bml[pos];

    if (msg->seq <= last_seq) {
      continue;
    }
    if (count >= out_cap) {
      return 0;
    }
    if (out != NULL) {
      out[count] = *msg;
    }
    count++;
  }
  *out_len = count;
  return 1;
}

int acc_genpp_memver(const acc_genpp_pp_t *pp,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx) {
  poly_q_vec_d lhs;
  poly_q_mat_d_d Atrap, Gblock, Mblock;
  poly_q tag_times_bexpi;
  int ok;

  poly_q_vec_d_init(lhs);
  poly_q_mat_d_d_init(Atrap);
  poly_q_mat_d_d_init(Gblock);
  poly_q_mat_d_d_init(Mblock);
  poly_q_init(tag_times_bexpi);

  poly_q_vec_d_set(lhs, sx->left);

  poly_q_mat_d_d_uniform(Atrap, pp->A.seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_d_muladd_vec_d(lhs, Atrap, sx->v12);

  coeff_q bexpi = ACC_GENPP_TRUNC_B_POW_L;
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_set(Mblock, pp->trunc_B[i]);
    poly_q_mul_scalar(tag_times_bexpi, sx->tag, bexpi);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_q_sub(Mblock->rows[j]->entries[j], Mblock->rows[j]->entries[j], tag_times_bexpi);
    }
    poly_q_mat_d_d_mulsub_vec_d(lhs, Mblock, sx->v2[i]);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    int xbit = acc_genpp_handle_bit(x, i);
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_set(Mblock, pp->Ba.bit[i].blocks[j]);
      if (xbit) {
        acc_genpp_gadget_block(Gblock, j);
        poly_q_mat_d_d_sub(Mblock, Mblock, Gblock);
      }
      poly_q_mat_d_d_muladd_vec_d(lhs, Mblock, sx->right[i][j]);
    }
  }

  ok = poly_q_vec_d_equal(lhs, A) && (acc_genpp_wit_norm2(sx) <= ACC_GENPP_WIT_BETA_SQ);
  poly_q_clear(tag_times_bexpi);
  poly_q_mat_d_d_clear(Mblock);
  poly_q_mat_d_d_clear(Gblock);
  poly_q_mat_d_d_clear(Atrap);
  poly_q_vec_d_clear(lhs);
  return ok;
}

int acc_genpp_memwitsync(acc_genpp_wit_t *sx,
                         const acc_genpp_pp_t *pp,
                         const acc_genpp_td_t *td,
                         const acc_genpp_state_t *st,
                         const poly_q_vec_d A,
                         uint32_t x) {
  poly_q_vec_d acc_i;
  acc_genpp_gmat_t B1y;
  poly_q_vec_d ginv_ua[ACC_GENPP_GADGET_DIM], delta;
  uint8_t ra[ACC_GENPP_SEED_BYTES];
  int pos = acc_genpp_dict_find(st, x);
  uint32_t add_index;

  if (pos < 0 || st->dict[pos].op != ACC_GENPP_OP_ADD) {
    return 0;
  }
  add_index = st->dict[pos].index;

  poly_q_vec_d_init(acc_i);
  poly_q_vec_d_init(delta);
  acc_genpp_gmat_init(&B1y);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_init(ginv_ua[i]);
  }

  poly_q_vec_d_set(acc_i, st->acc0);
  acc_genpp_ginv_vec(ginv_ua, pp->ua);
  for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
    if (st->dict[i].used && st->dict[i].op == ACC_GENPP_OP_DEL &&
        st->dict[i].index > 0 && st->dict[i].index <= add_index) {
      acc_genpp_evalf(&B1y, &pp->Ba, st->dict[i].handle);
      acc_genpp_gmat_mul_gvec(delta, &B1y, ginv_ua);
      poly_q_vec_d_add(acc_i, acc_i, delta);
    }
  }

  acc_genpp_Fk(ra, st->ka, x, acc_i);
  acc_genpp_sampleleft(sx, pp, td, acc_i, x, ra);
  for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
    if (st->dict[i].used && st->dict[i].op == ACC_GENPP_OP_DEL &&
        st->dict[i].index > add_index && st->dict[i].index <= st->ctr) {
      acc_genpp_update_t msg = {st->dict[i].index, ACC_GENPP_OP_DEL, st->dict[i].handle};
      acc_genpp_memwitup(sx, pp, x, &msg);
    }
  }

  (void)A;
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_clear(ginv_ua[i]);
  }
  acc_genpp_gmat_clear(&B1y);
  poly_q_vec_d_clear(delta);
  poly_q_vec_d_clear(acc_i);
  return 1;
}
