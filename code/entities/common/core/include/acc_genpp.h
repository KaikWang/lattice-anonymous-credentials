#ifndef ACC_GENPP_H
#define ACC_GENPP_H

#include <stdint.h>
#include <stddef.h>

#include "params.h"
#include "arith.h"
#include "sep.h"

/*
 * Prototype implementation of the AccGenpp.pdf accumulator algorithms.
 *
 * The PDF uses n = 1 and its own matrix dimensions.  This implementation
 * follows the existing code base and works over R_q^PARAM_D with gadget
 * dimension ceil(log2(q)).  Handles are fixed 32-bit strings as in the PDF.
 */
#define ACC_GENPP_HANDLE_BITS 32
#define ACC_GENPP_GADGET_DIM PARAM_Q_BITLEN
#define ACC_GENPP_MAX_DICT_ENTRIES 128
#define ACC_GENPP_MAX_DELETIONS ACC_GENPP_MAX_DICT_ENTRIES
#define ACC_GENPP_BML_CAPACITY 20
#define ACC_GENPP_SEED_BYTES 32
#define ACC_GENPP_TRUNC_K 5
#define ACC_GENPP_TRUNC_L 2
#define ACC_GENPP_TRUNC_K_L (ACC_GENPP_TRUNC_K - ACC_GENPP_TRUNC_L)
#define ACC_GENPP_TRUNC_B 14
#define ACC_GENPP_TRUNC_B_POW_L 196
#define ACC_GENPP_WIT_BETA_SQ UINT64_C(0x7fffffffffffffff)

typedef enum {
  ACC_GENPP_OP_EMPTY = 0,
  ACC_GENPP_OP_ADD = 1,
  ACC_GENPP_OP_DEL = 2
} acc_genpp_op_t;

typedef struct {
  poly_q_mat_d_d blocks[ACC_GENPP_GADGET_DIM];
} acc_genpp_gmat_t;

typedef struct {
  acc_genpp_gmat_t bit[ACC_GENPP_HANDLE_BITS];
} acc_genpp_Ba_t;

typedef struct {
  uint32_t handle;
  acc_genpp_op_t op;
  uint32_t index;
  int used;
} acc_genpp_dict_entry_t;

typedef struct {
  uint32_t seq;
  acc_genpp_op_t op;
  uint32_t handle;
} acc_genpp_update_t;

typedef struct {
  uint8_t ka[ACC_GENPP_SEED_BYTES];
  uint32_t ctr;
  poly_q_vec_d acc0;
  acc_genpp_dict_entry_t dict[ACC_GENPP_MAX_DICT_ENTRIES];
  acc_genpp_update_t bml[ACC_GENPP_BML_CAPACITY];
  size_t bml_start;
  size_t bml_len;
} acc_genpp_state_t;

typedef struct {
  sep_pk_t A;
  poly_q_mat_d_d trunc_B[ACC_GENPP_TRUNC_K_L];
  acc_genpp_Ba_t Ba;
  poly_q_vec_d ua;
} acc_genpp_pp_t;

typedef struct {
  sep_sk_t RA;
  poly_q_mat_d_d trunc_R[2][ACC_GENPP_TRUNC_K_L];
  poly_real_mat_2d_2d trunc_S;
} acc_genpp_td_t;

typedef struct {
  poly_q_vec_d left;
  poly_q_vec_d v12;
  poly_q_vec_d v2[ACC_GENPP_TRUNC_K_L];
  poly_q tag;
  poly_q_vec_d right[ACC_GENPP_HANDLE_BITS][ACC_GENPP_GADGET_DIM];
} acc_genpp_wit_t;

void acc_genpp_gmat_init(acc_genpp_gmat_t *m);
void acc_genpp_gmat_clear(acc_genpp_gmat_t *m);
void acc_genpp_Ba_init(acc_genpp_Ba_t *B);
void acc_genpp_Ba_clear(acc_genpp_Ba_t *B);
void acc_genpp_state_init(acc_genpp_state_t *st);
void acc_genpp_state_clear(acc_genpp_state_t *st);
void acc_genpp_pp_init(acc_genpp_pp_t *pp);
void acc_genpp_pp_clear(acc_genpp_pp_t *pp);
void acc_genpp_td_init(acc_genpp_td_t *td);
void acc_genpp_td_clear(acc_genpp_td_t *td);
void acc_genpp_wit_init(acc_genpp_wit_t *wit);
void acc_genpp_wit_clear(acc_genpp_wit_t *wit);
size_t acc_genpp_acc_size_bytes(void);
size_t acc_genpp_wit_compact_size_bytes(void);
size_t acc_genpp_wit_full_size_bytes(void);
uint64_t acc_genpp_wit_norm2(const acc_genpp_wit_t *wit);

void acc_genpp_Fk(uint8_t r[ACC_GENPP_SEED_BYTES],
                  const uint8_t ka[ACC_GENPP_SEED_BYTES],
                  uint32_t x,
                  const poly_q_vec_d acc);

void acc_genpp_evalf(acc_genpp_gmat_t *out,
                     const acc_genpp_Ba_t *B,
                     uint32_t selector);

void acc_genpp_evalfx_apply(acc_genpp_wit_t *delta,
                            const acc_genpp_Ba_t *B,
                            uint32_t selector,
                            uint32_t point,
                            const poly_q_vec_d u);

void acc_genpp_sampleleft(acc_genpp_wit_t *e,
                          const acc_genpp_pp_t *pp,
                          const acc_genpp_td_t *td,
                          const poly_q_vec_d u,
                          uint32_t x,
                          const uint8_t ra[ACC_GENPP_SEED_BYTES]);

void acc_genpp_gen(acc_genpp_pp_t *pp,
                   acc_genpp_td_t *td,
                   acc_genpp_state_t *st,
                   poly_q_vec_d A0);

int acc_genpp_add(poly_q_vec_d A_out,
                  acc_genpp_wit_t *sx,
                  acc_genpp_update_t *upmsg,
                  acc_genpp_state_t *st,
                  const acc_genpp_pp_t *pp,
                  const acc_genpp_td_t *td,
                  const poly_q_vec_d A,
                  uint32_t x);

int acc_genpp_delete(poly_q_vec_d A_out,
                     acc_genpp_update_t *upmsg,
                     acc_genpp_state_t *st,
                     const acc_genpp_pp_t *pp,
                     const acc_genpp_td_t *td,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx);

void acc_genpp_memwitup(acc_genpp_wit_t *sx,
                        const acc_genpp_pp_t *pp,
                        uint32_t x,
                        const acc_genpp_update_t *upmsg);

int acc_genpp_bml_get_missing(const acc_genpp_state_t *st,
                              uint32_t last_seq,
                              acc_genpp_update_t *out,
                              size_t out_cap,
                              size_t *out_len);

int acc_genpp_memver(const acc_genpp_pp_t *pp,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx);

int acc_genpp_memwitsync(acc_genpp_wit_t *sx,
                         const acc_genpp_pp_t *pp,
                         const acc_genpp_td_t *td,
                         const acc_genpp_state_t *st,
                         const poly_q_vec_d A,
                         uint32_t x);

#endif
