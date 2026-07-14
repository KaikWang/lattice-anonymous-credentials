#include "revocable_ac.h"

#include "fips202.h"

static void revac_absorb_u64(keccak_state *state, uint64_t x) {
  uint8_t buf[8];

  for (size_t i = 0; i < sizeof(buf); i++) {
    buf[i] = (uint8_t)(x >> (8 * i));
  }
  shake256_absorb(state, buf, sizeof(buf));
}

static void revac_absorb_poly_q(keccak_state *state, const poly_q p) {
  for (size_t i = 0; i < PARAM_N; i++) {
    revac_absorb_u64(state, (uint64_t)poly_q_get_coeff(p, i));
  }
}

uint32_t revac_handle_from_identity(const user_pk_t *pk, uint64_t id) {
  keccak_state state;
  uint8_t out[4];

  shake256_init(&state);
  revac_absorb_u64(&state, id);
  shake256_absorb(&state, pk->seed, sizeof(pk->seed));
  for (size_t i = 0; i < PARAM_D; i++) {
    revac_absorb_poly_q(&state, pk->t->entries[i]);
  }
  shake256_finalize(&state);
  shake256_squeeze(out, sizeof(out), &state);
  return ((uint32_t)out[0]) |
         ((uint32_t)out[1] << 8) |
         ((uint32_t)out[2] << 16) |
         ((uint32_t)out[3] << 24);
}
