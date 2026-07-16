#include "revac_proof_codec.h"

#include "huffman_compress.h"

#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define REVAC_DEFLATE_MAGIC "RVZD1"
#define REVAC_DEFLATE_MAGIC_LEN 5
#define REVAC_SIG_MAGIC "RVSGC1"
#define REVAC_SIG_MAGIC_LEN 6
#define REVAC_PROOF_MAGIC "RVACP1"
#define REVAC_PROOF_MAGIC_LEN 6

static void store_u32(uint8_t *out, uint32_t x) {
  out[0] = (uint8_t)x;
  out[1] = (uint8_t)(x >> 8);
  out[2] = (uint8_t)(x >> 16);
  out[3] = (uint8_t)(x >> 24);
}

static void store_u64(uint8_t *out, uint64_t x) {
  for (size_t i = 0; i < 8; i++) {
    out[i] = (uint8_t)(x >> (8 * i));
  }
}

static uint32_t load_u32(const uint8_t *in) {
  return ((uint32_t)in[0]) |
         ((uint32_t)in[1] << 8) |
         ((uint32_t)in[2] << 16) |
         ((uint32_t)in[3] << 24);
}

static uint64_t load_u64(const uint8_t *in) {
  uint64_t x = 0;
  for (size_t i = 0; i < 8; i++) {
    x |= ((uint64_t)in[i]) << (8 * i);
  }
  return x;
}

static int ensure_room(size_t off, size_t need, size_t len) {
  return off <= len && need <= len - off;
}

static void poly_to_uniform(int64_t *out, const poly_qshow p) {
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    out[i] = (int64_t)poly_qshow_get_coeff(p, i);
  }
}

static void poly_to_centered(int64_t *out, const poly_qshow p) {
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    out[i] = (int64_t)poly_qshow_get_coeff_centered(p, i);
  }
}

static void uniform_to_poly(poly_qshow p, const int64_t *in) {
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(p, i, (coeff_qshow)in[i]);
  }
}

static void centered_to_poly(poly_qshow p, const int64_t *in) {
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(p, i, (coeff_qshow)in[i]);
  }
}

static size_t pack_uniform_polys(uint8_t *out, const poly_qshow *polys,
                                 size_t npolys) {
  int64_t coeffs[PARAM_N_SHOW * PARAM_D_SHOW];
  size_t pos = 0;

  for (size_t i = 0; i < npolys; i++) {
    poly_to_uniform(coeffs + pos, polys[i]);
    pos += PARAM_N_SHOW;
  }
  return compact_pack_uniform(out, coeffs, pos);
}

static size_t unpack_uniform_polys(poly_qshow *polys, size_t npolys,
                                   const uint8_t *in, size_t in_len) {
  int64_t coeffs[PARAM_N_SHOW * PARAM_D_SHOW];
  size_t count = npolys * PARAM_N_SHOW;
  size_t pos = 0;

  compact_unpack_uniform(coeffs, in, in_len, count);
  for (size_t i = 0; i < npolys; i++) {
    uniform_to_poly(polys[i], coeffs + pos);
    pos += PARAM_N_SHOW;
  }
  return count;
}

static size_t pack_one_uniform(uint8_t *out, const poly_qshow p) {
  int64_t coeffs[PARAM_N_SHOW];

  poly_to_uniform(coeffs, p);
  return compact_pack_uniform(out, coeffs, PARAM_N_SHOW);
}

static void unpack_one_uniform(poly_qshow p, const uint8_t *in, size_t in_len) {
  int64_t coeffs[PARAM_N_SHOW];

  compact_unpack_uniform(coeffs, in, in_len, PARAM_N_SHOW);
  uniform_to_poly(p, coeffs);
}

static size_t pack_challenge(uint8_t *out, const poly_qshow p) {
  int64_t coeffs[PARAM_N_SHOW];

  poly_to_centered(coeffs, p);
  return compact_pack_challenge(out, coeffs, PARAM_N_SHOW);
}

static void unpack_challenge(poly_qshow p, const uint8_t *in, size_t in_len) {
  int64_t coeffs[PARAM_N_SHOW];

  compact_unpack_challenge(coeffs, in, in_len, PARAM_N_SHOW);
  centered_to_poly(p, coeffs);
}

static size_t pack_centered_vec_m1(uint8_t *out, const poly_qshow_vec_m1 v) {
  int64_t *coeffs = malloc(PARAM_M1_SHOW * PARAM_N_SHOW * sizeof(*coeffs));
  size_t off = 0;
  size_t got = 0;

  if (coeffs == NULL) {
    return 0;
  }
  for (size_t i = 0; i < PARAM_M1_SHOW; i++) {
    poly_to_centered(coeffs + off, v->entries[i]);
    off += PARAM_N_SHOW;
  }
  got = huffman_pack_z1(out, coeffs, off);
  free(coeffs);
  return got;
}

static size_t pack_centered_vec_m2(uint8_t *out, const poly_qshow_vec_m2 v) {
  int64_t *coeffs = malloc(PARAM_M2_SHOW * PARAM_N_SHOW * sizeof(*coeffs));
  size_t off = 0;
  size_t got = 0;

  if (coeffs == NULL) {
    return 0;
  }
  for (size_t i = 0; i < PARAM_M2_SHOW; i++) {
    poly_to_centered(coeffs + off, v->entries[i]);
    off += PARAM_N_SHOW;
  }
  got = huffman_pack_z2(out, coeffs, off);
  free(coeffs);
  return got;
}

static int unpack_centered_vec_m1(poly_qshow_vec_m1 v,
                                  const uint8_t *in, size_t in_len) {
  int64_t *coeffs = malloc(PARAM_M1_SHOW * PARAM_N_SHOW * sizeof(*coeffs));
  size_t off = 0;

  if (coeffs == NULL) {
    return 0;
  }
  huffman_unpack_z1(coeffs, in, in_len, PARAM_M1_SHOW * PARAM_N_SHOW);
  for (size_t i = 0; i < PARAM_M1_SHOW; i++) {
    centered_to_poly(v->entries[i], coeffs + off);
    off += PARAM_N_SHOW;
  }
  free(coeffs);
  return 1;
}

static int unpack_centered_vec_m2(poly_qshow_vec_m2 v,
                                  const uint8_t *in, size_t in_len) {
  int64_t *coeffs = malloc(PARAM_M2_SHOW * PARAM_N_SHOW * sizeof(*coeffs));
  size_t off = 0;

  if (coeffs == NULL) {
    return 0;
  }
  huffman_unpack_z2(coeffs, in, in_len, PARAM_M2_SHOW * PARAM_N_SHOW);
  for (size_t i = 0; i < PARAM_M2_SHOW; i++) {
    centered_to_poly(v->entries[i], coeffs + off);
    off += PARAM_N_SHOW;
  }
  free(coeffs);
  return 1;
}

static size_t pack_z3(uint8_t *out, const coeff_qshow z3[PARAM_ARP_SHOW]) {
  int64_t coeffs[PARAM_ARP_SHOW];

  for (size_t i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow z = z3[i];
    if (z > PARAM_Q_SHOW / 2) {
      z -= PARAM_Q_SHOW;
    }
    coeffs[i] = (int64_t)z;
  }
  return huffman_pack_z3(out, coeffs, PARAM_ARP_SHOW);
}

static void unpack_z3(coeff_qshow z3[PARAM_ARP_SHOW],
                      const uint8_t *in, size_t in_len) {
  int64_t coeffs[PARAM_ARP_SHOW];

  huffman_unpack_z3(coeffs, in, in_len, PARAM_ARP_SHOW);
  for (size_t i = 0; i < PARAM_ARP_SHOW; i++) {
    z3[i] = (coeff_qshow)coeffs[i];
  }
}

size_t revac_show_signature_raw_size_bytes(void) {
  return POLYQSHOW_VECD_PACKEDBYTES +
         POLYQSHOW_VEC256L_PACKEDBYTES +
         PARAM_ARP_SHOW * COEFFQSHOW_PACKEDBYTES +
         POLYQSHOW_VECL_PACKEDBYTES +
         3 * POLYQSHOW_PACKEDBYTES +
         POLYQSHOW_PACKEDBYTES +
         sizeof(uint32_t) +
         PARAM_M1_SHOW * POLYQSHOW_PACKEDBYTES +
         PARAM_M2_SHOW * POLYQSHOW_PACKEDBYTES;
}

size_t revac_show_signature_compact_bound_bytes(void) {
  return REVAC_SIG_MAGIC_LEN +
         COMPACT_TA_BYTES +
         COMPACT_TB_BYTES +
         HUFFMAN_Z3_COMPRESSED_BYTES_MAX +
         COMPACT_H_BYTES +
         3 * COMPACT_T1_BYTES +
         COMPACT_C_BYTES +
         COMPACT_CTRC_BYTES +
         HUFFMAN_Z1_COMPRESSED_BYTES_MAX +
         HUFFMAN_Z2_COMPRESSED_BYTES_MAX +
         64;
}

size_t revac_show_proof_wire_size_bytes(const revac_show_proof_t *proof) {
  return sizeof(proof->nonce) +
         revac_show_signature_raw_size_bytes() +
         sizeof(uint32_t) +
         sizeof(uint64_t) +
         proof->acc_zk_proof_len;
}

int revac_deflate_pack(uint8_t **out, size_t *out_len,
                       const uint8_t *in, size_t in_len) {
  uLongf cap;
  uint8_t *buf;
  int zret;

  *out = NULL;
  *out_len = 0;
  if (in == NULL && in_len != 0) {
    return 0;
  }
  cap = compressBound((uLong)in_len);
  buf = malloc(REVAC_DEFLATE_MAGIC_LEN + 8 + (size_t)cap);
  if (buf == NULL) {
    return 0;
  }
  memcpy(buf, REVAC_DEFLATE_MAGIC, REVAC_DEFLATE_MAGIC_LEN);
  store_u64(buf + REVAC_DEFLATE_MAGIC_LEN, (uint64_t)in_len);
  zret = compress2(buf + REVAC_DEFLATE_MAGIC_LEN + 8, &cap, in, (uLong)in_len,
                   Z_BEST_COMPRESSION);
  if (zret != Z_OK) {
    free(buf);
    return 0;
  }
  *out = buf;
  *out_len = REVAC_DEFLATE_MAGIC_LEN + 8 + (size_t)cap;
  return 1;
}

int revac_deflate_unpack(uint8_t **out, size_t *out_len,
                         const uint8_t *in, size_t in_len) {
  uint64_t raw_len64;
  uLongf raw_len;
  uint8_t *buf;
  int zret;

  *out = NULL;
  *out_len = 0;
  if (in_len < REVAC_DEFLATE_MAGIC_LEN + 8 ||
      memcmp(in, REVAC_DEFLATE_MAGIC, REVAC_DEFLATE_MAGIC_LEN) != 0) {
    return 0;
  }
  raw_len64 = load_u64(in + REVAC_DEFLATE_MAGIC_LEN);
  if (raw_len64 > (uint64_t)((size_t)-1)) {
    return 0;
  }
  raw_len = (uLongf)raw_len64;
  buf = malloc((size_t)raw_len == 0 ? 1 : (size_t)raw_len);
  if (buf == NULL) {
    return 0;
  }
  zret = uncompress(buf, &raw_len, in + REVAC_DEFLATE_MAGIC_LEN + 8,
                    (uLong)(in_len - REVAC_DEFLATE_MAGIC_LEN - 8));
  if (zret != Z_OK || raw_len != (uLongf)raw_len64) {
    free(buf);
    return 0;
  }
  *out = buf;
  *out_len = (size_t)raw_len;
  return 1;
}

int revac_show_signature_pack_compact(uint8_t **out, size_t *out_len,
                                      const show_proof_t *proof) {
  size_t cap = 2 * revac_show_signature_raw_size_bytes();
  uint8_t *buf = malloc(cap);
  size_t off = 0;
  size_t got;

  *out = NULL;
  *out_len = 0;
  if (buf == NULL) {
    return 0;
  }
  memcpy(buf, REVAC_SIG_MAGIC, REVAC_SIG_MAGIC_LEN);
  off += REVAC_SIG_MAGIC_LEN;

  got = pack_uniform_polys(buf + off, proof->tA->entries, PARAM_D_SHOW);
  off += got;
  got = pack_uniform_polys(buf + off, proof->tB->entries,
                           PARAM_ARP_DIV_N_L_SHOW);
  off += got;
  got = pack_z3(buf + off + 4, proof->z3);
  store_u32(buf + off, (uint32_t)got);
  off += 4 + got;
  got = pack_uniform_polys(buf + off, proof->h->entries, PARAM_L_SHOW);
  off += got;
  off += pack_one_uniform(buf + off, proof->t1);
  off += pack_one_uniform(buf + off, proof->t1_gadget);
  off += pack_one_uniform(buf + off, proof->t0_stored);
  off += pack_challenge(buf + off, proof->c);
  store_u32(buf + off, proof->ctr_c);
  off += 4;
  got = pack_centered_vec_m1(buf + off + 4, proof->z1);
  if (got == 0) {
    free(buf);
    return 0;
  }
  store_u32(buf + off, (uint32_t)got);
  off += 4 + got;
  got = pack_centered_vec_m2(buf + off + 4, proof->z2);
  if (got == 0) {
    free(buf);
    return 0;
  }
  store_u32(buf + off, (uint32_t)got);
  off += 4 + got;

  *out = buf;
  *out_len = off;
  return 1;
}

int revac_show_signature_unpack_compact(show_proof_t *proof,
                                        const uint8_t *in, size_t in_len) {
  size_t off = 0;
  uint32_t z3_len;
  uint32_t z1_len;
  uint32_t z2_len;

  if (in_len < REVAC_SIG_MAGIC_LEN ||
      memcmp(in, REVAC_SIG_MAGIC, REVAC_SIG_MAGIC_LEN) != 0) {
    return 0;
  }
  off += REVAC_SIG_MAGIC_LEN;
  if (!ensure_room(off, COMPACT_TA_BYTES, in_len)) {
    return 0;
  }
  unpack_uniform_polys(proof->tA->entries, PARAM_D_SHOW, in + off,
                       COMPACT_TA_BYTES);
  off += COMPACT_TA_BYTES;
  if (!ensure_room(off, COMPACT_TB_BYTES, in_len)) {
    return 0;
  }
  unpack_uniform_polys(proof->tB->entries, PARAM_ARP_DIV_N_L_SHOW, in + off,
                       COMPACT_TB_BYTES);
  off += COMPACT_TB_BYTES;
  if (!ensure_room(off, 4, in_len)) {
    return 0;
  }
  z3_len = load_u32(in + off);
  off += 4;
  if (!ensure_room(off, z3_len, in_len)) {
    return 0;
  }
  unpack_z3(proof->z3, in + off, z3_len);
  off += z3_len;
  if (!ensure_room(off, COMPACT_H_BYTES, in_len)) {
    return 0;
  }
  unpack_uniform_polys(proof->h->entries, PARAM_L_SHOW, in + off,
                       COMPACT_H_BYTES);
  off += COMPACT_H_BYTES;
  if (!ensure_room(off, 3 * COMPACT_T1_BYTES + COMPACT_C_BYTES + 4, in_len)) {
    return 0;
  }
  unpack_one_uniform(proof->t1, in + off, COMPACT_T1_BYTES);
  off += COMPACT_T1_BYTES;
  unpack_one_uniform(proof->t1_gadget, in + off, COMPACT_T1_BYTES);
  off += COMPACT_T1_BYTES;
  unpack_one_uniform(proof->t0_stored, in + off, COMPACT_T1_BYTES);
  off += COMPACT_T1_BYTES;
  unpack_challenge(proof->c, in + off, COMPACT_C_BYTES);
  off += COMPACT_C_BYTES;
  proof->ctr_c = load_u32(in + off);
  off += 4;
  if (!ensure_room(off, 4, in_len)) {
    return 0;
  }
  z1_len = load_u32(in + off);
  off += 4;
  if (!ensure_room(off, z1_len, in_len) ||
      !unpack_centered_vec_m1(proof->z1, in + off, z1_len)) {
    return 0;
  }
  off += z1_len;
  if (!ensure_room(off, 4, in_len)) {
    return 0;
  }
  z2_len = load_u32(in + off);
  off += 4;
  if (!ensure_room(off, z2_len, in_len) ||
      !unpack_centered_vec_m2(proof->z2, in + off, z2_len)) {
    return 0;
  }
  off += z2_len;
  return off == in_len;
}

int revac_show_proof_pack_compact(uint8_t **out, size_t *out_len,
                                  const revac_show_proof_t *proof,
                                  int deflate_acc_zk) {
  uint8_t *sig = NULL;
  uint8_t *acc = NULL;
  size_t sig_len = 0;
  size_t acc_len = proof->acc_zk_proof_len;
  revac_acc_zk_kind_t kind = proof->acc_zk_kind;
  size_t total;
  uint8_t *buf;
  size_t off = 0;

  *out = NULL;
  *out_len = 0;
  if (!revac_show_signature_pack_compact(&sig, &sig_len, &proof->sig_show)) {
    return 0;
  }
  if (deflate_acc_zk &&
      proof->acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER &&
      proof->acc_zk_proof != NULL &&
      revac_deflate_pack(&acc, &acc_len, proof->acc_zk_proof,
                         proof->acc_zk_proof_len) &&
      acc_len < proof->acc_zk_proof_len) {
    kind = REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE;
  } else {
    free(acc);
    acc = NULL;
    acc_len = proof->acc_zk_proof_len;
  }

  total = REVAC_PROOF_MAGIC_LEN + REVAC_SHOW_NONCE_BYTES + 4 + sig_len +
          4 + 8 + acc_len;
  buf = malloc(total);
  if (buf == NULL) {
    free(sig);
    free(acc);
    return 0;
  }
  memcpy(buf + off, REVAC_PROOF_MAGIC, REVAC_PROOF_MAGIC_LEN);
  off += REVAC_PROOF_MAGIC_LEN;
  memcpy(buf + off, proof->nonce, REVAC_SHOW_NONCE_BYTES);
  off += REVAC_SHOW_NONCE_BYTES;
  store_u32(buf + off, (uint32_t)sig_len);
  off += 4;
  memcpy(buf + off, sig, sig_len);
  off += sig_len;
  store_u32(buf + off, (uint32_t)kind);
  off += 4;
  store_u64(buf + off, (uint64_t)acc_len);
  off += 8;
  if (acc_len != 0) {
    memcpy(buf + off, acc != NULL ? acc : proof->acc_zk_proof, acc_len);
    off += acc_len;
  }

  free(sig);
  free(acc);
  *out = buf;
  *out_len = off;
  return 1;
}

int revac_show_proof_unpack_compact(revac_show_proof_t *proof,
                                    const uint8_t *in, size_t in_len) {
  size_t off = 0;
  uint32_t sig_len;
  uint32_t kind;
  uint64_t acc_len64;

  if (in_len < REVAC_PROOF_MAGIC_LEN ||
      memcmp(in, REVAC_PROOF_MAGIC, REVAC_PROOF_MAGIC_LEN) != 0) {
    return 0;
  }
  off += REVAC_PROOF_MAGIC_LEN;
  if (!ensure_room(off, REVAC_SHOW_NONCE_BYTES + 4, in_len)) {
    return 0;
  }
  memcpy(proof->nonce, in + off, REVAC_SHOW_NONCE_BYTES);
  off += REVAC_SHOW_NONCE_BYTES;
  sig_len = load_u32(in + off);
  off += 4;
  if (!ensure_room(off, sig_len, in_len) ||
      !revac_show_signature_unpack_compact(&proof->sig_show, in + off,
                                           sig_len)) {
    return 0;
  }
  off += sig_len;
  if (!ensure_room(off, 12, in_len)) {
    return 0;
  }
  kind = load_u32(in + off);
  off += 4;
  acc_len64 = load_u64(in + off);
  off += 8;
  if (acc_len64 > (uint64_t)((size_t)-1) ||
      !ensure_room(off, (size_t)acc_len64, in_len)) {
    return 0;
  }
  return revac_show_proof_set_acc_zk(proof, (revac_acc_zk_kind_t)kind,
                                     in + off, (size_t)acc_len64);
}
