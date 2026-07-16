#ifndef REVAC_PROOF_CODEC_H
#define REVAC_PROOF_CODEC_H

#include <stddef.h>
#include <stdint.h>

#include "revocable_ac.h"

size_t revac_show_signature_raw_size_bytes(void);
size_t revac_show_signature_compact_bound_bytes(void);
size_t revac_show_proof_wire_size_bytes(const revac_show_proof_t *proof);

int revac_deflate_pack(uint8_t **out, size_t *out_len,
                       const uint8_t *in, size_t in_len);
int revac_deflate_unpack(uint8_t **out, size_t *out_len,
                         const uint8_t *in, size_t in_len);

int revac_show_signature_pack_compact(uint8_t **out, size_t *out_len,
                                      const show_proof_t *proof);
int revac_show_signature_unpack_compact(show_proof_t *proof,
                                        const uint8_t *in, size_t in_len);

int revac_show_proof_pack_compact(uint8_t **out, size_t *out_len,
                                  const revac_show_proof_t *proof,
                                  int deflate_acc_zk);
int revac_show_proof_unpack_compact(revac_show_proof_t *proof,
                                    const uint8_t *in, size_t in_len);

#endif
