#include "revac_labrador_bridge.h"

#include "arith.h"
#include "poly_q_sampling.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#define REVAC_LABRADOR_DEFAULT_PYTHON "/home/wkk/miniconda3/envs/sage/bin/python"
#define REVAC_LABRADOR_DEFAULT_SDE "/home/wkk/tools/sde-external-10.8.0-2026-03-15-lin/sde64"
#define REVAC_LABRADOR_DEFAULT_SCRIPT "/home/wkk/pq_ac/lattice-anonymous-credentials/code/entities/zk_lazer/revac_acc_memver_import_labrador.py"

static const char *revac_env_or_default(const char *name, const char *fallback) {
  const char *value = getenv(name);
  return (value != NULL && value[0] != '\0') ? value : fallback;
}

int revac_labrador_online_enabled(void) {
  const char *enabled = getenv("REVAC_LABRADOR_ONLINE");
  return enabled != NULL &&
         (strcmp(enabled, "1") == 0 ||
          strcmp(enabled, "true") == 0 ||
          strcmp(enabled, "TRUE") == 0 ||
          strcmp(enabled, "yes") == 0);
}

static int revac_labrador_use_sde(void) {
  const char *use_sde = getenv("REVAC_LABRADOR_USE_SDE");
  return use_sde == NULL || strcmp(use_sde, "0") != 0;
}

static void json_poly(FILE *fp, const poly_q p) {
  fprintf(fp, "[");
  for (size_t i = 0; i < PARAM_N; i++) {
    if (i) {
      fprintf(fp, ",");
    }
    fprintf(fp, "%ld", (long)poly_q_get_coeff_centered(p, i));
  }
  fprintf(fp, "]");
}

static void json_vec_d(FILE *fp, const poly_q_vec_d v) {
  fprintf(fp, "[");
  for (size_t i = 0; i < PARAM_D; i++) {
    if (i) {
      fprintf(fp, ",");
    }
    json_poly(fp, v->entries[i]);
  }
  fprintf(fp, "]");
}

static void json_mat_d_d(FILE *fp, const poly_q_mat_d_d m) {
  fprintf(fp, "[");
  for (size_t r = 0; r < PARAM_D; r++) {
    if (r) {
      fprintf(fp, ",");
    }
    fprintf(fp, "[");
    for (size_t c = 0; c < PARAM_D; c++) {
      if (c) {
        fprintf(fp, ",");
      }
      json_poly(fp, m->rows[r]->entries[c]);
    }
    fprintf(fp, "]");
  }
  fprintf(fp, "]");
}

static int write_memver_json(FILE *fp,
                             const revac_ta_t *ta,
                             const poly_q_vec_d acc,
                             const revac_uav_t *uav) {
  poly_q_mat_d_d Atrap;
  int include_witness = (uav != NULL);

  poly_q_mat_d_d_init(Atrap);
  poly_q_mat_d_d_uniform(Atrap, ta->acc_pp.A.seed, DOMAIN_SEPARATOR_A, 0);

  fprintf(fp, "{\n");
  fprintf(fp,
          "\"params\":{\"n\":%d,\"d\":%d,\"q\":%ld,\"handle_bits\":%d,"
          "\"gadget_dim\":%d,\"trunc_k_l\":%d},\n",
          PARAM_N, PARAM_D, (long)PARAM_Q, ACC_GENPP_HANDLE_BITS,
          ACC_GENPP_GADGET_DIM, ACC_GENPP_TRUNC_K_L);

  if (include_witness) {
    fprintf(fp, "\"x\":%u,\n", uav->handle);
    fprintf(fp, "\"x_bits\":[");
    for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
      if (i) {
        fprintf(fp, ",");
      }
      fprintf(fp, "%u", (unsigned)((uav->handle >> i) & 1u));
    }
    fprintf(fp, "],\n");
  }

  fprintf(fp, "\"Atrap\":");
  json_mat_d_d(fp, Atrap);

  fprintf(fp, ",\n\"trunc_B\":[");
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (i) {
      fprintf(fp, ",");
    }
    json_mat_d_d(fp, ta->acc_pp.trunc_B[i]);
  }

  fprintf(fp, "],\n\"Ba\":[");
  for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
    if (h) {
      fprintf(fp, ",");
    }
    fprintf(fp, "[");
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (j) {
        fprintf(fp, ",");
      }
      json_mat_d_d(fp, ta->acc_pp.Ba.bit[h].blocks[j]);
    }
    fprintf(fp, "]");
  }

  fprintf(fp, "],\n\"acc\":");
  json_vec_d(fp, acc);

  if (include_witness) {
    fprintf(fp, ",\n\"witness\":{");
    fprintf(fp, "\"left\":");
    json_vec_d(fp, uav->witness.left);
    fprintf(fp, ",\"v12\":");
    json_vec_d(fp, uav->witness.v12);
    fprintf(fp, ",\"v2\":[");
    for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
      if (i) {
        fprintf(fp, ",");
      }
      json_vec_d(fp, uav->witness.v2[i]);
    }
    fprintf(fp, "],\"tag\":");
    json_poly(fp, uav->witness.tag);
    fprintf(fp, ",\"right\":[");
    for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
      if (h) {
        fprintf(fp, ",");
      }
      fprintf(fp, "[");
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        if (j) {
          fprintf(fp, ",");
        }
        json_vec_d(fp, uav->witness.right[h][j]);
      }
      fprintf(fp, "]");
    }
    fprintf(fp, "]}");
  }

  fprintf(fp, "\n}\n");
  poly_q_mat_d_d_clear(Atrap);
  return ferror(fp) == 0;
}

static int make_temp_path(char path[64]) {
  int fd;

  snprintf(path, 64, "/tmp/revac_labrador_%ld_XXXXXX", (long)getpid());
  fd = mkstemp(path);
  if (fd < 0) {
    return -1;
  }
  return fd;
}

static int write_json_temp(char path[64],
                           const revac_ta_t *ta,
                           const poly_q_vec_d acc,
                           const revac_uav_t *uav) {
  int fd = make_temp_path(path);
  FILE *fp;
  int ok;

  if (fd < 0) {
    return 0;
  }
  fp = fdopen(fd, "w");
  if (fp == NULL) {
    close(fd);
    unlink(path);
    return 0;
  }
  ok = write_memver_json(fp, ta, acc, uav);
  if (fclose(fp) != 0) {
    ok = 0;
  }
  if (!ok) {
    unlink(path);
  }
  return ok;
}

static int write_bytes_temp(char path[64], const uint8_t *bytes, size_t len) {
  int fd = make_temp_path(path);
  size_t off = 0;

  if (fd < 0) {
    return 0;
  }
  while (off < len) {
    ssize_t got = write(fd, bytes + off, len - off);
    if (got < 0) {
      if (errno == EINTR) {
        continue;
      }
      close(fd);
      unlink(path);
      return 0;
    }
    off += (size_t)got;
  }
  if (close(fd) != 0) {
    unlink(path);
    return 0;
  }
  return 1;
}

static int read_file(const char *path, uint8_t **out, size_t *out_len) {
  FILE *fp = fopen(path, "rb");
  long len;
  uint8_t *buf;
  int ok = 0;

  *out = NULL;
  *out_len = 0;
  if (fp == NULL) {
    return 0;
  }
  if (fseek(fp, 0, SEEK_END) != 0) {
    goto cleanup;
  }
  len = ftell(fp);
  if (len < 0 || fseek(fp, 0, SEEK_SET) != 0) {
    goto cleanup;
  }
  buf = malloc((size_t)len);
  if (buf == NULL && len != 0) {
    goto cleanup;
  }
  if (len != 0 && fread(buf, 1, (size_t)len, fp) != (size_t)len) {
    free(buf);
    goto cleanup;
  }
  *out = buf;
  *out_len = (size_t)len;
  ok = 1;

cleanup:
  fclose(fp);
  return ok;
}

static int spawn_wait(char *const argv[]) {
  pid_t pid = fork();
  int status;

  if (pid < 0) {
    return 0;
  }
  if (pid == 0) {
    execv(argv[0], argv);
    _exit(127);
  }
  do {
    if (waitpid(pid, &status, 0) < 0) {
      if (errno == EINTR) {
        continue;
      }
      return 0;
    }
    break;
  } while (1);
  return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

static int run_labrador(const char *mode,
                        const char *json_path,
                        const char *proof_flag,
                        const char *proof_path) {
  const char *python = revac_env_or_default("REVAC_LABRADOR_PYTHON",
                                            REVAC_LABRADOR_DEFAULT_PYTHON);
  const char *script = revac_env_or_default("REVAC_LABRADOR_SCRIPT",
                                            REVAC_LABRADOR_DEFAULT_SCRIPT);
  const char *sde = revac_env_or_default("REVAC_LABRADOR_SDE",
                                         REVAC_LABRADOR_DEFAULT_SDE);

  if (revac_labrador_use_sde()) {
    char *const argv[] = {
      (char *)sde, "-icx", "--", (char *)python, (char *)script,
      (char *)json_path, "--mode", (char *)mode,
      (char *)proof_flag, (char *)proof_path, NULL
    };
    return spawn_wait(argv);
  }

  {
    char *const argv[] = {
      (char *)python, (char *)script, (char *)json_path, "--mode",
      (char *)mode, (char *)proof_flag, (char *)proof_path, NULL
    };
    return spawn_wait(argv);
  }
}

int revac_labrador_memver_prove(revac_show_proof_t *proof,
                                const revac_ta_t *ta,
                                const revac_uav_t *uav) {
  char json_path[64] = {0};
  char proof_path[64] = {0};
  uint8_t *proof_bytes = NULL;
  size_t proof_len = 0;
  int proof_fd;
  int ok = 0;

  if (!write_json_temp(json_path, ta, ta->acc, uav)) {
    goto cleanup;
  }
  proof_fd = make_temp_path(proof_path);
  if (proof_fd < 0) {
    goto cleanup;
  }
  close(proof_fd);

  if (!run_labrador("prove", json_path, "--proof-out", proof_path)) {
    goto cleanup;
  }
  if (!read_file(proof_path, &proof_bytes, &proof_len)) {
    goto cleanup;
  }
  ok = revac_show_proof_set_acc_zk(proof, REVAC_ACC_ZK_LABRADOR_MEMVER,
                                   proof_bytes, proof_len);

cleanup:
  free(proof_bytes);
  if (json_path[0]) {
    unlink(json_path);
  }
  if (proof_path[0]) {
    unlink(proof_path);
  }
  return ok;
}

int revac_labrador_memver_verify(const revac_ta_t *ta,
                                 const revac_show_context_t *ctx,
                                 const revac_show_proof_t *proof) {
  char json_path[64] = {0};
  char proof_path[64] = {0};
  int ok = 0;

  if (proof->acc_zk_proof == NULL || proof->acc_zk_proof_len == 0) {
    return 0;
  }
  if (!write_json_temp(json_path, ta, ctx->acc, NULL)) {
    goto cleanup;
  }
  if (!write_bytes_temp(proof_path, proof->acc_zk_proof,
                        proof->acc_zk_proof_len)) {
    goto cleanup;
  }
  ok = run_labrador("verify", json_path, "--proof-in", proof_path);

cleanup:
  if (json_path[0]) {
    unlink(json_path);
  }
  if (proof_path[0]) {
    unlink(proof_path);
  }
  return ok;
}
