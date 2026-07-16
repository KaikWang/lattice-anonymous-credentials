# Revocable Anonymous Credential Entity Layout

This directory separates the revocable anonymous credential workflow by entity.

- `common/`: shared public types and helpers, including `handle = H(PK, ID)`.
- `ta/`: trusted authority logic: setup, UAV/GCS registration, credential issuance, revocation, and BML maintenance.
- `uav/`: UAV-side logic: key generation, update application, and missed-broadcast recovery.
- `gcs/`: ground control station logic: key generation and credential/witness verification.
- `zk_lazer/`: LaZer/Labrador prototype for the future revocable anonymous show proof.

For deployment experiments, copy the entity directory you need plus `common/`.
For example, the UAV-side portable package is:

- `uav/`: UAV-side wrapper and demo build file.
- `common/`: shared revocable-credential helpers.
- `common/core/include/`: copied cryptographic headers.
- `common/core/src/`: copied cryptographic core sources.
- `common/core/flint_headers/`: copied FLINT headers.
- `common/core/lib/libflint.so`: bundled FLINT dynamic library.

Build the UAV package after copying `uav/` and `common/` together:

```sh
cmake -S uav -B /tmp/revac_uav_build
cmake --build /tmp/revac_uav_build -j2
/tmp/revac_uav_build/uav_demo
```

Build the full standalone entity package after copying `entities/`:

```sh
cmake -S entities -B /tmp/revac_entities_build
cmake --build /tmp/revac_entities_build -j2
/tmp/revac_entities_build/revocable_ac_demo
```

The full demo in `test/revocable_ac_demo.c` follows `Plan_process.pdf`:
TA setup, UAV/GCS registration, handle derivation, UAV credential issuance,
GCS certificate gating, GCS verification, revocation by `Delete`, direct
broadcast witness update by `MemWitUp`, offline recovery through BML, fallback
resynchronization by `MemWitSync`, and final verification.

By default, the full demo keeps the expensive Labrador AccGenpp MemVer proof
disabled so the normal entity test remains quick.  In that mode the GCS verifies
the anonymous signature-show proof, then correctly rejects the full revocable
show because no accumulator-membership ZK proof bytes are attached.

To run the online Labrador-backed show path, enable:

```sh
cd entities/build
REVAC_LABRADOR_ONLINE=1 ./revocable_ac_demo
```

With `REVAC_LABRADOR_ONLINE=1`, the UAV serializes a Labrador proof into
`revac_show_proof_t.acc_zk_proof`, and `revac_gcs_verify_show()` verifies those
proof bytes through `zk_lazer/revac_acc_memver_import_labrador.py`.

The compact transport codec is in `common/revac_proof_codec.c`.  It does two
things:

- packs the signature show proof with field-aware bit packing plus the existing
  Huffman tables for Gaussian response vectors;
- optionally deflate-compresses the Labrador accumulator proof bytes and lets
  the GCS verifier inflate them before calling the normal Labrador verifier.

Run the codec demo with:

```sh
cd entities/build
./revac_proof_codec_demo
REVAC_LABRADOR_ONLINE=1 ./revac_proof_codec_demo
```

The online run verifies the compressed accumulator proof bytes end to end.  The
full `test/revocable_ac_demo.c` flow also uses this compact UAV-to-GCS wire
format through `revac_uav_show_prove_wire()` and
`revac_gcs_verify_show_wire()`, and prints the compressed show packet size in
KB during phase 3b.  The older struct-level APIs remain available for local
tests, but the entity demo models the compressed transport path.

The same hidden handle `x` is currently present in both the signature-show statement
and the accumulator MemVer statement, but the production C path still verifies
them as two proofs.  A single joint proof that cryptographically shares one
witness variable `x` across both statements still requires merging the show
protocol and Labrador relation rather than only changing serialization.

The UAV show packet intentionally does not carry `user_seed`, `acc_seq`,
`acc`, or `crs_seed`.  The GCS-side verifier receives those from a
`revac_show_context_t`, modeling data synchronized from TA before verification:
`user_seed` for rebuilding the public `Dx` matrix, and `acc_seq/acc` for the
current accumulator state.  The UAV packet contains only the nonce, signature
show proof, accumulator proof kind, and accumulator proof bytes.

The default paths are:

- `REVAC_LABRADOR_PYTHON=/home/wkk/miniconda3/envs/sage/bin/python`
- `REVAC_LABRADOR_SDE=/home/wkk/tools/sde-external-10.8.0-2026-03-15-lin/sde64`
- `REVAC_LABRADOR_SCRIPT=/home/wkk/pq_ac/lattice-anonymous-credentials/code/entities/zk_lazer/revac_acc_memver_import_labrador.py`

Set `REVAC_LABRADOR_USE_SDE=0` only on hardware that can run the LaZer AVX-512
code directly.  The `zk_lazer/` directory also keeps standalone proof-layer
demos and import checks:

```sh
cd entities/zk_lazer
make describe
make same-x-demo
```

The bundled FLINT library may still need compatible system GMP/MPFR runtime
libraries on the target. The current random backend uses AES-NI intrinsics, so
non-x86 UAV targets will need a portable random backend before deployment.
