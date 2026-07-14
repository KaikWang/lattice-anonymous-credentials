# Revocable Anonymous Credential Entity Layout

This directory separates the revocable anonymous credential workflow by entity.

- `common/`: shared public types and helpers, including `handle = H(PK, ID)`.
- `ta/`: trusted authority logic: setup, UAV/GCS registration, credential issuance, revocation, and BML maintenance.
- `uav/`: UAV-side logic: key generation, update application, and missed-broadcast recovery.
- `gcs/`: ground control station logic: key generation and credential/witness verification.

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

The bundled FLINT library may still need compatible system GMP/MPFR runtime
libraries on the target. The current random backend uses AES-NI intrinsics, so
non-x86 UAV targets will need a portable random backend before deployment.
