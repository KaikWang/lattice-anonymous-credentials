#!/usr/bin/env python3
"""Dimension layout for embedding the real AccGenpp MemVer relation in Labrador.

This is not a toy algebra description.  The dimensions mirror the current C
implementation in acc_genpp.h/c:

  left + Atrap*v12
  - sum_i (trunc_B_i - tag*b^e_i*I)*v2_i
  + sum_{h,j} (Ba[h][j] - x_h*G_j)*right[h][j] = acc

The companion exporter/importer can now feed real C coefficients into a
Labrador statement and lift the q=405769 equation to the proof modulus.
"""

from __future__ import annotations

PARAM_D = 4
PARAM_N = 256
PARAM_Q_BITLEN = 19
ACC_GENPP_HANDLE_BITS = 32
ACC_GENPP_GADGET_DIM = PARAM_Q_BITLEN
ACC_GENPP_TRUNC_K_L = 3


def main() -> int:
    base_vec_polys = PARAM_D
    right_vecs = ACC_GENPP_HANDLE_BITS * ACC_GENPP_GADGET_DIM
    right_polys = right_vecs * PARAM_D
    compact_polys = (
        PARAM_D +  # left
        PARAM_D +  # v12
        ACC_GENPP_TRUNC_K_L * PARAM_D +  # v2
        1 +  # tag
        ACC_GENPP_HANDLE_BITS  # x bits as degree-0 binary polynomials
    )
    total_witness_polys = compact_polys + right_polys

    # The vector equality has PARAM_D polynomial equations.  Each right vector
    # contributes one quadratic selector term per output row/gadget combination;
    # tag*v2 contributes one more quadratic term per truncation block and row.
    linear_terms = {
        "left": PARAM_D,
        "Atrap*v12": PARAM_D * PARAM_D,
        "trunc_B*v2": ACC_GENPP_TRUNC_K_L * PARAM_D * PARAM_D,
        "Ba*right": right_vecs * PARAM_D * PARAM_D,
    }
    quadratic_x_right_terms = right_vecs * PARAM_D
    quadratic_tag_v2_terms = ACC_GENPP_TRUNC_K_L * PARAM_D
    quadratic_terms = quadratic_x_right_terms + quadratic_tag_v2_terms

    full_witness_bytes = (total_witness_polys * PARAM_N * PARAM_Q_BITLEN + 7) // 8
    compact_witness_bytes = ((compact_polys - ACC_GENPP_HANDLE_BITS) * PARAM_N * PARAM_Q_BITLEN + 7) // 8

    print("Full AccGenpp MemVer -> Labrador layout")
    print("=======================================")
    print()
    print("Witness groups:")
    print(f"- x bits                  : {ACC_GENPP_HANDLE_BITS} degree-0 binary polynomials")
    print(f"- left                    : {PARAM_D} polynomials")
    print(f"- v12                     : {PARAM_D} polynomials")
    print(f"- v2                      : {ACC_GENPP_TRUNC_K_L * PARAM_D} polynomials")
    print("- tag                     : 1 polynomial")
    print(f"- right[h][j]             : {right_vecs} vec_d = {right_polys} polynomials")
    print(f"- total witness polys     : {total_witness_polys}")
    print()
    print("Constraints:")
    print(f"- polynomial equations    : {PARAM_D}")
    print(f"- binary constraints      : {ACC_GENPP_HANDLE_BITS}")
    print("- degree-0 constraints    : 32 for x bits, optionally 1 for tag if fixed")
    print("- norm constraints        : one aggregate ||s_x|| <= beta, or partitioned bounds")
    print()
    print("Linear term counts in MemVer:")
    for name, count in linear_terms.items():
        print(f"- {name:24s}: {count}")
    print(f"- quadratic x_h*right[h][j]: {quadratic_x_right_terms}")
    print(f"- quadratic tag*v2         : {quadratic_tag_v2_terms}")
    print(f"- quadratic total          : {quadratic_terms}")
    print()
    print("Payload scale:")
    print(f"- compact C witness payload: {compact_witness_bytes / 1024:.2f} KiB")
    print(f"- full C witness payload   : {full_witness_bytes / 1024:.2f} KiB")
    print()
    print("Remaining production integration:")
    print("- encode beta as a reviewed Labrador norm policy")
    print("- serialize/verify the resulting Labrador proof bytes in revac_gcs_verify_show")
    print("- replace the Python bridge with a C/C++ binding if the UAV/GCS path must be standalone")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
