#!/usr/bin/env python3
"""Import real AccGenpp MemVer coefficients and build the Labrador statement.

Input is produced by the C target `revac_acc_memver_export`.

The script first checks the exact integer lifting of the C relation

    left + Atrap*v12
    - sum_i trunc_B_i*v2_i
    + sum_i tag*b^e_i*v2_i
    + sum_{h,j} Ba[h][j]*right[h][j]
    - sum_{h,j} x_h*G_j*right[h][j]
    + q*lift
    = acc

over Z[X]/(X^256+1).  With `--mode labrador-simple`, it also appends the real
witnesses, public coefficients, quadratic terms, and lifting witnesses into a
Labrador `proof_statement` and runs `simple_verify`.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from io import BytesIO
from typing import Any


LAZER_PYTHON = pathlib.Path("/home/wkk/pq_ac/lazer/python")
sys.path.append(str(LAZER_PYTHON))


PROOF_MOD_38 = 2**38 - 107
LINEAR_DECOMP_BASE = 256
LINEAR_DECOMP_CHUNKS = 4
SMALL_CHUNK_NORM_SQ = 256 * 128 * 128
TAG_NORM_SQ = 1
V2_NORM_SQ = 256 * 1024 * 1024
RIGHT_NORM_SQ = 256
PROOF_MAGIC = b"REVAC-LABRADOR38-MEMVER\x02"


def zero(n: int) -> list[int]:
    return [0] * n


def add(a: list[int], b: list[int]) -> list[int]:
    return [x + y for x, y in zip(a, b)]


def sub(a: list[int], b: list[int]) -> list[int]:
    return [x - y for x, y in zip(a, b)]


def scalar_mul(c: int, a: list[int]) -> list[int]:
    return [c * x for x in a]


def negacyclic_mul(a: list[int], b: list[int]) -> list[int]:
    n = len(a)
    out = [0] * n
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            k = i + j
            if k >= n:
                out[k - n] -= ai * bj
            else:
                out[k] += ai * bj
    return out


def mat_vec_row(mat: list[list[list[int]]], row: int, vec: list[list[int]]) -> list[int]:
    acc = zero(len(vec[0]))
    for col, v in enumerate(vec):
        acc = add(acc, negacyclic_mul(mat[row][col], v))
    return acc


def const_poly(n: int, c: int) -> list[int]:
    out = [0] * n
    out[0] = c
    return out


def signed_decompose_poly(poly: list[int],
                          base: int = LINEAR_DECOMP_BASE,
                          chunks: int = LINEAR_DECOMP_CHUNKS) -> list[list[int]]:
    """Return small signed base chunks whose weighted sum is poly."""
    out: list[list[int]] = []
    vals = list(poly)
    for _ in range(chunks):
        digit_poly: list[int] = []
        next_vals: list[int] = []
        for v in vals:
            digit = ((v + base // 2) % base) - base // 2
            digit_poly.append(digit)
            next_vals.append((v - digit) // base)
        out.append(digit_poly)
        vals = next_vals
    if any(vals):
        raise ValueError("linear witness coefficient is too large for fixed decomposition")
    return out


def l2sq(poly: list[int]) -> int:
    return sum(c * c for c in poly)


def append_poly(ps: Any, poly_t: Any, ring: Any, polys: list[list[int]], norms: list[int], p: list[int], norm: int) -> int:
    polys.append(p)
    norms.append(norm)
    return len(polys) - 1


def compute_lift(data: dict[str, Any]) -> tuple[list[list[int]], list[list[int]]]:
    params = data["params"]
    n = params["n"]
    q = params["q"]
    handle_bits = params["handle_bits"]
    gadget_dim = params["gadget_dim"]
    trunc_k_l = params["trunc_k_l"]
    w = data["witness"]
    acc = data["acc"]
    x_bits = data["x_bits"]
    tag = w["tag"]

    expr_rows: list[list[int]] = []
    lifts: list[list[int]] = []
    for row in range(params["d"]):
        bexpi = 196
        expr = list(w["left"][row])
        expr = add(expr, mat_vec_row(data["Atrap"], row, w["v12"]))
        for i in range(trunc_k_l):
            expr = sub(expr, mat_vec_row(data["trunc_B"][i], row, w["v2"][i]))
            expr = add(expr, scalar_mul(bexpi, negacyclic_mul(tag, w["v2"][i][row])))
            bexpi *= 14
        for h in range(handle_bits):
            for j in range(gadget_dim):
                expr = add(expr, mat_vec_row(data["Ba"][h][j], row, w["right"][h][j]))
                expr = sub(expr, scalar_mul((1 << j) * x_bits[h], w["right"][h][j][row]))
        diff = sub(acc[row], expr)
        if "memver_lhs" in data:
            lhs_delta = sub(data["memver_lhs"][row], expr)
            lhs_bad = [c for c in lhs_delta if c % q != 0]
            if lhs_bad:
                raise ValueError(
                    f"row {row}: Python lhs disagrees with exported C lhs; "
                    f"first non-q-multiple delta {lhs_bad[0]}"
                )
        bad = [c for c in diff if c % q != 0]
        if bad:
            raise ValueError(f"row {row}: lifting is not divisible by q; first bad coefficient {bad[0]}")
        lifts.append([c // q for c in diff])
        expr_rows.append(expr)
    return expr_rows, lifts


def _write_u64(out: BytesIO, value: int) -> None:
    out.write(int(value).to_bytes(8, "little", signed=False))


def _read_u64(inp: BytesIO) -> int:
    raw = inp.read(8)
    if len(raw) != 8:
        raise ValueError("truncated proof bytes")
    return int.from_bytes(raw, "little", signed=False)


def _write_blob(out: BytesIO, blob: bytes) -> None:
    _write_u64(out, len(blob))
    out.write(blob)


def _read_blob(inp: BytesIO) -> bytes:
    size = _read_u64(inp)
    raw = inp.read(size)
    if len(raw) != size:
        raise ValueError("truncated proof bytes")
    return raw


def _lab_polz_len(pp: Any) -> int:
    if pp.compressed:
        return 4 * pp.kappa[2]
    if not pp.tail:
        return 2 * pp.kappa[1] + 3
    return pp.len[1] + pp.len[2] + pp.len[8] + pp.len[4]


def _set_lab_m_pointers(lab_proof: Any, pp: Any, arr: Any) -> None:
    lab_proof.m[0] = arr
    if pp.compressed:
        lab_proof.m[1] = arr + pp.kappa[2]
        lab_proof.m[2] = lab_proof.m[1] + pp.kappa[2]
        lab_proof.m[3] = lab_proof.m[2] + pp.kappa[2]
    elif not pp.tail:
        lab_proof.m[1] = None
        lab_proof.m[2] = arr + pp.kappa[1]
        lab_proof.m[3] = lab_proof.m[2] + 3
    else:
        lab_proof.m[1] = None
        lab_proof.m[2] = arr + pp.len[1] + pp.len[2]
        lab_proof.m[3] = lab_proof.m[2] + pp.len[8]


def _lnp_polz_len(pp: Any) -> int:
    return 4 * pp.kappa_linfmsis + 4


def _set_lnp_m_pointers(lnp_proof: Any, pp: Any, arr: Any) -> None:
    lnp_proof.m[0] = arr
    lnp_proof.m[1] = arr + pp.kappa_linfmsis
    lnp_proof.m[2] = arr + 2 * pp.kappa_linfmsis
    lnp_proof.m[3] = arr + 2 * pp.kappa_linfmsis + 4
    lnp_proof.m[4] = arr + 3 * pp.kappa_linfmsis + 4


def serialize_labrador38_proof(params: Any, proof: Any) -> bytes:
    from lazer import ffi, lib

    params_nbytes = lib.labrador38_py_params_nbytes(params)
    params_out = ffi.new("uint8_t[]", params_nbytes)
    err = lib.labrador38_py_params_serialize(params_out, params_nbytes, params)
    if err != 0:
        raise ValueError(f"Labrador params serialization failed: {err}")

    proof_nbytes = lib.labrador38_py_proof_nbytes(params, proof)
    proof_out = ffi.new("uint8_t[]", proof_nbytes)
    err = lib.labrador38_py_proof_serialize(proof_out, proof_nbytes, params, proof)
    if err != 0:
        raise ValueError(f"Labrador proof serialization failed: {err}")

    out = BytesIO()
    out.write(PROOF_MAGIC)
    _write_blob(out, bytes(ffi.buffer(params_out, params_nbytes)))
    _write_blob(out, bytes(ffi.buffer(proof_out, proof_nbytes)))
    return out.getvalue()


def deserialize_labrador38_proof(proof_bytes: bytes) -> tuple[Any, Any, list[Any]]:
    from lazer import ffi, lib

    inp = BytesIO(proof_bytes)
    magic = inp.read(len(PROOF_MAGIC))
    if magic != PROOF_MAGIC:
        raise ValueError("bad revocable Labrador proof envelope magic")
    params_bytes = _read_blob(inp)
    inner_proof_bytes = _read_blob(inp)
    if inp.read(1):
        raise ValueError("trailing data after revocable Labrador proof envelope")

    refs: list[Any] = []
    lib.labrador38_py_comkey_reset()

    params = ffi.new("labrador38_dch_pack_params")
    raw_params = ffi.new("uint8_t[]", params_bytes)
    refs.append(raw_params)
    err = lib.labrador38_py_params_deserialize(params, raw_params, len(params_bytes))
    if err != 0:
        raise ValueError(f"Labrador params deserialization failed: {err}")

    proof = ffi.new("labrador38_dch_pack_proof")
    raw = ffi.new("uint8_t[]", inner_proof_bytes)
    refs.append(raw)
    err = lib.labrador38_py_proof_deserialize(proof, params, raw, len(inner_proof_bytes))
    if err != 0:
        raise ValueError(f"Labrador proof deserialization failed: {err}")
    return params, proof, refs


def build_labrador(data: dict[str, Any],
                   lifts: list[list[int]] | None,
                   proof_out: pathlib.Path | None = None,
                   proof_in: pathlib.Path | None = None,
                   debug_hash: bool = False) -> bool:
    from lazer import lib, polyring_t, poly_t
    from labrados import pack_verify, proof_statement

    lib.labrador38_py_comkey_reset()

    params = data["params"]
    n = params["n"]
    q = params["q"]
    handle_bits = params["handle_bits"]
    gadget_dim = params["gadget_dim"]
    trunc_k_l = params["trunc_k_l"]
    ring = polyring_t(n, PROOF_MOD_38)
    one = poly_t(ring, const_poly(n, 1))

    have_witness = lifts is not None and "witness" in data
    witness_polys: list[list[int] | None] = []
    norms: list[int] = []
    w = data.get("witness", {})

    def append_slot(poly: list[int] | None, norm: int) -> int:
        witness_polys.append(poly)
        norms.append(norm)
        return len(witness_polys) - 1

    x_bits = data.get("x_bits", [0] * handle_bits)
    x_idx = [append_slot(const_poly(n, b) if have_witness else None, 0) for b in x_bits]
    tag_idx = append_slot(w["tag"] if have_witness else None, TAG_NORM_SQ)

    def append_decomposed(poly: list[int]) -> list[tuple[int, int]]:
        terms: list[tuple[int, int]] = []
        scale = 1
        for chunk in signed_decompose_poly(poly):
            idx = append_slot(chunk if have_witness else None, SMALL_CHUNK_NORM_SQ)
            terms.append((scale, idx))
            scale *= LINEAR_DECOMP_BASE
        return terms

    zero_poly = zero(n)
    left_idx = [append_decomposed(w["left"][i] if have_witness else zero_poly) for i in range(params["d"])]
    v12_idx = [append_decomposed(w["v12"][i] if have_witness else zero_poly) for i in range(params["d"])]
    v2_idx: list[list[int]] = []
    for i in range(trunc_k_l):
        v2_idx.append([
            append_slot(w["v2"][i][row] if have_witness else None, V2_NORM_SQ)
            for row in range(params["d"])
        ])
    right_idx: list[list[list[int]]] = []
    for h in range(handle_bits):
        h_idx = []
        for j in range(gadget_dim):
            h_idx.append([
                append_slot(w["right"][h][j][row] if have_witness else None, RIGHT_NORM_SQ)
                for row in range(params["d"])
            ])
        right_idx.append(h_idx)
    lift_idx = [append_decomposed(lifts[i] if have_witness else zero_poly) for i in range(params["d"])]

    ps = proof_statement([n] * len(witness_polys), [1] * len(witness_polys), norms,
                         num_constraints=params["d"], primesize="38",
                         num_deg0_constraints=handle_bits, zk=False)

    if have_witness:
        for p in witness_polys:
            assert p is not None
            ps.append_witness(poly_t(ring, p))
    for idx in x_idx:
        ps.append_deg0_statement(idx)

    for row in range(params["d"]):
        stat_lin = []
        wit_lin = []
        for scale, idx in left_idx[row]:
            stat_lin.append(poly_t(ring, const_poly(n, scale)))
            wit_lin.append(idx)
        for col in range(params["d"]):
            for scale, idx in v12_idx[col]:
                stat_lin.append(poly_t(ring, scalar_mul(scale, data["Atrap"][row][col])))
                wit_lin.append(idx)
        for i in range(trunc_k_l):
            for col in range(params["d"]):
                stat_lin.append(-poly_t(ring, data["trunc_B"][i][row][col]))
                wit_lin.append(v2_idx[i][col])
        for h in range(handle_bits):
            for j in range(gadget_dim):
                for col in range(params["d"]):
                    stat_lin.append(poly_t(ring, data["Ba"][h][j][row][col]))
                    wit_lin.append(right_idx[h][j][col])
        for scale, idx in lift_idx[row]:
            stat_lin.append(poly_t(ring, const_poly(n, q * scale)))
            wit_lin.append(idx)

        stat_quad: list[int] = []
        wit_quad1: list[int] = []
        wit_quad2: list[int] = []
        bexpi = 196
        for i in range(trunc_k_l):
            stat_quad.append(bexpi)
            wit_quad1.append(tag_idx)
            wit_quad2.append(v2_idx[i][row])
            bexpi *= 14
        for h in range(handle_bits):
            for j in range(gadget_dim):
                stat_quad.append(-(1 << j))
                wit_quad1.append(x_idx[h])
                wit_quad2.append(right_idx[h][j][row])

        ps.append_quadratic_statement(stat_lin, wit_lin, stat_quad, wit_quad1, wit_quad2,
                                      poly_t(ring, data["acc"][row]))

    if debug_hash:
        from lazer import ffi, lib
        st = ps.output_statement()
        h = ffi.new("uint8_t[16]")
        lib.labrador38_py_statement_hash(h, st)
        print("statement hash:", bytes(ffi.buffer(h, 16)).hex())
        if proof_in is None and proof_out is None:
            return True

    if proof_in is not None:
        proof_bytes = proof_in.read_bytes()
        params_from_proof, proof, refs = deserialize_labrador38_proof(proof_bytes)
        ok = pack_verify((params_from_proof, proof), ps.output_statement(), "38")
        _ = refs
        return ok

    if not have_witness:
        raise ValueError("proving requires witness data")

    proof = ps.pack_prove()
    if proof[0] != 0:
        print(f"Labrador pack_prove failed: {proof[0]}")
        return False
    proof_bytes = serialize_labrador38_proof(proof[1], proof[2])
    if proof_out is not None:
        proof_out.write_bytes(proof_bytes)
        print(f"Serialized Labrador proof bytes: {len(proof_bytes)}")
    lib.labrador38_py_comkey_reset()
    roundtrip_params, roundtrip, refs = deserialize_labrador38_proof(proof_bytes)
    ok = pack_verify((roundtrip_params, roundtrip), ps.output_statement(), "38")
    _ = refs
    return ok


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("json_path")
    parser.add_argument("--mode", choices=["check", "labrador-simple", "prove", "verify", "statement-hash"], default="check")
    parser.add_argument("--proof-out", type=pathlib.Path)
    parser.add_argument("--proof-in", type=pathlib.Path)
    parser.add_argument("--debug-hash", action="store_true")
    args = parser.parse_args()

    with open(args.json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    if args.mode == "verify":
        if args.proof_in is None:
            raise ValueError("--proof-in is required in verify mode")
        ok = build_labrador(data, None, proof_in=args.proof_in, debug_hash=args.debug_hash)
        print(f"Labrador proof-byte verify: {'passed' if ok else 'failed'}")
        return 0 if ok else 1
    if args.mode == "statement-hash":
        if "witness" in data:
            _, lifts = compute_lift(data)
        else:
            lifts = None
        ok = build_labrador(data, lifts, debug_hash=True)
        return 0 if ok else 1

    _, lifts = compute_lift(data)
    max_lift = max(abs(c) for row in lifts for c in row)
    print("real AccGenpp MemVer import")
    print("===========================")
    print(f"input json      : {args.json_path}")
    if "x" in data:
        print(f"x               : 0x{data['x']:08x}")
    print(f"q               : {data['params']['q']}")
    print(f"lift polynomials: {len(lifts)}")
    print(f"max |lift coeff|: {max_lift}")
    print("integer q-lifting check: passed")

    if args.mode in {"labrador-simple", "prove"}:
        if args.mode == "prove" and args.proof_out is None:
            raise ValueError("--proof-out is required in prove mode")
        ok = build_labrador(data, lifts, proof_out=args.proof_out, debug_hash=args.debug_hash)
        print(f"Labrador simple import/proof: {'passed' if ok else 'failed'}")
        return 0 if ok else 1
    print("Labrador import skipped; rerun with --mode labrador-simple under SDE.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
