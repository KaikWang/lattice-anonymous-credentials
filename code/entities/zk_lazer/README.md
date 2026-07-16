# LaZer prototype for revocable anonymous show

This directory is a prototype workspace for replacing the current hand-written
`show` proof with a LaZer/Labrador proof language that also covers revocation.

The current C implementation already has:

- credential issuance bound to a 32-bit handle `x`;
- accumulator `Add/Delete/MemWitUp/MemWitSync`;
- public `MemVer(acc, x, s_x)` checks.

It does not yet have a zero-knowledge proof that the same hidden `x` is used in
both the blind-signature relation and the accumulator membership relation.

## Why this is separate from the current C show proof

The current `show` proof has fixed dimensions in `params.h` and fixed vector
types under `arith/arith_qshow`.  Adding `x` and `s_x` changes the proof
language dimension and all norm bounds.  A local patch would create a proof that
compiles but does not have the intended soundness/zero-knowledge guarantees.

## Files

- `revac_show_relation.py`: pure-Python relation/dimension analysis that runs
  without LaZer.  It estimates how large a direct embedding of the current
  accumulator witness would be.
- `revac_show_labrador_skeleton.py`: executable skeleton for the future
  Labrador proof.  It explains the constraints and checks whether LaZer is
  available in the current environment.
- `revac_same_x_demo.py`: real LaZer proof-generation demo.  It proves a small
  relation where the same hidden handle lane `x` is used by both a
  credential/signature block and an accumulator-membership block.
- `revac_memver_quadratic_demo.py`: real Labrador proof-generation demo for
  the core revocation relation shape, including binary degree-0 handle bits and
  quadratic `x_i * right_i` selector terms from `MemVer`.
- `revac_acc_memver_full_layout.py`: dimension-accurate layout of the complete
  current C `AccGenpp` `MemVer` relation that still needs a coefficient
  exporter before it can become the production Labrador proof.
- `Makefile`: convenience targets.

## Commands

```sh
make describe
make check-lazer
make same-x-demo
make memver-quadratic-demo
make acc-memver-layout
```

The Makefile defaults to `/home/wkk/miniconda3/envs/sage/bin/python`, because
the built LaZer CFFI module in `/home/wkk/pq_ac/lazer/python` targets Python
3.11.  The quadratic Labrador target defaults to the local Intel SDE wrapper
because this machine does not expose native AVX-512:

```sh
make memver-quadratic-demo
make SDE= memver-quadratic-demo   # only on native AVX-512 hardware
```

`revac_same_x_demo.py` is the first runnable proof backend check.  It is not yet
the final revocable anonymous credential proof: the two linear blocks are still
stand-ins that must be replaced by the actual `Dx*x` signature relation and the
actual accumulator `MemVer(acc, x, s_x)` relation.

`revac_memver_quadratic_demo.py` is the next step: it proves the same hidden
binary handle bits in both a signature-side linear constraint and an
accumulator-side quadratic constraint of the form used by `MemVer`.  It still
uses toy public coefficients; the remaining work is to emit the full AccGenpp
matrix coefficients and modulus-lifting constraints.

On this build, the zero-knowledge packing works with `--primesize 38`.  The
smaller 32/36-bit Labrador parameter families verify the relation in non-ZK
debug mode, but do not have matching LNP parameters for this ZK statement.
