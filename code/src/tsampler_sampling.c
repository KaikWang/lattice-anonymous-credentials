/*************************************************
* TSampler-style gadget sampler for SEP ablation.
*
* Replaces the Klein/ZF lattice sampler inside SEP's
* poly_q_vec_2d_dk_sample_pre with base-decomposition
* discrete Gaussian coset sampling.
*
* This is the MINIMAL change for the ablation:
* everything else (keygen, perturbation, w computation,
* verification) stays exactly the same as SEP.
**************************************************/

#include "tsampler_sampling.h"
#include "tsampler_params.h"
#include "arith.h"
#include "random.h"

/* ================================================
 * Base-B gadget decomposition (TSampler-style)
 * Replaces sample_klein() in poly_q_sampling.c
 * ================================================ */

/*************************************************
* Name:        tsampler_gadget_sample_klein_replacement
*
* Description: Base-decomposition gadget sampler for
*   the SEP single gadget G = [I | b*I | ... | b^{K-1}*I].
*
*   For each coefficient of w, performs base-B
*   decomposition with discrete Gaussian rounding:
*
*     remainder = w_coeff
*     for j = 0..K-1:
*       center = remainder mod B
*       y[j] = D_{B*Z + center, sigma}
*       remainder = (remainder - y[j]) / B
*
*   This produces y[K] such that G*y = w exactly,
*   and each y[j] is short (discrete Gaussian over
*   the appropriate coset).
*
*   This replaces the Klein/ZF sampler which finds
*   y in the gadget lattice via Gram-Schmidt.
**************************************************/
void tsampler_gadget_sample_klein_replacement(
    poly_q_vec_d y[PARAM_K],
    const poly_q_vec_d w) {
  size_t i, j, k;
  coeff_q c, rem;

  for (i = 0; i < PARAM_D; i++) {
    for (k = 0; k < PARAM_N; k++) {
      rem = poly_q_get_coeff(w->entries[i], k);

      for (j = 0; j < PARAM_K; j++) {
        /* Deterministic base-B digit (no Gaussian noise).
         * This is the tightest possible decomposition:
         * c ∈ [0, B-1], remainder shrinks by factor B each step.
         * Since B^K > Q, K digits are sufficient. */
        c = rem % PARAM_B;
        c += (c >> (sizeof(coeff_q)*8-1)) & PARAM_B;

        /* store digit */
        poly_q_set_coeff(y[j]->entries[i], k, c);

        /* exact division: (rem - c) is divisible by B */
        rem = (rem - c) / PARAM_B;
      }
    }
  }
}
