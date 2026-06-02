#ifndef TSAMPLER_PARAMS_H
#define TSAMPLER_PARAMS_H

/*************************************************
* TSampler-style base-decomposition parameters
* for SEP ablation experiment.
*
* The only new parameter needed is the discrete
* Gaussian width for the coset sampling step:
*   s/B = standard deviation of D_{B*Z+center}
*
* Comparable to SEP's Klein widths (~3.43) and
* TSampler's sL/bL = sH/bH (~3.387).
**************************************************/

/* Gaussian width for base-B coset sampling.
 * Each coefficient of z_j is sampled from
 * D_{B*Z + center, s} where s = S_DIV_B * B.
 * The parameter S_DIV_B is the width passed
 * to SampleZ (i.e., the standard deviation
 * of the integer sample before base-B scaling).
 */
#define PARAM_S_DIV_B_TS 0.5  /* very narrow width for base-B coset sampling */

#endif /* TSAMPLER_PARAMS_H */
