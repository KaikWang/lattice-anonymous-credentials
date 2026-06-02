#ifndef TSAMPLER_SAMPLING_H
#define TSAMPLER_SAMPLING_H

#include <stdint.h>
#include "params.h"
#include "tsampler_params.h"
#include "arith.h"

/*************************************************
* TSampler-style base-decomposition gadget sampler.
*
* Replaces SEP's Klein/ZF sampler.  Given target w,
* finds y[PARAM_K] such that:
*   [I | b*I | b²*I | ... | b^{K-1}*I] * y = w
*
* Uses coefficient-wise base-B decomposition with
* discrete Gaussian coset rounding.
**************************************************/
void tsampler_gadget_sample_klein_replacement(
    poly_q_vec_d y[PARAM_K],
    const poly_q_vec_d w);

#endif /* TSAMPLER_SAMPLING_H */
