#ifndef SEP_TSAMPLER_H
#define SEP_TSAMPLER_H

#include <stdint.h>
#include "params.h"
#include "sep.h"  /* reuse sep_pk_t, sep_sk_t, sep_sig_t */

/*************************************************
* SEP signature scheme with TSampler-style
* gadget sampler (base decomposition instead
* of Klein/ZF sampling).
*
* Keygen and verify are identical to original SEP.
* Only the sign function's internal gadget sampler
* is replaced.
**************************************************/

/*************************************************
* Name:        sep_sign_tsampler
*
* Description: Sign using SEP with TSampler-style
*   base-decomposition gadget sampler.
*   Same interface as sep_sign().
**************************************************/
void sep_sign_tsampler(
    sep_sig_t      *sig,
    uint8_t         state[STATE_BYTES],
    const sep_sk_t *sk,
    const sep_pk_t *pk,
    const uint8_t   msg[PARAM_M * PARAM_N / 8]);

#endif /* SEP_TSAMPLER_H */
