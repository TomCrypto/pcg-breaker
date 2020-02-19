#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

uint64_t gen_random() {
    uint64_t a = (uint64_t)(rand()) % 65536;
    uint64_t b = ((uint64_t)(rand()) % 65536) << 16;
    uint64_t c = ((uint64_t)(rand()) % 65536) << 32;
    uint64_t d = ((uint64_t)(rand()) % 65536) << 48;

    return a | b | c | d;
}

void main() {
    srand(time(0));

    pcg32_random_t pcg = {
        .state = gen_random(),
        .inc   = gen_random()
    };

    fprintf(stderr, ">> PCG INITIAL STATE = %016" PRIx64 "\n", pcg.state);
    fprintf(stderr, ">> PCG INCREMENT     = %016" PRIx64 "\n", pcg.inc);

    while (1) {
        printf("0x%08x\n", pcg32_random_r(&pcg));
    }
}
