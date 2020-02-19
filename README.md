# PCG-Breaker

This program is able to observe the outputs of a [PCG-XSH-RR][1] pseudorandom number generator with a completely unknown state (both `state` and `inc` are unknown) and predict future outputs with near 100% accuracy. Optionally, given enough outputs, it can perform state recovery of the generator under observation. Note that state recovery will require several dozen million consecutive outputs; this is not a flaw and _any_ state recovery algorithm will need that many outputs to completely pin down the PCG generator's internal state, as there is significant overlap in output sequences for different internal states.

The algorithm basically works by performing an efficient initial filtering of the possible states and then applying basic linear congruential generator recovery techniques. This demonstrates that the PCG-XSH-RR variant is not at all challenging to predict. Note that the algorithm as-is doesn't currently extend to the 128-bit variants.

## Usage

First precompute the lookup table using the `gen-table` binary:

```
$ cargo run --release --bin gen-table
```

This will write a 1GiB file of precomputed data into the current working directory called `table.bin`, which will later be used by the `pcg-breaker` binary.

Then execute the `pcg-breaker` binary with the path to the table, piping/typing in outputs of a PCG generator into its standard input. After being given four outputs, the program will begin predicting the next output in the sequence produced by the generator. A typical session looks like this:

```text
[-] Starting clock.
[+] Loaded precomputed table.
[-] Reading 4 outputs to recover target (zeta, beta) parameters.
[+] Target parameters recovered after 0.50 seconds and 56 guess(es).

     (zeta, beta) = (-28992069, 89358047)

[-] Enumerating internal states...
[+] Enumeration complete.
[+] Output #5 will be one of the following:

    0x3F322013 (78.40% probability)
    0x78B3C9B9 (21.18% probability)
    0x78F3C9B9 (0.42% probability)

[-] Reading output #5 (with value 0x3F322013)
[+] Output #6 will be one of the following:

    0x124D5BF4 (100.00% probability)

[-] Reading output #6 (with value 0x124D5BF4)
[+] Output #7 will be one of the following:

    0x69CC5073 (100.00% probability)
```

If you are only interested in recovering the internal state of a generator rather than predicting future outputs, you can pass the `--recovery` flag to the program. Note that recovering the state requires a very large number of outputs, but it can be used for testing by e.g. piping the output of the provided `pcg.c` program. The `--recovery` flag is mostly a verbosity setting in disguise, as the algorithm always performs output prediction as part of state recovery. The only difference is that in recovery mode, we don't report predictions but simply consume outputs as quickly as possible to recover the PCG state, as we need a large number of outputs and displaying predictions would slow down the process unnecessarily.

After having recovered the generator's internal state, the program will display:

```text
[+] Generator internal state fully recovered after 29.57 seconds:

    pcg32_random_t state = {
        .state = 0xBD094A5E7A8A7587
        .inc   = 0x24E8930796B7B111
    };
```

Note that the `state` displayed will be the state used to produce the _very first output_ given to PCG-breaker; you can advance it yourself if needed. Also note that since the lowest bit of the increment in the PCG state is always masked to 1 (and is therefore irrelevant) the program will conventionally report the recovered increment `inc` with its lowest bit masked to 1 as well.


## Prediction Accuracy

PCG-Breaker will never make a mistake, and almost all of the time is able to identify exactly one future output which will be produced by the generator under observation. Rarely, one of multiple future outputs (usually four or less) may potentially occur; PCG-Breaker will report them all along with their respective probabilities. Once made aware of which output the generator produced, PCG-Breaker will get closer to recovering the internal state of the generator.

## Performance

The algorithm is currently single-threaded, but is able to begin predicting outputs after only four consecutive outputs and less than three seconds of computation on commodity hardware. Full state recovery will occur in less than a minute if enough outputs are available. When used for state recovery the program is currently bottlenecked by the stdin output parsing. It is possible to parallelize several parts of the program though it runs adequately fast already.

Currently there is no feature to skip unknown outputs from a generator, so all outputs must be consecutive.

[1]: https://www.pcg-random.org/download.html
